#!/bin/bash -l

OUTDIR=$1
PREFIX=$2
READS=$3

FEATURES_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/features.fasta
PLAM_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/pLAM.fasta
SEQ_4QB=/g/data/kr68/neysa/fshd_pipeline/inputs/4qB_rc.fa

# module load minimap2/2.22
# module load samtools/1.19
# module load bedtools/2.28.0

TEMP_BAM=$OUTDIR/tmp_bams
mkdir -p $TEMP_BAM

TMP_FA=$OUTDIR/tmp_read.fa
awk '/^>/ {if (seq) print seq > "${TMP_FA}"; print > "${TMP_FA}"; seq=""; next} {seq = seq $0} END {if (seq) print seq > "${TMP_FA}"}' "$READS"

while read -r header; do
    read_id=$(echo "$header" | cut -c2-)  # strip the leading ">"
    
    # Extract the single read into a temporary FASTA
    awk -v h="$header" 'BEGIN {p=0} $0==h {p=1; print; next} /^>/ {p=0} p {print}' "$READS" > ${TMP_FA}
    
    # Run minimap2 using the read as the reference
    minimap2 -x asm5 -t 1 -a ${TMP_FA} "$FEATURES_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}.bam"

    minimap2 -x asm5 -t 1 -a ${TMP_FA} "$SEQ_4QB" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}_4qB.bam"

    # Run minimap2 short for pLAM
    minimap2 -x sr -t 1 -a ${TMP_FA} "$PLAM_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}_short.bam"

    samtools merge -f "${TEMP_BAM}/${read_id}_merged.bam" \
        "${TEMP_BAM}/${read_id}.bam" \
        "${TEMP_BAM}/${read_id}_4qB.bam" \
        "${TEMP_BAM}/${read_id}_short.bam"

    rm "${TEMP_BAM}/${read_id}.bam"
    rm "${TEMP_BAM}/${read_id}_4qB.bam"
    rm "${TEMP_BAM}/${read_id}_short.bam"
    
    rm ${TMP_FA}
done < <(grep '^>' "$READS")

# Merge all into a single BAM
samtools merge -f "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" ${TEMP_BAM}/*.bam
samtools index "${OUTDIR}/${PREFIX}_mapped_features.merged.bam"

bedtools bamtobed -i "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" > "${OUTDIR}/${PREFIX}_mapped_features.merged.bed"

# convert bam -> sam -> psl
samtools view -h "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" > "${OUTDIR}/${PREFIX}_mapped_features.merged.sam"

python3 /g/data/kr68/neysa/fshd_pipeline/helper/sam2psl.py \
    -i "${OUTDIR}/${PREFIX}_mapped_features.merged.sam" \
    -o "${OUTDIR}/${PREFIX}_mapped_features.merged.psl"

cat "${OUTDIR}/${PREFIX}_mapped_features.merged.psl" > "${OUTDIR}/${PREFIX}_mapped_features.psl"

# Optional: clean up
rm -r $TEMP_BAM
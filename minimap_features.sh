#!/bin/bash -l

OUTDIR=$1
PREFIX=$2
READS=$3

REALPATH=$(dirname "$(readlink -f "$0")")

FEATURES_FASTA=${REALPATH}/inputs/features.fasta
PLAM_FASTA=${REALPATH}/inputs/short_features.fasta
SEQ_4QB=${REALPATH}/inputs/4qB_marker.fa

# module load minimap2/2.22
# module load samtools/1.19
# module load bedtools/2.28.0

TEMP_BAM=$OUTDIR/tmp_bams
mkdir -p "$TEMP_BAM"

TMP_FA=$OUTDIR/tmp_read.fa
awk -v out="$TMP_FA" '
    /^>/ {
        if (seq) print seq >> out;
        print >> out;
        seq = "";
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq) print seq >> out
    }
' "$READS"

while read -r header; do
    read_id=$(echo "$header" | cut -c2-)  # strip the leading ">"

    # Extract the single read into a temporary FASTA
    TMP_READ_FA="${OUTDIR}/tmp_read_${read_id}.fa"
    # awk -v h="$header" 'BEGIN {p=0} $0==h {p=1; print; next} /^>/ {p=0} p {print}' "$READS" > ${TMP_FA}
    awk -v h="$header" 'BEGIN {p=0} $0==h {p=1; print; next} /^>/ {p=0} p {print}' "$TMP_FA" > ${TMP_READ_FA}

    # Run minimap2 using the read as the reference
    minimap2 -x asm5 -t 1 -N 10 -a ${TMP_READ_FA} "$FEATURES_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}.bam"

    minimap2 -x asm5 -t 1 -a ${TMP_READ_FA} "$SEQ_4QB" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}_4qB.bam"

    # Run minimap2 short for pLAM
    minimap2 -x sr -t 1 -N 10 -a ${TMP_READ_FA} "$PLAM_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}_short.bam"

    samtools merge -f "${TEMP_BAM}/${read_id}_merged.bam" \
        "${TEMP_BAM}/${read_id}.bam" \
        "${TEMP_BAM}/${read_id}_4qB.bam" \
        "${TEMP_BAM}/${read_id}_short.bam"

    rm "${TEMP_BAM}/${read_id}.bam"
    rm "${TEMP_BAM}/${read_id}_4qB.bam"
    rm "${TEMP_BAM}/${read_id}_short.bam"

    rm ${TMP_READ_FA}
done < <(grep '^>' "$READS")

# Merge all into a single BAM
samtools merge -f "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" ${TEMP_BAM}/*.bam
samtools index "${OUTDIR}/${PREFIX}_mapped_features.merged.bam"

bedtools bamtobed -i "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" > "${OUTDIR}/${PREFIX}_mapped_features.merged.bed"

# convert bam -> sam -> psl
samtools view -h "${OUTDIR}/${PREFIX}_mapped_features.merged.bam" > "${OUTDIR}/${PREFIX}_mapped_features.merged.sam"

python3 ${REALPATH}/sam2psl.py \
    -i "${OUTDIR}/${PREFIX}_mapped_features.merged.sam" \
    -o "${OUTDIR}/${PREFIX}_mapped_features.merged.psl"

cat "${OUTDIR}/${PREFIX}_mapped_features.merged.psl" > "${OUTDIR}/${PREFIX}_mapped_features.psl"

# Optional: clean up
rm $TMP_FA
rm -r $TEMP_BAM

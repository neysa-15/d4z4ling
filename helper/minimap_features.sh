#!/bin/bash -l

OUTDIR=$1
PREFIX=$2
READS=$3

# module unload samtools
# module unload bedtools
# module unload blat
# module unload kentutils
# module unload minimap2

# module load samtools/1.19
# module load bedtools/2.28.0
# module load blat/37
# module load kentutils/0.0
# module load minimap2/2.22

PSLTOBED=/g/data/if89/apps/kentutils/0.0/bin/pslToBed

FEATURES_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/features.fasta
SHORT_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/short_features.fasta
PLAM_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/pLAM.fasta

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
    awk -v h="$header" 'BEGIN {p=0} $0==h {p=1; print; next} /^>/ {p=0} p {print}' "$TMP_FA" > ${TMP_READ_FA}
    
    # Run minimap2 using the read as the reference
    minimap2 -x asm5 -t 1 -N 10 -a ${TMP_READ_FA} "$FEATURES_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}.bam"

    # Run minimap2 short for p13-E11
    minimap2 -x sr -t 1 -N 10 -a ${TMP_READ_FA} "$SHORT_FASTA" | \
        samtools view -bS - | \
        samtools sort -o "${TEMP_BAM}/${read_id}_short.bam"

    samtools merge -f "${TEMP_BAM}/${read_id}_merged.bam" \
        "${TEMP_BAM}/${read_id}.bam" \
        "${TEMP_BAM}/${read_id}_short.bam"

    rm "${TEMP_BAM}/${read_id}.bam"
    rm "${TEMP_BAM}/${read_id}_short.bam"
    
    rm ${TMP_READ_FA}
done < <(grep '^>' "$READS")

# Merge all into a single BAM
samtools merge -f "${OUTDIR}/${PREFIX}_mapped_features.bam" ${TEMP_BAM}/*.bam
samtools index "${OUTDIR}/${PREFIX}_mapped_features.bam"

bedtools bamtobed -i "${OUTDIR}/${PREFIX}_mapped_features.bam" > "${OUTDIR}/${PREFIX}_mapped_features.bed"

# convert bam -> sam -> psl
samtools view -h "${OUTDIR}/${PREFIX}_mapped_features.bam" > "${OUTDIR}/${PREFIX}_mapped_features.merged.sam"

python3 /g/data/kr68/neysa/fshd_pipeline/helper/sam2psl.py \
    -i "${OUTDIR}/${PREFIX}_mapped_features.merged.sam" \
    -o "${OUTDIR}/${PREFIX}_mapped_features.merged.psl"

cat "${OUTDIR}/${PREFIX}_mapped_features.merged.psl" > "${OUTDIR}/${PREFIX}_mapped_features.psl"

# USE BLAT FOR PLAM DUE TO SIZE LIMITATION
blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 "$READS" "${PLAM_FASTA}" "${OUTDIR}/${PREFIX}_pLAM_matches.psl"

awk 'BEGIN{FS=OFS="\t"} 
     $1 ~ /^[0-9]+$/ { 
        split($19,a,","); 
        aln_size = $1 + $2 + $3; 
        if (aln_size > 0 && a[1] >= 100 && a[1]/aln_size >= 0.5) print 
     }' "${OUTDIR}/${PREFIX}_pLAM_matches.psl" \
> "${OUTDIR}/${PREFIX}_pLAM_matches_filtered.psl"

# Add to mapped features psl
cat "${OUTDIR}/${PREFIX}_pLAM_matches_filtered.psl" >> "${OUTDIR}/${PREFIX}_mapped_features.psl"

# Convert psl to bed
${PSLTOBED} "${OUTDIR}/${PREFIX}_pLAM_matches_filtered.psl" "${OUTDIR}/${PREFIX}_pLAM_matches.bed"

# merge with the rest of the features
cat "${OUTDIR}/${PREFIX}_pLAM_matches.bed" >> "${OUTDIR}/${PREFIX}_mapped_features.bed"

# Optional: clean up
rm $TMP_FA
rm -r $TEMP_BAM
rm ${OUTDIR}/${PREFIX}_pLAM_matches*
rm "${OUTDIR}/${PREFIX}_mapped_features.merged.psl"
rm "${OUTDIR}/${PREFIX}_mapped_features.merged.sam"
#!/bin/bash

INPUT_FASTA=""
REF=/g/data/kr68/genome/hs1.fa           # Reference genome
FEATURES_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/features.fasta
REGION_BED=/g/data/kr68/neysa/fshd_pipeline/inputs/d4z4_region.chm13.bed
PLAM_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/pLAM.fasta
SEQ_4QB=/g/data/kr68/neysa/fshd_pipeline/inputs/4qB_down.fa
PROBES=/g/data/kr68/neysa/fshd_pipeline/inputs/probes.fasta
PREFIX="SAMPLE"
OUTDIR="$PREFIX"
HAPLOTYPE_REFS=/g/data/kr68/neysa/fshd_pipeline/inputs/d4z4_repeats.fasta  # Fasta file containing haplotype-specific references
REPEATS_FASTA=/g/data/kr68/neysa/fshd_pipeline/inputs/dux4.gene_complete_genbank_20241127.reformatted.fasta

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix) PREFIX="$2"; OUTDIR="$PREFIX"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --input-fasta) INPUT_FASTA="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --region-bed) REGION_BED="$2"; shift 2;;
        --features-fasta) FEATURES_FASTA="$2"; shift 2;;
        --probes) PROBES="$2"; shift 2;;
        --haplotype-refs) HAPLOTYPE_REFS="$2"; shift 2;;
        --repeats-fasta) REPEATS_FASTA="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Paths to conversion tools
PSLTOBED=/g/data/if89/apps/kentutils/0.0/bin/pslToBed
BG2BW=/g/data/if89/apps/kentutils/0.0/bin/bedGraphToBigWig

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Winnowmap
echo "MERYL"
# Define meryl outputs relative to the reference genome directory
REF_NAME=$(basename $REF .fa)
MERYL_DB="/g/data/kr68/neysa/fshd_pipeline/inputs/merylDB_${REF_NAME}"
REPETITIVE_REGIONS="/g/data/kr68/neysa/fshd_pipeline/inputs/${REF_NAME}_repetitive_k15.txt"

# Step 1: Pre-compute k-mer frequency (only if meryl output does not exist)
if [ ! -d "$MERYL_DB" ]; then
    echo "Creating meryl database..."
    meryl count k=15 output "$MERYL_DB" "$REF"
else
    echo "Meryl database already exists at $MERYL_DB. Skipping this step."
fi

# Step 2: Extract repetitive regions (only if {refname}_repetitive_k15.txt does not exist)
if [ ! -f "$REPETITIVE_REGIONS" ]; then
    echo "Extracting repetitive regions..."
    meryl print greater-than distinct=0.9998 "$MERYL_DB" > "$REPETITIVE_REGIONS"
else
    echo "Repetitive regions file already exists at $REPETITIVE_REGIONS. Skipping this step."
fi

INPUT_FASTQ="${OUTDIR}/${PREFIX}_reads_of_interest.fastq"
seqtk seq -F 'I' "${INPUT_FASTA}" > "${INPUT_FASTQ}"

# Run winnowmap alignment
echo "WINNOWMAP" 
winnowmap -W ${REPETITIVE_REGIONS} -Y -y -ax $MODE "$REF" "$INPUT_FASTQ" | \
        samtools view -L "$REGION_BED" -q ${MAPQ} -Sb | \
        samtools sort -o "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

samtools index "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

samtools fasta "${OUTDIR}/${PREFIX}_reads_of_interest.bam" > "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"
samtools faidx "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

# Convert BAM to BED
echo "Converting BAM to BED"
bedtools bamtobed -i "${OUTDIR}/${PREFIX}_reads_of_interest.bam" > "${OUTDIR}/${PREFIX}_reads_of_interest.bed"

# map features
/g/data/kr68/neysa/fshd_pipeline/helper/minimap_features.sh "$OUTDIR" "$PREFIX" "$INPUT_FASTA"

echo "Converting PSL to BED"
 ${PSLTOBED} "${OUTDIR}/${PREFIX}_mapped_features.psl" "${OUTDIR}/${PREFIX}_mapped_features.bed"

# Map SSLP
seqkit amplicon --bed -F GGTGGAGTTCTGGTTTCAGC -R CCTGTGCTTCAGAGGCATTTG -m 2 "$INPUT_FASTA" > "${OUTDIR}/${PREFIX}_SSLP.bed"

# Step 9: Align reads to haplotype-specific references
echo "Aligning reads classified to specific haplotypes"
aligned_bams=()  # Array to collect BAM file paths for merging

for chrom in chr4 chr10; do

    if [[ "$chrom" == "chr4" ]]; then
        haplotype="4A"
    else
        haplotype="10A"
    fi
    # Extract haplotype-specific reference
    haplotype_ref="${OUTDIR}/${haplotype}_reference.fasta"
    grep -A1 "${haplotype}" "$HAPLOTYPE_REFS" > "$haplotype_ref"

    # Extract reads corresponding to the haplotype
    haplotype_reads="${OUTDIR}/${PREFIX}_${haplotype}_reads.fastq"
    # seqtk subseq <(samtools fastq "$INPUT_UBAM") <(grep ${chrom} "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" | cut -f1) > "$haplotype_reads"
    seqtk subseq "$INPUT_FASTQ" <(grep ${chrom} "${OUTDIR}/${PREFIX}_reads_of_interest.bed" | cut -f4 | sort | uniq) > "$haplotype_reads"

    # Align reads to the haplotype-specific reference
    aligned_bam="${OUTDIR}/${PREFIX}_${haplotype}_aligned.bam"
    minimap2 --secondary=no --MD -ax map-ont "$haplotype_ref" "$haplotype_reads" | samtools sort -o "$aligned_bam"
    samtools index "$aligned_bam"

    # Add the BAM file to the array
    aligned_bams+=("$aligned_bam")

    # Remove intermediate files
    rm -f "$haplotype_ref" "$haplotype_reads"
done
echo "Haplotype-specific alignment complete."

## Step 10: Merge haplotype-aligned BAM files
echo "Merging haplotype-aligned BAM files"
merged_bam="${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"

# Use samtools to merge all BAM files from the array
samtools merge -f "$merged_bam" "${aligned_bams[@]}"
samtools index "$merged_bam"

# Remove individual haplotype BAM files and their indexes
rm "${aligned_bams[@]}" "${aligned_bams[@]/%.bam/.bam.bai}"

# Step 8: Parse PSL to generate a summary table
echo "Parsing PSL to generate summary table"
python3 /g/data/kr68/neysa/fshd_pipeline/helper/read_classification.py \
    --psl "${OUTDIR}/${PREFIX}_mapped_features.psl" \
    --bed "${OUTDIR}/${PREFIX}_reads_of_interest.bed" \
    --output "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --fasta "$INPUT_FASTA" \
    --sslp "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --aligned-bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"

## Step 12: Generate ordered alignment sequences
echo "Generating ordered alignment sequences"
alignment_script="/g/data/kr68/neysa/fshd_pipeline/helper/alignment_order_clipping.py"

# Run the alignment order script
python3 "$alignment_script" \
    --bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam" \
    --fasta "$INPUT_FASTA" \
    --output_table "${OUTDIR}/${PREFIX}_d4z4_units.tsv" \
    --xapi_bed "${OUTDIR}/${PREFIX}_xapi_sites.bed" \
    --blni_bed "${OUTDIR}/${PREFIX}_blni_sites.bed" \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" 

# sort "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" | uniq | awk -v OFS='\t' '
# {
#     if ($1 != prev) {
#         count=1
#         prev=$1
#     } else {
#         count++
#     }
#     $4=$4""count
#     print
# }' > "${OUTDIR}/${PREFIX}_d4z4_repeats_annotated.bed"

# mv "${OUTDIR}/${PREFIX}_d4z4_repeats_annotated.bed" "${OUTDIR}/${PREFIX}_d4z4_repeats.bed"

# Unify features and repeats into a single bed and add colour
python3 /g/data/kr68/neysa/fshd_pipeline/helper/add_colors_to_bed.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --features_bed "${OUTDIR}/${PREFIX}_mapped_features.bed" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" \
    --sslp_bed "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --output_bed "${OUTDIR}/${PREFIX}_all_features.bed"


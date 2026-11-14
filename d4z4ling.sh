#!/bin/bash

set -e # <-- abort on any error

REALPATH=$(dirname "$(readlink -f "$0")")

# Default inputs
INPUT_UBAM=""                            # Unaligned BAM
INPUT_FASTQ=""
REF=${REALPATH}/genome/hs1.fa           # Reference genome
REGION_BED=${REALPATH}/inputs/d4z4_region.chm13.bed
FEATURES_FASTA=${REALPATH}/inputs/features.fasta
SHORT_FEATURES=${REALPATH}/inputs/short_features.fasta
PLAM=${REALPATH}/inputs/pLAM.fasta
PREFIX="SAMPLE"
OUTDIR="$PREFIX"
HAPLOTYPE_REFS=${REALPATH}/inputs/d4z4_repeats.fasta  # Fasta file containing haplotype-specific references
MMI=${REALPATH}/genome/hs1.mmi
REMOVE_INTERMEDIATE_FILES=false
FSHD1_STATUS=${REALPATH}/inputs/fshd1_status_template.tsv
THREADS=8

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix) PREFIX="$2"; OUTDIR="$PREFIX"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --input-ubam) INPUT_UBAM="$2"; shift 2;;
	    --input-fastq) INPUT_FASTQ="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --region-bed) REGION_BED="$2"; shift 2;;
        --features-fasta) FEATURES_FASTA="$2"; shift 2;;
        --short-features-fasta) SHORT_FEATURES="$2"; shift 2;;
        --plam) PLAM="$2"; shift 2;;
        --haplotype-refs) HAPLOTYPE_REFS="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --remove-intermediate-files) REMOVE_INTERMEDIATE_FILES="$2"; shift 2;;
        --minimap-ref-index) MMI="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Check if required arguments are set
if [[ (-z "$INPUT_UBAM" || ! -f "$INPUT_UBAM") && (-z "$INPUT_FASTQ" || ! -f "$INPUT_FASTQ") ]]; then
    echo "Error: Either --input-ubam or --input-fastq must be provided and must be a valid file."
    exit 1
fi

# Paths to conversion tools
[ -z ${PSLTOBED} ] && PSLTOBED=pslToBed
[ -z ${BG2BW} ] && BG2BW=bedGraphToBigWig

[ -z ${SAMTOOLS} ] && SAMTOOLS=samtools
[ -z ${MINIMAP2} ] && MINIMAP2=minimap2
[ -z ${FSHD_PYTHON} ] && FSHD_PYTHON=python3

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Keep the original file
# if have ubam, convert to fastq
# if have fastq, convert to ubam
echo "Convert ubam/fastq to fastq/ubam"
ubam_exist=false
fastq_exist=false
if [[  ( ! -z "$INPUT_UBAM" ) && ( -f "$INPUT_UBAM") ]]; then
    ubam_exist=true
    # INPUT_FASTQ=basename "$INPUT_UBAM" .bam
    UBAM_BASE=$(basename $INPUT_UBAM)
    $SAMTOOLS fastq -TMM,ML "$INPUT_UBAM" > "${OUTDIR}/${UBAM_BASE%.bam}.fastq"
    # samtools fastq -TMM,ML "$INPUT_UBAM" > "${INPUT_UBAM%.bam}.fastq"
    INPUT_FASTQ="${OUTDIR}/${UBAM_BASE%.bam}.fastq"
elif [[  ( ! -z "$INPUT_FASTQ" ) && ( -f "$INPUT_FASTQ") ]]; then
    fastq_exist=true
    FASTQ_BASE=$(basename $INPUT_FASTQ)
    $SAMTOOLS import -u -s "$INPUT_FASTQ" -o "${OUTDIR}/${FASTQ_BASE%.fastq}.bam"
    INPUT_UBAM="${OUTDIR}/${FASTQ_BASE%.fastq}.bam"
fi

echo "MERYL"
# Define meryl outputs relative to the reference genome directory
REF_NAME=$(basename $REF .fa)
MERYL_DB="${REALPATH}/inputs/merylDB_${REF_NAME}"
REPETITIVE_REGIONS="${REALPATH}/inputs/${REF_NAME}_repetitive_k15.txt"

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

# Run winnowmap alignment
# echo "WINNOWMAP"
# winnowmap -W ${REPETITIVE_REGIONS} -Y -y -ax $MODE "$REF" "$INPUT_FASTQ" -t ${PBS_NCPUS:-8} | \
#         $SAMTOOLS view -L "$REGION_BED" -Sb | \
#         $SAMTOOLS sort -o "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

echo "minimap winnowmap"
# MMI="${OUTDIR}/${REF_NAME}.mmi"
# $MINIMAP2 -x $MODE "$REF" "$INPUT_FASTQ" -t ${PBS_NCPUS:-8} | \
#  -r 500,100k --rmq=yes -f 0.001
$MINIMAP2 -x $MODE -a "$MMI" "$INPUT_FASTQ" -t ${PBS_NCPUS:-8} | \
        $SAMTOOLS view -L "$REGION_BED" -Sb | \
        $SAMTOOLS sort -o "${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam"

$SAMTOOLS index "${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam"

if [[ $($SAMTOOLS view -c "${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam") -eq 0 ]]; then
    echo "No reads mapped to the region of interest. Exiting."

    # Remove temp created fastq/ubam
    if [[ "$ubam_exist" == true ]]; then
        rm "${INPUT_FASTQ}"
    elif [[ "$fastq_exist" == true ]]; then
        rm "${INPUT_UBAM}"
    fi

    exit 0
fi

echo "Remapping with winnowmap to CHM13"
# Remap with winnowmap to CHM13
$SAMTOOLS view -N <($SAMTOOLS view "${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam" | cut -f1 | sort | uniq) -b ${INPUT_UBAM} | \
    $SAMTOOLS fastq -TMM,ML - | \
    winnowmap -W ${REPETITIVE_REGIONS} -Y -y -ax $MODE "$REF" - -t ${PBS_NCPUS:-8} | \
    $SAMTOOLS view -L "$REGION_BED" -Sb | \
    $SAMTOOLS sort -o "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

$SAMTOOLS index "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

# Generate FASTA from uBAM
echo "Generating FASTA from uBAM"
$SAMTOOLS view -N <($SAMTOOLS view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | cut -f1 | sort | uniq) -b "$INPUT_UBAM" | \
	$SAMTOOLS fasta > "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

# Convert BAM to BED
echo "Converting BAM to BED"
bedtools bamtobed -i "${OUTDIR}/${PREFIX}_reads_of_interest.bam" > "${OUTDIR}/${PREFIX}_reads_of_interest.bed"

#######################################
###    MAPPING FEATURES TO READS    ###
#######################################

# map features
${REALPATH}/helper/minimap_features.sh "$OUTDIR" "$PREFIX" "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "$FEATURES_FASTA" "$SHORT_FEATURES" "$PLAM"

# Map SSLP
seqkit amplicon --bed -F GGTGGAGTTCTGGTTTCAGC -R CCTGTGCTTCAGAGGCATTTG -m 2 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" > "${OUTDIR}/${PREFIX}_SSLP.bed"

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
    $MINIMAP2 --secondary=no --MD -ax map-ont "$haplotype_ref" "$haplotype_reads" | $SAMTOOLS sort -o "$aligned_bam"
    $SAMTOOLS index "$aligned_bam"

    # Add the BAM file to the array
    aligned_bams+=("$aligned_bam")

    # Remove intermediate files
    rm -f "$haplotype_ref" "$haplotype_reads"
done
echo "Haplotype-specific alignment complete."

## Step 10: Merge haplotype-aligned BAM files
echo "Merging haplotype-aligned BAM files"
merged_bam="${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"

# Use $SAMTOOLS to merge all BAM files from the array
$SAMTOOLS merge -f "$merged_bam" "${aligned_bams[@]}"
$SAMTOOLS index "$merged_bam"

# Remove individual haplotype BAM files and their indexes
rm "${aligned_bams[@]}" "${aligned_bams[@]/%.bam/.bam.bai}"

# Step 8: Parse PSL to generate a summary table
echo "Parsing PSL to generate summary table"
${FSHD_PYTHON} ${REALPATH}/helper/read_classification.py \
    --psl "${OUTDIR}/${PREFIX}_mapped_features.psl" \
    --bed "${OUTDIR}/${PREFIX}_reads_of_interest.bed" \
    --output "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" \
    --sslp "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --aligned-bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"

## Step 12: Generate ordered alignment sequences
echo "Generating ordered alignment sequences"
alignment_script="${REALPATH}/helper/alignment_order_clipping.py"

# Run the alignment order script
${FSHD_PYTHON} "$alignment_script" \
    --bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" \
    --output_table "${OUTDIR}/${PREFIX}_d4z4_units.tsv" \
    --xapi_bed "${OUTDIR}/${PREFIX}_xapi_sites.bed" \
    --blni_bed "${OUTDIR}/${PREFIX}_blni_sites.bed" \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed"

# Get potential hybrid reads by looking for supplementary alignments
$SAMTOOLS view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | awk '{
    primary_chr = $3
    for (i=12; i<=NF; i++) {
        if ($i ~ /^SA:Z:/) {
            sub(/^SA:Z:/, "", $i)
            n = split($i, sa_alignments, ";")
            for (j=1; j<=n; j++) {
                if (sa_alignments[j] == "") continue
                split(sa_alignments[j], sa_fields, ",")
                sa_chr = sa_fields[1]
                mapq = sa_fields[5]
                if (primary_chr != sa_chr && mapq >= 30) {
                    print $1, primary_chr, sa_alignments[j]
                }
            }
        }
    }
}' > "${OUTDIR}/${PREFIX}_potential_hybrid.tsv"

# Flag for gaps and overlaps
${FSHD_PYTHON} ${REALPATH}/helper/flag_repeats.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" \
    --hybrid "${OUTDIR}/${PREFIX}_potential_hybrid.tsv" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

# Unify features and repeats into a single bed and add colour
${FSHD_PYTHON} ${REALPATH}/helper/add_colors_to_bed.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --features_bed "${OUTDIR}/${PREFIX}_mapped_features.bed" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" \
    --sslp_bed "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --output_bed "${OUTDIR}/${PREFIX}_all_features.bed"

# reannotate distal haplotype when identified as the 4qA long haplotype (DUX4L)
${REALPATH}/helper/distal_haplotype_blast.sh "${OUTDIR}" "${PREFIX}"

# Check for overlapping repeats after distal haplotype re-annotation
ANNOTATED_BED="${OUTDIR}/${PREFIX}_all_features.bed"
python3 -c "
from helper import alignment_order_clipping
alignment_order_clipping.filter_d4z4_bed(
    '${ANNOTATED_BED}'
)
" 2>&1

$SAMTOOLS view -N <($SAMTOOLS view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | cut -f1 | sort | uniq) -b ${INPUT_UBAM} | \
    $SAMTOOLS fastq -TMM,ML - | $MINIMAP2 -t ${THREADS} -Y -y -x asm5 -a --secondary=no "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" - | \
    $SAMTOOLS view -F 2308 -b - | \
    $SAMTOOLS sort -@ ${THREADS} - > "${OUTDIR}/${PREFIX}_meth_reads.bam"

# Index the resulting BAM file
$SAMTOOLS index "${OUTDIR}/${PREFIX}_meth_reads.bam"

# use minimod to get meth probabilities for each CpG within a given read
minimod view -c m "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "${OUTDIR}/${PREFIX}_meth_reads.bam" | awk 'NR > 1 {print $4"\t"$5"\t"$5+1"\t.\t"$7}' > "${OUTDIR}/${PREFIX}_meth_reads.bed"

fshd1_read_status_tsv="${OUTDIR}/${PREFIX}_fshd1_status_counts.tsv"
${FSHD_PYTHON} ${REALPATH}/helper/methylation_summary.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --features_bed "${OUTDIR}/${PREFIX}_all_features.bed" \
    --meth "${OUTDIR}/${PREFIX}_meth_reads.bed" \
    --output "${OUTDIR}/${PREFIX}_methylation_summary.tsv" \
    --fshd1_read_status_tsv $fshd1_read_status_tsv \
    --status_template_tsv $FSHD1_STATUS \
    --updated_bed "${OUTDIR}/${PREFIX}_updated_features.bed"

mv "${OUTDIR}/${PREFIX}_methylation_summary.tsv" "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"
mv "${OUTDIR}/${PREFIX}_updated_features.bed" "${OUTDIR}/${PREFIX}_all_features.bed"

if [ ! -s "${OUTDIR}/${PREFIX}_meth_reads.bed" ]; then
    echo "No methylation calls found. Exiting."

    # Remove intermediate files
    if [[ "$REMOVE_INTERMEDIATE_FILES" == true ]]; then
        rm ${OUTDIR}/${PREFIX}_*blast*
        rm ${OUTDIR}/${PREFIX}_*top*
        rm ${OUTDIR}/${PREFIX}_*psl
        rm ${OUTDIR}/${PREFIX}_*d4z4*
    fi

    # Remove temp created fastq/ubam
    if [[ "$ubam_exist" == true ]]; then
        rm "${INPUT_FASTQ}"
    elif [[ "$fastq_exist" == true ]]; then
        rm "${INPUT_UBAM}"
    fi

    exit 0
fi

# Check if genome.chrom.sizes exists
if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" ]; then
    echo "Creating ${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes in $OUTDIR..."

    # Generate .fai index for the reference if it doesn't exist
    if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai" ]; then
        $SAMTOOLS faidx "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"
    fi

    # Create genome.chrom.sizes from the .fai index
    cut -f1,2 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai" > "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes"
else
    echo "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes already exists in $OUTDIR."
fi

# High probability (> 0.75)
awk '$5 > 0.75 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph"
if [ ! -s "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph" ]; then
    echo "No high probability methylation calls found. Skipping bigWig generation for high probability calls."
else
    ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bw"
fi

# Low probability (< 0.25)
awk '$5 < 0.25 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph"
if [ ! -s "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph" ]; then
    echo "No low probability methylation calls found. Skipping bigWig generation for low probability calls."
else
    ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bw"
fi

# Undefined (0.25 <= x <= 0.75)
awk '$5 >= 0.25 && $5 <= 0.75 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph"
if [ ! -s "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph" ]; then
    echo "No undefined probability methylation calls found. Skipping bigWig generation for undefined probability calls."
else
    ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.undefined.bw"
fi

# Remove intermediate files
rm "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph" "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph" "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph"

# Create plotly reports
${FSHD_PYTHON} ${REALPATH}/helper/report_plotly.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"

# Remove intermediate files
if [[ "$REMOVE_INTERMEDIATE_FILES" == true ]]; then
    rm ${OUTDIR}/${PREFIX}_*blast*
    rm ${OUTDIR}/${PREFIX}_*top*
    rm ${OUTDIR}/${PREFIX}_*psl
    rm ${OUTDIR}/${PREFIX}_*d4z4*
fi

# Remove temp created fastq/ubam
if [[ "$ubam_exist" == true ]]; then
    rm "${INPUT_FASTQ}"
elif [[ "$fastq_exist" == true ]]; then
    rm "${INPUT_UBAM}"
fi

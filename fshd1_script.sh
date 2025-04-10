#!/bin/bash
#PBS -P kr68
#PBS -q normalbw
#PBS -l walltime=01:15:00
#PBS -l ncpus=3
#PBS -l mem=40GB
#PBS -l jobfs=10GB
#PBS -l wd

# Set LMDB memory size for BLAST
export BLASTDB_LMDB_MAP_SIZE=200000000

# Default inputs
INPUT_UBAM=""                            # Unaligned BAM
INPUT_FASTQ=""
REF=/g/data/kr68/genome/hs1.fa           # Reference genome
REF_DIR=$(dirname "$REF")
REGION_BED=inputs/d4z4_region.chm13.bed
FEATURES_FASTA=inputs/features.v2.fasta
PROBES=inputs/probes.fasta
PREFIX="SAMPLE"
OUTDIR="$PREFIX"
HAPLOTYPE_REFS=inputs/d4z4_repeats.fasta  # Fasta file containing haplotype-specific references
REPEATS_FASTA=inputs/dux4.gene_complete_genbank_20241127.reformatted.fasta

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
        --probes) PROBES="$2"; shift 2;;
        --haplotype-refs) HAPLOTYPE_REFS="$2"; shift 2;;
        --repeats-fasta) REPEATS_FASTA="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Check if required arguments are set
if [[ (-z "$INPUT_UBAM" || ! -f "$INPUT_UBAM") && (-z "$INPUT_FASTQ" || ! -f "$INPUT_FASTQ") ]]; then
    echo "Error: Either --input-ubam or --input-fastq must be provided and must be a valid file."
    exit 1
fi

# Keep the original file
# if have ubam, convert to fastq
# if have fastq, convert to ubam
echo "Convert ubam/fastq to fastq/ubam"
ubam_exist=false
fastq_exist=false
if [[  ( ! -z "$INPUT_UBAM" ) && ( -f "$INPUT_UBAM") ]]; then
    ubam_exist=true
    samtools fastq -TMM,ML "$INPUT_UBAM" > "${INPUT_UBAM%.bam}.fastq"
    INPUT_FASTQ="${INPUT_UBAM%.bam}.fastq"
elif [[  ( ! -z "$INPUT_FASTQ" ) && ( -f "$INPUT_FASTQ") ]]; then
    fastq_exist=true
    samtools import -u -s "$INPUT_FASTQ" -o "${INPUT_FASTQ%.fastq}.bam"
    INPUT_UBAM="${INPUT_FASTQ%.fastq}.bam"
fi

# Paths to conversion tools
PSLTOBED=/g/data/if89/apps/kentutils/0.0/bin/pslToBed
BG2BW=/g/data/if89/apps/kentutils/0.0/bin/bedGraphToBigWig

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Identify 4qA & 4qB probes in raw data
echo "IDENTIFY PROBES"
samtools fasta ${INPUT_UBAM} > ${OUTDIR}/${PREFIX}_all_reads_uBAM.fasta
makeblastdb -in ${OUTDIR}/${PREFIX}_all_reads_uBAM.fasta -dbtype nucl -out ${OUTDIR}/${PREFIX}_db
blastn -query ${PROBES} -db ${OUTDIR}/${PREFIX}_db -out ${OUTDIR}/${PREFIX}_probes.blast.txt -outfmt 6

echo "MERYL"
# Define meryl outputs relative to the reference genome directory
REF_NAME=$(basename $REF .fa)
MERYL_DB="inputs/merylDB_${REF_NAME}"
REPETITIVE_REGIONS="inputs/${REF_NAME}_repetitive_k15.txt"

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
echo "WINNOWMAP" 
winnowmap -W ${REPETITIVE_REGIONS} -Y -y -ax $MODE "$REF" "$INPUT_FASTQ" | \
        samtools view -L "$REGION_BED" -q ${MAPQ} -Sb | \
        samtools sort -o "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

samtools index "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

# Generate FASTA from uBAM
echo "Generating FASTA from uBAM"
samtools view -N <(samtools view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | cut -f1 | sort | uniq) -b "$INPUT_UBAM" | \
	samtools fasta > "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

# Convert BAM to BED
echo "Converting BAM to BED"
bedtools bamtobed -i "${OUTDIR}/${PREFIX}_reads_of_interest.bam" > "${OUTDIR}/${PREFIX}_reads_of_interest.bed"

# Map features to reads using BLAT
echo "Mapping reads to features using BLAT"
blat -t=dna -q=dna -maxIntron=500 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "$FEATURES_FASTA" "${OUTDIR}/${PREFIX}_mapped_features.psl"

# Convert PSL to BED
echo "Converting PSL to BED"
 ${PSLTOBED} "${OUTDIR}/${PREFIX}_mapped_features.psl" "${OUTDIR}/${PREFIX}_mapped_features.bed"

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
python3 helper/read_classification.py \
    --psl "${OUTDIR}/${PREFIX}_mapped_features.psl" \
    --bed "${OUTDIR}/${PREFIX}_reads_of_interest.bed" \
    --output "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" \
    --sslp "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --aligned-bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"

## Step 11: get mapped estimated copies
# python3 helper/get_cne.py --outdir "${OUTDIR}" --prefix "${PREFIX}"

## Step 12: Generate ordered alignment sequences
echo "Generating ordered alignment sequences"
alignment_script="helper/alignment_order.clipping.py"

# Run the alignment order script
python3 "$alignment_script" \
    --bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" \
    --output_fasta "${OUTDIR}/${PREFIX}_d4z4_units.fasta" \
    --output_table "${OUTDIR}/${PREFIX}_d4z4_units.tsv" \
    --xapi_bed "${OUTDIR}/${PREFIX}_xapi_sites.bed" \
    --blni_bed "${OUTDIR}/${PREFIX}_blni_sites.bed" \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" 

# Map d4z4 repeats back to the reads of interest
# minimap2 --secondary=no --MD -ax asm5 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "${OUTDIR}/${PREFIX}_d4z4_units.fasta" | samtools sort -o "${OUTDIR}/${PREFIX}_d4z4_repeats.bam"
# samtools index "${OUTDIR}/${PREFIX}_d4z4_repeats.bam"

# bedtools bamtobed -i "${OUTDIR}/${PREFIX}_d4z4_repeats.bam" > "${OUTDIR}/${PREFIX}_d4z4_repeats.bed"

########## check if 4a haplotype is long ##############
./identify_long_repeats.sh "${OUTDIR}" "${PREFIX}"

# Unify features and repeats into a single bed and add colour
python3 helper/add_colors_to_bed.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --features_bed "${OUTDIR}/${PREFIX}_mapped_features.bed" \
    --repeats_bed "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" \
    --sslp_bed "${OUTDIR}/${PREFIX}_SSLP.bed" \
    --output_bed "${OUTDIR}/${PREFIX}_all_features.bed"

# Extract the read from the uBAM and convert it to FASTQ and FASTA formats abd align the extracted FASTQ back to the FASTA reference (self-mapping)
samtools view -N <(samtools view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | cut -f1 | sort | uniq) -b ${INPUT_UBAM} | \
    samtools fastq -TMM,ML - | minimap2 -t ${THREADS} -Y -y -x asm5 -a --secondary=no "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" - | \
    samtools sort -@ ${THREADS} - > "${OUTDIR}/${PREFIX}_meth_reads.bam"

# Index the resulting BAM file
samtools index "${OUTDIR}/${PREFIX}_meth_reads.bam"

# use minimod to get meth probabilities for each CpG within a given read
minimod view -c m "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "${OUTDIR}/${PREFIX}_meth_reads.bam" | awk 'NR > 1 {print $4"\t"$5"\t"$5+1"\t.\t"$7}' > "${OUTDIR}/${PREFIX}_meth_reads.bed"

python3 helper/methylation_summary.py \
    --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --features_bed "${OUTDIR}/${PREFIX}_all_features.bed" \
    --meth "${OUTDIR}/${PREFIX}_meth_reads.bed" \
    --output "${OUTDIR}/${PREFIX}_methylation_summary.tsv" \
    --updated_bed "${OUTDIR}/${PREFIX}_updated_features.bed"

mv "${OUTDIR}/${PREFIX}_methylation_summary.tsv" "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"
mv "${OUTDIR}/${PREFIX}_updated_features.bed" "${OUTDIR}/${PREFIX}_all_features.bed"

# Check if genome.chrom.sizes exists
if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" ]; then
    echo "Creating ${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes in $OUTDIR..."
    
    # Generate .fai index for the reference if it doesn't exist
    if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai" ]; then
        samtools faidx "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"
    fi

    # Create genome.chrom.sizes from the .fai index
    cut -f1,2 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai" > "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes"
else
    echo "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes already exists in $OUTDIR."
fi

# High probability (> 0.75)
awk '$5 > 0.75 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph"
 ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bw"

# Low probability (< 0.25)
awk '$5 < 0.25 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph"
 ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bw"

# Undefined (0.25 <= x <= 0.75)
awk '$5 >= 0.25 && $5 <= 0.75 {print $1, $2, $3, $5}' "${OUTDIR}/${PREFIX}_meth_reads.bed" | sort -k1,1 -k2,2n > "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph"
 ${BG2BW} "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph" "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" "${OUTDIR}/${PREFIX}_meth_reads.undefined.bw"

# Remove intermediate files
rm "${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph" "${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph" "${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph"

# Remove temp created fastq/ubam
if [[ $ubam_exist ]]; then
    rm "${INPUT_UBAM%.bam}.fastq"
elif [[ $fastq_exist ]]; then
    rm "${INPUT_FASTQ%.fastq}.bam"
fi

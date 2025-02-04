#!/bin/bash
#$ -cwd
#$ -N fshd_test
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=6G
#$ -l tmp_requested=10G
#$ -pe smp 8
#$ -l h='!epsilon*'
#$ -V

# Exit on errors
set -e

# Record start time
START_TIME=$(date +%s)

# Add relevant modules
export PATH=/directflow/KCCGGenometechTemp/projects/andmar/software/Winnowmap/bin:/home/iradev/scripts:$PATH
module load centos7.8/qiadu/minimap2/2.22
module load centos7.8/jamfer/samtools/gcc-11.1.0/samtools-1.18
module load centos6.10/joaach/bedtools/2.25.0
module load centos6.10/gi/blat/35
module load centos6.10/kseskv/seqtk/1.3
module load centos7.8/aletan/gcc/11.2.0

MINIMOD=/directflow/KCCGGenometechTemp/projects/iradev/clinical/fshd_test/november_24/JURA89_RF.fshd_out/minimod/minimod
BG2BW=/directflow/KCCGGenometechTemp/projects/iradev/tandem_repeats/scripts/kentUtils/bin/bedGraphToBigWig

# Set LMDB memory size for BLAST
export BLASTDB_LMDB_MAP_SIZE=200000000

# Inputs
INPUT_UBAM=${2}      # Unaligned BAM
REF=${3}             # Reference genome
REF_DIR=$(dirname "$REF")
REGION_BED=d4z4_region.chm13.bed
FEATURES_FASTA=features.v2.fasta
PROBES=probes.fasta
PREFIX=${1}
OUTDIR=${1}
HAPLOTYPE_REFS=d4z4_repeats.fasta  # Fasta file containing haplotype-specific references
REPEATS_FASTA=dux4.gene_complete_genbank_20241127.reformatted.fasta

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Paths for PSL to BED conversion tool
PSLTOBED=/directflow/KCCGGenometechTemp/projects/iradev/tandem_repeats/scripts/kentUtils/bin/pslToBed

# Identify 4qA & 4qB probes in raw data
samtools fasta ${INPUT_UBAM} > ${OUTDIR}/${PREFIX}_all_reads_uBAM.fasta
makeblastdb -in ${OUTDIR}/${PREFIX}_all_reads_uBAM.fasta -dbtype nucl -out ${OUTDIR}/${PREFIX}_db
blastn -query ${PROBES} -db ${OUTDIR}/${PREFIX}_db -out ${OUTDIR}/${PREFIX}_probes.blast.txt -outfmt 6

# Define meryl outputs relative to the reference genome directory
MERYL_DB="${REF_DIR}/merylDB"
REPETITIVE_REGIONS="${REF_DIR}/repetitive_k15.txt"

# Step 1: Pre-compute k-mer frequency (only if meryl output does not exist)
if [ ! -d "$MERYL_DB" ]; then
    echo "Creating meryl database..."
    meryl count k=15 output "$MERYL_DB" "$REF"
else
    echo "Meryl database already exists at $MERYL_DB. Skipping this step."
fi

# Step 2: Extract repetitive regions (only if repetitive_k15.txt does not exist)
if [ ! -f "$REPETITIVE_REGIONS" ]; then
    echo "Extracting repetitive regions..."
    meryl print greater-than distinct=0.9998 "$MERYL_DB" > "$REPETITIVE_REGIONS"
else
    echo "Repetitive regions file already exists at $REPETITIVE_REGIONS. Skipping this step."
fi

# Run Winnowmap
winnowmap -W ${REPETITIVE_REGIONS} -Y -y -ax $MODE "$REF" <(samtools fastq -TMM,ML "$INPUT_UBAM") | \
	samtools view -L "$REGION_BED" -q ${MAPQ} -Sb | \
	samtools sort -o "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

samtools index "${OUTDIR}/${PREFIX}_reads_of_interest.bam"

# Generate FASTA from uBAM
echo "Generating FASTA from uBAM"
samtools view -N <(samtools view "${OUTDIR}/${PREFIX}_reads_of_interest.bam" | cut -f1 | sort | uniq) -b "$INPUT_UBAM" | \
	samtools fastq -@ ${THREADS} - | seqtk seq -A  > "${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

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

# Step 8: Parse PSL to generate a summary table
echo "Parsing PSL to generate summary table"
python3 parse_psl_output.v5.py \
    --psl "${OUTDIR}/${PREFIX}_mapped_features.psl" \
    --bed "${OUTDIR}/${PREFIX}_reads_of_interest.bed" \
    --output "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" \
    --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" \
    --sslp "${OUTDIR}/${PREFIX}_SSLP.bed"

# Step 9: Align reads to haplotype-specific references
#echo "Aligning reads classified to specific haplotypes"
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
    seqtk subseq <(samtools fastq "$INPUT_UBAM") <(grep ${chrom} "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" | cut -f1) > "$haplotype_reads"
 
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

# Step 11: Estimate the number of copies
mapped_copies_file="${OUTDIR}/${PREFIX}_mapped_copies.tsv"
echo -e "ReadID\tMappedEstimatedCopies\tMappedEstimatedCopiesPlus\tMappedEstimatedCopiesMinus" > "$mapped_copies_file"

while read -r line; do
    readID=$(echo "$line" | awk '{print $1}')
    chrom=$(echo "$line" | awk '{split($2, a, ":"); print a[1]}')  # Extract chromosome from column 2

    # Determine haplotype based on chromosome
    if [[ "$chrom" == "chr4" ]]; then
        haplotype="4A"
    elif [[ "$chrom" == "chr10" ]]; then
        haplotype="10A"
    else
        haplotype="NA"
    fi

    # Set the repeat_region based on haplotype
    if [[ "$haplotype" == "4A" ]]; then
        repeat_region="HM190196.1_dux4_4A166:1-3199"
    elif [[ "$haplotype" == "4B" ]]; then
        repeat_region="HM190164.1_dux4_4B168:1-3121"
    elif [[ "$haplotype" == "10A" ]]; then
        repeat_region="HM190192.1_dux4_10A180T:1-3204"
    elif [[ "$haplotype" == "10B" ]]; then
        repeat_region="HM190165.1_dux4_10B161T:1-3125"
    else
        repeat_region="NA"
    fi

    # Execute the main commands only if haplotype is not NA
    if [[ "$haplotype" != "NA" ]]; then
        # For the plus strand
        samtools view -N <(echo ${readID}) -h -b -F 16 "$merged_bam" > "${OUTDIR}/temp_read_plus.bam"
        samtools index "${OUTDIR}/temp_read_plus.bam"
        number_copies_plus=$(samtools depth -r "${repeat_region}" "${OUTDIR}/temp_read_plus.bam" | awk '{sum += $3; count++} END {if (count > 0) printf "%.2f", sum / count; else print 0}')

        # For the minus strand
        samtools view -N <(echo ${readID}) -h -b -f 16 "$merged_bam" > "${OUTDIR}/temp_read_minus.bam"
        samtools index "${OUTDIR}/temp_read_minus.bam"
        number_copies_minus=$(samtools depth -r "${repeat_region}" "${OUTDIR}/temp_read_minus.bam" | awk '{sum += $3; count++} END {if (count > 0) printf "%.2f", sum / count; else print 0}')

        rm "${OUTDIR}/temp_read_plus.bam"*
        rm "${OUTDIR}/temp_read_minus.bam"*

        # Compare the two numbers and set number_copies to the greater value
        if (( $(echo "$number_copies_plus > $number_copies_minus" | bc -l) )); then
            number_copies=$number_copies_plus
        else
            number_copies=$number_copies_minus
        fi

    else
        number_copies="NA"
    fi  
    echo -e "${readID}\t${number_copies}\t${number_copies_plus}\t${number_copies_minus}" >> "$mapped_copies_file"

done < <(grep -v "ReadID" "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv")

# Add the MappedEstimatedCopies column using paste
paste "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv" <(cut -f2-4 "$mapped_copies_file") > "${OUTDIR}/${PREFIX}_updated_features_summary.tsv"

# Replace the original summary file with the updated file
mv "${OUTDIR}/${PREFIX}_updated_features_summary.tsv" "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"

# Cleanup
rm "$mapped_copies_file"

# Step 12: Generate ordered alignment sequences
echo "Generating ordered alignment sequences"
alignment_script="alignment_order.clipping.v3.py"

# Run the alignment order script
 python3 "$alignment_script" --bam "${OUTDIR}/${PREFIX}_aligned_haplotypes.bam" --fasta "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" --output_fasta "${OUTDIR}/${PREFIX}_d4z4_units.fasta" --output_table "${OUTDIR}/${PREFIX}_d4z4_units.tsv" --xapi_bed "${OUTDIR}/${PREFIX}_xapi_sites.bed" --blni_bed "${OUTDIR}/${PREFIX}_blni_sites.bed" --main_tsv "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"

# Map d4z4 repeats back to the reads of interest
minimap2 --secondary=no --MD -ax asm5 "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "${OUTDIR}/${PREFIX}_d4z4_units.fasta" | samtools sort -o "${OUTDIR}/${PREFIX}_d4z4_repeats.bam"
samtools index "${OUTDIR}/${PREFIX}_d4z4_repeats.bam"

bedtools bamtobed -i "${OUTDIR}/${PREFIX}_d4z4_repeats.bam" > "${OUTDIR}/${PREFIX}_d4z4_repeats.bed"

# Unify features and repeats into a single bed and add colour
 python3 add_colors_to_bed.py \
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
 ${MINIMOD} view -c m "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" "${OUTDIR}/${PREFIX}_meth_reads.bam" | awk 'NR > 1 {print $4"\t"$5"\t"$5+1"\t.\t"$7}' > "${OUTDIR}/${PREFIX}_meth_reads.bed"

 python3 methylation_summary.v5.py \
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

# Record end time
END_TIME=$(date +%s)

# Calculate elapsed time
ELAPSED_TIME=$((END_TIME - START_TIME))

# Format elapsed time as hours:minutes:seconds
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

# Print the formatted time
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"

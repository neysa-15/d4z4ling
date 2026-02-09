#!/bin/bash -l
#PBS -P kr68
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68
#PBS -l jobfs=10GB
#PBS -l wd

# set -e

module unload python3
module unload blast
module unload seqkit
module unload bedtools

module load python3/3.12.1
module load blast/2.14.1
module load seqkit/2.9.0
module load bedtools/2.31.0

OUTDIR=$1
PREFIX=$2

repeat_threshold=2640 # 80% of 3300 bp

# Take reads of interest that is labeled "Complete" or "Partial distal"
echo "Take read id of reads with complete and partial distal"
TSV="${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"
FEATURES_BED="${OUTDIR}/${PREFIX}_all_features.bed"
ORIGINAL_SEQ="${OUTDIR}/${PREFIX}_reads_of_interest.fasta"

read_id_list=$(egrep -i 'Complete 4qA|Partial distal 4qA' $TSV | cut -f1 | sort | uniq)

# Extract sequence of distal copy until pLAM to a fasta file
EXTRACTED_SEQUENCE="${OUTDIR}/${PREFIX}_read_distal_seq.fasta"
DISTAL_COORDS_BED="${OUTDIR}/${PREFIX}_distal_coords.bed"

# Empty output files first
> "$EXTRACTED_SEQUENCE"
> "$DISTAL_COORDS_BED"

echo "get coordinates and sequence of distal copy and the pLAM"
# For each read ID
for read_id in $read_id_list; do
    echo $read_id

    # Extract pLAM line for this read
    tmp_read_plam="${OUTDIR}/${PREFIX}_tmp_${read_id}_pLAM.bed"
    grep $read_id $FEATURES_BED | grep "pLAM" > ${tmp_read_plam}

    [[ ! -s $tmp_read_plam ]] && continue  # skip if none

    # Extract repeats for this read and take the one with size more than 80% of 3300 bp
    tmp_read_repeats="${OUTDIR}/${PREFIX}_tmp_${read_id}_repeats.bed"
    grep "$read_id" "$FEATURES_BED" | grep "d4z4_" | grep -v "proximal" | awk -v out="$tmp_read_repeats" '
    {
        size = $3 - $2
        if (size > 2640) {
            print $0 >> out
            has_large = 1
        }
        if (size > max_size) {
            max_size = size
            max_line = $0
        }
    }
    END {
        if (!has_large) {
            print max_line >> out
        }
    }'

    [[ ! -s $tmp_read_repeats ]] && continue

    sort -k1,1 -k2,2n $tmp_read_repeats > "${OUTDIR}/${PREFIX}_tmp_${read_id}_repeats.sorted.bed"
    mv "${OUTDIR}/${PREFIX}_tmp_${read_id}_repeats.sorted.bed" $tmp_read_repeats

    # Find the closest repeats to pLAM
    tmp_distal_bed="${OUTDIR}/${PREFIX}_tmp_${read_id}_distal.bed"
    # bedtools closest -s -a $tmp_read_plam -b $tmp_read_repeats > $tmp_distal_bed
    bedtools closest -s -io -a "$tmp_read_plam" -b "$tmp_read_repeats" | uniq > $tmp_distal_bed

    # Take the lowest and highest coordinate, strand and repeat number of distal copy
    while IFS=$'\t' read -r chrA startA endA nameA scoreA strandA thickStartA thickEndA colorA blockCountA blockSizesA blockStartsA \
                      chrB startB endB nameB scoreB strandB thickStartB thickEndB colorB blockCountB blockSizesB blockStartsB; do

        # Compute max coordinates
        start_coord=$(( startA < startB ? startA : startB ))
        end_coord=$(( endA > endB ? endA : endB ))
        strand=$strandA
        repeat_number="${nameB}"  # Extract number from d4z4_* label

        [[ -z "$start_coord" || -z "$end_coord" ]] && continue

        if [[ "$start_coord" -le 0 ]]; then
            start_coord=1
        fi

        size=$(( end_coord - start_coord ))

        echo -e "${read_id}\t${start_coord}\t${end_coord}\t${repeat_number}\t.\t${strand}\t${size}," >> "$DISTAL_COORDS_BED"
        echo "Extracting $read_id : ${start_coord}-${end_coord} on strand $strand" >&2
        seqkit subseq --chr "$read_id" -r "${start_coord}:${end_coord}" "$ORIGINAL_SEQ" >> "$EXTRACTED_SEQUENCE"

    done < "$tmp_distal_bed"

    rm $tmp_read_plam
    rm $tmp_read_repeats
    rm $tmp_distal_bed
done

echo "✅ Extracted sequences written to $EXTRACTED_SEQUENCE"

# Clean up
# rm ${OUTDIR}/${PREFIX}_tmp_*

# Set LMDB memory size for BLAST
export BLASTDB_LMDB_MAP_SIZE=200000000

# BLAST distal copy to DUX4 sequence
DUX4_SEQUENCES=inputs/dux4_genbank.fasta
DUX4_DB=${OUTDIR}/${PREFIX}_db
BLAST_RESULTS=${OUTDIR}/${PREFIX}_dux4.blast.tsv
makeblastdb -in "${DUX4_SEQUENCES}" -dbtype nucl -out ${DUX4_DB}
blastn -query ${EXTRACTED_SEQUENCE} -db ${DUX4_DB} -out ${BLAST_RESULTS} -outfmt "6 qseqid sseqid pident length slen evalue bitscore qcovs qcovhsp"

# Set output file names
blast_with_score="${OUTDIR}/${PREFIX}_blast_with_score.tsv"
top1_read="${OUTDIR}/${PREFIX}_top1_blast.tsv"
# top3_read="${OUTDIR}/${PREFIX}_top3_blast.tsv"

# Add header to full results file
header="ReadID\tDUX4_genbank\tpercent_identity\talignment_length\tsubject_length\tevalue\tbitscore\tcoverage\tcoverage_hsp\treciprocal_coverage\tfinal_score"
echo -e "$header" > "$blast_with_score"

# Compute reciprocal coverage + final weighted score
awk -F'\t' -v OFS='\t' '{
    recip_cov = ($4/$5)*100
    evalue_score = ($6 > 0) ? -log($6) : 100  # transform evalue
    final_score = (recip_cov * 40) + ($8 * 40) + ($3 * 15) + (evalue_score * 5)
    print $0, recip_cov, final_score
}' "$BLAST_RESULTS" >> "$blast_with_score"

# Write header to top1 and top3 files
head -n1 "$blast_with_score" > "$top1_read"

# Sort by: final_score desc, then coverage desc, then evalue asc, then percent_identity desc
awk 'NR>1' "$blast_with_score" | \
sort -t$'\t' -k11,11gr -k8,8gr -k6,6g -k3,3gr | \
awk -F'\t' -v top1="$top1_read" '{
    count[$1]++
    if (count[$1]==1) print >> top1
}'

# Separate coordinates from read id
awk -F'\t' -v OFS='\t' '
NR==1 {
    print $0, "start_coords", "end_coords"
    next
}
{
    n = split($1, arr, "_")
    coords = arr[n]
    read_id = $1
    sub("_" coords "$", "", read_id)
    split(coords, pos, "-")
    $1 = read_id
    print $0, pos[1], pos[2]
}' "$top1_read" > "${OUTDIR}/${PREFIX}_top1_read.edited.tsv"

mv "${OUTDIR}/${PREFIX}_top1_read.edited.tsv" "$top1_read"

cut -f2 "$top1_read" | awk '{
    if ($0 ~ /haplotype_[^_]+/) {
        match($0, /haplotype_[^_]+/, arr)
        hap = arr[0]
        count[hap]++
    } else {
        count[$0]++
    }
}
END {
    for (h in count) print count[h], h
}' | sort -nr | grep -v "DUX4_genbank" > "${OUTDIR}/${PREFIX}_top1_match_counts.txt"

# Clean up
rm ${OUTDIR}/${PREFIX}_db*

# Count percentage of DUX4L in 4qA reads
percentage_longvar=$(awk '
{
    count[$2] = $1
    total += $1
}
END {
    dux4_count = count["MF693913.1_Homo_sapiens_DUX4_(DUX4L)_gene"]
    percentage = (dux4_count / total) * 100
    printf "%.2f", percentage
}' "${OUTDIR}/${PREFIX}_top1_match_counts.txt")

echo "$percentage_longvar"

# Reannotate reads with long var
output_bed="${OUTDIR}/${PREFIX}_updated_long_var.bed"

if (( $(echo "$percentage_longvar > 70" | bc -l) )); then
    echo "$PREFIX have long var"
    match_read_tmp=${OUTDIR}/${PREFIX}_tmp_matching_reads.txt
    tmp_update_bed=${OUTDIR}/${PREFIX}_tmp_updates.bed

    # Step 1: Get ReadIDs with DUX4L match
    awk -F'\t' '$2=="MF693913.1_Homo_sapiens_DUX4_(DUX4L)_gene" {print $1}' "$top1_read" > "$match_read_tmp"

    # Step 2: grep updated coords for matching reads
    fgrep -f  "$match_read_tmp" "$DISTAL_COORDS_BED" > "$tmp_update_bed"

    # Step 3: Apply updates to all_features.bed
    awk 'BEGIN{FS=OFS="\t"} 
    NR==FNR {
        key=$1"\t"$4
        updates[key]=$2 OFS $3 OFS $7
        next
    } 
    {
        key=$1"\t"$4
        if (key in updates) {
            split(updates[key], coords, "\t")
            $2=coords[1]
            $3=coords[2]
            $7=coords[1]
            $8=coords[2]
            $11=coords[3]
        }
        print
    }' "$tmp_update_bed" "$FEATURES_BED" > "$output_bed"

    mv "$output_bed" "$FEATURES_BED"

    # Clean up
    rm $match_read_tmp $tmp_update_bed
fi
 
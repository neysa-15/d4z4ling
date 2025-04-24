#!/bin/bash
#PBS -P kr68
#PBS -q normalbw
#PBS -l walltime=00:15:00
#PBS -l ncpus=5
#PBS -l mem=40GB
#PBS -l jobfs=10GB
#PBS -l wd
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/

module unload bedtools
module unload kentutils
module unload blat
module unload seqkit

module load bedtools/2.28.0
module load seqkit/2.9.0
module load blat/37
module load kentutils/0.0
PSLTOBED=/g/data/if89/apps/kentutils/0.0/bin/pslToBed
BLAT=/g/data/kr68/andre/software/blat

# Long D4z4 repeats fasta
LONG_D4Z4="/g/data/kr68/neysa/fshd_pipeline/inputs/4A161L.fa"

OUTDIR="AS2603_pos" # variable
PREFIX="AS2603" # variable

# Extract just pLAMs
grep -E 'pLAM' "${OUTDIR}/${PREFIX}_mapped_features.bed" > "${OUTDIR}/${PREFIX}_pLAM.bed"

# Take reads that have NON-intersecting d4z4 and pLAM as potential reads with 4A161l
bedtools intersect -v -a "${OUTDIR}/${PREFIX}_pLAM.bed" -b "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" > "${OUTDIR}/${PREFIX}_potential_long.bed"

# grep potential reads from the reads_of_interest fasta
cut -f1 "${OUTDIR}/${PREFIX}_potential_long.bed" | seqkit grep -f - "${OUTDIR}/${PREFIX}_reads_of_interest.fasta" > "${OUTDIR}/${PREFIX}_long.fasta"

# Align potential long d4z4 repeats with long d4z4 sequence
blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 "${OUTDIR}/${PREFIX}_long.fasta" "${LONG_D4Z4}" "${OUTDIR}/${PREFIX}_4a_long_blat_output.psl"

# Convert psl to bed
${PSLTOBED} "${OUTDIR}/${PREFIX}_4a_long_blat_output.psl" "${OUTDIR}/${PREFIX}_4a_long_blat_output.bed"

# Take reads if the long d4z4 intersect with both the distal repeat and the pLAM
bedtools intersect -wa -wb -a "${OUTDIR}/${PREFIX}_4a_long_blat_output.bed" -b "${OUTDIR}/${PREFIX}_pLAM.bed" | cut -f1 > "${OUTDIR}/long_intersect_pLAM.txt"
bedtools intersect -wa -wb -a "${OUTDIR}/${PREFIX}_4a_long_blat_output.bed" -b "${OUTDIR}/${PREFIX}_d4z4_repeats.bed" | cut -f1 > "${OUTDIR}/long_intersect_d4z4.txt"
comm -12 <(sort "${OUTDIR}/long_intersect_pLAM.txt" | uniq) <(sort "${OUTDIR}/long_intersect_d4z4.txt" | uniq) > "${OUTDIR}/${PREFIX}_long_d4z4_readIDs.txt"

grep -F -f "${OUTDIR}/${PREFIX}_long_d4z4_readIDs.txt" "${OUTDIR}/${PREFIX}_4a_long_blat_output.bed" > "${OUTDIR}/${PREFIX}_4a_long_repeats_final_blat37_webparameter.bed"

rm "${OUTDIR}/long_intersect_pLAM.txt" "${OUTDIR}/long_intersect_d4z4.txt"
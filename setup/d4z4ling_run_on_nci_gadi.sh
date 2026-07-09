#!/bin/bash -l
#PBS -P kr68
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68+gdata/ox63+scratch/ox63
#PBS -l jobfs=10GB
#PBS -l wd

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/

module unload Miniconda3
module load Miniconda3/4.12.0

conda activate /path/to/d4z4ling/conda_env_d4z4ling

if [[ $UBAM == *.fasta ]] || [[ $UBAM == *.fa ]]; then
    echo "Processing fasta file: $SAMPLE $UBAM"
    ./d4z4ling_no_methylation.sh --prefix ${SAMPLE} --outdir ${OUTDIR} --input-fasta ${UBAM}
else
    echo "Processing uBAM file: $SAMPLE $UBAM"
    ./d4z4ling.sh --prefix ${SAMPLE} --outdir ${OUTDIR} --input-ubam ${UBAM} --remove-intermediate-files true
fi

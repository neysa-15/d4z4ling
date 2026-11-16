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

# unload modules before loading to ensure correct version
module unload python3
module unload minimap2
module unload samtools
module unload bedtools
module unload seqtk
module unload blast
module unload gcc
module unload seqkit
module unload minimod
module unload meryl
module unload Winnowmap
module unload kentutils
module unload blat

# load modules
module load python3/3.12.1
module load minimap2/2.28
module load samtools/1.19
module load bedtools/2.28.0
module load seqtk/1.3
module load blast/2.11.0
module load gcc/12.2.0
module load seqkit/2.9.0
module load minimod/0.4.0
module load meryl/1.4.1
module load Winnowmap/2.03
module load kentutils/0.0
module load blat/37

REALPATH=${REALPATH:-$PBS_O_WORKDIR}
echo "REALPATH: $REALPATH"

source ${REALPATH}/d4z4ling_venv/bin/activate

PSLTOBED=/g/data/if89/apps/kentutils/0.0/bin/pslToBed
BG2BW=/g/data/if89/apps/kentutils/0.0/bin/bedGraphToBigWig

if [[ $UBAM == *.fasta ]] || [[ $UBAM == *.fa ]]; then
    echo "Processing fasta file: $SAMPLE $UBAM"
    ${REALPATH}/d4z4ling_no_methylation.sh --prefix ${SAMPLE} --outdir ${OUTDIR} --input-fasta ${UBAM}
else
    echo "Processing uBAM file: $SAMPLE $UBAM"
    ${REALPATH}/d4z4ling.sh --prefix ${SAMPLE} --outdir ${OUTDIR} --input-ubam ${UBAM} \
		--remove-intermediate-files true \
		--minimap-ref-index ${REALPATH}/genome/hs1.mmi \
		--ref ${REALPATH}/genome/hs1.fa
fi

#!/bin/bash -l
#PBS -P kr68
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68
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
module unload blat

# load modules
module load python3/3.12.1
module load minimap2/2.22
module load samtools/1.19
module load bedtools/2.28.0
module load seqtk/1.3
module load blast/2.11.0
module load gcc/12.2.0
module load blat/37
module load seqkit/2.9.0
module load minimod/0.2.0
module load meryl/1.4.1
module load Winnowmap/2.03
module load kentutils/0.0

./fshd1_new.sh --prefix AS2603 --outdir AS2603_6-3 --input-ubam /g/data/kr68/neysa/uBAM/QGXHXX240243_supmeth_unaligned.bam

#!/bin/bash -l
#PBS -P {YOUR_PROJECT}
#PBS -q copyq
#PBS -l walltime=01:30:00
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l storage=gdata/if89+gdata/{YOUR_PROJECT}+scratch/{YOUR_PROJECT}
#PBS -l jobfs=10GB
#PBS -l wd

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/

module load Miniconda3/4.12.0

conda config --set channel_priority flexible
conda config --add pkgs_dirs /scratch/{YOUR_PROJECT}/{USER}/conda_pkgs

conda env create -f /path/to/d4z4ling/environment.yml --prefix /path/to/d4z4ling/conda_env_d4z4ling

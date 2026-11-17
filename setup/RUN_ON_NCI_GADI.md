# d4z4ling set up for NCI user
- [d4z4ling set up for NCI user](#d4z4ling-set-up-for-nci-user)
  - [What is NCI?](#what-is-nci)
  - [NCI Pre-requisite](#nci-pre-requisite)
  - [Reference Genome](#reference-genome)
    - [Generate meryl database for repetitive region mapping](#generate-meryl-database-for-repetitive-region-mapping)
  - [OPTION1: WITHOUT CONDA, with python venv](#option1-without-conda-with-python-venv)
    - [How to run](#how-to-run)
  - [OPTION2: WITH CONDA](#option2-with-conda)
    - [Setup conda environment](#setup-conda-environment)
    - [Run](#run)
  - [Available default input file](#available-default-input-file)

## What is NCI?
NCI is Australia's National Computational Infrastructure, a supercomputing facility used by researchers across Australia. For more information: [nci.org.au](https://nci.org.au/about-us/who-we-are).

## NCI Pre-requisite
Request to be in if89 project on NCI website 

Add this to your `~/.bash_profile` once approved to the project
```
module use -a /g/data/if89/apps/modulefiles
```

## Reference Genome
Load module
```
module load samtools
module load minimap2
```

Get a copy of the hs1 (T2T-CHM13) reference genome and create indexes.
```
bash ./genome/download_ref_chm13.sh
```

### Generate meryl database for repetitive region mapping
```
bash ./genome/meryl.sh
```

## OPTION1: WITHOUT CONDA, with python venv
Load python3:
```
module unload python3
module load python3/3.12.1
```

Open your python3 virtual environment and install requirements

```
python3 -m venv d4z4ling_venv
source d4z4ling_venv/bin/activate
pip install --upgrade pip
pip install -r setup/requirements.txt
```

### How to run

Change at least these two config in `nci_run.sh`:
```
#PBS -P {YOUR PROJECT}
#PBS -l storage=gdata/{YOUR PROJECT}+gdata/if89+scratch/{YOUR PROJECT}
```

Last line of `nci_run.sh` should be custom edited depending on the samples  
Required field:

```
./d4z4ling.sh --prefix {Prefix of output file} --outdir {Output directory name} --input-ubam {Path to uBAM input}
```

Optional field would be adding these flags (which have default input if not added):

```
--input-fastq {Path to fastq input} (can replace --input-bam, either one have to be provided)
--ref {Path to reference genome fasta}
--region-bed {Path to d4z4 region bed file}
--features-fasta {Path to features fasta file}
--short-features-fasta {Path to features fasta file with sequence under 1000bp}
--plam {plam sequence fasta}
--haplotype-refs {Path to D4Z4 haplotype references}
--remove-intermediate-files {true/false}
--minimap-ref-index {Path to minimap index file of reference}
--threads {Number of threads}
--analysis-mode {true/false to generate bigwig and plotly visualisation files}
```

Job script:

```
qsub nci_run.sh
```

Job script with variable through job submission:
```
qsub nci_run.sh -v SAMPLE=${sample},OUTDIR=${outdir},UBAM=${ubam},REALPATH=/path/to/d4z4ling
```

## OPTION2: WITH CONDA
### Setup conda environment
Load conda:
```
module load Miniconda3/4.12.0
```

Change project to your project on gadi in `setup/nci_gadi_conda_create_venv.sh` and run the script to setup conda environment
```
qsub setup/nci_conda_create_venv.sh
```

### Run
Change these two config in `setup/d4z4ling_run_on_nci_gadi.sh`:
```
#PBS -P {YOUR PROJECT}
#PBS -l storage=gdata/{YOUR PROJECT}+gdata/if89+scratch/{YOUR PROJECT}
```

Last line of `setup/d4z4ling_run_on_nci_gadi.sh` should be custom edited depending on the samples  
Required field:

```
./d4z4ling.sh --prefix {Prefix of output file} --outdir {Output directory name} --input-ubam {Path to uBAM input}
```

Optional field would be adding these flags (which have default input if not added):

```
--input-fastq {Path to fastq input} (can replace --input-ubam, either one have to be provided)
--ref {Path to reference genome fasta}
--region-bed {Path to d4z4 region bed file}
--features-fasta {Path to features fasta file}
--short-features-fasta {Path to features fasta file with sequence under 1000bp}
--plam {plam sequence fasta}
--haplotype-refs {Path to D4Z4 haplotype references}
--remove-intermediate-files {true/false}
--minimap-ref-index {Path to minimap index file of reference}
--threads {Number of threads}
--analysis-mode {true/false to generate bigwig and plotly visualisation files}
```

Job script:

```
qsub setup/d4z4ling_run_on_nci_gadi.sh
```

Job script with variable through job submission:
```
qsub setup/d4z4ling_run_on_nci_gadi.sh -v SAMPLE=${sample},OUTDIR=${outdir},UBAM=${ubam},REALPATH=/path/to/d4z4ling
```

## Available default input file

1\. `d4z4_region.chm13.bed`  
2\. `d4z4_repeats.fasta`  
3\. `dux4.gene_complete_genbank_20241127.reformatted.fasta`  
4\. `features.fasta`  
5\. `short_features.fasta`  
6\. `pLAM.fasta`  
7\. `Fshd1_status_template.tsv`  


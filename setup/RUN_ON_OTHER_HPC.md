# d4z4ling set up other HPC user
- [d4z4ling set up other HPC user](#d4z4ling-set-up-other-hpc-user)
  - [Pre-requisite](#pre-requisite)
  - [Setup conda environment](#setup-conda-environment)
  - [Reference Genome](#reference-genome)
  - [Run](#run)
  - [Available default input file](#available-default-input-file)

## Pre-requisite
Have conda/miniconda installed. Install [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Setup conda environment
```
conda config --set channel_priority flexible
conda env create -f setup/environment.yml --prefix /path/to/conda_env_d4z4ling
```

Activate:
```
conda activate /path/to/conda_env_d4z4ling
```

## Reference Genome
Get a copy of the hs1 (T2T-CHM13) reference genome
```
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz inputs/.
gunzip inputs/hs1.fa.gz
samtools faidx inputs/hs1.fa
```

Run minimap2 reference indexing for script input (For efficiency purposes when running for multiple samples to the same reference)
```
minimap2 -d inputs/hs1.mmi {Path to chm13 fasta}
```

## Run

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
--PLAM {plam sequence fasta}
--haplotype-refs {Path to D4Z4 haplotype references}
--remove-intermediate-files {true/false}
--minimap-ref-index {Path to minimap index file of reference}
```

## Available default input file

1\. `d4z4_region.chm13.bed`  
2\. `d4z4_repeats.fasta`  
3\. `dux4.gene_complete_genbank_20241127.reformatted.fasta`  
4\. `features.fasta`  
5\. `short_features.fasta`  
6\. `pLAM.fasta`  
7\. `Fshd1_status_template.tsv`  


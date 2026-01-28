# d4z4ling set up for NCI user

## What is NCI?
NCI is Australia's National Computational Infrastructure. For more information: (https://nci.org.au/about-us/who-we-are).

## NCI Pre-requisite
Request to be in if89 project on NCI website 

Add this to your \~/.bash\_profile once approved to project
```
module use -a /g/data/if89/apps/modulefiles
```

## Installation
Load python3 (NCI users):
```
module unload python3
module load python3/3.12.1
```

Open your python3 virtual environment and install requirements

```
python3 -m venv fshd1_script
source fshd1_script/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## How to run

Last line of `nci_run.sh` should be custom edited depending on the samples  
Required field:

```
./fshd1_script.sh --prefix {Prefix of output file} --outdir {Output directory name} --input-ubam {Path to uBAM input}
```

Optional field would be adding these flags (which have default input if not added):

```
--input-fastq {Path to fastq input} (can replace --input-bam, either one have to be provided)
--ref {Path to reference genome fasta}
--region-bed {Path to d4z4 region bed file}
--features-fasta {Path to features fasta file}
--short-features-fasta {Path to features fasta file with sequence under 1000bp}
--PLAM {plam sequence fasta}
--haplotype-refs {Path to D4Z4 haplotype references}
--remove-intermediate-files {true/false}
--minimap-ref-index {Path to minimap index file of reference}
```

Job script:

```
qsub nci_run.sh
```

# fshd_pipeline

# Run script on NCI
Load python3:
```
module unload python3
module load python3/3.12.1
```

Open your python3 virtual environment
```
python3 -m venv fshd1_script
source fshd1_script/bin/activate
```

Install python requirements
```
pip install --upgrade pip
pip install -r requirements.txt
```

Request to be in IF89 project on NCI website

Add this to your ~/.bash_profile once approved to project
```
module use -a /g/data/if89/apps/modulefiles
```

Job script:
```
qsub nci_run.sh
```

Last line of nci_run.sh should be customed edited depending on the samples
Required field:
```
./fshd1_script.sh --prefix {PREFIX OF OUTPUT FILES} --outdir {OUTPUT DIRECTORY NAME} --input-ubam {PATH-TO-UBAM-INPUT}
```

Optional field would be adding these flags (which have default input if not added):
```
--input-fastq {PATH-TO-FASTQ-INPUT} (can replace --input-bam, either one have to be provided)
--ref {PATH-TO-REFERENCE}
--region-bed {PATH-TO-D4Z4-REGION-BED-FILE}
--features-fasta {PATH-TO-FEATURES-FASTA-FILE}
--probes {PATH-TO-PROBES-FASTA-FILE}
--haplotype-refs {PATH-TO-D4Z4-HAPLOTYPE-REFS}
--repeats-fasta {PATH-TO-REPEATS-FASTA-FILE}
```

# Available default input file
1. d4z4_region.chm13.bed
2. d4z4_repeats.fasta
3. dux4.gene_complete_genbank_20241127.reformatted.fasta
4. features.v2.fasta
5. probes.fasta

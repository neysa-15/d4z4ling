# d4z4ling demo with HG01811 sample from ONT1000G project

# Demo data source

HG01811 is a sample from ONT1000g project public data [https://s3.amazonaws.com/1000g-ont/index.html](https://s3.amazonaws.com/1000g-ont/index.html).  
For demonstration purposes, the data has been subsetted to only a few reads that maps to the D4Z4 region, with the majority having diagnostic value to FSHD diagnosis.

# Run setup

```shell
SAMPLE=HG01811
OUTDIR=demo/result_HG01811
UBAM=demo/HG01811.ONT.R10.LSK114.dorado090_sup_5mCG_5hmCG_v500.subsetted.unmapped.bam
REF=${your path to hs1.fasta file}
MMI=${path to hs1.mmi minimap index file}

./fshd1_script.sh --prefix ${SAMPLE} --outdir ${OUTDIR} --input-ubam ${UBAM} --remove-intermediate-files true --ref ${REF} --minimap-ref-index ${MMI}
```

# Expected output

Expected output will show on `demo/result_HG01811` directory, but the main output to look at for base check are:

1. Visualisation report in `demo/result_HG01811/HG01811_report.html`  
2. Read by read detail in `demo/result_HG01811/HG01811_mapped_features_summary.tsv`  
3. Summary counts of reads diagnostic value in `demo/result_HG01811/HG01811_fshd1_status_counts.tsv`

| FSHD1\_status | Haplotype | Read\_type | Copy\_number | Methylation | PolyA | Count |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| Diagnostic/Positive | 4qA | Complete | \<=10 | \<50% | Permissive | 0 |
| Diagnostic/Negative | 4qA | Complete/Partial distal | \>10 | \>=50% | Any | 2 |
| Diagnostic/Negative | 4qB | Complete/Partial distal | \>3 | Any | Non-permissive / absent | 11 |
| Diagnostic/Conflicting | 4qA | Complete | \<=10 | \>=50% | Non-permissive / absent | 0 |
| Support/Positive | 4qA | Partial distal | 3-Oct | \<50% | Permissive | 0 |
| Support/Negative | 4qA | Partial distal | 3-Oct | \>=50% | Any | 2 |
| Non-diagnostic | Any | Any | Any | Any | Any | 6 |

# Time

1. Cloning \+ installation time \< 5 minutes  
2. Run time of `demo/HG01811.ONT.R10.LSK114.dorado090\_sup\_5mCG\_5hmCG\_v500.subsetted.unmapped.bam` on Australia's National Computational Infrastructure:


```shell
Resource requested:
#PBS -l walltime=05:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB

Resource used:
   Service Units: 7.89
   CPU Time Used: 00:06:35        
   Memory Used  : 17.23GB         
   Walltime Used: 00:04:56  
```


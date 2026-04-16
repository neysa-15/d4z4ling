# d4z4ling run time benchmark on gadi

Unaligned BAM (uBAM) files produced by ONT readfish for our myopathy cohort from our paper were used for this benchmark. The tests were done on the Australia’s National Computer Infrastructure (NCI) Gadi supercomputer, with varying uBAM size.  
Although it’s not an exhaustive benchmark, the following conclusions can be drawn:

* Our largest uBAM file with 61GB takes 1 hour and 40.5 minutes to run, meaning that files up to this size are expected to finish in under 2 hours.  
* Processing time increases roughly proportionally with UBAM size and number of reads.

# NCI resource request

```
#PBS -l walltime=05:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
```

# Software versions

| Tool | Version |
| :---- | :---- |
| python3 | 3.12.1 |
| minimap2 | 2.28 |
| samtools | 1.19 |
| bedtools | 2.28.0 |
| seqtk | 1.3 |
| blast | 2.11.0 |
| gcc | 12.2.0 |
| seqkit | 2.9.0 |
| minimod | 0.2.0 |
| meryl | 1.4.1 |
| Winnowmap | 2.03 |
| kentutils | 0.0 |
| blat | 37 |

# Python library versions

| Library | Version |
| :---- | :---- |
| Bio | 1.8.1 |
| biopython | 1.81 |
| numpy | 2.3.5 |
| pandas | 2.3.3 |
| plotly | 5.18.0 |
| pysam | 0.22.1 |

# Pre-computation run time

| Indexing steps | Run time (hh:mm:ss) |
| :---- | :---- |
| meryl kmer frequencies and repetitive regions extraction | 00:02:50 |
| Minimap2 reference indexing to chm13 | 00:01:15 |

# Run time with different size of uBAM files

| uBAM size (GB) | Number of reads in uBAM | Minimap2 reference alignment of all reads (sec) | Winnowmap time (sec) | Run time (hh:mm:ss) | Reference alignment time (%) |
| :---- | :---- | :---- | :---- | :---- | :---- |
| 2.9 | 242109 | 150.879 | 140.928 | 00:10:29 | 46.39 |
| 3.8 | 327805 | 162.779 | 144.905 | 00:10:39 | 48.15 |
| 15 | 12444528 | 643.900 | 221.308 | 00:30:42 | 46.97 |
| 19 | 14928162 | 732.305 | 236.005 | 00:32:31 | 49.63 |
| 37 | 23786170 | 1008.987 | 349.437 | 00:49:22 | 45.86 |
| 42 | 18516022 | 1136.588 | 384.134 | 01:12:20 | 35.04 |
| 42 | 22202688 | 1157.658 | 385.416 | 01:13:08 | 35.17 |
| 48 | 25589601 | 2233.813 | 425.889 | 01:24:56 | 52.19 |
| 50 | 22266720 | 1363.484 | 425.927 | 01:12:28 | 41.15 |
| 54 | 27308552 | 2449.523 | 484.331 | 02:00:11 | 40.69 |
| 56 | 25793121 | 1565.560 | 476.853 | 01:35:42 | 35.57 |
| 61 | 38581010 | 2486.657 | 505.358 | 01:39:32 | 50.1 |

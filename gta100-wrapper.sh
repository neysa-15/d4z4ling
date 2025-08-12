#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")


export PATH=$PATH:/install/merqury/meryl-1.4.1/bin/ #only for creatiob

export PATH=$PATH:${REALPATH}/install/winnowmap-v2.03/bin
export PATH=$PATH:${REALPATH}/install/seqkit-v2.3.0/
export PATH=$PATH:${REALPATH}/install/seqtk-1.5
export PATH=$PATH:${REALPATH}/install/ncbi-blast-2.17.0+/bin/
export PATH=$PATH:${REALPATH}/install/minimod-v0.4.0/
export PATH=$PATH:${REALPATH}/install/kentutils/
#samtools
#bedtools
#slow5-dorado
#minimap2
# python

# check for all /data/hasindu paths and edit accordingly
#check for all the relative paths and edit accordnly

source ${REALPATH}/venv3/bin/activate
# pip3 install pandas Bio pysam plotly scipy

#/usr/bin/time -v /install/slow5-dorado-0.8.3/bin/slow5-dorado basecaller --models-directory /install/dorado-0.8.3/models /install/dorado-0.8.3/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 /data/hasindu/realfshd/test/data/PGXXXX230339_reads_20k.blow5  --modified-bases 5mCG_5hmCG -x cuda:all > reads.bam
# 0:44.51

#/usr/bin/time -v ${REALPATH}/fshd1_script.sh --input-ubam reads.bam
/usr/bin/time -v ${REALPATH}/fshd1_script.sh --input-ubam  /data/hasindu/realfshd/test/data/positive/AS2603/QGXHXX240243_supmeth_unaligned.bam
# 1:47.76

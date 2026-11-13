#!/bin/bash

# Install hs1 (T2T-CHM13) reference genome and create minimap2 index

die() {
    echo "$1" >&2
    exit 1
}

[ -z ${SAMTOOLS} ] && SAMTOOLS=samtools
[ -z ${MINIMAP2} ] && MINIMAP2=minimap2

# FASTA sourced from hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
# then minimap2 indexed

rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz genome/. || die "Failed to download genome data"
gunzip genome/hs1.fa.gz || die "Failed to extract genome data"
samtools faidx genome/hs1.fa || die "Failed to index genome data"

minimap2 -d genome/hs1.mmi genome/hs1.fa || die "Failed to create minimap2 index for genome data"

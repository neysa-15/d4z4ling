#!/bin/bash

module load samtools

samtools view ${1} | awk '{
    primary_chr = $3
    for (i=12; i<=NF; i++) {
        if ($i ~ /^SA:Z:/) {
            sub(/^SA:Z:/, "", $i)
            n = split($i, sa_alignments, ";")
            for (j=1; j<=n; j++) {
                if (sa_alignments[j] == "") continue
                split(sa_alignments[j], sa_fields, ",")
                sa_chr = sa_fields[1]
                mapq = sa_fields[5]
                if (primary_chr != sa_chr && mapq >= 30) {
                    print $1, primary_chr, sa_alignments[j]
                }
            }
        }
    }
}'


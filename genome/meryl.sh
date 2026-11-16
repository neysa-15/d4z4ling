#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")

REF=${REALPATH}/hs1.fa
REF_NAME=$(basename $REF .fa)
MERYL_DB="${REALPATH}/merylDB_${REF_NAME}"
REPETITIVE_REGIONS="${REALPATH}/${REF_NAME}_repetitive_k15.txt"

# MERYL=${REALPATH}/../../install/winnowmap-v2.03/bin/meryl
MERYL=meryl

# Step 1: Pre-compute k-mer frequency (only if meryl output does not exist)
if [ ! -d "$MERYL_DB" ]; then
    echo "Creating meryl database..."
    ${MERYL} count k=15 output "$MERYL_DB" "$REF"
else
    echo "Meryl database already exists at $MERYL_DB. Skipping this step."
fi

# Step 2: Extract repetitive regions (only if {refname}_repetitive_k15.txt does not exist)
if [ ! -f "$REPETITIVE_REGIONS" ]; then
    echo "Extracting repetitive regions..."
    ${MERYL} print greater-than distinct=0.9998 "$MERYL_DB" > "$REPETITIVE_REGIONS"
else
    echo "Repetitive regions file already exists at $REPETITIVE_REGIONS. Skipping this step."
fi
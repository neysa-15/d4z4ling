#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file.tsv [target_len] [panel_len_eff] [conf]"
    exit 1
fi

TSV="$1"
TARGET_LEN="${2:-33000}"
PANEL_LEN_EFF="${3:-110214}"
CONF="${4:-0.95}"

# Extract sample name from filename
# Expecting something like SAMPLE_mapped_features_summary.tsv
sample=$(basename "$TSV" | sed 's/_mapped_features_summary\.tsv$//')

# ---- spanning ruleout function ----
spanning_ruleout() {
    local target_len="$1"
    local panel_len_eff="$2"
    local conf="$3"

    awk -v target_len="$target_len" \
        -v panel_len_eff="$panel_len_eff" \
        -v conf="$conf" '
    {
        len = $1
        diff = len - target_len
        if (diff > 0) {
            p = diff / panel_len_eff
            log_nohit += log(1 - p)
        }
    }
    END {
        if (NR == 0) {
            prob = 0
        } else {
            prob = 1 - exp(log_nohit)
        }
        meets = (prob >= conf ? "TRUE" : "FALSE")
        printf "%s\t%.10g\t%s\n", sample, prob, meets
    }' sample="$sample"
}

# ---- Extract lengths and compute ----
awk -F'\t' 'NR > 1 && $8 == "4qA" { print $3 }' "$TSV" \
    | spanning_ruleout "$TARGET_LEN" "$PANEL_LEN_EFF" "$CONF"

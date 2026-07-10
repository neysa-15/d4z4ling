#!/usr/bin/env python3
"""
Compare two d4z4ling mapped_features_summary TSVs and output a filtered TSV
of changed rows with diffs, change_status (critical/medium/low), and notes.

Use case:
- Compare results to sanity check code installation and environment.
  - Run: d4z4ling script on HG1811 subsetted data
  - Check with: python3 check_summary_tsv.py demo/result_HG01811/HG01811_mapped_features_summary.tsv {YOUR OUTPUT DIR}/HG01811_mapped_features_summary.tsv {YOUR OUTPUT DIR}/HG01811_mapped_features_summary_diffs.tsv
- Testing on different versions of tools when updating the d4z4ling pipeline.

Usage:
    python check_summary_tsv.py <control_tsv> <test_tsv> <output_tsv>
"""

import pandas as pd
import numpy as np
import sys
import re

# --------------------------------------------------
# Thresholds
# --------------------------------------------------
NUMERIC_PCT_THRESHOLD = 10  # % diff to escalate from low -> medium
COORD_BP_THRESHOLD = 10     # bp diff to escalate from low -> medium

# --------------------------------------------------
# Column type definitions
# --------------------------------------------------
ID_COLS = ['ReadID', 'strand']

COL_INT = [
    'AlignmentLength', 'MAPQ', 'ReadLength', 'XapI_Sensitive_Repeats',
    'BlnI_Sensitive_Repeats', 'gap_distance', 'd4z4_CpG_Total',
    'd4z4_CpG_Methylated', 'pLAM_CpG_Total', 'pLAM_CpG_Methylated'
]
COL_FLOAT = [
    'd4z4_chr4_proximal_score', 'd4z4_chr4_proximal_completeness',
    'p13-E11_score', 'p13-E11_completeness', 'pLAM_score', 'pLAM_completeness',
    '4qA_marker_score', '4qA_marker_completeness', '4qA_marker_percent_identity',
    '4qB_marker_score', '4qB_marker_completeness', '4qB_marker_percent_identity',
    'SSLP_length', 'MappedEstimatedCopies', 'MappedEstimatedCopiesPlus',
    'MappedEstimatedCopiesMinus', 'XapI_RE_Ratio_(%)', 'BlnI_RE_Ratio_(%)',
    'd4z4_Methylation_Percentage', 'pLAM_Methylation_Percentage'
]
COL_STRING = [
    'ReadLabel', 'Haplotype', 'duplex', 'optimal_duplex_strand',
    'd4z4_chr4_proximal_mapped', 'p13-E11_mapped', 'pLAM_mapped',
    '4qA_marker_mapped', '4qB_marker_mapped', 'pLAM_contains_polyA',
    'pLAM_polyA_signal', 'gaps', 'overlaps', 'FSHD1_status'
]
COL_COORDS = [
    'GenomeCoords', 'd4z4_chr4_proximal_coords', 'p13-E11_coords',
    'pLAM_coords', '4qA_marker_coords', '4qB_marker_coords',
    'pLAM_polyA_coords', 'SSLP_coords'
]
COL_SETS = ['overlapping_repeats', 'overlapping_repeats_coords', 'misclassification_flags']

# String cols that trigger critical if changed
CRITICAL_STRING_COLS = {'FSHD1_status', 'Haplotype'}


# --------------------------------------------------
# Diff helpers
# --------------------------------------------------

def set_diff(a, b, sep=','):
    """Return a string describing symmetric difference between two comma-sets."""
    sa = set(str(a).split(sep)) - {'NA', '', 'nan'} if pd.notna(a) else set()
    sb = set(str(b).split(sep)) - {'NA', '', 'nan'} if pd.notna(b) else set()
    added = sorted(sb - sa)
    removed = sorted(sa - sb)
    return f"added={added};removed={removed}" if (added or removed) else None


def assign_status_and_notes(row):
    """Return (change_status, notes) for a single changed row."""
    feature = row['feature']
    diffs = row.get('Diffs')
    pct = row.get('Diffs_pct')

    # --- string cols ---
    if feature in CRITICAL_STRING_COLS:
        return 'critical', f"{feature} changed: {diffs}"
    if feature in COL_STRING:
        return 'medium', f"{feature} changed: {diffs}"

    # --- numeric cols ---
    if feature in COL_INT + COL_FLOAT:
        if pd.isna(pct):
            return 'medium', f"{feature} diff={diffs} (pct unavailable, control=0)"
        status = 'low' if abs(pct) <= NUMERIC_PCT_THRESHOLD else 'medium'
        return status, f"{feature} diff={diffs} ({pct}%)"

    # --- coord cols ---
    if feature in COL_COORDS:
        if diffs is None:
            return 'low', f"{feature}: one value was NA"
        if isinstance(diffs, tuple) and len(diffs) == 2:
            ds, de = diffs
            status = 'low' if max(abs(ds), abs(de)) <= COORD_BP_THRESHOLD else 'medium'
            return status, f"{feature} start_diff={ds}, end_diff={de}"
        if isinstance(diffs, list) and len(diffs) == 3:  # genome coords
            chr_status, ds, de = diffs
            is_chr_change = chr_status == 'diff chr'
            bp_large = max(abs(ds), abs(de)) > COORD_BP_THRESHOLD
            status = 'critical' if is_chr_change else ('medium' if bp_large else 'low')
            return status, f"{feature} {chr_status}, start_diff={ds}, end_diff={de}"

    # --- set cols ---
    if feature in COL_SETS:
        return 'medium', f"{feature} set changed: {diffs}"

    return 'low', str(diffs)


# --------------------------------------------------
# Main
# --------------------------------------------------

def main(control_tsv, test_tsv, output_tsv):
    # --- Load ---
    df_control = pd.read_csv(control_tsv, sep='\t').fillna("NA")
    df_test = pd.read_csv(test_tsv, sep='\t').fillna("NA")

    # --- Validate columns ---
    control_vars = [c for c in df_control.columns if c not in ID_COLS]
    test_vars = [c for c in df_test.columns if c not in ID_COLS]
    only_control = set(control_vars) - set(test_vars)
    only_test = set(test_vars) - set(control_vars)
    if only_control:
        print(f"Warning: columns only in control: {only_control}")
    if only_test:
        print(f"Warning: columns only in test: {only_test}")

    n_categorised = len(COL_INT + COL_FLOAT + COL_STRING + COL_COORDS + COL_SETS)
    assert len(control_vars) == n_categorised, \
        f"Column count mismatch: {len(control_vars)} features vs {n_categorised} categorised"

    # --- Melt to long form ---
    shared_vars = [c for c in control_vars if c in test_vars]

    df_control_long = pd.melt(
        df_control, id_vars=ID_COLS, value_vars=shared_vars,
        var_name='feature', value_name='control_values'
    )
    df_test_long = pd.melt(
        df_test, id_vars=ID_COLS, value_vars=shared_vars,
        var_name='feature', value_name='test_values'
    )

    # --- Outer join ---
    merged = pd.merge(df_control_long, df_test_long, on=ID_COLS + ['feature'], how='outer')

    # Flag reads present in only one file
    merged['row_status'] = 'both'
    merged.loc[merged['control_values'].isna(), 'row_status'] = 'test_only'
    merged.loc[merged['test_values'].isna(), 'row_status'] = 'control_only'

    merged['same_results'] = merged['control_values'] == merged['test_values']
    merged['Diffs'] = None
    merged['Diffs_pct'] = np.nan

    # --------------------------------------------------
    # Compute Diffs per column type
    # --------------------------------------------------

    # --- Numeric ---
    mask_num = (merged['feature'].isin(COL_INT + COL_FLOAT)) & ~merged['same_results']
    test_num = pd.to_numeric(merged.loc[mask_num, 'test_values'], errors='coerce')
    ctrl_num = pd.to_numeric(merged.loc[mask_num, 'control_values'], errors='coerce')
    abs_diff = (test_num - ctrl_num).round(4)
    pct_diff = (abs_diff / ctrl_num.replace(0, np.nan) * 100).round(2)
    merged.loc[mask_num, 'Diffs'] = abs_diff
    merged.loc[mask_num, 'Diffs_pct'] = pct_diff

    # --- Coords ---
    mask_coords = (merged['feature'].isin(COL_COORDS)) & ~merged['same_results']
    temp = merged.loc[mask_coords].copy()

    # Genome coords (chr4/chr10)
    is_genome = temp['control_values'].str.contains(r'^chr(?:4|10):', na=False)

    # Simple coords (start-end only), excluding genome and NAs
    is_simple = ~is_genome & temp['control_values'].notna() & temp['test_values'].notna() \
                & (temp['control_values'] != 'NA') & (temp['test_values'] != 'NA')

    if is_simple.any():
        c = temp.loc[is_simple, 'control_values'].str.extract(r'(\d+)-(\d+)').astype(int)
        t = temp.loc[is_simple, 'test_values'].str.extract(r'(\d+)-(\d+)').astype(int)
        results = list(zip(t[0] - c[0], t[1] - c[1]))
        temp.loc[is_simple, 'Diffs'] = pd.Series(results, index=temp.index[is_simple])

    if is_genome.any():
        # Further filter to rows where both values are non-NA
        is_genome_valid = is_genome \
            & temp['control_values'].notna() & (temp['control_values'] != 'NA') \
            & temp['test_values'].notna()    & (temp['test_values']    != 'NA')

        if is_genome_valid.any():
            c_gen = temp.loc[is_genome_valid, 'control_values'].str.extract(r'(chr\w+):(\d+)-(\d+)')
            t_gen = temp.loc[is_genome_valid, 'test_values'].str.extract(r'(chr\w+):(\d+)-(\d+)')
            chr_label = np.where(c_gen[0] == t_gen[0], 'same chr', 'diff chr')
            ds = t_gen[1].astype(int) - c_gen[1].astype(int)
            de = t_gen[2].astype(int) - c_gen[2].astype(int)
            results = [list(x) for x in zip(chr_label, ds, de)]
            temp.loc[is_genome_valid, 'Diffs'] = pd.Series(results, index=temp.index[is_genome_valid])

    merged.loc[mask_coords, 'Diffs'] = temp['Diffs']

    # --- Strings ---
    mask_str = (merged['feature'].isin(COL_STRING)) & ~merged['same_results']
    merged.loc[mask_str, 'Diffs'] = (
        merged.loc[mask_str, 'control_values'].astype(str)
        + ' → '
        + merged.loc[mask_str, 'test_values'].astype(str)
    )

    # --- Sets ---
    mask_sets = (merged['feature'].isin(COL_SETS)) & ~merged['same_results']
    merged.loc[mask_sets, 'Diffs'] = merged.loc[mask_sets].apply(
        lambda r: set_diff(r['control_values'], r['test_values']), axis=1
    )

    # --------------------------------------------------
    # Assign change_status and notes (changed rows only)
    # --------------------------------------------------
    changed = merged[~merged['same_results']].copy()

    statuses, notes = zip(*changed.apply(assign_status_and_notes, axis=1)) \
        if len(changed) else ([], [])
    changed['change_status'] = statuses
    changed['notes'] = notes

    # Sort by severity then ReadID
    status_order = {'critical': 0, 'medium': 1, 'low': 2}
    changed['_sort'] = changed['change_status'].map(status_order)
    changed = changed.sort_values(['_sort', 'ReadID', 'strand']).drop(columns='_sort')

    # --------------------------------------------------
    # Output
    # --------------------------------------------------
    col_order = [
        'ReadID', 'strand', 'feature', 'row_status',
        'control_values', 'test_values',
        'Diffs', 'Diffs_pct', 'change_status', 'notes'
    ]
    changed[col_order].to_csv(output_tsv, sep='\t', index=False)

    # Check total unique ReadID that changed
    unique_reads = changed['ReadID'].nunique()
    print(f"Total unique ReadID with changes: {unique_reads}")
    print(f"Wrote {len(changed)} changed rows -> {output_tsv}")
    # only print counts if there are any changed rows
    if len(changed):
        print(changed['change_status'].value_counts().to_string())


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
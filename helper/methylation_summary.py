import pandas as pd
import argparse

# Methylation threshold
methylation_threshold = 0.75
repeat_threshold = 0.8 * 3300  # 80% of 3300 bp

# Function to get the start and end based on strand
def get_distal_d4z4_coords(df):
    if df.empty:
        return None, None

    # Extract numeric values safely as a Series
    df = df.copy()  # Avoid modifying the original dataframe
    df["d4z4_num"] = df["Feature"].str.extract(r"d4z4_(\d+)$")[0].astype(float)

    strand = df["Strand"].iloc[0]  # Assuming all rows have the same strand

    if strand == "+":  
        # Get the row with max d4z4_num
        distal_d4z4 = df.loc[df["d4z4_num"].idxmax()]

        # skip small d4z4 copies
        size = distal_d4z4["BlockSizes"].split(",")[0]
        if int(size) < repeat_threshold:
            # If the first block size is 0, use the second largest
            distal_copy = distal_d4z4["d4z4_num"] - 1
            distal_d4z4 = df[df["d4z4_num"] == distal_copy].iloc[0] if (df["d4z4_num"] == distal_copy).any() else df.loc[df["d4z4_num"].idxmax()]
    else:  
        # Get d4z4_1, or fallback to min if missing
        distal_d4z4 = df[df["d4z4_num"] == 1].iloc[0] if (df["d4z4_num"] == 1).any() else df.loc[df["d4z4_num"].idxmin()]

        size = distal_d4z4["BlockSizes"].split(",")[0]
        if int(size) < repeat_threshold:
            # If the first block size is 0, use the second largest
            distal_copy = distal_d4z4["d4z4_num"] + 1
            distal_d4z4 = df[df["d4z4_num"] == distal_copy].iloc[0] if (df["d4z4_num"] == distal_copy).any() else df.loc[df["d4z4_num"].idxmin()]

    return distal_d4z4["Start"], distal_d4z4["End"]

def calculate_distal_copy_methylation(idx, summary_df, read_methylation, d4z4_start, d4z4_end):
    # Calculate CpG statistics for the identified (overlapping or extended) d4z4 repeat
    overlapping_cpgs_d4z4 = read_methylation[(read_methylation["Start"] >= d4z4_start) & (read_methylation["End"] <= d4z4_end)]
    plam_total_count = len(overlapping_cpgs_d4z4)
    plam_methylated_count = (overlapping_cpgs_d4z4["Methylation"] >= methylation_threshold).sum()

    # Update pLAM columns in the DataFrame
    summary_df.loc[idx, "pLAM_CpG_Total"] = plam_total_count
    summary_df.loc[idx, "pLAM_CpG_Methylated"] = plam_methylated_count
    summary_df.loc[idx, "pLAM_Methylation_Percentage"] = round((plam_methylated_count / plam_total_count * 100), 2) if plam_total_count > 0 else 0.0

    return summary_df

def find_closest_coordinates(features_df, read_id, read_features, plam_start, plam_end):
    # If no d4z4 repeat overlaps, find the closest d4z4 repeat
    if not read_features.empty:
        read_features["distance"] = read_features.apply(
            lambda feature_row: min(abs(feature_row["Start"] - plam_end), abs(feature_row["End"] - plam_start)),
            axis=1
        )
        closest_d4z4 = read_features.loc[read_features["distance"].idxmin()]
        d4z4_start, d4z4_end = closest_d4z4["Start"], closest_d4z4["End"]

        # Extend the closest d4z4 repeat to include the target coordinates
        new_d4z4_start = min(d4z4_start, plam_start)
        new_d4z4_end = max(d4z4_end, plam_end)

        # Update the features DataFrame
        # features_df.loc[
        #     (features_df["ReadID"] == read_id) &
        #     (features_df["Start"] == d4z4_start) &
        #     (features_df["End"] == d4z4_end),
        #     ["Start", "End"]
        # ] = [new_d4z4_start, new_d4z4_end]

        # Update the d4z4 coordinates for CpG calculation
        return new_d4z4_start, new_d4z4_end

def methylation_summary(summary_file, features_file, methylation_file, output_file, updated_bed_file):
    # Read the files into dataframes
    summary_df = pd.read_csv(summary_file, sep="\t")
    features_df = pd.read_csv(features_file, sep="\t", header=None, names=[
        "ReadID", "Start", "End", "Feature", "Score", "Strand", "ThickStart", "ThickEnd", "RGB", "BlockCount", "BlockSizes", "BlockStarts"
    ])
    methylation_df = pd.read_csv(methylation_file, sep="\t", header=None, names=["ReadID", "Start", "End", "Dot", "Methylation"])

    # Filter features for d4z4 repeats
    # d4z4_features = features_df[features_df["Feature"].str.startswith(r"d4z4_")]
    d4z4_features = features_df[features_df["Feature"].str.match(r"d4z4_\d+$")]

    # Initialize new columns
    summary_df["d4z4_CpG_Total"] = 0
    summary_df["d4z4_CpG_Methylated"] = 0
    summary_df["d4z4_Methylation_Percentage"] = 0.0
    summary_df["pLAM_CpG_Total"] = 0
    summary_df["pLAM_CpG_Methylated"] = 0
    summary_df["pLAM_Methylation_Percentage"] = 0.0

    # Iterate through each row in the main TSV
    for idx, row in summary_df.iterrows():
        read_id = row["ReadID"]

        # Filter methylation data for this read
        read_methylation = methylation_df[methylation_df["ReadID"] == read_id]

        # Get d4z4 repeat features for this read
        read_features = d4z4_features[d4z4_features["ReadID"] == read_id].copy()

        if read_features.empty:
            continue

        # Count CpG sites overlapping d4z4 repeats
        methylated_d4z4_count = 0
        total_d4z4_cpg_count = 0

        for _, feature_row in read_features.iterrows():
            d4z4_start, d4z4_end = feature_row["Start"], feature_row["End"]

            # Find CpG sites overlapping this d4z4 repeat
            overlapping_cpgs = read_methylation[(read_methylation["Start"] >= d4z4_start) & (read_methylation["End"] <= d4z4_end)]

            # Filter CpG sites with probability > 0.75 (methylated) or < 0.25 (unmethylated)
            valid_cpgs = overlapping_cpgs[(overlapping_cpgs["Methylation"] > 0.75) | (overlapping_cpgs["Methylation"] < 0.25)]

            # Count total and methylated CpG sites
            total_d4z4_cpg_count += len(valid_cpgs)
            methylated_d4z4_count += (valid_cpgs["Methylation"] > 0.75).sum()

        # Update d4z4 columns in the DataFrame
        summary_df.loc[idx, "d4z4_CpG_Total"] = total_d4z4_cpg_count
        summary_df.loc[idx, "d4z4_CpG_Methylated"] = methylated_d4z4_count
        summary_df.loc[idx, "d4z4_Methylation_Percentage"] = round((methylated_d4z4_count / total_d4z4_cpg_count * 100), 2) if total_d4z4_cpg_count > 0 else 0.0

        # if read_id == "ea060016-69ea-44a9-9ac8-c2ce34195111":
        #     print(read_features)

        # Check if pLAM is available for this read
        if row["pLAM_mapped"]:
            plam_coords = row["pLAM_coords"]
            if pd.notna(plam_coords):
                plam_start, plam_end = map(int, plam_coords.split("-"))

                # Find the d4z4 repeat that overlaps with the pLAM
                overlapping_d4z4 = read_features[(read_features["Start"] <= plam_end) & (read_features["End"] >= plam_start)]

                if not overlapping_d4z4.empty:
                    # If a d4z4 repeat overlaps the pLAM, calculate CpG statistics for it
                    d4z4_start, d4z4_end = overlapping_d4z4.iloc[0]["Start"], overlapping_d4z4.iloc[0]["End"]

                else:
                    # If no d4z4 repeat overlaps, find the closest d4z4 repeat
                    d4z4_start, d4z4_end = find_closest_coordinates(features_df, read_id, read_features, plam_start, plam_end)

                summary_df = calculate_distal_copy_methylation(idx, summary_df, read_methylation, d4z4_start, d4z4_end)

        else:
            ###################################################################################
            # Check if it is partial distal but pLAM not mapping
            if "distal" in row["ReadLabel"] or row["ReadLabel"] == "Complete 4qB":
                # choose strands, important for duplex reads
                curr_strand = row["strand"]
                strand_read_features = read_features[read_features["Strand"] == curr_strand]

                # Then take the last copy of the d4z4 repeat
                # d4z4_{max} if strand is + and d4z4_1 if strand is -
                distal_start, distal_end = get_distal_d4z4_coords(strand_read_features)

                if distal_start is not None and distal_end is not None:
                    summary_df = calculate_distal_copy_methylation(idx, summary_df, read_methylation, distal_start, distal_end)

    # Save the updated features DataFrame to a new BED file
    # features_df.to_csv(updated_bed_file, sep="\t", index=False, header=False)
    return summary_df

def fshd1_individual_read_status(df, output_file):
    df.astype({'MappedEstimatedCopies': 'int', 'pLAM_Methylation_Percentage': 'float'}).dtypes
    
    # ----------- FSHD1_diagnostic_positive -----------
    df.loc[(df['ReadLabel'] == 'Complete 4qA') & (df['MappedEstimatedCopies'] < 11) & (df['pLAM_Methylation_Percentage'] < 50.0) & ((df['pLAM_polyA_signal'] == 'ATTAAA') | (df['pLAM_polyA_signal'] == 'TTTAAT')), 'FSHD1_status'] = 'FSHD1_diagnostic_positive'

    # ----------- FSHD1_diagnostic_negative -----------
    # Type 1: Uncontracted 4qA
    df.loc[(df['ReadLabel'] == 'Complete 4qA') & (df['MappedEstimatedCopies'] >= 11) & (df['pLAM_Methylation_Percentage'] >= 50.0) & ((df['pLAM_polyA_signal'] == 'ATTAAA') | (df['pLAM_polyA_signal'] == 'TTTAAT')), 'FSHD1_status'] = 'FSHD1_diagnostic_negative'
    # Type 2: 4qB
    df.loc[(df['ReadLabel'] == 'Complete 4qB') | ((df['ReadLabel'] == 'Partial distal 4qB') & (df['MappedEstimatedCopies'] >= 3)), 'FSHD1_status'] = 'FSHD1_diagnostic_negative'

    # -----------  FSHD1_partial_consistent_support -----------
    # Read label: Partial distal 4qA
    # 3 <= MappedEstimatedCopies <= 10
    # Distal methylation <50%
    # Permissive polyA signal (ATTAAA, TTTAAT)
    df.loc[(df['ReadLabel'] == 'Partial distal 4qA') & ((df['MappedEstimatedCopies'] >= 3) & (df['MappedEstimatedCopies'] < 11)) & (df['pLAM_Methylation_Percentage'] < 50.0) & ((df['pLAM_polyA_signal'] == 'ATTAAA') | (df['pLAM_polyA_signal'] == 'TTTAAT')), 'FSHD1_status'] = 'FSHD1_partial_consistent_support'

    # -----------  FSHD1_complete_conflicting  -----------
    # Read label: Complete 4qA
    # MappedEstimatedCopies <= 10
    # (Distal methylation >= 50%) OR (PolyA non-permissive or absent)
    df.loc[(df['ReadLabel'] == 'Complete 4qA') & (df['MappedEstimatedCopies'] < 11) & ((df['pLAM_Methylation_Percentage'] >= 50.0) | ((df['pLAM_polyA_signal'] != 'ATTAAA') & (df['pLAM_polyA_signal'] != 'TTTAAT'))), 'FSHD1_status'] = 'FSHD1_complete_conflicting'

    # -----------  FSHD1_partial_conflicting_support  -----------
    # Read label: Partial distal 4qA
    # MappedEstimatedCopies <= 10
    # (Distal methylation >= 50%) OR (PolyA non-permissive or absent)
    df.loc[(df['ReadLabel'] == 'Partial distal 4qA') & ((df['MappedEstimatedCopies'] >= 3) & (df['MappedEstimatedCopies'] < 11)) & ((df['pLAM_Methylation_Percentage'] >= 50.0) | ((df['pLAM_polyA_signal'] != 'ATTAAA') & (df['pLAM_polyA_signal'] != 'TTTAAT'))), 'FSHD1_status'] = 'FSHD1_partial_conflicting_support'

    # If didn't get any status from above then FSHD1_status = FSHD1_non_diagnostic
    df['FSHD1_status'] = df['FSHD1_status'].fillna('FSHD1_non_diagnostic')
    
    return df

def fshd2_individual_read_status(df, output_file):
    df.astype({'MappedEstimatedCopies': 'int', 'pLAM_Methylation_Percentage': 'float'}).dtypes

    # ----------- FSHD2_positive_support -----------
    # first: no 4qA contraction but hypomethylated
    df.loc[(df['ReadLabel'] == 'Complete 4qA') & (df['MappedEstimatedCopies'] >= 11) & (df['pLAM_Methylation_Percentage'] < 50.0) & ((df['pLAM_polyA_signal'] == 'ATTAAA') | (df['pLAM_polyA_signal'] == 'TTTAAT')), 'FSHD2_status'] = 'FSHD2_positive_support'

    # second: the other haplotype also hypomethylated
    df.loc[((df['Haplotype'] == '4qB') | (df['Haplotype'] == '10qA') | (df['Haplotype'] == '10qB')) & (df['MappedEstimatedCopies'] >= 3) & (df['pLAM_Methylation_Percentage'] < 50.0), 'FSHD2_status'] = 'FSHD2_positive_support'

    # If didn't get any status from above then FSHD2_status = FSHD2_non_diagnostic
    df['FSHD2_status'] = df['FSHD2_status'].fillna('FSHD2_non_diagnostic')

    # Save the updated main TSV file, replacing empty values with "NA"
    df.fillna("NA").to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process methylation data and update the main TSV file with CpG statistics.")
    parser.add_argument("--main_tsv", required=True, help="Path to the main summary TSV file (e.g., mapped_features_summary.tsv).")
    parser.add_argument("--features_bed", required=True, help="Path to the features BED file (e.g., all_features.bed).")
    parser.add_argument("--meth", required=True, help="Path to the methylation BED file (e.g., meth_reads.bed).")
    parser.add_argument("--output", required=True, help="Path to save the updated main TSV file.")
    parser.add_argument("--updated_bed", required=True, help="Path to save the updated features BED file.")

    # Parse arguments
    args = parser.parse_args()

    # File paths
    summary_file = args.main_tsv
    features_file = args.features_bed
    methylation_file = args.meth
    output_file = args.output
    updated_bed_file = args.updated_bed

    summary_df = methylation_summary(summary_file, features_file, methylation_file, output_file, updated_bed_file)
    summary_df = fshd1_individual_read_status(summary_df, output_file)
    fshd2_individual_read_status(summary_df, output_file)

    print(f"Updated main TSV file saved to {output_file}")
    print(f"Updated features BED file saved to {updated_bed_file}")
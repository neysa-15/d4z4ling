import pandas as pd
from Bio import SeqIO
import argparse


# Ensure all reads from the BED file are included in results
def ensure_all_reads_in_results(results, bed_df, features, lengths_dict):
    for read_id in bed_df["read_id"].unique():
        if read_id not in results:
            # Initialize default values for missing reads
            results[read_id] = {f"{feature}_mapped": False for feature in features}
            results[read_id].update({f"{feature}_coords": None for feature in features})
            results[read_id].update({f"{feature}_score": None for feature in features})
            results[read_id].update({f"{feature}_completeness": None for feature in features})
            results[read_id]["pLAM_contains_polyA"] = False
            results[read_id]["pLAM_polyA_coords"] = None
            results[read_id]["pLAM_polyA_signal"] = None
            results[read_id]["ReadLength"] = lengths_dict.get(read_id, "NA")
    return results

def load_fasta_to_dict(fasta_file):
    """
    Loads a FASTA file into a dictionary for efficient access.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with read IDs as keys and sequences as values.
        dict: A dictionary with read IDs as keys and sequence lengths as values.
    """
    sequences = {}
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
        lengths[record.id] = len(record.seq)
    return sequences, lengths

def read_bed_file(bed_file):
    """
    Reads a BED file and returns a DataFrame.
    
    Args:
        bed_file (str): Path to the BED file.
        
    Returns:
        pd.DataFrame: DataFrame containing the BED file data.
    """
    bed_columns = ["chrom", "start", "end", "read_id", "MAPQ", "strand"]
    return pd.read_csv(bed_file, sep="\t", names=bed_columns)

def add_bed_info_to_results(results_df, bed_df):
    """
    Adds GenomeCoords, strand, and MAPQ to the results DataFrame from the BED file.
    
    Args:
        results_df (pd.DataFrame): Results DataFrame.
        bed_df (pd.DataFrame): BED file DataFrame.
    
    Returns:
        pd.DataFrame: Updated results DataFrame.
    """
    def get_bed_info(row, column):
        bed_entry = bed_df[bed_df["read_id"] == row["ReadID"]]
        if not bed_entry.empty:
            if column == "GenomeCoords":
                chrom = bed_entry.iloc[0]["chrom"]
                start = bed_entry.iloc[0]["start"]
                end = bed_entry.iloc[0]["end"]
                return f"{chrom}:{start}-{end}"
            elif column == "AlignmentLength":
                start = bed_entry.iloc[0]["start"]
                end = bed_entry.iloc[0]["end"]
                return end - start
            return bed_entry.iloc[0][column]
        return "NA"

    results_df["GenomeCoords"] = results_df.apply(get_bed_info, axis=1, column="GenomeCoords")
    results_df["AlignmentLength"] = results_df.apply(get_bed_info, axis=1, column="AlignmentLength")
    results_df["strand"] = results_df.apply(get_bed_info, axis=1, column="strand")
    results_df["MAPQ"] = results_df.apply(get_bed_info, axis=1, column="MAPQ")
    return results_df


def calculate_estimated_copies(row):
    """
    Calculate the estimated number of copies for a read based on its label and orientation (strand),
    and round the result to 2 decimal places. Estimates are not calculated for "Complete 4qB" and "Partial distal 4qB".
    """
    try:
        # Ensure valid coordinates before splitting
        if row["pLAM_coords"] and row["pLAM_coords"] != "NA":
            pLAM_start, pLAM_end = map(int, row["pLAM_coords"].split("-"))
        else:
            pLAM_start, pLAM_end = None, None

        if row["p13-E11_coords"] and row["p13-E11_coords"] != "NA":
            p13_start, p13_end = map(int, row["p13-E11_coords"].split("-"))
        else:
            p13_start, p13_end = None, None

        read_length = int(row["ReadLength"]) if row["ReadLength"] != "NA" else None

        # Skip calculation for specific ReadLabels
        if row["ReadLabel"] in ["Complete 4qB", "Partial distal 4qB"]:
            return "NA"

        # Initialize copies as None
        copies = None

        # Determine estimated copies based on ReadLabel
        if row["ReadLabel"] in ["Complete 4qA", "Complete 10qA"]:
            if row["strand"] == "+":
                # Direct orientation
                copies = (pLAM_start - p13_end) / 3300 if pLAM_start and p13_end else None
            elif row["strand"] == "-":
                # Complement orientation
                copies = (p13_start - pLAM_end) / 3300 if p13_start and pLAM_end else None

        elif row["ReadLabel"] in ["Partial distal 4qA", "Partial distal 10qA"]:
            if row["strand"] == "+":
                # Direct orientation
                copies = pLAM_start / 3300 if pLAM_start else None
            elif row["strand"] == "-":
                # Complement orientation
                copies = (read_length - pLAM_end) / 3300 if read_length and pLAM_end else None

        elif row["ReadLabel"] == "Partial proximal Unclassified":
            if row["strand"] == "+":
                # Direct orientation
                copies = (read_length - p13_end) / 3300 if read_length and p13_end else None
            elif row["strand"] == "-":
                # Complement orientation
                copies = p13_start / 3300 if p13_start else None

        # Return rounded copies or "NA" if not calculable
        return round(copies, 2) if copies is not None else "NA"

    except Exception as e:
        print(f"Error calculating EstimatedCopies for row {row.get('ReadID', 'Unknown')}: {e}")
        return "NA"

def parse_sslp_bed(sslp_file):
    """
    Parse the SSLP BED file to extract coordinates and sequence lengths.

    Args:
        sslp_file (str): Path to the SSLP BED file.

    Returns:
        dict: A dictionary mapping read_id to (coords, length).
    """
    sslp_dict = {}
    with open(sslp_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            read_id = cols[0]
            start, end = int(cols[1]), int(cols[2])
            sequence = cols[6]
            coords = f"{start}-{end}"
            length = len(sequence)
            sslp_dict[read_id] = (coords, length)
    return sslp_dict


def find_polyA_coords_and_signal(sequence, start, end):
    """
    Find the polyA signal within a sequence, prioritising ATTAAA and its reverse complement TTTAAT.
    If not found, search for ATCAAA and ATTTAA as alternative signals.
    """
    signals = [
        ("ATTAAA", "TTTAAT"),  # Primary signals
        ("ATCAAA", "TTTGAT"),  # Alternative signal 1
        ("ATTTAA", "TTAAAT")   # Alternative signal 2
    ]

    aligned_seq = sequence[start:end]

    # Helper function to search for direct and reverse signals
    def search_signals(direct, reverse):
        direct_index = aligned_seq.find(direct)
        if direct_index != -1:
            polyA_start = start + direct_index
            polyA_end = polyA_start + len(direct)
            return f"{polyA_start}-{polyA_end}", direct

        reverse_index = aligned_seq.find(reverse)
        if reverse_index != -1:
            polyA_start = start + reverse_index
            polyA_end = polyA_start + len(reverse)
            return f"{polyA_start}-{polyA_end}", reverse

        return None, None

    # Search for signals in order of priority
    for direct, reverse in signals:
        coords, signal = search_signals(direct, reverse)
        if coords:
            return coords, signal

    # Return None if no signals are found
    return None, None



def parse_psl_to_table(psl_file, output_table, bed_file, sslp_file, fasta_file=None):
    """
    Parse a PSL file and create a table summarising alignments for each feature per read.
    
    Args:
        psl_file (str): Path to the PSL file.
        output_table (str): Path to save the output table (TSV format).
        fasta_file (str): Path to the FASTA file (optional) for pLAM sequence validation.
        
    Returns:
        None
    """
 
    sequences_dict, lengths_dict = load_fasta_to_dict(fasta_file)

    # Column names for PSL format
    psl_columns = [
        "matches", "mismatches", "rep_matches", "n_count", "q_gap_count",
        "q_gap_bases", "t_gap_count", "t_gap_bases", "strand", "q_name",
        "q_size", "q_start", "q_end", "t_name", "t_size", "t_start", "t_end",
        "block_count", "block_sizes", "q_starts", "t_starts"
    ]

    # Read the PSL file into a DataFrame
    psl_df = pd.read_csv(psl_file, sep="\t", names=psl_columns, comment="#", dtype=str)

    # Convert relevant columns to numeric types
    numeric_columns = [
        "matches", "mismatches", "rep_matches", "n_count", "q_gap_count",
        "q_gap_bases", "t_gap_count", "t_gap_bases", "q_size", "q_start",
        "q_end", "t_size", "t_start", "t_end"
    ]
    for col in numeric_columns:
        psl_df[col] = pd.to_numeric(psl_df[col], errors="coerce")

    # List of features to check
    features = ["d4z4_chr4_proximal", "p13-E11", "pLAM","4qA_probe","4qB_probe"]

    # Prepare results dictionary
    results = {}

    for _, row in psl_df.iterrows():
        # Assign target (read) and query (feature)
        read_name = row["t_name"]  # Target is the read
        feature_name = row["q_name"]  # Query is the feature
        # Check if t_start or t_end is NaN and handle gracefully
        if pd.isna(row["t_start"]) or pd.isna(row["t_end"]):
            #print(f"Skipping row due to missing t_start or t_end for ReadID: {read_name}, feature: {feature_name}")
            continue
        t_start, t_end = int(row["t_start"]), int(row["t_end"])
        alignment_score = row["matches"] - row["mismatches"] - row["q_gap_bases"] - row["t_gap_bases"]
        completeness = round((row["matches"] / row["q_size"]) * 100, 2)  # Feature completeness in percentage
        read_length = int(row["t_size"])

        # Initialize results for this read if not already present
        if read_name not in results:
            results[read_name] = {f"{feature}_mapped": False for feature in features}
            results[read_name].update({f"{feature}_coords": None for feature in features})
            results[read_name].update({f"{feature}_score": None for feature in features})
            results[read_name].update({f"{feature}_completeness": None for feature in features})
            results[read_name]["pLAM_contains_polyA"] = None
            results[read_name]["pLAM_polyA_coords"] = None
            results[read_name]["pLAM_polyA_signal"] = None
            results[read_name]["ReadLength"] = read_length
            results[read_name]["4qA_probe_percent_identity"] = "NA"
            results[read_name]["4qB_probe_percent_identity"] = "NA"

        # Populate details for the current feature
        if feature_name in features:
            # Skip this alignment if completeness is below the threshold
            COMPLETENESS_THRESHOLD = 10
            if completeness < COMPLETENESS_THRESHOLD:
                continue
            # Check if the feature is already filled; if so, skip to avoid overwriting
            if results[read_name][f"{feature_name}_mapped"]:
                continue
            results[read_name][f"{feature_name}_mapped"] = True
            results[read_name][f"{feature_name}_coords"] = f"{t_start}-{t_end}"
            results[read_name][f"{feature_name}_score"] = alignment_score
            results[read_name][f"{feature_name}_completeness"] = completeness

            # Check for polyadenylation signal if feature is pLAM
            if feature_name == "pLAM" and fasta_file:
                sequence = sequences_dict.get(read_name, "")
                if sequence:
                    #print(read_name)
                    #print(sequence)
                    polyA_coords, polyA_signal = find_polyA_coords_and_signal(sequence, t_start, t_end)
                    if polyA_signal in ["ATTAAA", "TTTAAT"]:
                        results[read_name]["pLAM_contains_polyA"] = True
                    else:
                        results[read_name]["pLAM_contains_polyA"] = False
                    results[read_name]["pLAM_polyA_coords"] = polyA_coords
                    results[read_name]["pLAM_polyA_signal"] = polyA_signal
                else:
                    results[read_name]["pLAM_contains_polyA"] = False
                    results[read_name]["pLAM_polyA_coords"] = None
                    results[read_name]["pLAM_polyA_signal"] = None
            if feature_name in ["4qA_probe", "4qB_probe"]:
                matches = row["matches"]
                mismatches = row["mismatches"]
                rep_matches = row["rep_matches"]

                # Calculate PercentIdentity
                percent_identity = (matches / (matches + mismatches + rep_matches)) * 100

                # Add to results dictionary
                if read_name not in results:
                    results[read_name] = {}
                results[read_name][f"{feature_name}_percent_identity"] = round(percent_identity, 2)

    # Ensure all BED reads are in results
    bed_df = read_bed_file(bed_file)
    results = ensure_all_reads_in_results(results, bed_df, features, lengths_dict)

    print(results)

    # Convert results dictionary to a DataFrame
    results_df = pd.DataFrame.from_dict(results, orient="index").reset_index()
    results_df.rename(columns={"index": "ReadID"}, inplace=True)

    print(results_df)

    def assign_read_label_and_haplotype(row):
        p13_mapped = row["p13-E11_mapped"]
        pLAM_mapped = row["pLAM_mapped"]
        q4b_mapped = row["4qB_probe_mapped"]
        genome_coords = row.get("GenomeCoords", "NA")  # Default to "NA" if not available
    
        # Check PercentIdentity for q4B
        q4b_identity = row.get("4qB_probe_percent_identity", 0)  # Default to 0 if not available
        q4b_high_identity = q4b_mapped and q4b_identity > 95
    
        # Check chromosomes
        starts_with_chr4 = genome_coords.startswith("chr4")
        starts_with_chr10 = genome_coords.startswith("chr10")
        #print(p13_mapped,pLAM_mapped,q4b_mapped)
        #print(starts_with_chr4,starts_with_chr10)
    
        # Determine label and haplotype
        if p13_mapped and pLAM_mapped and starts_with_chr4:
            return "Complete 4qA", "4qA"
        elif p13_mapped and q4b_mapped and q4b_high_identity and starts_with_chr4:
            return "Complete 4qB", "4qB"
        elif pLAM_mapped and starts_with_chr4:
            return "Partial distal 4qA", "4qA"
        elif q4b_mapped and q4b_high_identity and starts_with_chr4:
            return "Partial distal 4qB", "4qB"
        elif p13_mapped and starts_with_chr4:
            return "Partial proximal Unclassified", "NA"  # Haplotype is NA for this label
        elif p13_mapped and pLAM_mapped and starts_with_chr10:
            return "Complete 10qA", "10qA"
        elif pLAM_mapped and starts_with_chr10:
            return "Partial distal 10qA", "10qA"
        else:
            return "no_features", "NA"


    # Add GenomeCoords, strand, and MAPQ to results DataFrame
    results_df = add_bed_info_to_results(results_df, bed_df)

    # Apply the updated function to assign both ReadLabel and Haplotype
    print(results_df.apply(lambda row: pd.Series(assign_read_label_and_haplotype(row)), axis=1).head())
    results_df[["ReadLabel", "Haplotype"]] = results_df.apply(
        lambda row: pd.Series(assign_read_label_and_haplotype(row)), axis=1
    )
    print("ERROR?")

    # results_df["EstimatedCopies"] = results_df.apply(calculate_estimated_copies, axis=1)

    # Parse the SSLP BED file if provided
    sslp_data = {}
    if sslp_file:
        print(f"Parsing SSLP BED file: {sslp_file}")
        sslp_data = parse_sslp_bed(sslp_file)

    # Add SSLP_coords and SSLP_length columns
    if sslp_data:
        results_df["SSLP_coords"] = results_df["ReadID"].map(lambda x: sslp_data.get(x, ("NA", "NA"))[0])
        results_df["SSLP_length"] = results_df["ReadID"].map(lambda x: sslp_data.get(x, ("NA", "NA"))[1])
 
    # Rearrange columns in the desired order
    columns = ["ReadID","GenomeCoords", "AlignmentLength", "strand", "MAPQ","ReadLength", "ReadLabel", "Haplotype"]
    for feature in features:
        columns.extend([
            f"{feature}_mapped",
            f"{feature}_coords",
            f"{feature}_score",
            f"{feature}_completeness"
        ])
        # Include PercentIdentity for specific features
        if feature in ["4qA_probe", "4qB_probe"]:
            columns.append(f"{feature}_percent_identity")
    columns.extend(["pLAM_contains_polyA", "pLAM_polyA_coords", "pLAM_polyA_signal"])

    if sslp_data:
        columns.extend(["SSLP_coords", "SSLP_length"])

    results_df = results_df[columns]

    # Replace empty values with NA
    results_df.fillna("NA", inplace=True)

    # Exclude rows with NA in the ReadLabel column
    results_df["ReadLabel"].fillna("NA", inplace=True)

    # Add a temporary column for custom chromosome sorting
    results_df["ChromosomeOrder"] = results_df["GenomeCoords"].apply(lambda x: 0 if x.startswith("chr4") else 1 if x.startswith("chr10") else 2)

    # Define the custom order for ReadLabel
    read_label_order = [
        "Complete 4qA",
        "Partial distal 4qA",
        "Complete 4qB",
        "Partial distal 4qB",
        "Partial proximal Unclassified",
        "Complete 10qA",
        "Partial distal 10qA",
        "no_features"
    ]

    # Convert ReadLabel column to a categorical type with the specified order
    results_df["ReadLabel"] = pd.Categorical(
        results_df["ReadLabel"],
        categories=read_label_order,
        ordered=True
    )

    # Sort rows by ReadLabel and Haplotype
    results_df = results_df.sort_values(by=["ChromosomeOrder","ReadLabel", "Haplotype"])
 
    # Drop the temporary sorting column
    results_df = results_df.drop(columns=["ChromosomeOrder"])

    # Save to TSV
    results_df.to_csv(output_table, index=False, sep="\t")
    print(f"Results saved to {output_table}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse PSL file and generate a summary table.")
    parser.add_argument("--psl", required=True, help="Path to the PSL file.")
    parser.add_argument("--bed", required=True, help="Path to the BED file.")
    parser.add_argument("--output", required=True, help="Path to the output TSV file.")
    parser.add_argument("--fasta", required=False, help="Path to the FASTA file for polyA signal checking.")
    parser.add_argument("--sslp", help="Path to SSLP BED file.", default=None)
    
    args = parser.parse_args()

    # Parse PSL to table
    parse_psl_to_table(
        psl_file=args.psl,
        bed_file=args.bed,
        output_table=args.output,
        fasta_file=args.fasta,
        sslp_file=args.sslp
    )

import pandas as pd
from Bio import SeqIO
import argparse

from copy_number_estimation import get_cne

# Ensure all reads from the BED file are included in results
def initialise_results(results, bed_df, features, lengths_dict):
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
            results[read_id]["4qA_probe_percent_identity"] = "NA"
            results[read_id]["4qB_probe_percent_identity"] = "NA"
            results[read_id]["duplex"] = False
            results[read_id]["optimal_duplex_strand"] = None
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

def read_features_bed_file(features_bed):
    """
    Reads a BED file containing features and returns a DataFrame.
    
    Args:
        features_bed (str): Path to the features BED file.
        
    Returns:
        pd.DataFrame: DataFrame containing the features BED file data.
    """
    features_columns = ["read_id", "start", "end", "feature_name", "MAPQ", "strand"]
    return pd.read_csv(features_bed, sep="\t", names=features_columns)

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

def check_duplex(read_name, psl_df, results_df, COMPLETENESS_THRESHOLD=10):
    """
    Checking if a read is duplex

    Args:
        read_name (str): t_name of psl
        psl_df (df): df of features

    Returns:
        boolean
    """
    strands = set()
    for _, f_row in psl_df.iterrows():
        completeness = round((f_row["matches"] / f_row["q_size"]) * 100, 2)
        if completeness < COMPLETENESS_THRESHOLD:
            continue

        strands.add(f_row["strand"])

    read_row = results_df[results_df["ReadID"] == read_name]
    mapped_estimated_copies_plus = read_row["MappedEstimatedCopiesPlus"].values[0]
    mapped_estimated_copies_minus = read_row["MappedEstimatedCopiesMinus"].values[0]
    if (mapped_estimated_copies_plus > 0.00) and (mapped_estimated_copies_minus > 0.00):
        return True

    return len(strands) > 1

def process_duplex_row(row, psl_df, results_df, features, COMPLETENESS_THRESHOLD=10):
    """
    Process a row: determine if it's duplex and duplicate if needed.

    Args:
        row (pd.Series): A row from results_df
        psl_df (pd.DataFrame): Feature alignment dataframe
        results_df (pd.DataFrame): The original results dataframe
        features (list): List of features
        COMPLETENESS_THRESHOLD (int): Minimum completeness percentage

    Returns:
        list of pd.Series: Either one or two rows, depending on duplex status
    """
    read_name = row["ReadID"]
    filtered_df = psl_df.loc[psl_df["t_name"] == read_name]

    # if filtered_df.empty:
    #     row["duplex"] = False
    #     row["optimal duplex strand"] = None
    #     return [row]

    is_duplex = check_duplex(read_name, filtered_df, results_df)

    if not is_duplex:
        row["duplex"] = False
        # row["optimal duplex strand"] = filtered_df["strand"].iloc[0]
        row["optimal_duplex_strand"] = row["strand"]
        return [row]

    # If duplex, determine the optimal strand
    pos_feature, neg_feature = set(), set()

    for _, f_row in filtered_df.iterrows():
        completeness = round((f_row["matches"] / f_row["q_size"]) * 100, 2)
        if completeness < COMPLETENESS_THRESHOLD:
            continue

        if f_row["q_name"] in features:
            if f_row["strand"] == '+':
                pos_feature.add(f_row["q_name"])
            else:
                neg_feature.add(f_row["q_name"])

    if read_name == 'c95c3baa-460b-4c09-a80c-4e2a229603d0':
        print(f"POS {pos_feature}")
        print(f"NEG {neg_feature}")

    if len(pos_feature) == len(neg_feature):
        read_row = results_df.loc[results_df["ReadID"] == read_name]
        optimal_strand = '+' if float(read_row["MappedEstimatedCopiesPlus"].values[0]) > \
                               float(read_row["MappedEstimatedCopiesMinus"].values[0]) else '-'
    elif len(pos_feature) > len(neg_feature):
        optimal_strand = '+'
    else:
        optimal_strand = '-'

    # Create two rows (one for each strand)
    row["duplex"] = True
    row["optimal_duplex_strand"] = optimal_strand
    row["strand"] = optimal_strand

    row_dupe = row.copy()
    row_dupe["strand"] = '+' if optimal_strand == '-' else '-'

    print(row["ReadID"], row["strand"], row["duplex"], row["optimal_duplex_strand"])

    return [row, row_dupe]  # Returning a list of rows

# def populate_features(psl_df, features_df, results_df, features, fasta_file, sequences_dict, COMPLETENESS_THRESHOLD=10):
#     # Set MultiIndex on ReadID and strand for uniqueness
#     results_df.set_index(["ReadID", "strand"], inplace=True)

#     for _, row in features_df.iterrows():
#         read_name = row["read_id"]
#         feature_name = row["feature_name"]
#         strand = row["strand"]
#         start = row["start"]
#         end = row["end"]

#         # Skip row if 't_start' or 't_end' is NaN
#         if pd.isna(start) or pd.isna(end):
#             continue

#         if (read_name, strand) not in results_df.index:
#             continue

#         #TODO
#         if feature_name == "pLAM":
#             continue

#         # Skip if feature_name is not in the list
#         if feature_name not in features:
#             continue

#         # if results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True:
#         if (results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True) and (feature_name != "d4z4_chr4_proximal"):
#             continue

#         # Update feature mapping details
#         results_df.at[(read_name, strand), f"{feature_name}_mapped"] = True
#         results_df.at[(read_name, strand), f"{feature_name}_coords"] = f"{start}-{end}"

#         # Handle percent identity for 4qA_probe and 4qB_probe
#         # if feature_name in ["4qA_probe", "4qB_probe"]:
#         #     matches = row["matches"]
#         #     mismatches = row["mismatches"]
#         #     rep_matches = row["rep_matches"]
#         #     percent_identity = (matches / (matches + mismatches + rep_matches)) * 100
#         #     results_df.at[(read_name, strand), f"{feature_name}_percent_identity"] = round(percent_identity, 2)

#     # OLD LOOPS
#     # MAYBE CHANGE TO JUST PROCESSING PLAM DF
#     for _, row in psl_df.iterrows():
#         read_name = row["t_name"]  # Target is the read
#         feature_name = row["q_name"]  # Query is the feature

#         # Skip row if 't_start' or 't_end' is NaN
#         if pd.isna(row["t_start"]) or pd.isna(row["t_end"]):
#             continue

#         t_start, t_end = int(row["t_start"]), int(row["t_end"])
#         alignment_score = row["matches"] - row["mismatches"] - row["q_gap_bases"] - row["t_gap_bases"]
#         completeness = round((row["matches"] / row["q_size"]) * 100, 2)  # Feature completeness in percentage
#         read_length = int(row["t_size"])
#         strand = row["strand"]

#         if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a' and feature_name == "d4z4_chr4_proximal":
#             print("f30f0390-4265-4c0a-a90c-956928fa701a")
#             print(f"completeness {completeness}")
#             # print(feature_name)

#         # Ensure the read-strand combination exists in results_df
#         if (read_name, strand) not in results_df.index:
#             continue

#         # Update ReadLength if not set
#         if pd.isna(results_df.at[(read_name, strand), "ReadLength"]):
#             results_df.at[(read_name, strand), "ReadLength"] = read_length

#         # Skip if feature_name is not in the list
#         if feature_name not in features:
#             continue

#         # Skip this alignment if completeness is below threshold
#         if completeness < COMPLETENESS_THRESHOLD:
#             continue

#         # # If the feature is already mapped, skip to avoid overwriting
#         # print("FEATURE")
#         # print(results_df.at[(read_name, strand), f"{feature_name}_mapped"])
#         # print("----")

#         # if results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True:
#         if (results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True) and (feature_name != "d4z4_chr4_proximal"):
#             continue

#         # To combine fragmented d4z4_chr4_proximal, taking the coordinates that covers the most
#         # For now it's currently only changing the coordinate (not the other components)
#         if (results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True) and (feature_name == "d4z4_chr4_proximal") and (results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"] is not None):
#             prev_tstart = results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"].split("-")[0]
#             prev_tend = results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"].split("-")[1]

#             prev_tstart = int(prev_tstart)
#             prev_tend = int(prev_tend)

#             if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a':
#                 print("f30f0390-4265-4c0a-a90c-956928fa701a")
#                 print(f"prev {prev_tstart}-{prev_tend}")
#                 print(f"curr {t_start}-{t_end}")

#             if prev_tstart < t_start:
#                 t_start = prev_tstart
#             if prev_tend > t_end:
#                 t_end = prev_tend

#             if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a':
#                 print(f"Final {t_start}-{t_end}")

#             results_df.at[(read_name, strand), f"{feature_name}_coords"] = f"{t_start}-{t_end}"
            
#             continue

#         # Update feature mapping details
#         results_df.at[(read_name, strand), f"{feature_name}_mapped"] = True
#         results_df.at[(read_name, strand), f"{feature_name}_coords"] = f"{t_start}-{t_end}"
#         results_df.at[(read_name, strand), f"{feature_name}_score"] = alignment_score
#         results_df.at[(read_name, strand), f"{feature_name}_completeness"] = completeness

#         # Handle polyA signal for pLAM
#         if feature_name == "pLAM" and fasta_file:
#             sequence = sequences_dict.get(read_name, "")
#             if sequence:
#                 polyA_coords, polyA_signal = find_polyA_coords_and_signal(sequence, t_start, t_end)
#                 results_df.at[(read_name, strand), "pLAM_contains_polyA"] = polyA_signal in ["ATTAAA", "TTTAAT"]
#                 results_df.at[(read_name, strand), "pLAM_polyA_coords"] = polyA_coords
#                 results_df.at[(read_name, strand), "pLAM_polyA_signal"] = polyA_signal
#             else:
#                 results_df.at[(read_name, strand), "pLAM_contains_polyA"] = False

#         # Handle percent identity for 4qA_probe and 4qB_probe
#         if feature_name in ["4qA_probe", "4qB_probe"]:
#             matches = row["matches"]
#             mismatches = row["mismatches"]
#             rep_matches = row["rep_matches"]
#             percent_identity = (matches / (matches + mismatches + rep_matches)) * 100
#             results_df.at[(read_name, strand), f"{feature_name}_percent_identity"] = round(percent_identity, 2)

#     # Reset index before returning
#     results_df.reset_index(inplace=True)
#     return results_df

def populate_features(psl_df, results_df, features, fasta_file, sequences_dict, COMPLETENESS_THRESHOLD=10):
    # Set MultiIndex on ReadID and strand for uniqueness
    results_df.set_index(["ReadID", "strand"], inplace=True)

    for _, row in psl_df.iterrows():
        read_name = row["t_name"]  # Target is the read
        feature_name = row["q_name"]  # Query is the feature

        # Skip row if 't_start' or 't_end' is NaN
        if pd.isna(row["t_start"]) or pd.isna(row["t_end"]):
            continue

        t_start, t_end = int(row["t_start"]), int(row["t_end"])
        alignment_score = row["matches"] - row["mismatches"] - row["q_gap_bases"] - row["t_gap_bases"]
        completeness = round((row["matches"] / row["q_size"]) * 100, 2)  # Feature completeness in percentage
        read_length = int(row["t_size"])
        strand = row["strand"]

        if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a' and feature_name == "d4z4_chr4_proximal":
            print("f30f0390-4265-4c0a-a90c-956928fa701a")
            print(f"completeness {completeness}")
            # print(feature_name)

        # Ensure the read-strand combination exists in results_df
        if (read_name, strand) not in results_df.index:
            continue

        # Update ReadLength if not set
        if pd.isna(results_df.at[(read_name, strand), "ReadLength"]):
            results_df.at[(read_name, strand), "ReadLength"] = read_length

        # Skip if feature_name is not in the list
        if feature_name not in features:
            continue

        # Skip this alignment if completeness is below threshold
        if completeness < COMPLETENESS_THRESHOLD:
            continue

        # # If the feature is already mapped, skip to avoid overwriting
        # print("FEATURE")
        # print(results_df.at[(read_name, strand), f"{feature_name}_mapped"])
        # print("----")

        # if results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True:
        if (results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True) and (feature_name != "d4z4_chr4_proximal"):
            continue

        # To combine fragmented d4z4_chr4_proximal, taking the coordinates that covers the most
        # For now it's currently only changing the coordinate (not the other components)
        if (results_df.at[(read_name, strand), f"{feature_name}_mapped"] == True) and (feature_name == "d4z4_chr4_proximal") and (results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"] is not None):
            prev_tstart = results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"].split("-")[0]
            prev_tend = results_df.at[(read_name, strand), "d4z4_chr4_proximal_coords"].split("-")[1]

            prev_tstart = int(prev_tstart)
            prev_tend = int(prev_tend)

            if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a':
                print("f30f0390-4265-4c0a-a90c-956928fa701a")
                print(f"prev {prev_tstart}-{prev_tend}")
                print(f"curr {t_start}-{t_end}")

            if prev_tstart < t_start:
                t_start = prev_tstart
            if prev_tend > t_end:
                t_end = prev_tend

            if read_name == 'f30f0390-4265-4c0a-a90c-956928fa701a':
                print(f"Final {t_start}-{t_end}")

            results_df.at[(read_name, strand), f"{feature_name}_coords"] = f"{t_start}-{t_end}"
            
            continue

        # Update feature mapping details
        results_df.at[(read_name, strand), f"{feature_name}_mapped"] = True
        results_df.at[(read_name, strand), f"{feature_name}_coords"] = f"{t_start}-{t_end}"
        results_df.at[(read_name, strand), f"{feature_name}_score"] = alignment_score
        results_df.at[(read_name, strand), f"{feature_name}_completeness"] = completeness

        # Handle polyA signal for pLAM
        if feature_name == "pLAM" and fasta_file:
            sequence = sequences_dict.get(read_name, "")
            if sequence:
                polyA_coords, polyA_signal = find_polyA_coords_and_signal(sequence, t_start, t_end)
                results_df.at[(read_name, strand), "pLAM_contains_polyA"] = polyA_signal in ["ATTAAA", "TTTAAT"]
                results_df.at[(read_name, strand), "pLAM_polyA_coords"] = polyA_coords
                results_df.at[(read_name, strand), "pLAM_polyA_signal"] = polyA_signal
            else:
                results_df.at[(read_name, strand), "pLAM_contains_polyA"] = False

        # Handle percent identity for 4qA_probe and 4qB_probe
        if feature_name in ["4qA_probe", "4qB_probe"]:
            matches = row["matches"]
            mismatches = row["mismatches"]
            rep_matches = row["rep_matches"]
            percent_identity = (matches / (matches + mismatches + rep_matches)) * 100
            results_df.at[(read_name, strand), f"{feature_name}_percent_identity"] = round(percent_identity, 2)

    # Reset index before returning
    results_df.reset_index(inplace=True)
    return results_df

def read_psl(psl_file):
    """
    Read a PSL file

    Args:
        psl_file (str): Path to the PSL file

    Returns:
        psl_df (pd.DataFrame): DataFrame containing the PSL data
    """
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

    return psl_df

def read_classification(psl_file, output_table, bed_file, sslp_file, bam_file, fasta_file=None):
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

    psl_df = read_psl(psl_file)

    # List of features to check
    features = ["d4z4_chr4_proximal", "p13-E11", "pLAM", "4qA_probe", "4qB_probe"]

    # Prepare results dictionary
    results = {}

    # Ensure all BED reads are in results
    bed_df = read_bed_file(bed_file)
    results = initialise_results(results, bed_df, features, lengths_dict)

    # Read features BED file
    # features_df = read_features_bed_file(features_bed)

    # Convert results dictionary to a DataFrame
    results_df = pd.DataFrame.from_dict(results, orient="index").reset_index()
    results_df.rename(columns={"index": "ReadID"}, inplace=True)

    # Add GenomeCoords, strand, and MAPQ to results DataFrame
    results_df = add_bed_info_to_results(results_df, bed_df)

    results_df = get_cne(results_df, bam_file)

    mapped_estimated_copies_column = ['MappedEstimatedCopies', 'MappedEstimatedCopiesPlus', 'MappedEstimatedCopiesMinus']
    for col in mapped_estimated_copies_column:
        results_df[col] = pd.to_numeric(results_df[col], errors="coerce")

    # results_df[["duplex", "optimal_duplex_strand"]] = results_df.apply(
    #     lambda row: select_strand_if_duplex(row['ReadID'], psl_df, results_df, features), axis=1
    # )

    # **Applying process_duplex_row correctly**
    expanded_rows = []
    for _, row in results_df.iterrows():
        expanded_rows.extend(process_duplex_row(row, psl_df, results_df, features))

    # Create a new DataFrame
    results_df = pd.DataFrame(expanded_rows).reset_index(drop=True)

    # print(results_df)
    # results_df.to_csv(output_table, index=False, sep="\t")
    # return
    
    # Populate features from psl file
    results_df = populate_features(psl_df, results_df, features, fasta_file, sequences_dict, COMPLETENESS_THRESHOLD=10)

    def assign_read_label_and_haplotype(row):
        p13_mapped = row["p13-E11_mapped"]
        pLAM_mapped = row["pLAM_mapped"]
        q4b_mapped = row["4qB_probe_mapped"]
        duplex = row["duplex"]
        genome_coords = row.get("GenomeCoords", "NA")  # Default to "NA" if not available
    
        # Check PercentIdentity for q4B
        q4b_identity = row.get("4qB_probe_percent_identity", 0)  # Default to 0 if not available
        q4b_high_identity = q4b_mapped and q4b_identity > 95
    
        # Check chromosomes
        starts_with_chr4 = genome_coords.startswith("chr4")
        starts_with_chr10 = genome_coords.startswith("chr10")

        # Determine label and haplotype
        if p13_mapped and pLAM_mapped and starts_with_chr4:
            return "Complete 4qA", "4qA"
        elif p13_mapped and q4b_mapped and q4b_high_identity and starts_with_chr4:  
            return "Complete 4qB", "4qB"
        elif p13_mapped and pLAM_mapped and starts_with_chr10:
            return "Complete 10qA", "10qA"
        elif pLAM_mapped and starts_with_chr4:
            return "Partial distal 4qA", "4qA"
        elif q4b_mapped and q4b_high_identity and starts_with_chr4:  
            return "Partial distal 4qB", "4qB"
        elif pLAM_mapped and starts_with_chr10:
            return "Partial distal 10qA", "10qA"
        elif q4b_mapped  and q4b_high_identity and starts_with_chr10:  
            return "Partial distal 10qB", "10qB"
        elif p13_mapped and (starts_with_chr4 or starts_with_chr10):
            return "Partial proximal Unclassified", "NA"  # Haplotype is NA for this label
        else:
            # Classify no feature reads
            try:
                # Calculate expected read length based on MappedEstimatedCopies
                mapped_estimated_copies = float(row.get("MappedEstimatedCopies", 0))  # Default to 0 if missing
                expected_length_kb = mapped_estimated_copies * 3.3  # Convert to kb by multiplying by 3.3
                actual_length_kb = float(row.get("ReadLength", 0)) / 1000  # Convert ReadLength to kb

                # Compare expected length with actual length (+/- 3.3)
                if (actual_length_kb - 3.3 <= expected_length_kb <= actual_length_kb + 3.3) or duplex:
                    return "Partial internal Unclassified", "NA"
                else:
                    # ADD if it's partial distal unclassified, if GenomeCoords starts with chr4 then 4qB or chr10 then 10qB
                    if row["GenomeCoords"].startswith("chr4"):
                        return "Partial distal 4qB", "4qB"
                    elif row["GenomeCoords"].startswith("chr10"):
                        return "Partial distal 10qB", "10qB"
                    # print(row["ReadID"], row["ReadLabel"])

            except ValueError as e:
                print(f"Error processing row {row.get('ReadID', 'Unknown')}: {e}")

    # Apply the updated function to assign both ReadLabel and Haplotype
    # print(results_df.apply(lambda row: pd.Series(assign_read_label_and_haplotype(row)), axis=1).head())
    results_df[["ReadLabel", "Haplotype"]] = results_df.apply(
        lambda row: pd.Series(assign_read_label_and_haplotype(row)), axis=1
    )

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
    columns = ["ReadID","GenomeCoords", "AlignmentLength", "strand", "MAPQ","ReadLength", "ReadLabel", "Haplotype", "duplex", "optimal_duplex_strand"]
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

    columns.extend(["MappedEstimatedCopies", "MappedEstimatedCopiesPlus", "MappedEstimatedCopiesMinus"])
    # print(results_df)

    results_df = results_df[columns]

    # Replace empty values with NA
    results_df.fillna("NA", inplace=True)

    # Exclude rows with NA in the ReadLabel column
    # results_df["ReadLabel"].fillna("NA", inplace=True)
    results_df["ReadLabel"] = results_df["ReadLabel"].fillna("NA")

    # Add a temporary column for custom chromosome sorting
    results_df["ChromosomeOrder"] = results_df["GenomeCoords"].apply(lambda x: 0 if x.startswith("chr4") else 1 if x.startswith("chr10") else 2)

    # Define the custom order for ReadLabel
    read_label_order = [
        "Complete 4qA",
        "Partial distal 4qA",
        "Complete 4qB",
        "Partial distal 4qB",
        "Partial proximal Unclassified",
        "Partial internal Unclassified",
        "Complete 10qA",
        "Partial distal 10qA",
        "Complete 10qB",
        "Partial distal 10qB",
        # "no_features"
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
    parser.add_argument("--bed", required=True, help="Path to the reads of interest BED file.")
    parser.add_argument("--output", required=True, help="Path to the output TSV file.")
    parser.add_argument("--fasta", required=False, help="Path to the FASTA file for polyA signal checking.")
    parser.add_argument("--sslp", help="Path to SSLP BED file.", default=None)
    parser.add_argument("--aligned-bam", help="Path to merged haplotype bam file.", default=None)
    
    args = parser.parse_args()

    # Parse PSL to table
    read_classification(
        psl_file=args.psl,
        bed_file=args.bed,
        output_table=args.output,
        fasta_file=args.fasta,
        sslp_file=args.sslp,
        bam_file=args.aligned_bam
    )

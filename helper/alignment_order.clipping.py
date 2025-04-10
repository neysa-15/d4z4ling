import pysam
import argparse
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re

def read_tsv(file_path):
    """
    Read a TSV file and return a list of dictionaries.
    Each dictionary represents a row, with keys as column headers.
    """
    with open(file_path, "r") as file:
        header = file.readline().strip().split("\t")
        return [dict(zip(header, line.strip().split("\t"))) for line in file]

def count_sensitive_repeats(d4z4_table):
    """
    Count the number of repeat units with XapI_Sites and BlnI_Sites per read.
    """
    xapi_sensitive_counts = {}
    blni_sensitive_counts = {}

    for row in d4z4_table:
        read_id = row["ReadID"]
        if row["XapI_Sites"] != "NA":
            xapi_sensitive_counts[read_id] = xapi_sensitive_counts.get(read_id, 0) + 1
        if row["BlnI_Sites"] != "NA":
            blni_sensitive_counts[read_id] = blni_sensitive_counts.get(read_id, 0) + 1

    # xapi_ratio = round(xapi_sensitive_counts.get(read_id, 0) / (xapi_sensitive_counts.get(read_id, 0) + blni_sensitive_counts.get(read_id, 0)), 2)

    return xapi_sensitive_counts, blni_sensitive_counts #, xapi_ratio, (100 - xapi_ratio)

def update_main_tsv(main_tsv, xapi_counts, blni_counts):
    """
    Update the main TSV with counts of sensitive repeats.
    """
    updated_data = []
    with open(main_tsv, "r") as file:
        header = file.readline().strip().split("\t")
        rows = [dict(zip(header, line.strip().split("\t"))) for line in file]
        for row in rows:
            read_id = row["ReadID"]
            row["XapI_Sensitive_Repeats"] = xapi_counts.get(read_id, 0)
            row["BlnI_Sensitive_Repeats"] = blni_counts.get(read_id, 0)

            # Get ratio
            if (row["XapI_Sensitive_Repeats"] == 0) and (row["BlnI_Sensitive_Repeats"] == 0):
                row["XapI_RE_Ratio_(%)"] = 0
                row["BlnI_RE_Ratio_(%)"] = 0
            else:
                row["XapI_RE_Ratio_(%)"] = round((row["XapI_Sensitive_Repeats"] / (row["XapI_Sensitive_Repeats"] + row["BlnI_Sensitive_Repeats"])) * 100, 2)
                row["BlnI_RE_Ratio_(%)"] = 100 - row["XapI_RE_Ratio_(%)"]
            updated_data.append(row)

    # Add new columns to the header
    if "XapI_Sensitive_Repeats" not in header:
        header.append("XapI_Sensitive_Repeats")
    if "BlnI_Sensitive_Repeats" not in header:
        header.append("BlnI_Sensitive_Repeats")
    if "XapI_RE_Ratio_(%)" not in header:
        header.append("XapI_RE_Ratio_(%)")
    if "BlnI_RE_Ratio_(%)" not in header:
        header.append("BlnI_RE_Ratio_(%)")

    # Write the updated table
    with open(main_tsv, "w") as file:
        file.write("\t".join(header) + "\n")
        for row in updated_data:
            file.write("\t".join(str(row.get(col, "")) for col in header) + "\n")

def detect_restriction_sites(sequence):
    """
    Detect restriction sites for XapI (R^AATTY) and BlnI (CCTAGG) in both strands of a genome.
    
    Args:
        sequence (str): The genome sequence to search.
 
    Returns:
        dict: A dictionary with enzyme names as keys and lists of positions (start, end, strand) as values.
    """
    # Define patterns for the enzymes
    patterns = {
        "XapI": r"[AG]AATT[CT]",       # Forward pattern
        "BlnI": r"CCTAGG"              # Palindromic, same for reverse
    }
 
    # Reverse complement patterns for non-palindromic enzymes
    reverse_patterns = {
        "XapI": r"[AG]TTAAT[CT]",      # Reverse complement of XapI
        "BlnI": r"CCTAGG"              # Same as forward for BlnI
    }
 
    # Dictionary to store results
    restriction_sites = {enzyme: [] for enzyme in patterns}
 
    for enzyme, pattern in patterns.items():
        # Search for the forward pattern
        for match in re.finditer(pattern, sequence):
            start = match.start()
            end = match.end()
            restriction_sites[enzyme].append((start, end, "+"))
 
        # Search for the reverse complement pattern in the forward strand
        for match in re.finditer(reverse_patterns[enzyme], sequence):
            start = match.start()
            end = match.end()
            restriction_sites[enzyme].append((start, end, "-"))
 
    return restriction_sites

def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence.
    """
    return str(Seq(sequence).reverse_complement())

def get_clipping(cigar, strand):
    """
    Rearrange the clipping for orientation based on the strand.
    """
    if cigar == "*" or not isinstance(cigar, str) or cigar == "":
        return [0, 0]  # Return default values for unmapped or invalid CIGAR strings

    left_match = re.match(r'^(\d+)([SH])', cigar)
    right_match = re.search(r'(\d+)([SH])$', cigar)

    # Extract values or default to 0 if not found
    left = int(left_match.group(1)) if left_match else 0
    right = int(right_match.group(1)) if right_match else 0

    return [left, right] if strand == '+' else [right, left]

def calculate_alignment_length(cigar):
    """
    Calculate alignment length from the CIGAR string.
    """
    matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    length = sum(int(length) for length, op in matches if op in "M=X")  # Count only matching positions
    return length

def extract_aligned_sequence(read_sequence, cigar, query_start):
    """
    Extract the sequence from the read based on the CIGAR string and alignment start position.
    """
    aligned_sequence = []
    query_start = 0
    cursor = query_start

    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        length, op = int(match.group(1)), match.group(2)
        if op in "M=X":  # Matching segment
            aligned_sequence.append(read_sequence[cursor:cursor + length])
            cursor += length
        elif op in "I":  # Insertion in the read
            aligned_sequence.append(read_sequence[cursor:cursor + length])
            cursor += length
        elif op in "S":  # Soft clipping
            cursor += length
        elif op in "DNP":  # Deletion, padding, or skipped reference
            continue
    return "".join(aligned_sequence)

def get_repeat_coords(repeat_sequence, original_sequence):
    """
    Get coordinates of repeats in the original sequence
    """
    if repeat_sequence in original_sequence:
        result = re.search(repeat_sequence, original_sequence)
        start_coords, end_coords = result.span()
    elif repeat_sequence in str(Seq(original_sequence).reverse_complement()):
        final_coord = len(original_sequence)
        result = re.search(repeat_sequence, str(Seq(original_sequence).reverse_complement()))
        start, end = result.span()
        start_coords, end_coords = final_coord - end, final_coord - start
    else:
        return "NA", "NA"
    
    return start_coords, end_coords

def get_order(read, original_sequence):
    """
    Determine the order of supplementary alignments for the given read, sorted by clipping values.
    """
    if read.cigarstring == "*" or not read.cigarstring:
        print(f"Warning: Read {read.query_name} is unmapped or lacks a valid CIGAR string.")
        return []  # Skip processing for unmapped reads

    sa_tag = read.get_tag("SA") if read.has_tag("SA") else None
    alignment_details = []

    # Add the primary alignment
    primary_strand = "+" if not read.is_reverse else "-"
    cigar = read.cigarstring
    clipping = get_clipping(cigar, primary_strand)

    primary = {
        "chrom": read.reference_name,
        "start": read.reference_start,
        "end": read.reference_end,
        "alignment_length": read.reference_length,
        "strand": primary_strand,
        "cigar": cigar,
        "clipping": clipping,
    }
    alignment_details.append(primary)

    # Track strands across all alignments
    all_strands = {primary_strand}

    if sa_tag:
        # Parse SA tag: Format -> "chr,start,strand,CIGAR,mapQ,NM;..."
        supplementary_alignments = [sa for sa in sa_tag.split(";") if sa]

        for alignment in supplementary_alignments:
                           
            alignment_info = alignment.split(",")
            chrom = alignment_info[0]
            start = int(alignment_info[1]) - 1  # Convert to 0-based
            strand = alignment_info[2]
            cigar = alignment_info[3]
            all_strands.add(strand)

            alignment_details.append({
                "chrom": chrom,
                "start": start,
                "end": start + calculate_alignment_length(cigar),
                "alignment_length": calculate_alignment_length(cigar),
                "strand": strand,
                "cigar": cigar,
                "clipping": get_clipping(cigar, strand),
            })

    # Reverse complement if alignments exist on both strands
    both_strands = len(all_strands) > 1

    # Extract sequences for each alignment
    for alignment in alignment_details:
        strand = alignment["strand"]
        sequence_to_use = read.query_sequence
        if both_strands and strand == "+" and str(Seq(read.query_sequence).reverse_complement()) == original_sequence:
            sequence_to_use = reverse_complement(read.query_sequence)
        elif both_strands and strand == "-" and read.query_sequence == original_sequence:
            sequence_to_use = reverse_complement(read.query_sequence)
        alignment["sequence"] = extract_aligned_sequence(sequence_to_use, alignment["cigar"], 0)

    # Sort alignments by clipping values
    alignment_details.sort(key=lambda x: (x["clipping"][0], x["clipping"][1]))

    if read.query_name == "6ee9eeb3-3ae4-4176-bc64-d6c76a6c7849":
        for idx, alignment in enumerate(alignment_details, start=1):
            print(f"  Alignment {idx}:")
            print(f"    Strand: {alignment['strand']}")
            print(f"    CIGAR: {alignment['cigar']}")
            print(f"    Clipping: {alignment['clipping']}")
            print(f"    Sequence: {alignment['sequence']}")

    return alignment_details


def process_bam(bam_file, read_id, fasta_file, table_file, fasta_out_file, xapi_bed, blni_bed, repeats_bed):
    """
    Process a BAM file and determine alignment details for a specific read ID,
    outputting aligned sequences to a FASTA file and alignment details to a table file.
    """
    # Parse the input FASTA file to get the original read sequences
    original_reads = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    if read_id not in original_reads:
        print(f"Warning: Read {read_id} not found in input FASTA.")
        return

    original_sequence = original_reads[read_id]

    # Detect restriction sites
    restriction_sites = detect_restriction_sites(original_sequence)

    with pysam.AlignmentFile(bam_file, "rb") as bam, open(fasta_out_file, "a") as fasta_out, open(table_file, "a") as table_out, open(xapi_bed, "a") as xapi_bed_out, open(blni_bed, "a") as blni_bed_out, open(repeats_bed, "a") as repeats_bed_out:

        # Write restriction sites to BED files
        for start, end, orientation in restriction_sites["XapI"]:
            xapi_bed_out.write(f"{read_id}\t{start}\t{end}\t{orientation}\n")
        for start, end, orientation in restriction_sites["BlnI"]:
            blni_bed_out.write(f"{read_id}\t{start}\t{end}\t{orientation}\n")

        for read in bam.fetch(until_eof=True):
            if read.query_name == read_id and not read.is_supplementary:
                alignments = get_order(read, original_sequence)
                for idx, alignment in enumerate(alignments, 1):
                    aligned_seq = alignment["sequence"]
                    strand = alignment["strand"]

                    # get d4z4 repeat coordinates
                    start_pos, end_pos = get_repeat_coords(aligned_seq, original_sequence)

                    # Write coordinates of repeat to output BED file
                    repeats_bed_out.write(
                        f"{read.query_name}\t"
                        f"{start_pos}\t"
                        f"{end_pos}\t"
                        f"{read.query_name}:d4z4_{idx}\t"
                        f"{strand}\n"
                    )

                    # Intersect restriction sites with the alignment
                    xapi_intersections = [
                        f"{s}-{e}" for s, e, o in restriction_sites["XapI"]
                        if start_pos != "NA" and start_pos <= s < end_pos
                    ]
                    blni_intersections = [
                        f"{s}-{e}" for s, e, o in restriction_sites["BlnI"]
                        if start_pos != "NA" and start_pos <= s < end_pos
                    ]

                    table_out.write(
                        f"{read_id}\t"
                        f"{start_pos}-{end_pos}\t"
                        f"{strand}\t"
                        f"d4z4_{idx}\t"
                        f"{','.join(sorted(set(xapi_intersections))) if xapi_intersections else 'NA'}\t"
                        f"{','.join(sorted(set(blni_intersections))) if blni_intersections else 'NA'}\n"
                    )


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract alignment details and sequences from a BAM file.")
    parser.add_argument("--bam", required=True, help="Input BAM file.")
    parser.add_argument("--fasta", required=True, help="Input FASTA file with full reads.")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file with extracted alignments.")
    parser.add_argument("--output_table", required=True, help="Output table file with alignment details.")
    parser.add_argument("--xapi_bed", required=True, help="Output BED file for XapI restriction sites.")
    parser.add_argument("--blni_bed", required=True, help="Output BED file for BlnI restriction sites.")
    parser.add_argument("--main_tsv", required=True, help="Main output TSV file to update with sensitive repeat counts.")
    parser.add_argument("--repeats_bed", required=True, help="Output BED file for d4z4 repeats.")

    args = parser.parse_args()

    # Check for existing output files
    for file_path in [args.output_fasta, args.output_table, args.repeats_bed]:
        if os.path.exists(file_path):
            print(f"Output file {file_path} already exists. Deleting it to start fresh.")
            os.remove(file_path)

    # Collect all unique read IDs from the BAM file
    unique_read_ids = set()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            unique_read_ids.add(read.query_name)

    # Open the table file and write the header
    with open(args.output_table, "w") as table_out:
        table_out.write("ReadID\tPosition\tStrand\tD4Z4_Order\tXapI_Sites\tBlnI_Sites\n")

    # Process each unique read ID and generate outputs
    for read_id in unique_read_ids:
        #print(f"Processing read ID: {read_id}")
        process_bam(args.bam, read_id, args.fasta, args.output_table, args.output_fasta,args.xapi_bed, args.blni_bed, args.repeats_bed)

    print(f"FASTA file generated: {args.output_fasta}")
    print(f"Alignment details table generated: {args.output_table}")

    # Count sensitive repeats
    d4z4_table = read_tsv(args.output_table)
    xapi_counts, blni_counts = count_sensitive_repeats(d4z4_table)

    # Update the main TSV file
    update_main_tsv(args.main_tsv, xapi_counts, blni_counts)
    # refine_read_label(args.main_tsv)
    print(f"Main TSV file updated with sensitive repeat counts: {args.main_tsv}")

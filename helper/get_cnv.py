import pysam
import pandas as pd
import multiprocessing
import os
import argparse

parser = argparse.ArgumentParser(description="Calculate mapped estimated copies")
parser.add_argument("--outdir", required=True, help="output directory")
parser.add_argument("--prefix", required=True, help="sample prefix")

# Parse arguments
args = parser.parse_args()

# Paths
OUTDIR = os.path.abspath(args.outdir)
PREFIX = args.prefix
merged_bam_path = f"{OUTDIR}/merged.bam"
input_summary_path = f"{OUTDIR}/{PREFIX}_mapped_features_summary.tsv"
output_summary_path = f"{OUTDIR}/{PREFIX}_updated_features_summary.tsv"

# Haplotype regions
repeat_regions = {
    "4A": "HM190196.1_dux4_4A166:1-3199",
    "4B": "HM190164.1_dux4_4B168:1-3121",
    "10A": "HM190192.1_dux4_10A180T:1-3204",
    "10B": "HM190165.1_dux4_10B161T:1-3125"
}

# Function to determine haplotype
def get_haplotype(chrom):
    if chrom == "chr4":
        return "4A"
    elif chrom == "chr10":
        return "10A"
    return "NA"

# Function to estimate number of copies
def estimate_copies(read_id, chrom):
    haplotype = get_haplotype(chrom)
    if haplotype == "NA":
        return read_id, "NA", "NA", "NA"

    repeat_region = repeat_regions.get(haplotype, "NA")
    if repeat_region == "NA":
        return read_id, "NA", "NA", "NA"

    with pysam.AlignmentFile(merged_bam_path, "rb") as bam:
        # Filter reads for given ReadID
        reads_plus = [read for read in bam.fetch() if read.query_name == read_id and not read.is_reverse]
        reads_minus = [read for read in bam.fetch() if read.query_name == read_id and read.is_reverse]

    # Calculate coverage using pysam.depth
    number_copies_plus = sum(r.query_length for r in reads_plus) / len(reads_plus) if reads_plus else 0
    number_copies_minus = sum(r.query_length for r in reads_minus) / len(reads_minus) if reads_minus else 0
    number_copies = max(number_copies_plus, number_copies_minus)

    return read_id, number_copies, number_copies_plus, number_copies_minus

# Read input data
df = pd.read_csv(input_summary_path, sep="\t")
df_filtered = df[df.columns[:2]]  # Keep only ReadID and Chromosome columns

# Run in parallel
with multiprocessing.Pool(processes=4) as pool:
    results = pool.starmap(estimate_copies, df_filtered.values)

# Convert results to DataFrame
df_results = pd.DataFrame(results, columns=["ReadID", "MappedEstimatedCopies", "MappedEstimatedCopiesPlus", "MappedEstimatedCopiesMinus"])

# Merge with original summary file
df_final = pd.merge(df, df_results, on="ReadID", how="left")

# Save output
df_final.to_csv(input_summary_path, sep="\t", index=False)

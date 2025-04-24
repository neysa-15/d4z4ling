import pysam
import pandas as pd
import multiprocessing
import os
import subprocess
import argparse

# parser = argparse.ArgumentParser(description="Calculate mapped estimated copies")
# parser.add_argument("--outdir", required=True, help="Output directory")
# parser.add_argument("--prefix", required=True, help="Sample prefix")

# args = parser.parse_args()

# OUTDIR = os.path.abspath(args.outdir)
# PREFIX = args.prefix
# merged_bam_path = f"{OUTDIR}/{PREFIX}_aligned_haplotypes.bam"
# input_summary_path = f"{OUTDIR}/{PREFIX}_mapped_features_summary.tsv"
# output_summary_path = f"{OUTDIR}/{PREFIX}_updated_features_summary.tsv"

repeat_regions = {
    "4A": "HM190196.1_dux4_4A166:1-3199",
    "4B": "HM190164.1_dux4_4B168:1-3121",
    "10A": "HM190192.1_dux4_10A180T:1-3204",
    "10B": "HM190165.1_dux4_10B161T:1-3125"
}

def get_haplotype(chrom):
    if "chr4" in chrom:
        return "4A"
    elif "chr10" in chrom:
        return "10A"
    return "NA"

def compute_depth(region, bam_file, read_id, strand_flag):
    # temp_bam_path = f"{OUTDIR}/{read_id}_{os.getpid()}_temp.bam"
    # read_id_file = f"{OUTDIR}/{read_id}.txt"
    temp_bam_path = f"{read_id}_{os.getpid()}_temp.bam"
    read_id_file = f"{read_id}.txt"

    # Write Read ID to file
    with open(read_id_file, "w") as f:
        f.write(read_id + "\n")

    # Run samtools view with -N (file input)
    cmd = f"samtools view -N {read_id_file} -h -b {bam_file} {strand_flag} -o {temp_bam_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Remove temp read ID file
    os.remove(read_id_file)

    if result.returncode != 0:
        print(f"Error creating {temp_bam_path}: {result.stderr}")
        return 0

    # Index and calculate depth
    index_cmd = f"samtools index {temp_bam_path}"
    subprocess.run(index_cmd, shell=True, capture_output=True, text=True)

    depth_cmd = f"samtools depth -r {region} {temp_bam_path}"
    output = subprocess.run(depth_cmd, shell=True, capture_output=True, text=True).stdout.strip().split("\n")

    os.remove(temp_bam_path)
    os.remove(f"{temp_bam_path}.bai")

    depths = [int(line.split("\t")[2]) for line in output if len(line.split("\t")) == 3]
    return round(sum(depths) / len(depths), 2) if depths else 0


def estimate_copies(read_id, chrom, merged_bam_path):
    haplotype = get_haplotype(chrom)
    if haplotype == "NA":
        return read_id, "NA", "NA", "NA"

    repeat_region = repeat_regions.get(haplotype, "NA")
    if repeat_region == "NA":
        return read_id, "NA", "NA", "NA"

    number_copies_plus = compute_depth(repeat_region, merged_bam_path, read_id, "-F 16")
    number_copies_minus = compute_depth(repeat_region, merged_bam_path, read_id, "-f 16")
    number_copies = max(number_copies_plus, number_copies_minus)

    return read_id, number_copies, number_copies_plus, number_copies_minus

def get_cne(df, merged_bam_path):
    filtered_df = df[["ReadID", "GenomeCoords"]]

    print(filtered_df)

    with multiprocessing.Pool(processes=56) as pool:
        # results = pool.starmap(estimate_copies, filtered_df.values)
        results = pool.starmap(estimate_copies, [(read_id, chrom, merged_bam_path) for read_id, chrom in filtered_df.values])

    # Convert results to DataFrame
    df_results = pd.DataFrame(results, columns=["ReadID", "MappedEstimatedCopies", "MappedEstimatedCopiesPlus", "MappedEstimatedCopiesMinus"])

    # Merge with original summary file
    df_final = pd.merge(df, df_results, on="ReadID", how="left")

    return df_final

# Read input data
# df = pd.read_csv(input_summary_path, sep="\t")
# df_filtered = df[df.columns[:2]]

# # Run in parallel
# with multiprocessing.Pool(processes=56) as pool:
#     results = pool.starmap(estimate_copies, df_filtered.values)

# Convert results to DataFrame
# df_results = pd.DataFrame(results, columns=["ReadID", "MappedEstimatedCopies", "MappedEstimatedCopiesPlus", "MappedEstimatedCopiesMinus"])

# # Merge with original summary file
# df_final = pd.merge(df, df_results, on="ReadID", how="left")

# # Save output
# df_final.to_csv(input_summary_path, sep="\t", index=False)

# print("Processing complete. Output saved to", input_summary_path)
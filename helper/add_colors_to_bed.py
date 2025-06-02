import argparse
import pandas as pd

def color_d4z4_chr4_proximal(main_tsv, outfile):
    df = pd.read_csv(main_tsv, sep="\t")

    for _, row in df.iterrows():
        if not row["d4z4_chr4_proximal_mapped"]:
            continue
        
        chrom = row["ReadID"]
        start = row["d4z4_chr4_proximal_coords"].split("-")[0]
        end = row["d4z4_chr4_proximal_coords"].split("-")[1]
        score = row["d4z4_chr4_proximal_score"]
        strand = row["strand"]
        thick_start, thick_end = start, end
        color = "128,0,0"
        block_count, block_sizes, block_starts = 1, f"{int(end) - int(start)},", "0,"

        outfile.write("\t".join([chrom, start, end, "d4z4_chr4_proximal", str(score), strand,
            thick_start, thick_end, color, str(block_count),
            block_sizes, block_starts]) + "\n")

def process_beds(main_tsv, features_bed, repeats_bed, sslp_bed, output_bed):
    # Define color mapping
    even_color = "255,100,0"  # Darker orange
    odd_color = "255,200,100"  # Lighter orange
    feature_color_map = {
        "d4z4_chr4_proximal": "255,0,0",  # Red
        "p13-E11": "0,255,0",            # Green
        "pLAM": "0,0,255",                # Blue
        "4qA_marker": "122,48,160",       # Dark purple
        "4qB_marker": "0,100,0",       # Dark green
    }
    sslp_color = "255,0,255"  # Magenta for SSLP

    with open(main_tsv, "r") as tsv, open(features_bed, "r") as fbed, open(repeats_bed, "r") as rbed, open(output_bed, "w") as outfile, open(sslp_bed, "r") as sbed:

        # Process d4z4_chr4_proximal
        color_d4z4_chr4_proximal(main_tsv, outfile)
        
        # Process feature BED entries
        for line in fbed:
            if line.strip():  # Skip empty lines
                fields = line.strip().split("\t")
                chrom, start, end, feature_name = fields[:4]

                # if feature_name == "d4z4_chr4_proximal":
                #     continue

                score = fields[4] if len(fields) > 4 else "0"
                strand = fields[5] if len(fields) > 5 else "+"
                color = feature_color_map.get(feature_name, "0,0,0")  # Default to black
                thick_start, thick_end = start, end
                block_count, block_sizes, block_starts = 1, f"{int(end) - int(start)},", "0,"
                outfile.write("\t".join([chrom, start, end, feature_name, score, strand,
                                         thick_start, thick_end, color, str(block_count),
                                         block_sizes, block_starts]) + "\n")

        # Process repeat BED entries
        for line in rbed:
            if line.strip():  # Skip empty lines
                fields = line.strip().split("\t")
                chrom, start, end, name = fields[:4]
                # score = fields[4] if len(fields) > 4 else "0"
                # strand = fields[5] if len(fields) > 5 else "+"
                strand = fields[4] if len(fields) > 4 else "+"

                # Extract read ID and repeat index
                if ":" in name:
                    read_id, repeat_name = name.split(":")
                    if chrom != read_id:
                        print(f"Warning: Read ID {read_id} does not match chromosome {chrom}.")
                        print(fields)
                else:
                    repeat_name = name
                    read_id = None

                # Ensure read ID in name matches mapped read ID
                if read_id and read_id == chrom:
                    # Remove read ID and leave just the repeat name
                    color = even_color if strand == "+" else odd_color

                    # Generate BED12 fields
                    thick_start, thick_end = start, end
                    block_count, block_sizes, block_starts = 1, f"{int(end) - int(start)},", "0,"
                    outfile.write("\t".join([chrom, start, end, repeat_name, ".", strand,
                                             thick_start, thick_end, color, str(block_count),
                                             block_sizes, block_starts]) + "\n")

        # Get the length for each read
        read_length_dict = {}
        for line in tsv:
            if line.strip() and not line.startswith("ReadID"):
                fields = line.strip().split("\t")
                read_length = fields[5]

                # Skip entries where the length is "NA" or invalid
                if read_length == "NA":
                    print(f"Warning: Skipping read {fields[0]} due to missing read length.")
                    continue

                try:
                    read_length_dict[fields[0]] = int(read_length)
                except ValueError:
                    print(f"Error: Invalid read length for read {fields[0]}: {read_length}")
                    continue

        # Process SSLP BED entries
        for line in sbed:
            if line.strip():  # Skip empty lines
                fields = line.strip().split("\t")
                chrom, start, end, name = fields[:4]
                score = fields[4] if len(fields) > 4 else "0"
                strand = fields[5] if len(fields) > 5 else "+"

                # Assign magenta color for SSLP entries
                thick_start, thick_end = start, end
                block_count, block_sizes, block_starts = 1, f"{int(end) - int(start)},", "0,"
                outfile.write("\t".join([chrom, start, end, "SSLP", score, strand,
                                         thick_start, thick_end, sslp_color, str(block_count),
                                         block_sizes, block_starts]) + "\n")

    print(f"Unified BED file written to {output_bed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process features and repeats BED files and create a unified color-coded BED12 file.")
    parser.add_argument("--main_tsv", required=True, help="Path to the main TSV file.")
    parser.add_argument("--features_bed", required=True, help="Path to the input features BED file.")
    parser.add_argument("--repeats_bed", required=True, help="Path to the input repeats BED file.")
    parser.add_argument("--sslp_bed", required=True, help="Path to the input SSLP BED file.")
    parser.add_argument("--output_bed", required=True, help="Path to the output unified BED file.")

    args = parser.parse_args()
    print(args.main_tsv)

    process_beds(args.main_tsv, args.features_bed, args.repeats_bed, args.sslp_bed, args.output_bed)


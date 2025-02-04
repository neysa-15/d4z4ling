import argparse

def parse_blast_results(blast_file):
    """
    Parse BLAST results and return a dictionary of grouped results by ReadID and D4Z4_Order.
    """
    results = {}
    with open(blast_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            query = fields[0]
            subject = fields[1]
            identity = float(fields[2])
            bit_score = float(fields[11])

            read_id, d4z4_order = query.split(":")
            key = (read_id, d4z4_order)

            if key not in results:
                results[key] = []
            results[key].append((fields, bit_score, identity))
    return results

def parse_d4z4_units(d4z4_file):
    """
    Parse D4Z4 units into a list of tuples (ReadID, D4Z4_Order, other_fields).
    """
    d4z4_units = []
    with open(d4z4_file, "r") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            read_id = fields[0]
            d4z4_order = fields[3]
            d4z4_units.append((read_id, d4z4_order, fields))
    return header, d4z4_units

def find_best_blast_results(d4z4_file, blast_file, output_file):
    # Parse input files
    blast_results = parse_blast_results(blast_file)
    header, d4z4_units = parse_d4z4_units(d4z4_file)

    # Add BLAST columns to header
    header.extend(["Subject", "E_Value", "Bit_Score", "Identity"])

    # Open output file and write header
    with open(output_file, "w") as out:
        out.write("\t".join(header) + "\n")

        for read_id, d4z4_order, d4z4_fields in d4z4_units:
            key = (read_id, d4z4_order)
            if key in blast_results:
                # Find the best result by Bit_Score and use Identity as a tie-breaker
                best_result = max(
                    blast_results[key],
                    key=lambda x: (x[1], x[2])  # Sort by Bit_Score, then Identity
                )
                best_fields = best_result[0]
                out.write("\t".join(d4z4_fields + [best_fields[1], best_fields[10], best_fields[11], best_fields[2]]) + "\n")
            else:
                # No BLAST results, output NA for missing fields
                out.write("\t".join(d4z4_fields + ["NA", "NA", "NA", "NA"]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the best BLAST result for each D4Z4 unit.")
    parser.add_argument("--d4z4", required=True, help="Input TSV file with D4Z4 units.")
    parser.add_argument("--blast", required=True, help="Input BLAST results file (outfmt 6).")
    parser.add_argument("--output", required=True, help="Output file for the best BLAST results.")

    args = parser.parse_args()

    find_best_blast_results(args.d4z4, args.blast, args.output)


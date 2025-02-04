import re

# Input and output file paths
input_fasta = "dux4.gene_complete_genbank_20241119.fasta"
output_fasta = "dux4.gene_complete_genbank_20241119.reformatted.fasta"

# Open the input and output files
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    write_sequence = False  # Flag to track whether to write the sequence
    for line in infile:
        if line.startswith(">"):  # Identify header lines
            # Extract the accession and haplotype information
            match = re.search(r"^(>[^ ]+).*haplotype ([^ ]+)", line)
            if match:
                accession = match.group(1)
                haplotype = match.group(2)
                simplified_header = f"{accession}_dux4_{haplotype}"
                outfile.write(simplified_header + "\n")
                write_sequence = True  # Enable writing the following sequence
            else:
                write_sequence = False  # Disable writing if no match
        elif write_sequence:
            outfile.write(line)  # Write sequence lines only if header matched

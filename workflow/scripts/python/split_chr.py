#!/usr/bin/python3
from Bio import SeqIO

def split_and_overlap(input_fasta, output_prefix, num_files):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    records = records[0]

    # Calculate the length of each split
    #print(len(records.seq))
    split_length = len(records.seq) // num_files
    #print(split_length)
    initial_chr = records.id

    # Overlapping by 189 nt
    overlap = 189

    for i in range(num_files):
        start = i * split_length
        end = min((i + 1) * split_length + overlap, len(records))
        records.id = f'{initial_chr}_{i+1}'
        print(records.id)
        #print(start, end, end - start + 1)

        # Create a new file for each split
        output_file = f"{output_prefix}_{i + 1}.fa"
        SeqIO.write(records[start:end], output_file, "fasta")

if __name__ == "__main__":
    input_fasta = "data/references/genome_fa/drosophila_melanogaster/2R.fa"  # Replace with your input FASTA file
    output_prefix = "2R"  # Adjust the output file prefix
    num_files = 10

    split_and_overlap(input_fasta, output_prefix, num_files)


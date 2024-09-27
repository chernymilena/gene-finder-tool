import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence):
    """Find all Open Reading Frames in the given DNA sequence."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    seq_len = len(sequence)

    for frame in range(3):  # For each reading frame
        for i in range(frame, seq_len - 2, 3):
            if sequence[i:i + 3] == start_codon:  # Found a start codon
                for j in range(i + 3, seq_len - 2, 3):
                    if sequence[j:j + 3] in stop_codons:  # Found a stop codon
                        orfs.append(sequence[i:j + 3])  # Capture the ORF
                        break
    return orfs

def main():
    # Get the input files from command line arguments
    parser = argparse.ArgumentParser(description='Find Open Reading Frames in multiple genomes.')
    parser.add_argument('input_files', nargs='+', help='Paths to the input FASTA (.fna) files')
    args = parser.parse_args()

    # Iterate over all input files
    for input_file in args.input_files:
        if os.path.exists(input_file):
            with open(input_file, 'r') as file:
                sequence = str(SeqIO.read(file, "fasta").seq)

            orfs = find_orfs(sequence)
            print(f"Open Reading Frames in {input_file}:")
            for orf in orfs:
                print(orf)  # Print each ORF
            print()  # Print a blank line between files
        else:
            print(f"File {input_file} does not exist.")

if __name__ == "__main__":
    main()

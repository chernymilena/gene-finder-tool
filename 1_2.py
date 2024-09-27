import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def read_fasta_file(file_path):
    """Read a FASTA file and return the sequence."""
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq) 

def find_open_reading_frames(sequence):
    """Find all open reading frames (ORFs) in a DNA sequence."""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    # Search in three reading frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            if sequence[i:i+3] == start_codon:  # Check for start codon
                for j in range(i + 3, len(sequence) - 2, 3):
                    if sequence[j:j+3] in stop_codons:  # Check for stop codons
                        orfs.append(sequence[i:j+3])  # Capture the ORF
                        break

    return orfs

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Gene Finder Tool')
    parser.add_argument('input_file', help='Path to the input FASTA (.fna) file')
    args = parser.parse_args()

    # Read the FASTA file
    original_sequence = read_fasta_file(args.input_file)

    # Find ORFs in the original sequence
    forward_orfs = find_open_reading_frames(original_sequence)

    # Find ORFs in the reverse complement
    reverse_sequence = str(Seq(original_sequence).reverse_complement())
    reverse_orfs = find_open_reading_frames(reverse_sequence)

    # Combine results
    all_orfs = forward_orfs + reverse_orfs

    # Output found ORFs
    for orf in all_orfs:
        print(orf)

if __name__ == "__main__":
    main()

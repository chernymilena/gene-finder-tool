import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Function to reverse complement a DNA sequence
def reverse_complement(dna_sequence):
    return str(Seq(dna_sequence).reverse_complement())

# Function to translate a DNA sequence into a protein
def translate_dna_to_protein(dna_sequence):
    return str(Seq(dna_sequence).translate())

# Function to find all ORFs in a sequence
def find_orfs(dna_sequence, min_length):
    orfs = []
    seq_length = len(dna_sequence)
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    for frame in range(3):  # Reading frames 0, 1, and 2
        for i in range(frame, seq_length - 2):
            if dna_sequence[i:i+3] == start_codon:  # Found a start codon
                for j in range(i + 3, seq_length - 2, 3):
                    if dna_sequence[j:j+3] in stop_codons:  # Found a stop codon
                        orf = dna_sequence[i:j+3]
                        if len(orf) >= min_length * 3:  # Filtering by codon length
                            orfs.append(orf)
                        break
    return orfs

# Main function
def main():
    parser = argparse.ArgumentParser(description="Filter ORFs by length and discard short ORFs")
    parser.add_argument("input_file", help="Path to the input file with DNA sequence (output from previous question)")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum number of codons for an ORF to be considered a gene")
    args = parser.parse_args()

    # Read DNA sequence from the FASTA file
    with open(args.input_file, 'r') as file:
        record = SeqIO.read(file, "fasta")
        dna_sequence = str(record.seq)

    # Find ORFs in the original and reverse complement sequences
    original_orfs = find_orfs(dna_sequence, args.min_length)
    reverse_sequence = reverse_complement(dna_sequence)
    reverse_orfs = find_orfs(reverse_sequence, args.min_length)

    # Translate ORFs to protein strings and collect distinct proteins
    protein_strings = set()
    
    for orf in original_orfs:
        protein = translate_dna_to_protein(orf)
        if protein:
            protein_strings.add(protein)

    for orf in reverse_orfs:
        protein = translate_dna_to_protein(orf)
        if protein:
            protein_strings.add(protein)

    # Write filtered results to output file
    with open('output_5.txt', 'w') as output_file:
        for protein in protein_strings:
            output_file.write(f"{protein}\n")

if __name__ == "__main__":
    main()

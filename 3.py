import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(dna_sequence):
    """Calculate the reverse complement of a DNA sequence."""
    return str(Seq(dna_sequence).reverse_complement())

def translate_dna_to_protein(dna_sequence):
    """Translate a DNA sequence into a protein sequence."""
    protein_seq = Seq(dna_sequence).translate(to_stop=True)
    return str(protein_seq)

def find_orfs(dna_sequence):
    """Find all ORFs in a DNA sequence."""
    orfs = []
    seq_length = len(dna_sequence)
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    for frame in range(3):  # Reading frames 0, 1, and 2
        for i in range(frame, seq_length - 2, 3):
            codon = dna_sequence[i:i+3]
            if codon == start_codon:
                for j in range(i, seq_length - 2, 3):
                    stop_codon = dna_sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = dna_sequence[i:j+3]
                        orfs.append(orf)
                        break

    return orfs

def main():
    parser = argparse.ArgumentParser(description="Find protein strings from ORFs in both strands of a DNA sequence")
    parser.add_argument("input_file", help="Path to the input DNA sequence file (FASTA format)")
    args = parser.parse_args()

    # Read DNA sequence from the FASTA file
    with open(args.input_file, 'r') as file:
        record = SeqIO.read(file, "fasta")
        dna_sequence = str(record.seq)

    # Find ORFs in the original and reverse complement sequences
    original_orfs = find_orfs(dna_sequence)
    reverse_sequence = reverse_complement(dna_sequence)
    reverse_orfs = find_orfs(reverse_sequence)

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

    # Write results to output file
    with open('output_proteins.txt', 'w') as output_file:
        for protein in protein_strings:
            output_file.write(f"{protein}\n")

if __name__ == "__main__":
    main()

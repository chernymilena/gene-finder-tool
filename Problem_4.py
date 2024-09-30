import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def read_fasta_file(file_path):
    """Read a FASTA file and return the sequence."""
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)

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
    dna_sequence = read_fasta_file(args.input_file)

    # Find ORFs in the original sequence
    original_orfs = find_orfs(dna_sequence)

    # Find ORFs in the reverse complement sequence
    reverse_sequence = reverse_complement(dna_sequence)
    reverse_orfs = find_orfs(reverse_sequence)

    # Combine ORFs from both original and reverse complement sequences
    all_orfs = original_orfs + reverse_orfs

    # Translate ORFs to protein strings and collect distinct proteins
    protein_strings = set()
    
    for orf in all_orfs:
        protein = translate_dna_to_protein(orf)
        if protein:
            protein_strings.add(protein)
    for protein in protein_strings:
        print(f"{protein}")

if __name__ == "__main__":
    main()

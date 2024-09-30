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

def find_orfs(dna_sequence, min_length):
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
                        if len(orf) >= min_length * 3:  # Filtering by codon length
                            orfs.append(orf)
                        break

    return orfs

def contains_rbs(orf, rbs_sequences, upstream_length):
    """Check if the ORF contains any RBS within the upstream length."""
    start_codon_index = orf.find('ATG')  # Finding the start codon
    if start_codon_index != -1:
        # Check the sequence upstream of the start codon
        upstream_sequence = orf[max(0, start_codon_index - upstream_length):start_codon_index]
        return any(rbs in upstream_sequence for rbs in rbs_sequences)
    return False

def filter_orfs_with_rbs(orfs, rbs_sequences, upstream_length):
    """Filter ORFs that contain any of the RBS sequences."""
    filtered_orfs = []
    found_rbs_count = 0
    
    for orf in orfs:
        if contains_rbs(orf, rbs_sequences, upstream_length):
            filtered_orfs.append(orf)
            found_rbs_count += 1

    return filtered_orfs, found_rbs_count


def main():
    parser = argparse.ArgumentParser(description="Find protein strings from ORFs in both strands of a DNA sequence")
    parser.add_argument("input_file", help="Path to the input DNA sequence file (FASTA format)")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum number of codons for an ORF to be considered a gene")
    parser.add_argument('--rbs_sequences', nargs='+', help='Ribosome Binding Site sequences', default=['AGGAGG'])  # Common RBS sequences
    parser.add_argument('--upstream_length', type=int, help='Length to scan upstream for RBS', default=20)
    args = parser.parse_args()

    # Read DNA sequence from the FASTA file
    dna_sequence = read_fasta_file(args.input_file)

    # Find ORFs in the original sequence
    original_orfs = find_orfs(dna_sequence, args.min_length)

    # Find ORFs in the reverse complement sequence
    reverse_sequence = reverse_complement(dna_sequence)
    reverse_orfs = find_orfs(reverse_sequence, args.min_length)

    # Combine ORFs from both original and reverse complement sequences
    all_orfs = original_orfs + reverse_orfs

    # Translate ORFs to protein strings and collect distinct proteins
    protein_strings = set()
    
    filtered_orfs, found_rbs_count = filter_orfs_with_rbs(all_orfs, args.rbs_sequences, args.upstream_length)
    
    print(f"Found {found_rbs_count} ORFs with RBS.")
    for orf in filtered_orfs:
        protein = translate_dna_to_protein(orf)
        if protein:
            protein_strings.add(protein)
    for protein in protein_strings:
        print(f"{protein}")

if __name__ == "__main__":
    main()

from Bio import SeqIO
import argparse

def extract_scaffold_lengths(fasta_file):
    scaffold_lengths = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        scaffold_id = record.id
        scaffold_length = len(record.seq)
        scaffold_lengths[scaffold_id] = scaffold_length

    return scaffold_lengths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract scaffold lengths from a .fa file.")
    parser.add_argument('fasta_file', type=str, help="Path to the .fa file.")
    args = parser.parse_args()

    scaffold_lengths = extract_scaffold_lengths(args.fasta_file)

    # Save to a file without sorting
    with open("scaffold_lengths.txt", "w") as out_file:
        for scaffold, length in scaffold_lengths.items():
            out_file.write(f"{scaffold}\t{length}\n")

    # Save to another file with only lengths without sorting
    with open("only_lengths.txt", "w") as out_file:
        for _, length in scaffold_lengths.items():
            out_file.write(f"{length}\n")

    print(f"Extracted {len(scaffold_lengths)} scaffold lengths and saved to 'scaffold_lengths.txt'")
    print(f"Extracted {len(scaffold_lengths)} scaffold lengths without scaffold name (unsorted) and saved to 'only_lengths.txt'")

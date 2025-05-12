#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sequences from a fasta file using gene IDs. Made by Gabriel Dall'Alba."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input gene IDs, either as a .txt file or a single gene ID.",
    )
    parser.add_argument(
        "-f",
        "--fasta_database_file",
        required=True,
        help="Fasta file containing sequences to be extracted.",
    )
    return parser.parse_args()


def read_gene_ids(input_arg):
    if input_arg.endswith(".txt"):
        with open(input_arg, "r") as f:
            gene_ids = [line.strip() for line in f.readlines()]
    else:
        gene_ids = [input_arg]
    return gene_ids


def extract_sequences(gene_ids, fasta_database_file):
    sequences = SeqIO.index(fasta_database_file, "fasta")

    for gene_id in gene_ids:
        if gene_id in sequences:
            output_file = f"{gene_id}.fasta"
            SeqIO.write(sequences[gene_id], output_file, "fasta")
            print(f"Sequence for '{gene_id}' saved to '{output_file}'")
        else:
            print(f"Gene ID '{gene_id}' not found in the fasta file")


def main():
    args = parse_args()
    gene_ids = read_gene_ids(args.input)
    extract_sequences(gene_ids, args.fasta_database_file)


if __name__ == "__main__":
    main()

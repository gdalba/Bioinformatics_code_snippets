#!/usr/bin/env python3

import sys

def reformat_fasta_header(input_file, output_file):
    """
    Reformats headers in a FASTA file.
    """    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:        
        for line in infile:           
            if line.startswith('>'):
                print("Processing header:", line.strip())
                header_parts = line.strip().split(' ')
                new_header = header_parts[0]
                print("Writing new header:", new_header)
                outfile.write(new_header + '\n')            
            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python reformat_fasta_header.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    reformat_fasta_header(input_fasta, output_fasta)
import argparse
import numpy as np

def extract_gene_lengths_from_gff3(gff3_file, chromosome):
    data = []

    with open(gff3_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):  # Skip header lines
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == "gene" and parts[0] == chromosome: 
                    # Assuming columns are: seqid, source, type, start, end, score, strand, phase, attributes
                    start, end = int(parts[3]), int(parts[4])
                    length = end - start + 1
                    data.append((start, end, length))

    return np.array(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene lengths from a .gff3 file for a specified chromosome.")
    parser.add_argument('gff3_file', type=str, help="Path to the .gff3 file.")
    parser.add_argument('chromosome', type=str, help="Chromosome identifier to extract gene lengths for.")
    args = parser.parse_args()

    
    gene_data = extract_gene_lengths_from_gff3(args.gff3_file, args.chromosome)

    header = "gene start\tgene end\tgene length"
    np.savetxt(f"gene_lengths_{args.chromosome}.txt", gene_data, fmt='%d', delimiter='\t', comments='')
    print(f"Extracted {len(gene_data)} gene lengths for chromosome {args.chromosome} and saved to 'gene_lengths_{args.chromosome}.txt'")
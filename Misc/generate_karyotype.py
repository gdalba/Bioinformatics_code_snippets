import argparse
import os


def main():
    welcome_header = """
******************************************************
*                                                    *
*           Welcome to generate_karyotype.py         *
*                                                    *
******************************************************
"""

    script_name = "generate_karyotype.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-03-13"
    description = "This script reads a FASTA index file, extracts chromosome information, and saves the data to a txt file in Circos format."
    usage = "python generate_karyotype.py -i <input_file> -o <output_file> -c <colour>"

    print(welcome_header)

    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    parser = argparse.ArgumentParser(description='Convert GFF3 file to Circos format for gene annotations.')
    parser.add_argument('-i', '--input', required=True, help='Input GFF3 file path.')
    parser.add_argument('-o', '--output', required=True, help='Output csv file path.')
    parser.add_argument('-c', '--colours', type=str , required=False, help='Colour file path. Default: Black')

    args = parser.parse_args()

    generate_karyotype(args.input, args.output, args.colours)
    

def generate_karyotype(fai_file, output_file, colours):
    # Define a color for the chromosomes.
    if colours != None:
        color = colours
    else:       
        color = "black"
    
    with open(fai_file, 'r') as fai, open(output_file, 'w') as out:
        for line in fai:
            fields = line.strip().split('\t')
            chrom_id = fields[0]
            chrom_length = fields[1]
            # Format: chr - ID LABEL START END COLOR
            karyotype_line = f"chr - {chrom_id} {chrom_id} 0 {chrom_length} {color}\n"
            out.write(karyotype_line)


if __name__ == "__main__":
    main()

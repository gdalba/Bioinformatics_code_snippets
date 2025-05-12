import pandas as pd
import argparse

print('-----------------------------------')
print('GFF Cleaner V1 - 07/28/2023')
print('Made by Gabriel Dall\'Alba')
print('Please provide 2 arguments:')
print('-i or --input for input GFF file,')
print('-o or --output for output GFF file')
print('-----------------------------------')

parser = argparse.ArgumentParser(description='Clean GFF files.')
parser.add_argument('-i', '--input', required=True, help='Input GFF file')
parser.add_argument('-o', '--output', required=True, help='Output directory for cleaned GFF file')
args = parser.parse_args()

input_gff = args.input
output_gff = args.output

gff_data = pd.read_csv(input_gff, sep='\t', header=None)

# Shift the values in the 8th column to the 9th column
gff_data[8] = gff_data[7]
# Fill the 8th column with "."
gff_data[7] = "."

# convert the end column from float to int by first removing the period
gff_data[4] = gff_data[4].apply(lambda x: int(float(str(x).rstrip('.'))))

gff_data.to_csv(output_gff, sep='\t', header=False, index=False)

print("Processing of "+ input_gff +" done")

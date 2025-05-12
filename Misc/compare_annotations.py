import os
import pandas
import argparse
import pprint

def main():

    # Welcome header
    welcome_header = """
    ******************************************************
    *                                                    *
    *      Welcome to compare_annotations.py!            *
    *                                                    *
    ******************************************************
    """

    # Script information
    script_name = "compare_annotations.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-03-20"
    description = "This script will compare the two Mle annotations and output a table with matches."
    usage = "python compare_annotations.py -g1 <reference.gff3> -g2 <query.gff3> -o1 <matches.csv> -o2 <partials.csv> -o3 <unified.csv>"

    # Print welcome header
    print(welcome_header)

    # Print script information
    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Compare two Mle annotations and output a table with matches.")
    parser.add_argument("-g1", "--gff3_1", help="Path to the first GFF3 file. Treat it as reference.")
    parser.add_argument("-g2", "--gff3_2", help="Path to the second GFF3 file. Treat it as query.")
    parser.add_argument("-o1", "--output_1", help="Path to the output file with matches, will be .csv.")
    parser.add_argument("-o2", "--output_2", help="Path to the output file with partial-matches, will be .csv.")
    parser.add_argument("-o3", "--output_3", help="Path to the output file with ALL entries, will be .csv.")

    #parser.add_argument("-o", "--output", help="Path to the output file")
    args = parser.parse_args()

    reference_dict = read_gff3_to_dict(args.gff3_1)
    print(f"Reading GFF3 file: {args.gff3_1}")

    first_key, first_value = list(reference_dict.items())[0]
    first_item_key, first_item_value = list(first_value.items())[0]
    print(f"First item of the first entry in {first_key}: {first_item_key}: {first_item_value}")



    query_dict = read_gff3_to_dict(args.gff3_2)
    print(f"Reading GFF3 file: {args.gff3_2}")

    first_key, first_value = list(query_dict.items())[0]
    first_item_key, first_item_value = list(first_value.items())[0]
    print(f"First item of the first entry in {first_key}: {first_item_key}: {first_item_value}")

    print("Comparing genes...")
    matches, non_matches, partial_matches  = compare_genes_by_coordinates(reference_dict, query_dict)
    #print non-matches one line
    print(f"Matches: {len(matches)}, Partial-Matches: {len(partial_matches)}")


    output_path_1 = args.output_1
    output_path_2 = args.output_2
    output_unified = args.output_3

    print(f"Writing output of matches to {output_path_1}...")
    write_output(matches, output_path_1)
    print("Done!")
    
    print(f"Writing output of partial-matches to {output_path_2}...")
    write_output_of_partials(partial_matches, output_path_2)
    print("Done!")

    print(f"Writing output of all entries to {output_unified}...")
    write_output_to_single_file(matches, partial_matches, output_unified)
    print("Done!")


# function to compare the two dictionaries for a given chromosome, we will then iterate each gene of query to all genes of reference as long as
# we are comparing the SAME CHROMOSOME of query gene, for instance, gene 1 of chromosome 1 in query will be compared to all genes of chromosome 1 in reference
# we consider coordinates, we do not care about the gene name, we just want to know if the gene is the same in both annotations


'''
reference --------------- -------------- ------------- -> g1 g2 g3
query        ------       --------------  ---    ---- -> q1 q2 q3 q4
matches: q2-g2
partial-matches: q1-g1, q3-g3, q4-g3



def compare_genes_by_coordinates(reference_dict, query_dict):
    matches = []
    non_matches = []
    for chromosome in query_dict:
        if chromosome in reference_dict:
            print(f"Comparing genes on {chromosome}")
            for query_gene_id, query_gene_info in query_dict[chromosome].items():
                query_start = int(query_gene_info['start'])
                query_end = int(query_gene_info['end'])
                for ref_gene_id, ref_gene_info in reference_dict[chromosome].items():
                    ref_start = int(ref_gene_info['start'])
                    ref_end = int(ref_gene_info['end'])
                    if query_start == ref_start and query_end == ref_end:
                        matches.append((chromosome, query_gene_id, ref_gene_id, query_start, query_end, ref_start, ref_end))
                    else:
                        non_matches.append((chromosome, query_gene_id, ref_gene_id, query_start, query_end, ref_start, ref_end))
    return matches, non_matches
'''

def compare_genes_by_coordinates(reference_dict, query_dict):
    matches = []
    partial_matches = []
    non_matches = []

    for chromosome in query_dict:
        if chromosome in reference_dict:
            print(f"Comparing genes on {chromosome}")
            for query_gene_id, query_gene_info in query_dict[chromosome].items():
                query_start = int(query_gene_info['start'])
                query_end = int(query_gene_info['end'])
                query_strand = query_gene_info['strand']

                found_match = False  # track match or partial match

                for ref_gene_id, ref_gene_info in reference_dict[chromosome].items():
                    ref_start = int(ref_gene_info['start'])
                    ref_end = int(ref_gene_info['end'])
                    ref_strand = ref_gene_info['strand']

                    if query_start == ref_start and query_end == ref_end:
                        matches.append((chromosome, query_gene_id, ref_gene_id, query_start, query_end, ref_start, ref_end, ref_strand, query_strand))
                        found_match = True
                        break  # no need to keep checking once exact match is found
                    elif query_start >= ref_start and query_end <= ref_end:
                        partial_matches.append((chromosome, query_gene_id, ref_gene_id, query_start, query_end, ref_start, ref_end, ref_strand, query_strand))
                        found_match = True  # other better matches may exist

                if not found_match:
                    non_matches.append((chromosome, query_gene_id, query_start, query_end, query_strand))

    return matches, non_matches, partial_matches


#function to compare the non-matches from function compare_genes_by_coordinates, we want to check if a query gene start and end are 
#at least contained within a reference gene start and end, if so, we will consider it a partial match
#the idea is the following, a given gene of reference might be, say, from reference 100 to 800, and a query gene might be from 200 to 600
#we will consider it a partial match, as the query gene is contained within the reference gene, there will be situations where
#multiple query genes are contained within a single reference gene, we want to know these cases
#HOWEVER, we do not want to consider a query gene that is contained within a reference gene if it is already a match, so we will check for that

def find_partial_matches(matches, non_matches):
    partial_matches = []
    for non_match in non_matches:
        query_chromosome = non_match[0]
        query_start = non_match[3]
        query_end = non_match[4]
        ref_chromosome = non_match[2]
        ref_start = non_match[5]
        ref_end = non_match[6]
        if query_chromosome != ref_chromosome:
            continue
        print(f"Comparing genes on {query_chromosome}")
        if any(match[3] == query_start and match[4] == query_end for match in matches):
            continue
        if query_start >= ref_start and query_end <= ref_end:
            partial_matches.append(non_match)
    return partial_matches


#function to write output with two columns in csv: reference gene id and query gene id, for now, also write chromosome
#we also want to write the start and stop of query and reference genes, so we can check if the coordinates are the same
#we also want to write the strand of each query and ref gene

def write_output(matches, output_file):
    with open(output_file, 'w') as file:
        file.write("chromosome,query_gene_id,reference_gene_id,query_start, query_end, reference_start, reference_end, ref_strand, query_strand, same_gene\n")
        for match in matches:
            file.write(f"{match[0]},{match[1]},{match[2]},{match[3]},{match[4]},{match[5]},{match[6]}, {match[7]}, {match[8]}, 'yes'\n")

def write_output_of_partials(partial_matches, output_file): 
    with open(output_file, 'w') as file:
        file.write("chromosome,query_gene_id,reference_gene_id,query_start, query_end, reference_start, reference_end, ref_strand, query_strand, same_gene\n")
        for match in partial_matches:
            file.write(f"{match[0]},{match[1]},{match[2]},{match[3]},{match[4]},{match[5]},{match[6]}, {match[7]}, {match[8]}, 'no'\n")


#function to write output to a single file

def write_output_to_single_file(matches, partial_matches, output_file):
    with open(output_file, 'w') as file:
        file.write("chromosome,query_gene_id,reference_gene_id,query_start, query_end, reference_start, reference_end, ref_strand, query_strand, same_gene\n")
        for match in matches:
            file.write(f"{match[0]},{match[1]},{match[2]},{match[3]},{match[4]},{match[5]},{match[6]}, {match[7]}, {match[8]}, 'yes'\n")
        for match in partial_matches:
            file.write(f"{match[0]},{match[1]},{match[2]},{match[3]},{match[4]},{match[5]},{match[6]}, {match[7]}, {match[8]}, 'no'\n")

def parse_gff3(gff3_file):
    chromosomes = {}
    with open(gff3_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chromosome = line[0]
            if chromosome not in chromosomes:
                chromosomes[chromosome] = {}
            if line[2] == 'gene':
                gene_id = line[8].split(';')[0].split('=')[1]
                chromosomes[chromosome][gene_id] = []
            elif line[2] in ['exon', 'CDS']:
                chromosomes[chromosome][gene_id].append((int(line[3]), int(line[4])))
                print(f"Parsed chromosomes from {gff3_file}: {chromosomes}") 

    return chromosomes

def compare_genes(chromosomes1, chromosomes2):
    matches = []
    non_matches = []
    for chromosome, genes1 in chromosomes1.items():
        if chromosome in chromosomes2:
            genes2 = chromosomes2[chromosome]
            for gene_id, gene_info in genes1.items():
                if gene_id in genes2:
                    if sorted(genes2[gene_id]) == sorted(gene_info):
                        matches.append((chromosome, gene_id))
                    else:
                        non_matches.append((chromosome, gene_id))
    return matches, non_matches

def read_gff3_to_dict(gff3_file_path):
    """
    Reads a GFF3 file and processes it into a nested dictionary. 
    The outer dictionary has chromosomes as keys, 
    and each value is another dictionary with gene IDs as keys.
    The inner dictionaries store information about each gene, 
    including its start position, end position, strand, and attributes.
    
    """
    chromosome_dict = {}
    with open(gff3_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip headers and empty lines
            
            parts = line.strip().split('\t')
            if parts[2].lower() == 'gene':  # Focus on gene features
                chromosome = parts[0]
                gene_id = parts[8].split('ID=')[1].split(';')[0]
                gene_info = {
                    'start': parts[3],
                    'end': parts[4],
                    'strand': parts[6],
                    'attributes': parts[8]
                }
                #print the gene that was parsed
                
                
                if chromosome not in chromosome_dict:
                    chromosome_dict[chromosome] = {}
                chromosome_dict[chromosome][gene_id] = gene_info
                    
    return chromosome_dict

if __name__ == "__main__":
    main()



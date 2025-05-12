#!/usr/bin/env python3

"""
Author: Gabriel Dall'Alba
Creation Date: 31-07-2023
Last Modification: 31-07-2023
Required Conda Environment: gffutils
Libraries included in the Conda environment: gffutils
Usage: python script.py -i input.gff3 -o output.csv
Or:
Usage: python script.py -d input.db -o output.csv
"""

import argparse
import csv
import gffutils
import numpy as np

def print_header():
    header = """
    =====================================
          M. leidyi annotation stats
                 GFF3 PROCESSOR         
    Author: Gabriel Dall'Alba
    Date: 31-07-2023
    Environment: gffutils
    Libs: gffutils
    =====================================
    """
    print(header)

def read_gff3(filename):
    print("STATUS: Reading .gff3 file. Generating .db file.")
    db = gffutils.create_db(filename, dbfn='temp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    print("STATUS: .db file created.")
    for gene in db.features_of_type('gene'):
        yield gene

def read_db(dbname):
    print("STATUS: Reading .db file.")
    db = gffutils.FeatureDB(dbname)
    for gene in db.features_of_type('gene'):
        yield gene
    print(".")

def check_nested(genes):
    host_genes = set()
    nested_genes = set()
    partially_nested_genes = set()
    for gene1 in genes:
        for gene2 in genes:
            if gene1.seqid == gene2.seqid and gene1.id != gene2.id:
                if (gene2.start >= gene1.start and gene2.end <= gene1.end):
                    nested_genes.add(gene2.id)
                    host_genes.add(gene1.id)
                elif gene2.start > gene1.end or gene2.end < gene1.start:
                    continue
                elif (gene2.start < gene1.start and gene2.end < gene1.end) or (gene2.end > gene1.end and gene2.start > gene1.start):
                    partially_nested_genes.add(gene2.id)
                    partially_nested_genes.add(gene1.id)
    return host_genes, nested_genes, partially_nested_genes


def get_exon_info(gene, db):
    exons = list(db.children(gene, featuretype='exon'))
    num_exons = len(exons)
    avg_exon_length = np.mean([exon.end - exon.start + 1 for exon in exons]) if exons else 0
    return num_exons, avg_exon_length

def get_intron_info(gene, db):
    exons = list(db.children(gene, featuretype='exon'))
    total_exon_length = sum([exon.end - exon.start + 1 for exon in exons])
    total_gene_length = gene.end - gene.start + 1
    total_intron_length = total_gene_length - total_exon_length
    num_introns = len(exons) - 1 if len(exons) > 1 else 0
    return num_introns, total_intron_length

def get_intron_boundaries(gene, db):
    exons = sorted(list(db.children(gene, featuretype='exon')), key=lambda x: x.start)
    intron_boundaries = [(exons[i].end + 1, exons[i + 1].start - 1) for i in range(len(exons) - 1)]
    return intron_boundaries

def check_transcripts_in_introns(genes, db):
    genes_with_transcripts_in_introns = set()
    for gene1 in genes:
        gene1_introns = get_intron_boundaries(gene1, db)
        for gene2 in genes:
            if gene1.id != gene2.id:
                for transcript in db.children(gene2, featuretype='transcript'):
                    for intron in gene1_introns:
                        if intron[0] <= transcript.start and transcript.end <= intron[1]:
                            genes_with_transcripts_in_introns.add(gene1.id)
                            break
    return genes_with_transcripts_in_introns



def analyze_genes(genes, db):
    exon_lengths = []
    exon_counts = []
    intron_lengths=[]
    intron_counts=[]
    single_exon_genes = 0
    genes_with_transcripts_in_introns = check_transcripts_in_introns(genes, db)

    for gene in genes:
        num_exons, avg_exon_length = get_exon_info(gene, db)
        num_introns, total_intron_length = get_intron_info(gene, db)
        exon_lengths.append(avg_exon_length)
        exon_counts.append(num_exons)
        intron_lengths.append(total_intron_length)
        intron_counts.append(num_introns)
        if num_exons == 1:
            single_exon_genes += 1

    avg_exon_length = np.mean(exon_lengths)
    avg_exon_count = np.mean(exon_counts)
    avg_intron_length = np.mean(intron_lengths)
    avg_intron_count = np.mean(intron_counts)

    host_genes_completely, completely_nested_genes, partially_nested_genes = check_nested(genes)

    return avg_exon_length, avg_exon_count, avg_intron_count, avg_intron_length,single_exon_genes, host_genes_completely, completely_nested_genes, partially_nested_genes, genes_with_transcripts_in_introns

def write_to_csv(genes, db, host_genes_completely, completely_nested_genes, partially_nested_genes, genes_with_transcripts_in_introns, output_file):
    print("STATUS: Writing to output.")
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Gene ID', 'Scaffold', 'Start', 'End', 'Is Host', 'Is Nested', 'Is Overlapping', 'Has Transcript in Intron', 'Strand', 'Number of Exons', 'Average Exon Length', 'Number of Introns', 'Average Intron Length'])
        for gene in genes:
            num_exons, avg_exon_length = get_exon_info(gene, db)
            num_introns, avg_intron_length = get_intron_info(gene,db)
            is_host = 'YES' if gene.id in host_genes_completely else 'NO'
            is_nested = 'YES' if gene.id in completely_nested_genes else 'NO'
            is_overlapping = 'YES' if gene.id in partially_nested_genes else 'NO'
            has_transcripts_in_intron = 'YES' if gene.id in genes_with_transcripts_in_introns else 'NO'        
            strand = gene.strand
            writer.writerow([gene.id, gene.seqid, gene.start, gene.end, is_host, is_nested, is_overlapping, has_transcripts_in_intron, strand, num_exons, avg_exon_length, num_introns, avg_intron_length])
    print("STATUS: Output file processed.")

def main():
    print_header()
    parser = argparse.ArgumentParser(description='Process GFF3 file.')
    parser.add_argument('-i', '--input', type=str, help='Input GFF3 file.')
    parser.add_argument('-d', '--db', type=str, help='Input database file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output CSV file.')

    args = parser.parse_args()

    if args.db:
        db = gffutils.FeatureDB(args.db)
        genes = list(db.all_features(featuretype='gene'))
    elif args.input:
        db = gffutils.create_db(args.input, dbfn='temp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        genes = list(db.all_features(featuretype='gene'))
    else:
        raise ValueError('Either --input or --db must be given.')
    
    print("STATUS: Computing statistics.")
    total_genes = len(genes)
    avg_exon_length, avg_exon_count, avg_intron_count, avg_intron_length, single_exon_genes, host_genes_completely, completely_nested_genes, partially_nested_genes, genes_with_transcripts_in_introns = analyze_genes(genes, db)
    #host_genes = host_genes_completely.union(host_genes_partially)

    print(f"STATUS: Number of genes analyzed: {len(genes)}")
    print(f"STATUS: Average exon length: {avg_exon_length}")
    print(f"STATUS: Average exon count per gene: {avg_exon_count}")
    print(f"STATUS: Average intron length: {avg_intron_length}")
    print(f"STATUS: Average intron count per gene: {avg_intron_count}")
    print(f"STATUS: Number of genes with transcripts in introns: {len(genes_with_transcripts_in_introns)}")
    print(f"STATUS: Number of single-exon genes: {single_exon_genes} ({(single_exon_genes/total_genes)*100:.2f}%)")
    #print(f"STATUS: Total number of genes involved in nesting (old function): {len(nested_genes)} ({(len(nested_genes)/total_genes)*100:.2f}%)")
    print(f"STATUS: Total number of genes involved in nesting (new function): {len(completely_nested_genes)+len(partially_nested_genes)} ({((len(completely_nested_genes)+len(partially_nested_genes))/total_genes)*100:.2f}%)")
    print(f"STATUS: Number of completely nested genes: {len(completely_nested_genes)} ({(len(completely_nested_genes)/total_genes)*100:.2f}%)")
    print(f"STATUS: Number of partially nested (overlapping) genes: {len(partially_nested_genes)} ({(len(partially_nested_genes)/total_genes)*100:.2f}%)")

    write_to_csv(genes, db, host_genes_completely, completely_nested_genes, partially_nested_genes, args.output)


if __name__ == "__main__":
    main()

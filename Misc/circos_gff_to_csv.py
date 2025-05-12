import gffutils
import pandas as pd
import argparse
import os


def main():
    welcome_header = """
******************************************************
*                                                    *
*           Welcome to circos_gff_to_csv.py          *
*                                                    *
******************************************************
"""

    script_name = "circos_gff_to_csv.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-03-12"
    description = "This script reads a GFF3/GTF file, extracts gene annotations, and saves the data to a CSV file in Circos format."
    usage = "python circos_gff_to_csv.py -i <input_file> -o <output_file>"

    print(welcome_header)

    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    parser = argparse.ArgumentParser(description='Convert GFF3 file to Circos format for gene annotations.')
    parser.add_argument('-i', '--input', required=True, help='Input GFF3 file path.')
    parser.add_argument('-o', '--output', required=True, help='Output csv file path.')
    parser.add_argument('-d', '--db', required=True, help='Database file path.') #can this be optional?
    parser.add_argument('-f', '--force', required=False, help='Force database creation.Boolean value (0 or 1). Default = 0.')
    parser.add_argument('-n', '--db_name', required=False, help='Database name. Default = "Annotation"')
    args = parser.parse_args()

    if args.force == None:
        args.force = "0"

    if args.force == "0":
        args.force = False
    elif args.force == "1":
        args.force = True

    if args.db_name == "None":
        args.db_name = "Annotation"


    #handle if no db is provided
    if args.db == "None":
        args.db = f"{args.input}.db"

    #check if db exists
    if os.path.exists(args.db) and args.force == False:
        db = gffutils.FeatureDB(args.db)
    else:
        db = gffutils.create_db(args.input, dbfn=args.db, force=args.force, keep_order=True, merge_strategy="merge", sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=True) 
    
    genes_df = extract_genes(db)
    genes_df.to_csv(args.output, sep='\t', index=False, header=False)

    #check if file has correct number of genes per
    print(f"Number of genes: {len(genes_df)}")
    
# Fetch gene annotations and convert to Circos format -- should be .gtf format
def extract_genes(db):
    genes = []
    for gene in db.features_of_type('gene'):
        # Extract ID, chromosome (seqid), start, end, strand
        gene_id = gene.id
        chromosome = gene.seqid
        start = gene.start
        end = gene.end
        strand = gene.strand
        genes.append([chromosome, start, end, gene_id, strand])
    return pd.DataFrame(genes, columns=['chromosome', 'start', 'end', 'gene_id', 'strand'])

if __name__ == "__main__":
    main()
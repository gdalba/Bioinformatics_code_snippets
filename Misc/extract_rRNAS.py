import sys
import csv

def gff3_to_tsv(gff3_file_path, tsv_file_path):
    """
    Convert GFF3 file to TSV file with specified columns.
    
    Parameters:
    - gff3_file_path (str): Path to the input GFF3 file.
    - tsv_file_path (str): Path to the output TSV file.
    """
    
    with open(gff3_file_path, 'r') as gff3_file:
        
        with open(tsv_file_path, 'w', newline='') as tsv_file:
            
            tsv_writer = csv.writer(tsv_file, delimiter='\t')
            tsv_writer.writerow(['chromosome', 'class_of_rRNA', 'start', 'end', 'sense/antisense'])
            
            for line in gff3_file:
                
                if line.startswith("##"):
                    continue
                
                fields = line.strip().split("\t")
                
                chromosome = fields[0]
                class_of_rRNA = fields[8].split(";")[0].split("=")[1]
                start = fields[3]
                end = fields[4]
                sense_antisense = "+" if fields[6] == "+" else "-"
                
                tsv_writer.writerow([chromosome, class_of_rRNA, start, end, sense_antisense])

if __name__ == "__main__":
    input_gff3 = sys.argv[1]  
    output_tsv = sys.argv[2]  
    gff3_to_tsv(input_gff3, output_tsv)

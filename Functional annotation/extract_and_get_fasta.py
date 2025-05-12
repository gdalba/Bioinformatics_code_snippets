import pandas as pd
import argparse
import sys

def main():

    welcome_header = """
    ******************************************************
    *                                                    *
    *            Welcome to extract_and_get_fasta.py     *
    *                                                    *
    ******************************************************
    """

    script_name = "extract_and_get_fasta.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-02-28"
    description = "This script reads an Excel file and extracts IDs based on empty cells in a given column. It then recovers FASTA sequences based on the extracted IDs and saves them to a new file."
    usage = "python extract_and_get_fasta.py -i <file> -c <column> -o <output_file> -f <fasta_file>"

    print(welcome_header)

    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    parser = argparse.ArgumentParser(description="Read an Excel file and extract IDs based on empty cells in a given column")
    parser.add_argument("-i", "--input", help="Path to the Excel file")
    parser.add_argument("-c", "--col", help="Column to check for empty cells (1-based index)", type=int)
    parser.add_argument("-o", "--output", help="Path to the output file")
    parser.add_argument("-f", "--fasta", help="Path to the fasta file")

    args = parser.parse_args() 
    ids = extract_ids_from_excel(args.input, args.col)

    recover_and_save_fasta_sequences(ids, args.fasta, args.output)


def extract_ids_from_excel(excel_path, column_to_check):
    # Load the Excel file
    df = pd.read_excel(excel_path)
    # Check for empty cells in the specified column
    empty_cells = df[df.columns[column_to_check - 1]].isnull()
    #print(empty_cells)
    # Extract IDs (assuming the IDs are in the second column)
    ids = df.loc[empty_cells, df.columns[1]].tolist()  # Adjusted to 0-based index for ID column
    print(len(ids))
    return ids


def recover_and_save_fasta_sequences(ids, fasta_path, output_path):
    # Initialize a flag to track when a matching ID is found
    write_sequence = False
    # Open the output file
    with open(output_path, 'w') as output_file:
        # Read the input FASTA file
        # with open(fasta_path, 'r') as fasta_file:
#            for line in fasta_file:
                        
            # check if id in the list of ids match any part of a given id in fasta file
            # for id in ids:
            #     print("Looking for: ", id)

            #     for line in fasta_file:
            #         if id in line:
            #             print("Found id: ", id)
            #             print("Found in .fasta file: ", line)
            #             write_sequence = True
            #             print("State of write_sequence: ", write_sequence)
            #             output_file.write(line)
            #             break
            #print specific line in fasta_file by the number of line
            # print(fasta_file.readlines()[130308 - 1])
            
            # ids1 = [ids[0], ids[1]]
            for id in ids:
                #print("Looking for: ", id)
                
                with open(fasta_path, 'r') as fasta_file:
                    for line in fasta_file:
                        # print("trying in all lines for: ", id)
                        if id in line:
                            #print("Found id: ", id)
                            #print("Found in .fasta file: ", line)
                            write_sequence = True
                            #print("State of write_sequence: ", write_sequence)
                            output_file.write(line)
                            # break
                        elif line.startswith('>'):
                            write_sequence = False
                        elif write_sequence:
                            output_file.write(line)
                # print"not found in all lines for: ", id)


                    # elif line.startswith('>'):
                    #     write_sequence = False
                    # elif write_sequence:
                    #     output_file.write(line)
                    
'''                    
                    if fasta_id in ids:
                        write_sequence = True
                        
                        # Write the header line to the output file
                        output_file.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    # Write the sequence line to the output file
                    output_file.write(line)
'''


if __name__ == "__main__":
    main()

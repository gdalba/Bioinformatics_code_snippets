import openpyxl
import requests
import argparse
import psutil
import os
import time
import threading
from bioservices import KEGG


def main():

    welcome_header = """
    ******************************************************
    *                                                    *
    *                Welcome to kegg_IDs.py              *
    *                                                    *
    ******************************************************
    """

    script_name = "kegg_IDs.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-02-27"
    description = "This script imports kegg IDs from an Excel file and retrieves their names using the API from kegg."
    usage = "python kegg_IDs.py /path/to/excel_file.xlsx"

    print(welcome_header)

    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    parser = argparse.ArgumentParser(description="Import kegg terms from an Excel file")
    parser.add_argument("file", help="Path to the Excel file")
    args = parser.parse_args()

    workbook = openpyxl.load_workbook(args.file)
    sheet = workbook.active
    
    kegg_terms_dict = {}
    

    # Store the locations of each kegg term
    kegg_term_locations = {}

    for row in sheet.iter_rows():
        for cell in row:
            if cell.value and cell.value.startswith("ko:"):
                kegg_term = cell.value.split("ko:")[1]  # Split the term and get the value after "ko:" otherwise it will return the whole string
                if kegg_term not in kegg_term_locations:
                    kegg_term_locations[kegg_term] = []
                kegg_term_locations[kegg_term].append(cell)

    # Replace each kegg term
    for kegg_term, locations in kegg_term_locations.items():
        # Get the name of the kegg term
        if kegg_term in kegg_terms_dict:
            kegg_term_name = kegg_terms_dict[kegg_term]
        else:
            #search for kegg through bioservices KEGG API
            k_db = KEGG() 
            #use try to keep the script running even if a potential timeout occurs
            try:
                result = k_db.get(kegg_term)
            except requests.Timeout:
                print(f"Timeout occurred while retrieving kegg term {kegg_term}")
                continue


            #break the string into a list and retrieve entry symbol name from the list
            #if result is not int object            
            if not isinstance(result, int):
                kegg_term_name = result
                kegg_term_name = kegg_term_name.split("\n")            
                kegg_line_one = kegg_term_name[0] # entry
                kegg_line_one = kegg_line_one.split() 
                kegg_line_two = kegg_term_name[1] # symbol
                kegg_line_two = kegg_line_two.split()
                kegg_line_three = kegg_term_name[2] #name
                kegg_line_three = kegg_line_three.split()
                
                #merge elements 1 to the end of the list
                kegg_line_three = " ".join(kegg_line_three[1:])
                
                # Output format NAME (SYMBOL) (ENTRY)
                #line three is name
                #line two is symbol
                #line one is entry

                #print(kegg_line_one[1], kegg_line_two[1], kegg_line_three[1])

                # Merge the elements to be format name (symbol) (entry)
                kegg_term_name = f"{kegg_line_three} ({kegg_line_two[1]}) ({kegg_line_one[1]})"

                print(f"Retrieved kegg term {kegg_term} with name {kegg_term_name}")
                kegg_terms_dict[kegg_term] = kegg_term_name
            else:
                kegg_term_name = kegg_term

        # Replace the kegg term in each location
        for cell in locations:
            cell.value = f"{kegg_term_name}"

        print(f"Substituted {len(locations)} occurrences of kegg term {kegg_term} with {kegg_term_name}")

    # Save the changes and close the Excel file
    workbook.save(args.file)
    workbook.close()

if __name__ == "__main__":
    main()
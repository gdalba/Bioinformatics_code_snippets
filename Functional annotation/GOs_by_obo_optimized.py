import openpyxl
import requests
import argparse
import psutil
import os
import time
import threading

def main():

    welcome_header = """
    ******************************************************
    *                                                    *
    *      Welcome to GOs_by_obo_optimized.py            *
    *                                                    *
    ******************************************************
    """

    script_name = "GOs_by_obo_optimized.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-02-26"
    description = "This script imports GO terms from an Excel file and retrieves their names using the .obo database."
    usage = "python GOs_by_obo_optimized.py <file> <obo_file>"

    print(welcome_header)

    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    parser = argparse.ArgumentParser(description="Import GO terms from an Excel file")
    parser.add_argument("file", help="Path to the Excel file")
    parser.add_argument("obo_file", help="Path to the .obo file")
    args = parser.parse_args()

    workbook = openpyxl.load_workbook(args.file)
    sheet = workbook.active
    
    go_terms_dict = {}
    obo_go_terms_dict = parse_obo_file(args.obo_file)

    go_term_locations = {}

    for row in sheet.iter_rows():
        for cell in row:
            if cell.value and cell.value.startswith("GO"):
                go_term = cell.value
                if go_term not in go_term_locations:
                    go_term_locations[go_term] = []
                go_term_locations[go_term].append(cell)

    # Replace each GO term
    for go_term, locations in go_term_locations.items():
        if go_term in go_terms_dict:
            go_term_name = go_terms_dict[go_term]
        else:
            go_term_name = obo_go_terms_dict.get(go_term)
            go_terms_dict[go_term] = go_term_name

        for cell in locations:
            cell.value = f"{go_term_name} ({go_term})"

        print(f"Substituted {len(locations)} occurrences of GO term {go_term} with {go_term_name} ({go_term})")

    workbook.save(args.file)
    workbook.close()

def get_go_term_name(go_term, go_terms_dict):
    go_term_name = go_terms_dict.get(go_term)
    if go_term_name:
        return go_term_name
    else:
        print(f"GO term {go_term} not found in local .obo file.")
        return None


def parse_obo_file(obo_file_path):
    go_terms = {}
    with open(obo_file_path, 'r') as file:
        current_id = None
        current_name = None
        for line in file:
            line = line.strip()
            if line.startswith('id: GO:'):
                current_id = line.split("id: ")[1]
            elif line.startswith('name: '):
                current_name = line.split("name: ")[1]
            elif line == '[Term]' or line == '':
                if current_id and current_name:
                    go_terms[current_id] = current_name
                current_id, current_name = None, None
    return go_terms

if __name__ == "__main__":
    main()
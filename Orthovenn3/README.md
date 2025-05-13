# Orthovenn3 Orthogroup Filter

This project is designed to filter entries from an Excel (.xlsx) file based on identifiers provided in a text (.txt) file. The application reads the identifiers, searches for them in the Excel file, and outputs a new Excel file containing only the matching entries. 

## Project Structure

```
excel-filter-project
├── src
│   ├── main.py                # Entry point of the application
│   └── utils
│       ├── __init__.py        # Initializes the utils package
│       ├── excel_processor.py  # Functions for processing Excel files
│       └── file_handler.py     # Functions for handling file operations
├── environment.yml            # Conda environment configuration
├── setup.py                   # Packaging information
├── requirements.txt           # Python package dependencies
└── README.md                  # Project documentation
```

## Installation

To set up the project, you need to create a conda environment. You can do this by running the following command in your terminal:

```bash
conda env create -f environment.yml
```

This will install all the necessary dependencies specified in the `environment.yml` file.

## Usage

1. Prepare a text file containing the identifiers you want to filter, with one identifier per line.
2. Place the .xlsx file you want to filter in an accessible location.
3. Run the script using the following command:

```bash
python src/main.py <path_to_identifiers.txt> <path_to_input.xlsx> <path_to_output.xlsx>
```

Replace `<path_to_identifiers.txt>`, `<path_to_input.xlsx>`, and `<path_to_output.xlsx>` with the actual file paths.

## Example

Given a text file `identifiers.txt` with the following content:

```
ID1
ID2
ID3
```

And an Excel file `data.xlsx`, running the script will produce a new Excel file containing only the rows that match the identifiers from `identifiers.txt`.
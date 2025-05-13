import sys
import os
from utils.excel_processor import filter_excel, save_filtered_data
from utils.file_handler import read_clusters

def main():    
    if len(sys.argv) != 3:
        print("Usage: python main.py <cluster_file.txt> <input_excel.xlsx>")
        print("Output files will be created as <input>_<cluster>_filtered.xlsx")
        sys.exit(1)    
    
    cluster_file = sys.argv[1]
    input_excel_file = sys.argv[2]
    clusters = read_clusters(cluster_file)
    
    
    for cluster_name, identifiers in clusters.items():
        print(f"\nProcessing {cluster_name} with {len(identifiers)} identifiers...")
        
        output_dir = os.path.dirname(input_excel_file)
        base_name = os.path.basename(input_excel_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = os.path.join(output_dir, f"{name_without_ext}_{cluster_name}_filtered.xlsx")
        
        save_filtered_data(filtered_data, output_file)
        
        print(f"Filtered data for {cluster_name} saved to {output_file}")

if __name__ == "__main__":
    main()
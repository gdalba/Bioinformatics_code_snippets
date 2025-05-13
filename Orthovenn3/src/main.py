import sys
import os
from utils.excel_processor import filter_excel, save_filtered_data
from utils.file_handler import read_clusters

def main():    
    if len(sys.argv) != 4:
        print("Usage: python main.py <cluster_file.txt> <input_excel.xlsx> <species_id> <- Ensure to use the same name you used when uploading that species file to orthovenn3.")
        print("Output files will be created as <input>_<cluster>_filtered.xlsx")
        sys.exit(1)    
    
    cluster_file = sys.argv[1]
    input_excel_file = sys.argv[2]
    clusters = read_clusters(cluster_file)
    species_id = sys.argv[3]
    
    
    for cluster_name, identifiers in clusters.items():
        print(f"\nProcessing {cluster_name} with {len(identifiers)} identifiers...")
        # Filter Excel data based on identifiers for this cluster
        filtered_data = filter_excel(identifiers, input_excel_file, species_id)
        
        output_dir = os.path.dirname(input_excel_file)
        base_name = os.path.basename(input_excel_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = os.path.join(output_dir, f"{name_without_ext}_{cluster_name}_filtered.xlsx")
        
        save_filtered_data(filtered_data, output_file)
        
        print(f"Filtered data for {cluster_name} saved to {output_file}")

if __name__ == "__main__":
    main()
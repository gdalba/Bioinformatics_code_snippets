import sys
import re

'''
def read_identifiers(file_path):
    """
    Read identifiers from a text file.
    Each line in the file is considered as one identifier.
    """
    identifiers = []
    try:
        with open(file_path, 'r') as file:
            identifiers = [line.strip() for line in file if line.strip()]
        print(f"Read {len(identifiers)} identifiers from {file_path}")
        return identifiers
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        sys.exit(1)

'''

def read_clusters(file_path):
    """
    Read cluster information from a text file.
    The file format has lines like:
    cluster<N>    <count>    <identifier_1>;<identifier_2>;...    
    
    Returns:
        Dictionary mapping cluster names to lists of identifiers
    """
    clusters = {}
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                
                # Split by tabs or multiple spaces
                parts = re.split(r'\s+', line, maxsplit=2)
                if len(parts) < 3:
                    print(f"Skipping line with insufficient parts: {line}")
                    continue
                
                # First part is the cluster name
                if not parts[0].startswith('cluster'):
                    print(f"Skipping line that doesn't start with 'cluster': {line}")
                    continue
                    
                cluster_name = parts[0]
                count = int(parts[1])
                
                # Third part contains semicolon-separated identifiers
                identifiers_text = parts[2]
                identifiers = identifiers_text.split(';')
                
                # Clean up any extra identifiers that might be tab-separated at the end
                if len(parts) > 3:
                    for extra in parts[3:]:
                        if extra.startswith('unique_'):
                            identifiers.append(extra)
                
                clusters[cluster_name] = identifiers
                print(f"Found {cluster_name} with {count} identifiers, parsed {len(identifiers)} identifiers")
        
        if not clusters:
            print("No clusters found in the file. Check if the file format is correct.")
            sys.exit(1)
            
        total_identifiers = sum(len(ids) for ids in clusters.values())
        print(f"Read {len(clusters)} clusters with a total of {total_identifiers} identifiers")
        
        return clusters
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        sys.exit(1)
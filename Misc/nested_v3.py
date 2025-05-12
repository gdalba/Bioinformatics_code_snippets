import numpy as np
import argparse

# Define the genemerge function
def genemerge(G):
    fully_nested_genes = []
    partly_nested_genes = []
    Nfullynested = 0
    Npartlynested = 0

    k = 0
    while k < G.shape[1]:

        m = k + 1
        while m < G.shape[1]:

            # If the starting position of gene 'm' is before the ending position of gene 'k'
            if G[0, m] < G[1, k]:

                # If gene 'm' starts and ends inside gene 'k' (fully nested)
                if G[1, m] < G[1, k]:
                    fully_nested_genes.append((G[0, m], G[1, m]))
                    G = np.delete(G, m, axis=1)
                    Nfullynested += 1
                    continue

                # If gene 'm' starts inside gene 'k' but ends outside (partly nested)
                else:
                    # Check if 'm' and 'k' are the same
                    if not (G[0, m] == G[0, k] or G[1, m] == G[1, k]):
                        partly_nested_genes.append((G[0, m], G[1, m]))
                        G[:, k] = [G[0, k], G[1, m]]
                        G = np.delete(G, m, axis=1)
                        Npartlynested += 1
                        continue
                    else:
                        m += 1

            else:
                m += 1

        k += 1

    return G, Nfullynested, Npartlynested, fully_nested_genes, partly_nested_genes


# Main script starts here
# Read chromosome lengths
# Add argparse logic
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute nested genes from provided data files.")
    
    parser.add_argument('-NN', dest='NN_file', type=str, required=True,
                        help="Path to the chromosome lengths file. Expected format: one-length-per-line, with the target chromosome length at index 12 (13th line).")
    
    parser.add_argument('-Lidat', dest='Lidat_file', type=str, required=True,
                        help="Path to the gene lengths and positions file. Expected format: three columns, with the first two columns representing gene start and end positions respectively, and the third column representing the gene length.")
    args = parser.parse_args()

    # Read chromosome lengths
    NN = np.loadtxt(args.NN_file)
    N = NN[12]  # Python uses 0-based indexing

    # Read gene lengths and positions
    Lidat = np.loadtxt(args.Lidat_file)

    LidatL = Lidat[:, 2]  # 3rd column contains lengths
    Li = np.sort(LidatL)[::-1]  # sort in descending order

    # Sort based on the 1st column
    Liseqdat = Lidat[np.argsort(Lidat[:, 0])]

    # For experimental nesting:
    GE = Liseqdat[:, 0:2].T
    GGE, NfulE, NpartE, fully_nested_E, partly_nested_E = genemerge(GE)
   
    # Extract chromosome identifier from filename

    chromosome_identifier = args.Lidat_file.split("_")[-1].replace(".txt", "")
    print(f"Chr {chromosome_identifier}: NfullyNested = {NfulE} | NpartlyNested = {NpartE}")
#    print("Fully Nested Genes:", fully_nested_E)
#    print("Partly Nested Genes:", partly_nested_E)
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math 

scaffold_sizes = {
    "scaffold_1": 21919503,
    "scaffold_2": 19571729,
    "scaffold_3": 18510050,
    "scaffold_4": 16784600,
    "scaffold_5": 16175346,
    "scaffold_6": 15697342,
    "scaffold_7": 15457567,
    "scaffold_8": 14934027,
    "scaffold_9": 13864370,
    "scaffold_10": 11962182,
    "scaffold_11": 11619698,
    "scaffold_12": 11554503,
    "scaffold_13": 10559445
}

def parse_gff3(input_file):
    data = []
    with open(input_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                attr_dict = {field.split('=')[0]: field.split('=')[1] for field in columns[8].split(';')}
                classification = attr_dict.get('Classification', None)
                identity = attr_dict.get('Identity', None)
                method = attr_dict.get('Method', None)
                kimura = float(attr_dict.get('Kimura', '0'))
                data.append([columns[0], columns[2], int(columns[3]), int(columns[4]), columns[6],
                             classification, identity, method, kimura])
    df = pd.DataFrame(data, columns=['seqid', 'type', 'start', 'end', 'strand', 'classification', 'identity', 'method', 'kimura'])
    return df

def calculate_kimura_distance(te_df):
    kimura_distances = te_df['kimura']
    mean_kimura = np.mean(kimura_distances)
    median_kimura = np.median(kimura_distances)
    
    print(f"Mean Kimura Distance for Repeats: {mean_kimura}")
    print(f"Median Kimura Distance for Repeats: {median_kimura}")

    plt.figure(figsize=(10,5))
    sns.histplot(te_df['kimura'], kde=True, color='blue')
    plt.title('Distribution of Kimura distances')
    plt.xlabel('Kimura Distance')
    plt.ylabel('Frequency')
    plt.savefig('kimura_distance_plot.pdf', format='pdf')

def main(input_folder):
    te_gff3_file = os.path.join(input_folder, 'Mnemiopsis_leidyi_genome.fa.mod.EDTA.TEanno.gff3')
    te_df = parse_gff3(te_gff3_file)
    calculate_kimura_distance(te_df)

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2 and (sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        print('+' + '-'*67 + '+')
        print('| PLOT EDTA OUTPUT + genes V2 - MADE BY GABRIEL DALL\'ALBA - 06/12/2023|')
        print('+' + '-'*67 + '+')
        print('| Usage: python script_name.py input_folder                            |')
        print('| Requirements: Python 3, pandas, seaborn, matplotlib                  |')
        print('| Environments ready to run: edta_visualization                        |')
        print('+' + '-'*67 + '+')
    elif len(sys.argv) != 2:
        print('Usage: python script_name.py input_folder')
    else:
        input_folder = sys.argv[1]
        main(input_folder)

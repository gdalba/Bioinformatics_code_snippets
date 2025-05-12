import os
import pandas
import argparse
import pprint
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.stats import chi2_contingency
from sklearn.metrics import r2_score


def main():

    
    welcome_header = """
    ******************************************************
    *                                                    *
    *      Welcome to TE_log_scatterplot.py script!      *
    *                                                    *
    ******************************************************
    """

    
    script_name = "TE_log_scatterplot.py"
    creator = "Gabriel Dall'Alba"
    date = "2024-03-22"
    description = " Reads gff3 files of the three ctenophores and outputs a log scatterplot of the TEs by class."
    usage = "python TE_log_scatterplot.py -g1 <gff3_1> -g2 <gff3_2> -g3 <gff3_3> -o <output>"

    
    print(welcome_header)

    
    print(f"Script Name: {script_name}")
    print(f"Creator: {creator}")
    print(f"Date: {date}")
    print(f"Description: {description}")
    print(f"Usage: {usage}")

    
    parser = argparse.ArgumentParser(description=" Reads gff3 files of the three ctenophores and outputs a log scatterplot of the TEs by class.")
    parser.add_argument("-g1", "--gff3_1", help="Path to the first GFF3 file.")
    parser.add_argument("-g2", "--gff3_2", help="Path to the second GFF3 file.")
    parser.add_argument("-g3", "--gff3_3", help="Path to the third GFF3 file.")
    parser.add_argument("-o", "--output", help="Path to the output file.")
    args = parser.parse_args()

    print("Parsing GFF3 files...")
    print(f"First GFF3 file: {args.gff3_1}")
    print("Parsing...")
    te_classes_1 = parse_gff3(args.gff3_1)
    print(f"Second GFF3 file: {args.gff3_2}")
    print("Parsing...")
    te_classes_2 = parse_gff3(args.gff3_2)
    print(f"Third GFF3 file: {args.gff3_3}")
    print("Parsing...")
    te_classes_3 = parse_gff3(args.gff3_3)
    print("Parsing complete!")

    print("Counts for the first GFF3 file:")
    print(te_classes_1)
    print("__________________________________________________________________________________________________________________________")
    print("Counts for the second GFF3 file:")
    print(te_classes_2)
    print("__________________________________________________________________________________________________________________________")
    print("Counts for the third GFF3 file:")
    print(te_classes_3)
    print("__________________________________________________________________________________________________________________________")

    print("Converting dictionaries to dataframes...")
    # Convert dictionaries to pandas DataFrames
    df1 = pd.DataFrame(list(te_classes_1.items()), columns=['TE Class', 'Count'])
    df2 = pd.DataFrame(list(te_classes_2.items()), columns=['TE Class', 'Count'])
    df3 = pd.DataFrame(list(te_classes_3.items()), columns=['TE Class', 'Count'])
    
    print("Dataframes created!")
    
    print("Dataframe 1:")
    print(df1)

    # plots to do:
    # genome 1 x genome 2
    # genome 1 x genome 3
    # genome 3 x genome 2
    
    # I'm giving hormiphora as 1, mle as 2, bol as 3 so
    # hormiphora x mnemiopsis
    # hormiphora x bolinopsis
    # bolinopsis x mnemiopsis


    df1['GFF3 File'] = os.path.basename(args.gff3_1) #will be the red squares -- this line is adding a column with identifier of the file
    df2['GFF3 File'] = os.path.basename(args.gff3_2) #will be the green triangles
    df3['GFF3 File'] = os.path.basename(args.gff3_3) #will be the blue circles

    
    #genome 1 x genome 2 -- hormiphora x mnemiopsis

    # Merge the two DataFrames on 'TE Class'
    print("Merging dataframes...")	
    merged_df = pd.merge(df1, df2, on='TE Class', suffixes=('_genome1', '_genome2'))

    # Log transformation
    merged_df['Log Count_genome1'] = np.log(merged_df['Count_genome1'] + 1) # Adding 1 to avoid log(0)
    merged_df['Log Count_genome2'] = np.log(merged_df['Count_genome2'] + 1)

    print("Plotting scatterplot...")
    
    plt.figure(figsize=(10, 6))

    te_classes = merged_df['TE Class'].unique()

    markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|', '_', '1', '2', '3', '4']
    colors = plt.cm.rainbow(np.linspace(0, 1, len(merged_df['TE Class'].unique())))

    # Iterate over each TE Class and filter data for the current TE Class
    for te_class, (marker, color) in zip(te_classes, zip(markers, colors)):    
        data = merged_df[merged_df['TE Class'] == te_class]
    
    
        plt.scatter(data['Log Count_genome1'], data['Log Count_genome2'], label=te_class, marker=marker, color=color, edgecolors='black')
    
    # Regression line for all data
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_df['Log Count_genome1'], merged_df['Log Count_genome2'])
    line = slope * merged_df['Log Count_genome1'] + intercept
    plt.plot(merged_df['Log Count_genome1'], line, color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, $R^2$={r_value**2:.2f}')

    plt.xlabel(f'Log Count Genome 1  - {args.gff3_1}')
    plt.ylabel(f'Log Count Genome 2  - {args.gff3_2}')
    plt.title('Scatterplot of TE Counts between Two Genomes by TE Class')
    plt.legend(title='TE Class', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    plt.savefig(args.output)
    print(f"Scatterplot saved to {args.output}")

    # genome 1 x genome 3 -- hormiphora x bolinopsis

    merged_df = pd.merge(df1, df3, on='TE Class', suffixes=('_genome1', '_genome2'))

    # Log transformation
    merged_df['Log Count_genome1'] = np.log(merged_df['Count_genome1'] + 1) # Adding 1 to avoid log(0)
    merged_df['Log Count_genome2'] = np.log(merged_df['Count_genome2'] + 1)

    print("Plotting scatterplot...")
    
    plt.figure(figsize=(10, 6))

    te_classes = merged_df['TE Class'].unique()

    markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|', '_', '1', '2', '3', '4']
    colors = plt.cm.rainbow(np.linspace(0, 1, len(merged_df['TE Class'].unique())))

    # Iterate over each TE Class and filter data for the current TE Class
    for te_class, (marker, color) in zip(te_classes, zip(markers, colors)):    
        data = merged_df[merged_df['TE Class'] == te_class]
    
    
        plt.scatter(data['Log Count_genome1'], data['Log Count_genome2'], label=te_class, marker=marker, color=color, edgecolors='black')
    
    # Regression line for all data
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_df['Log Count_genome1'], merged_df['Log Count_genome2'])
    line = slope * merged_df['Log Count_genome1'] + intercept
    plt.plot(merged_df['Log Count_genome1'], line, color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, $R^2$={r_value**2:.2f}')

    plt.xlabel(f'Log Count Genome 1  - {args.gff3_1}')
    plt.ylabel(f'Log Count Genome 2  - {args.gff3_3}')
    plt.title('Scatterplot of TE Counts between Two Genomes by TE Class')
    plt.legend(title='TE Class', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    plt.savefig('hor_bol.pdf')
    print(f"Scatterplot saved to hor_bol.pdf")

    # genome 2 x genome 3 -- mnemiopsis x bolinopsis

    merged_df = pd.merge(df3, df2, on='TE Class', suffixes=('_genome1', '_genome2'))

    # Log transformation
    merged_df['Log Count_genome1'] = np.log(merged_df['Count_genome1'] + 1) # Adding 1 to avoid log(0)
    merged_df['Log Count_genome2'] = np.log(merged_df['Count_genome2'] + 1)

    print("Plotting scatterplot...")
    
    plt.figure(figsize=(10, 6))

    te_classes = merged_df['TE Class'].unique()

    markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|', '_', '1', '2', '3', '4']
    colors = plt.cm.rainbow(np.linspace(0, 1, len(merged_df['TE Class'].unique())))

    # Iterate over each TE Class and filter data for the current TE Class
    for te_class, (marker, color) in zip(te_classes, zip(markers, colors)): 
        data = merged_df[merged_df['TE Class'] == te_class]
    
    
        plt.scatter(data['Log Count_genome1'], data['Log Count_genome2'], label=te_class, marker=marker, color=color, edgecolors='black')
    
    # Regression line for all data
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_df['Log Count_genome1'], merged_df['Log Count_genome2'])
    line = slope * merged_df['Log Count_genome1'] + intercept
    plt.plot(merged_df['Log Count_genome1'], line, color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, $R^2$={r_value**2:.2f}')

    plt.xlabel(f'Log Count Genome 1  - {args.gff3_3}')
    plt.ylabel(f'Log Count Genome 2  - {args.gff3_2}')
    plt.title('Scatterplot of TE Counts between Two Genomes by TE Class')
    plt.legend(title='TE Class', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    plt.savefig('scatterplot_log_sns_bol_mle.pdf')
    print(f"Scatterplot saved to scatterplot_log_sns_bol_mle.pdf")

def parse_gff3(gff3_file):
    te_classes = {}
    with open(gff3_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            attributes = line[8].split(';')
            for attribute in attributes:
                if attribute.startswith('Classification='):
                    classification = attribute.split('=')[1]
                    if classification not in te_classes:
                        te_classes[classification] = 1
                    else:
                        te_classes[classification] += 1
    return te_classes

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

def generate_heatmap(input_file, output_file):
    try:
        # Step 1: Load the data
        data = pd.read_csv(input_file, sep="\t")

        # Step 2: Preprocess the data
        data.set_index('clade_name', inplace=True)  # Set 'clade_name' as index
        data = data.drop('UNCLASSIFIED', errors='ignore')  # Remove 'UNCLASSIFIED' row if it exists

        # Step 3: Create the heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(data.astype(float), cmap="viridis", cbar=True)

        # Customize heatmap
        plt.title("Heatmap of Metaphlan Abundance Data", fontsize=16)
        plt.xlabel("Samples", fontsize=12)
        plt.ylabel("Clades", fontsize=12)

        # Save the plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        print(f"Heatmap saved as {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    # Ensure the script is run with the correct arguments
    if len(sys.argv) != 3:
        print("Usage: python generate_heatmap.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' does not exist.")
        else:
            generate_heatmap(input_file, output_file)

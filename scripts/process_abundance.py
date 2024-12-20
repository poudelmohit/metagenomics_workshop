#!/usr/bin/env python3
import pandas as pd

# Load the CSV file into a DataFrame
file_path = 'modified_abundance_table.txt'  # Replace with your file path
df = pd.read_csv(file_path, sep="\t")  # Adjust the delimiter if necessary

df['clade_name'].str.split('|', expand = True)
df[['kingdom', 'phylum', 'class', 'order' , 'family', 'genus', 'species', 'type']] = df['clade_name'].str.split('|', expand=True)

df = df.drop('clade_name', axis=1)

# Apply the transformation only to the last 8 columns
last_columns = df.columns[-8:]
df[last_columns] = df[last_columns].apply(lambda col: col.str.replace(r'^.*?__', '', regex=True))

# Remove leading and trailing spaces from the last 8 columns
df[last_columns] = df[last_columns].apply(lambda col: col.str.strip())

# Reorder columns so that the last 8 columns come to the first
df = pd.concat([df[last_columns], df.drop(last_columns, axis=1)], axis=1)

# Save the new DataFrame with updated columns (if needed)
df.to_csv('updated_abundance_table.csv', sep="\t", index=False)

# Print the updated DataFrame
print(df)

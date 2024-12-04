import os
import shutil
import pandas as pd

# Define the directory containing the FASTA files and the path to the CSV file
fasta_dir = './fasta'
csv_file_path = 'attributes_merged_300k.csv'  # Assuming the file is in the current directory
filtered_dir = './filtered'

# Define the completeness percentage threshold
threshold = 100.0

# Create the filtered directory if it does not exist
if not os.path.exists(filtered_dir):
    os.makedirs(filtered_dir)

# Read the CSV file into a DataFrame
try:
    df = pd.read_csv(csv_file_path, low_memory=False)
except ImportError as e:
    print(f"ImportError: {e}")
    print("Make sure you have installed pandas library.")
    print("Use 'conda install pandas' to install it.")

# Print the DataFrame to verify it was read correctly
# print(df)

# List all FASTA files in the directory
fasta_files_in_dir = [f for f in os.listdir(fasta_dir) if f.endswith('.fna.gz')]

# Filter the FASTA files based on the completeness percentage
for fasta_file in fasta_files_in_dir:
    # Check if the FASTA file is in the CSV file
    if fasta_file in df['filePath'].values:
        # Get the completeness percentage
        completeness_percentage = df[df['filePath'] == fasta_file]['checkm_completeness'].values[0]
        # Check if it meets the threshold
        if completeness_percentage >= threshold:
            # Move the file to the filtered directory
            shutil.move(os.path.join(fasta_dir, fasta_file), os.path.join(filtered_dir, fasta_file))



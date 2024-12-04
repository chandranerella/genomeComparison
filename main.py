import itertools
import pandas as pd
import subprocess
import glob
import os
import time
from scipy.stats import kendalltau
from scipy.stats import spearmanr



def log_message(message, log_file="log_file.log"):
    with open(log_file, "a") as f:
        f.write(message + "\n")
        

def run_pyani(input_dir, output_dir, log_file):
    command = [
        "average_nucleotide_identity.py",
        "-i", input_dir,
        "-o", output_dir,
        "-m", "ANIm"
    ]
    
    try:
        start_time = time.time()
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        log_message(f"pyANI elapsed time: {elapsed_time:.2f} seconds", log_file)
        df = pd.read_csv(f'{output_dir}/ANIm_percentage_identity.tab', sep='\t', index_col=0)
        # Get the list of genomes
        genomes = df.index.tolist() 
        
        # Create a list to store the pairwise comparisons
        pairwise_comparisons = []
        
        # Iterate through all unique pairs of genomes
        for genome1, genome2 in itertools.combinations(genomes, 2):
            ani_value = df.loc[genome1, genome2]
            pairwise_comparisons.append({
                'Pair': f"{genome1} vs {genome2}",
                'pyANI': ani_value * 100  # Convert to percentage if necessary
            })
        result_df = pd.DataFrame(pairwise_comparisons)
        result_df.set_index('Pair', inplace=True)
        result_df.index.name = None
        genome_pairs = list(itertools.combinations(genomes, 2))
        formatted_genomes = [f"{g1} vs {g2}" for g1, g2 in genome_pairs]
        return result_df, formatted_genomes
    except subprocess.CalledProcessError as e:
        print("Error running pyANI command:", e)
        print("Error output:", e.stderr)

def run_fastani(reference_list, log_file):
    # Implementation for fastANI
    original_dir = os.getcwd()
    build_path = "../fastani/FastANI/build/"
    os.chdir(build_path)
    print("fastANI has started.")
    try:
        command = ["./fastANI", "--ql", f"{original_dir}/{reference_list}", "--rl", f"{original_dir}/{reference_list}", "-o", f"{original_dir}/fastani_output.txt"]
        start_time = time.time()
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"FastANI elapsed time: {elapsed_time:.2f} seconds")
        with open("example.txt", "w") as file:
            file.write(f"FastANI elapsed time: {elapsed_time:.2f} seconds")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        print(f"Error output: {e.stderr}")
    finally:
    # Change back to the original directory
        os.chdir(original_dir)
    
        # Read the file
    with open('fastani_output.txt', 'r') as file:
        lines = file.readlines()

    # Process each line
    data = []
    seen_pairs = set()
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            genome1 = extract_genome_name(parts[0])
            genome2 = extract_genome_name(parts[1])
            ani = float(parts[2])
            
            # Skip self-comparisons
            if genome1 == genome2:
                continue
            
            # Ensure unique pairs
            pair = tuple(sorted([genome1, genome2]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            
            data.append({
                'Genome1': genome1,
                'Genome2': genome2,
                'fastANI': ani
            })

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Set the index to be "Genome1 vs Genome2"
    df['Pair'] = df.apply(lambda row: f"{row['Genome1']} vs {row['Genome2']}", axis=1)
    df.set_index('Pair', inplace=True)
    df.index.name = None
    # Keep only the ANI column
    df = df[['fastANI']]

    print(df)
    return df
    
    
def extract_genome_name(path):
    return os.path.basename(path).split('.fna')[0]

def run_sourmash(repository, final_df, log_file):
    # Implementation for sourmash
    
    try:
    # Run the command
        fna_files = glob.glob(os.path.join(repository, "*.fna"))
        output_dir = "signature_output/"

        command = [
            "sourmash", "sketch", "dna", 
            "-p", "scaled=1000,k=10"
        ] + fna_files + ["--outdir", output_dir]
        # command = [
        #     "sourmash", "sketch", "dna",
        #     "-p", "scaled=1000,k=31",
        #     f"{repository}/*",
        #     "--outdir", "signature_output/"
        # ]
        start_time = time.time()
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        end_time = time.time()
        elapsed_time1 = end_time - start_time
        log_message(f"Sourmash sketch elapsed time: {end_time - start_time:.2f} seconds", log_file)
        # Print the output
        print("Command output:")
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        print(f"Error output: {e.stderr}")  
        
    #Now comparison of signatures
    # Define the directory containing your signature files
    # Define the directory containing the .sig files
    signature_dir = "signature_output/*.sig"
    output_npy = "comparison_output.npy"
    output_csv = "comparison_output.csv"
    
    # Construct the list of .sig files in the signature directory
    sig_files = glob.glob(signature_dir)

    # Construct the command
    command = [
        "sourmash", "compare"
    ] + sig_files + ["-o", output_npy, "--csv", output_csv]

    try:
        start_time = time.time()
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        end_time = time.time()
        elapsed_time2 = end_time - start_time
        log_message(f"Sourmash compare elapsed time: {end_time - start_time:.2f} seconds", log_file)
        log_message(f"Sourmash total elapsed time: {elapsed_time1 + elapsed_time2:.2f} seconds", log_file)
        print("Sourmash compare command executed successfully.")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        print(f"Error output: {e.stderr}")
    # Load the Sourmash comparison results
    comparison_df = pd.read_csv("comparison_output.csv", header=None)

    # Exclude the first row and column with file names
    data = comparison_df.iloc[1:, 0:].astype(float)

    # Initialize a list to store the extracted values
    upper_triangular_values = []

    # Iterate over the rows
    for i in range(0, data.shape[0] + 1):  # Iterate through rows
        for j in range(i + 1, data.shape[1]):  # Start from the column after the diagonal
            upper_triangular_values.append(round(data.iloc[i, j], 4))  # Round to 4 decimal places

    # Check if the length of the values matches the expected number of pairwise comparisons
    final_df['sourmash'] = upper_triangular_values
            

def run_bindash(input_file, final_df, log_file):
    # Implementation for Bindash 
    output_file = "bindash_sketches"

    # Ensure the input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found")

    # Construct the command
    command = [
        "bindash",
        "sketch",
        f"--listfname={input_file}",
        f"--outfname={output_file}",
        "--kmerlen=10"
    ]
    try:
        # Run the command
        start_time = time.time()
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        end_time = time.time()
        elapsed_time1 = end_time - start_time
        log_message(f"BinDash sketch elapsed time: {end_time - start_time:.2f} seconds", log_file)
        print("BinDash sketch command executed successfully.")
        print("Output:")
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running BinDash: {e}")
        print("Error output:")
        print(e.stderr)

    sketch_file = output_file    
    output_file = "bindash_distances.txt"
    command = [
        "bindash",
        "dist",
        f"--outfname={output_file}",
        f"{sketch_file}"
    ]
    
    try:
        # Run the command
        start_time = time.time()
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        end_time = time.time()
        elapsed_time2 = end_time - start_time
        log_message(f"BinDash dist elapsed time: {end_time - start_time:.2f} seconds", log_file)
        total_time = elapsed_time1 + elapsed_time2
        log_message(f"BinDash total elapsed time: {elapsed_time1 + elapsed_time2:.2f} seconds", log_file)
        print("BinDash distance calculation executed successfully.")
        print("Output:")
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running BinDash: {e}")
        print("Error output:")
        print(e.stderr)
    
    # Load the Bindash distances file, focusing on the last column for decimal extraction
    bindash_df = pd.read_csv("bindash_distances.txt", sep='\t', header=None, usecols=[4])

    # Extract decimal values from the last column in bindash_df
    bindash_values = bindash_df[4].str.split('/').apply(lambda x: round(int(x[0]) / int(x[1]), 4))

    # Map the extracted values to the Bindash column in final_df
    final_df["Bindash"] = bindash_values.values
        
def main():
    # List of input genomes
    
    # Initialize results dictionary
    # Perform comparisons
    log_file = "log_file.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    
    # pyani_df, formatted_genomes = run_pyani("input_filtered", "output_pyani", log_file=log_file)
    # final_df = pd.DataFrame(index=formatted_genomes, columns=['pyANI', 'fastANI', 'sourmash', 'Bindash'])
    # final_df.update(pyani_df) 
    # First, ensure that the 'Pair' column in both dataframes match exactly
    fastani_df = run_fastani("reference_fastani.txt", log_file)
    # final_df.update(fastani_df)
    
    
    #function for sourmash comparison
    # run_sourmash("input_filtered", final_df, log_file)
    
    # run_bindash("reference_fastani.txt", final_df, log_file)
  
    final_df.to_csv("final_df.csv", index=True)
    
    # Create DataFrame and save to CSV
    ranked_df = final_df.copy()
    for col in ranked_df.columns[0:]:  # Skip the first column, "Genome Pair"
        ranked_df[col] = ranked_df[col].rank(ascending=False, method='min').astype(int)
        
    # Extract rank values for each tool from rank_df
    pyANI_ranks, fastANI_ranks, sourmash_ranks, binDash_ranks = [
        ranked_df[col].values for col in ranked_df.columns
    ]
    # Calculate Kendall's Tau for each pair
    ranked_df.to_csv("ranked_df.csv", index=True)
    tool_pairs_kendall = {
        pair: kendalltau(*rank_pair)
        for pair, rank_pair in [
            ("pyANI vs FastANI", (pyANI_ranks, fastANI_ranks)),
            ("pyANI vs Sourmash", (pyANI_ranks, sourmash_ranks)),
            ("pyANI vs BinDash", (pyANI_ranks, binDash_ranks)),
            ("FastANI vs Sourmash", (fastANI_ranks, sourmash_ranks)),
            ("FastANI vs BinDash", (fastANI_ranks, binDash_ranks)),
            ("Sourmash vs BinDash", (sourmash_ranks, binDash_ranks))
        ]
    }
    for pair, (tau, p_value) in tool_pairs_kendall.items():
        log_message(f"Kendall's Tau for {pair}: {tau:.4f}, p-value: {p_value:.4f}", log_file)


    tool_pairs_spearman = {
        pair: spearmanr(*rank_pair)
        for pair, rank_pair in [
            ("pyANI vs FastANI", (pyANI_ranks, fastANI_ranks)),
            ("pyANI vs Sourmash", (pyANI_ranks, sourmash_ranks)),
            ("pyANI vs BinDash", (pyANI_ranks, binDash_ranks)),
            ("FastANI vs Sourmash", (fastANI_ranks, sourmash_ranks)),
            ("FastANI vs BinDash", (fastANI_ranks, binDash_ranks)),
            ("Sourmash vs BinDash", (sourmash_ranks, binDash_ranks))
        ]
    }
    for pair, (rho, p_value) in tool_pairs_spearman.items():
        log_message(f"Spearman's rho for {pair}: {rho:.4f}, p-value: {p_value:.4f}", log_file)



if __name__ == "__main__":
    main()
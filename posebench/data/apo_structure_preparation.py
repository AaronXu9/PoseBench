import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import shutil
import re

def update_paths(input_csv, output_csv, base_protein_dir, base_ligand_dir):
    # Read input CSV file
    with open(input_csv, 'r') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        
        # Open output CSV file for writing
        with open(output_csv, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for row in reader:
                # Extract the complex name from the current row
                complex_name = row['complex_name']
                
                # Modify the protein path and ligand path
                protein_path = os.path.join(base_protein_dir, f"{complex_name}", f"{complex_name}_protein.pdb")
                ligand_path = os.path.join(base_ligand_dir, f"{complex_name}", f"{complex_name}_ligand.sdf")
                # check if the protein and ligand paths exist
                if not os.path.exists(protein_path):
                    print(f"Protein path {protein_path} does not exist.")
                if not os.path.exists(ligand_path):
                    print(f"Ligand path {ligand_path} does not exist.")
                # Convert the ligand SDF file to SMILES
                ligand_smile = Chem.MolToSmiles(Chem.MolFromMolFile(ligand_path))
                # Update the row with the new paths
                row['protein_path'] = protein_path
                row['ligand_description'] = ligand_smile  # Assuming you want to change this
                # Write the updated row to the new file
                writer.writerow(row)

def copy_protein_files(source_dir, target_dir, file_extension=".pdb"):
    """
    Copies protein files with a given extension from subdirectories of the source directory
    into a target directory.

    Args:
        source_dir (str): Path to the source directory containing multiple subdirectories.
        target_dir (str): Path to the target directory where files will be copied.
        file_extension (str): The file extension to filter for (default is '.pdb').
    """
    # Create the target directory if it doesn't exist
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    
    # Traverse through the source directory
    for subdir in os.listdir(source_dir):
        subdir_path = os.path.join(source_dir, subdir)
        
        # Ensure we are processing only directories
        if os.path.isdir(subdir_path) and re.match(r'^[A-Za-z0-9]{4}_[A-Za-z0-9]{3}$', subdir):
            # Look for files with the specified extension
            for file_name in os.listdir(subdir_path):
                if file_name.endswith(file_extension):
                    source_file = os.path.join(subdir_path, file_name)
                    target_file = os.path.join(target_dir, file_name)
                    
                    # Copy the file to the target directory
                    shutil.copy(source_file, target_file)
                    print(f"Copied {source_file} to {target_file}")

if __name__ == "__main__":
    method = "diffdock"  # or "dynamic_bind"
    if method == "diffdock":
        # Example usage
        root_dir = '/home/aoxu/projects/PoseBench/'
        dataset = 'posebusters_benchmark_set' # or 'astex_diverse_set'
        input_csv_name = 'posebusters_benchmark' # or 'astex_diverse'
        input_csv = root_dir + f'forks/DiffDock/inference/diffdock_{input_csv_name}_inputs.csv'
        output_csv = root_dir + f'forks/DiffDock/inference/diffdock_{input_csv_name}_input_orig.csv'
        base_protein_dir = f'/home/aoxu/projects/PoseBench/data/{dataset}'
        base_ligand_dir = f'/home/aoxu/projects/PoseBench/data/{dataset}'

        update_paths(input_csv, output_csv, base_protein_dir, base_ligand_dir)

    elif method == "dynamic_bind":
        # Example usage
        source_directory = "/home/aoxu/projects/PoseBench/data/astex_diverse_set"
        target_directory = "/home/aoxu/projects/PoseBench/data/astex_diverse_set/astex_diverse_orig"

        copy_protein_files(source_directory, target_directory)

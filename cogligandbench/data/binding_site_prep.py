import rootutils
import pandas as pd
import os
rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
from cogligandbench.data.binding_site_crop_preparation import get_binding_site_residue_indices, crop_protein_binding_site

def save_cropped_binding_site(csv_path: str, output_dir: str):
    # Path to the CSV file
    csv_path = "forks/GNINA/inference/gnina_plinder_benchmark_inputs.csv"

    # Output directory for cropped proteins
    output_dir = "/mnt/katritch_lab2/aoxu/data/plinder/cropped_binding_sites"
    os.makedirs(output_dir, exist_ok=True)

    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    for idx, row in df.iterrows():
        complex_name = row['complex_name']
        protein_path = row['protein_path']
        ligand_path = row['ligand_path']
        
        print(f"Processing {complex_name}...")
        
        try:
            # Get binding site residue indices
            binding_site_indices = get_binding_site_residue_indices(
                protein_filepath=protein_path,
                ligand_filepath=ligand_path,
                protein_ligand_distance_threshold=4.0,
                num_buffer_residues=7
            )
            
            # Crop and save the binding site
            crop_protein_binding_site(
                protein_filepath=protein_path,
                binding_site_residue_indices=binding_site_indices,
                output_dir=output_dir,
                pdb_id=complex_name,
                filename_midfix="",
                filename_suffix=""
            )
            
            print(f"  Successfully cropped binding site with {len(binding_site_indices)} residues")
            
        except Exception as e:
            print(f"  Error processing {complex_name}: {str(e)}")

    print(f"Processed {len(df)} complexes. Results saved to {output_dir}")

def prep_cropped_input_csv(original_csv_path: str, cropped_dir: str, new_csv_path: str):
    # Read the original CSV file
    df = pd.read_csv(original_csv_path)

    # Create a new column with paths to cropped proteins
    df['cropped_protein_path'] = df['complex_name'].apply(
        lambda x: os.path.join(cropped_dir, f"{x}_protein.pdb")
    )

    # Verify all files exist
    missing_files = df[~df['cropped_protein_path'].apply(os.path.exists)]
    if len(missing_files) > 0:
        print(f"Warning: {len(missing_files)} cropped protein files not found:")
        for idx, row in missing_files.iterrows():
            print(f"  - {row['complex_name']}: {row['cropped_protein_path']}")

    # Create a new DataFrame with the updated protein paths
    new_df = df.copy()
    new_df['protein_path'] = new_df['cropped_protein_path']
    new_df = new_df.drop(columns=['cropped_protein_path'])

    # Save the new CSV file
    new_df.to_csv(new_csv_path, index=False)

    print(f"Created new CSV file with updated protein paths: {new_csv_path}")
    print(f"Total entries: {len(new_df)}")

if __name__ == '__main__':
    csv_path = "forks/GNINA/inference/gnina_plinder_benchmark_inputs.csv"
    output_dir = "/mnt/katritch_lab2/aoxu/data/plinder/cropped_binding_sites"
    new_csv_path = "forks/DiffDock/inference/diffdock_pocket_only_plinder_benchmark_inputs.csv"
    # save_cropped_binding_site(csv_path, output_dir)
    prep_cropped_input_csv(csv_path, output_dir, new_csv_path)

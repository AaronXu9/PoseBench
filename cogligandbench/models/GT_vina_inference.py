import os
import subprocess
import sys
import math
from typing import List
import re 
import gzip 
import time 
import logging
from omegaconf import OmegaConf
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
# If you have custom classes for parsing .pdbqt, import them
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
import rootutils

rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
from cogligandbench.utils.log import setup_base_logging, get_custom_logger

# # Configure logging: logs will be written to 'docking.log'
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s [%(levelname)s] %(message)s',
#     handlers=[
#         logging.FileHandler('gnina_timing.log'),
#         logging.StreamHandler()  # optional: this sends log output to the console as well
#     ]
# )

# def get_custom_logger(logger_name: str, log_filename: str) -> logging.Logger:
#     """
#     Returns a logger configured with a FileHandler that writes to log_filename.
#     Any existing handlers will be removed.
#     """
#     logger = logging.getLogger(logger_name)
#     # Remove any existing handlers attached to this logger
#     if logger.hasHandlers():
#         logger.handlers.clear()
    
#     logger.setLevel(logging.INFO)
    
#     # Create file handler with the custom filename
#     file_handler = logging.FileHandler(log_filename)
#     file_handler.setLevel(logging.INFO)
    
#     # Define a formatter and set it for the file handler
#     formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
#     file_handler.setFormatter(formatter)
    
#     logger.addHandler(file_handler)
    
#     # Optionally, also log to the console by adding a StreamHandler
#     console_handler = logging.StreamHandler()
#     console_handler.setLevel(logging.INFO)
#     console_handler.setFormatter(formatter)
#     logger.addHandler(console_handler)
    
#     return logger


def compute_ligand_center_and_size(ligand_sdf):
    """
    Parses ligand coordinates from an SDF file.
    Returns the center (center_x, center_y, center_z)
    and size (size_x, size_y, size_z) suitable for defining
    the AutoDock Vina box.
    """
    x_coords, y_coords, z_coords = [], [], []
    with open(ligand_sdf, 'r') as f:
        lines = f.readlines()

    atom_block = False
    for line in lines:
        if 'END' in line:
            break
        parts = line.strip().split()
        if len(parts) >= 4:  # Very loose check
            try:
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
                atom_block = True
            except ValueError:
                if atom_block:
                    break

    if not x_coords:
        raise ValueError(f"No atomic coordinates found in {ligand_sdf}. Check file format.")

    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)
    min_z, max_z = min(z_coords), max(z_coords)

    center_x = (max_x + min_x) / 2.0
    center_y = (max_y + min_y) / 2.0
    center_z = (max_z + min_z) / 2.0

    size_x = (max_x - min_x) + 5.0
    size_y = (max_y - min_y) + 5.0
    size_z = (max_z - min_z) + 5.0

    return (center_x, center_y, center_z), (size_x, size_y, size_z)


def prepare_and_run_vina(protein_pdb, ligand_sdf, out_dir, exhaustiveness=8):
    """
    Converts the protein PDB and ligand SDF to PDBQT, computes the
    box, and runs AutoDock Vina. Docked pose is written to out_dir/vina_out.pdbqt.
    Returns the path to the Vina output file for further processing.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    protein_pdbqt = os.path.join(out_dir, "protein.pdbqt")
    ligand_pdbqt  = os.path.join(out_dir, "ligand.pdbqt")
    vina_out      = os.path.join(out_dir, "vina_out.pdbqt")

    # Convert protein from PDB -> PDBQT
    cmd_prot = f"obabel -i pdb {protein_pdb} -o pdbqt -O {protein_pdbqt} --partialcharge gasteiger -xh"
    subprocess.run(cmd_prot, shell=True, check=True)

    # Convert ligand from SDF -> PDBQT
    cmd_lig = f"obabel -i sdf {ligand_sdf} -o pdbqt -O {ligand_pdbqt} --partialcharge gasteiger -xh"
    subprocess.run(cmd_lig, shell=True, check=True)

    # Compute the docking box from ligand
    center, size = compute_ligand_center_and_size(ligand_sdf)
    cx, cy, cz = center
    sx, sy, sz = size

    # Run AutoDock Vina
    cmd_vina = (
        f"vina --receptor {protein_pdbqt} --ligand {ligand_pdbqt} "
        f"--center_x {cx:.3f} --center_y {cy:.3f} --center_z {cz:.3f} "
        f"--size_x {sx:.3f} --size_y {sy:.3f} --size_z {sz:.3f} "
        f"--exhaustiveness {exhaustiveness} --out {vina_out}"
    )
    subprocess.run(cmd_vina, shell=True, check=True)
    
    return vina_out


def extract_and_write_top_poses(vina_pdbqt_file: str, out_dir: str, prefix: str = "docked", 
                                remove_hs=True, top_n=10):
    """
    Parse the Vina output (PDBQT) which may contain multiple poses.
    Write the top N poses to individual .sdf files in out_dir.

    :param vina_pdbqt_file: Path to the .pdbqt file produced by Vina, containing multiple poses.
    :param out_dir: Where to write the .sdf files.
    :param prefix: A prefix to use in naming the .sdf files.
    :param remove_hs: Whether to remove hydrogens in final .sdf files.
    :param top_n: Number of top conformations to extract.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # -- Example parse of Vina scores from REMARK lines
    scores = []
    with open(vina_pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('REMARK VINA RESULT:'):
                parts = line.strip().split()
                if len(parts) >= 4:
                    score = float(parts[3])
                    scores.append(score)

    # scores[i] corresponds to pose i (the order they appear in the file)
    # We assume we have at least as many poses as lines with "REMARK VINA RESULT".
    # In practice, you would parse the actual poses out carefully.
    
    # -- Here we assume you have some code to convert each conformation from the
    #    multi-POSE PDBQT into an RDKit molecule. For example:
    #    1. Use a custom parser that returns a list of RDKit molecules (one per pose).
    #    2. Or run "obabel -ipdbqt vina_out.pdbqt -osdf multi.sdf --uniq" and then read them with RDKit.
    # For demonstration, let's assume you have a function get_rdkit_poses(vina_pdbqt_file)
    # that returns (list_of_molecules, list_of_scores). We'll show a placeholder:

    # Placeholder: you MUST implement something like this to extract actual poses
    # from your .pdbqt or from a converted .sdf. We’ll keep this super-simplified:
    def get_rdkit_poses(pdbqt_path: str):
        """
        Parse a multi-POSE PDBQT file into a list of RDKit Mol objects.

        Returns:
            (all_mols, all_scores_placeholder)
            - all_mols: A list of RDKit Mol objects (one for each pose or ligand).
            Each Mol may have multiple conformers if the PDBQT had multiple poses.
            - all_scores_placeholder: Return an empty list or None here,
            since we parse scores separately from REMARK lines.
        """
        # 1) Parse the multi-pose PDBQT file into a PDBQTMolecule object
        pdbqt_mol = PDBQTMolecule.from_file(pdbqt_path, skip_typing=True)
        
        # 2) Convert the PDBQTMolecule into RDKit Mol objects
        #    Typically returns a list of Mols. 
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        
        # 3) You can do any additional filtering or rearranging as needed.
        #    For example, skip any None entries.
        all_mols = []
        for mol in rdkit_mols:
            if mol is not None:
                all_mols.append(mol)
        
        # 4) We return an empty list for 'all_scores' because 
        #    in the main code, we parse the scores from the REMARK lines.
        return all_mols, []
    
    # Actually retrieve the poses. In real usage, each pose might be an RDKit Mol.
    # For demonstration, assume it returns a list of length = len(scores).
    all_poses, all_scores = get_rdkit_poses(vina_pdbqt_file)

    # If your 'get_rdkit_poses' doesn’t already parse scores, use the 'scores' list we collected from remarks
    # to keep them consistent. Otherwise, you can rely on all_scores if it’s accurate.
    if not all_scores and len(scores) == len(all_poses):
        all_scores = scores

    if not all_poses:
        print(f"[WARNING] No poses found in {vina_pdbqt_file}. Skipping extraction.")
        return

    # Sort the poses by score if needed (lowest score is best for Vina).
    # We'll create an index-based sort:
    # pairs like (pose_index, score)
    indexed_scores = list(enumerate(all_scores))
    # Sort by the score value
    indexed_scores.sort(key=lambda x: x[1])

    # Now keep only top_n
    top_indices = indexed_scores[:top_n]

    for rank, (pose_idx, pose_score) in enumerate(top_indices, start=1):
        rdkit_mol = all_poses[pose_idx]

        if remove_hs:
            rdkit_mol = Chem.RemoveHs(rdkit_mol)

        # Write out as .sdf
        out_sdf = os.path.join(out_dir, f"{prefix}_pose{rank}_score{pose_score:.2f}.sdf")
        # RDKit needs the conformer. If your mol has multiple conformers, ensure
        # you select the correct one. For example:
        # conformer_id = 0
        # tmp_mol = Chem.Mol(rdkit_mol)
        # tmp_mol.RemoveAllConformers()
        # tmp_mol.AddConformer(rdkit_mol.GetConformer(conformer_id))

        # For now, just write it directly if each RDKit Mol is a single conformer:
        Chem.MolToMolFile(rdkit_mol, out_sdf)

    print(f"[INFO] Wrote top-{top_n} poses to SDF in {out_dir}")


def run_gnina_docking(protein_pdb, ligand_sdf, out_dir, exhaustiveness=8):
    gnina_command = f"forks/GNINA/gnina --receptor {protein_pdb} --ligand {ligand_sdf} --autobox_ligand {ligand_sdf} --out {os.path.join(out_dir, 'docked.sdf.gz')}"
    subprocess.run(gnina_command, shell=True, check=True)
    return 

def decompress_file(file_path):
    decompressed_file_path = file_path[:-3]  # Remove '.gz' extension
    
    with gzip.open(file_path, 'rb') as compressed_file:
        with open(decompressed_file_path, 'wb') as decompressed_file:
            decompressed_file.write(compressed_file.read())
    
    return decompressed_file_path


def load_sdf(file_path):
    suppl = Chem.ForwardSDMolSupplier(file_path, removeHs=False)
    mols = [mol for mol in suppl if mol is not None]
    return mols


def extract_scores(mols):
    results = []
    for mol in mols:
        base_name = mol.GetProp('_Name')
        
        if mol.HasProp('CNNscore'):
            cnn_score = float(mol.GetProp('CNNscore'))
        else:
            cnn_score = None
        
        results.append({
            'base_name': base_name,
            'cnn_score': cnn_score,
            'mol': mol
        })
    
    return results


def rank_and_save_poses(results, output_dir):
    results.sort(key=lambda x: x['cnn_score'], reverse=True)
    
    for i, result in enumerate(results, start=1):
        mol = result['mol']
        output_name = f"rank{i}_score{result['cnn_score']:.2f}.sdf"
        output_path = os.path.join(output_dir, output_name)
        
        writer = Chem.SDWriter(output_path)
        writer.write(mol)
        writer.close()
    
    print(f"Ranked poses saved to directory: {output_dir}")


def main(posebusters_dir="posebusters_dataset", results_dir="forks/Vina/inference/GT_pocket_vina_posebusters_benchmark_outputs_2", top_n=10, method='gnina'):
    """
    Iterates through each target directory in 'posebusters_dataset'.
    For each, finds the protein.pdb and ligand.sdf, runs docking, and extracts top poses.
    """
    # Filter directories
    # protein_dir_pattern = re.compile(r'^[0-9]+[A-Z]+_[A-Z0-9]+$')

    protein_dirs = [
        d for d in os.listdir(posebusters_dir)
        if os.path.isdir(os.path.join(posebusters_dir, d))
    ]

    logger = get_custom_logger(f"{method}", f"{method}_timing.log")

    if not protein_dirs:
        print(f"[ERROR] No subdirectories found in {posebusters_dir}.")
        sys.exit(1)

    for protein_dir in protein_dirs:
        target_path = os.path.join(posebusters_dir, protein_dir)
        protein_pdb = os.path.join(target_path, f"{protein_dir}_protein.pdb")
        ligand_sdf  = os.path.join(target_path, f"{protein_dir}_ligand.sdf")

        if not (os.path.exists(protein_pdb) and os.path.exists(ligand_sdf)):
            print(f"[WARNING] Missing protein.pdb or ligand.sdf in {target_path}. Skipping.")
            continue

        print(f"\n[INFO] Processing target: {protein_dir}")
        out_dir = os.path.join(results_dir, protein_dir)
        os.makedirs(out_dir, exist_ok=True)
        
        start_time = time.time()
        try:
            if method == 'vina':
                vina_out = prepare_and_run_vina(protein_pdb, ligand_sdf, out_dir)
                # Now parse the top N conformations from the Vina output
                extract_and_write_top_poses(
                    vina_pdbqt_file=vina_out, 
                    out_dir=out_dir,
                    prefix=f"{protein_dir}",
                    remove_hs=True,
                    top_n=top_n
                )
            elif method == 'gnina':
                run_gnina_docking(protein_pdb, ligand_sdf, out_dir)
        except Exception as e:
            print(f"[ERROR] Docking failed for {protein_dir}. Reason: {str(e)}")
        elapsed_time = time.time() - start_time
        print(f"[INFO] Docking completed in {elapsed_time:.2f} seconds.")
        logger.info(f"{protein_dir},{elapsed_time:.2f}")


def run(config: dict):
    """
    Iterates through each target directory in 'posebusters_dataset'.
    For each, finds the protein.pdb and ligand.sdf, runs docking, and extracts top poses.
    """
    log_path = os.path.join(config['log_dir'], f"{config['method']}_timing_{config['repeat_index']}.log")
    logger = get_custom_logger(f"{config['method']}", config, f"{config['method']}_timing_{config['repeat_index']}.log") 
    inputs_df = pd.read_csv(config['inputs_csv'])

    for index, row in inputs_df.iterrows():
        protein_dir = row['complex_name']
        protein_pdb = row['protein_path']
        ligand_sdf  = row['ligand_path']

        if not (os.path.exists(protein_pdb) and os.path.exists(ligand_sdf)):
            print(f"[WARNING] Missing protein.pdb or ligand.sdf in {protein_dir}. Skipping.")
            continue

        print(f"\n[INFO] Processing target: {protein_dir}")
        out_dir = os.path.join(config['output_dir'], protein_dir)
        os.makedirs(out_dir, exist_ok=True)
        
        start_time = time.time()
        try:
            if config['method'] == 'vina':
                vina_out = prepare_and_run_vina(protein_pdb, ligand_sdf, out_dir)
                # Now parse the top N conformations from the Vina output
                extract_and_write_top_poses(
                    vina_pdbqt_file=vina_out, 
                    out_dir=out_dir,
                    prefix=f"{protein_dir}",
                    remove_hs=True,
                    top_n=config['top_n']
                )
            elif config['method'] == 'gnina':
                run_gnina_docking(protein_pdb, ligand_sdf, out_dir)
        except Exception as e:
            print(f"[ERROR] Docking failed for {protein_dir}. Reason: {str(e)}")
        elapsed_time = time.time() - start_time
        print(f"[INFO] Docking completed in {elapsed_time:.2f} seconds.")
        logger.info(f"{protein_dir},{elapsed_time:.2f}")


def parse_config_file_with_omegaconf(config_path):
    """
    Parse a configuration file using OmegaConf and resolve variables.
    
    Args:
        config_path (str): Path to the YAML configuration file
        
    Returns:
        dict: Resolved configuration as a dictionary
    """
    # Make sure the file exists
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    # Set up environment variables for resolution
    # In a real application, these should be set before calling this function
    os.environ["PROJECT_ROOT"] = "/home/aoxu/projects/PoseBench"  # Example value
    
    # Load the configuration from file
    cfg = OmegaConf.load(config_path)
    
    # Resolve interpolation
    resolved_cfg = OmegaConf.to_container(cfg, resolve=True)
    
    return resolved_cfg


if __name__ == "__main__":
    # Optionally pass the dataset directory or number of top poses from command line
    # Usage:
    #   python vina_posebusters.py /path/to/posebusters_dataset 10
    import yaml
    config = parse_config_file_with_omegaconf(sys.argv[1])
    run(config)
    exit()

    if len(sys.argv) >= 4:
        main(posebusters_dir=sys.argv[1], results_dir=sys.argv[2], method=sys.argv[3], top_n=int(sys.argv[4]))
    if len(sys.argv) >= 3:
        main(posebusters_dir=sys.argv[1], top_n=int(sys.argv[2]))
    elif len(sys.argv) == 2:
        main(posebusters_dir=sys.argv[1])
    
    else:
        base_dir = "forks/GNINA/inference/GT_pocket_gnina_posebusters_benchmark_outputs_1"
        for protein_dir in os.listdir(base_dir):
            print(f"processing protien {protein_dir}")
            output_directory = os.path.join(base_dir, protein_dir)
            sdf_file_path = os.path.join(output_directory, "docked.sdf.gz")
            if not os.path.exists(sdf_file_path):
                print(f"{protein_dir} does not exist")
                continue
            # Create the output directory if it doesn't exist
            os.makedirs(output_directory, exist_ok=True)

            # Decompress the SDF file
            if not os.path.exists(sdf_file_path.replace("docked.sdf.gz", "docked.sdf")):
                decompressed_file_path = decompress_file(sdf_file_path)

            # Load the compressed SDF file using RDKit
            mols = load_sdf(sdf_file_path.replace("docked.sdf.gz", "docked.sdf"))

            # Extract scores from the loaded molecules
            results = extract_scores(mols)

            # Rank the poses and save to a single SDF file
            rank_and_save_poses(results, output_directory)


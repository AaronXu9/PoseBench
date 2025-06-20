#!/usr/bin/env python3

import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdMolAlign
import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def debug_coordinate_transformation():
    """Debug the coordinate transformation issue step by step"""
    
    # Use the system we know works (1afb)
    system_id = '1afb__1__1.A__1.D_1.F'
    pdb_id = '1afb'
    
    logger.info(f"=== DEBUGGING COORDINATE TRANSFORMATION FOR {system_id} ===")
    
    approach = BoltzApproach()
    
    # 1. Get the prediction file
    protein_dir = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}"
    pred_files = approach.list_top_n_files(protein_dir, 1)
    
    if not pred_files:
        logger.error("No prediction files found")
        return
        
    pred_path = pred_files[0]
    logger.info(f"Prediction file: {pred_path}")
    
    # 2. Get the original PDB file (not the processed SDF)
    logger.info("=== STEP 1: ORIGINAL LIGAND COORDINATES ===")
    original_pdb_path = pred_path.replace('_ligand_protein_aligned.sdf', '.pdb')
    logger.info(f"Original PDB file: {original_pdb_path}")
    
    original_mol = Chem.MolFromPDBFile(original_pdb_path, removeHs=False)
    if original_mol is None:
        original_mol = Chem.MolFromPDBFile(original_pdb_path, removeHs=True)
    
    if original_mol:
        original_conf = original_mol.GetConformer()
        original_coords = []
        for i in range(original_mol.GetNumAtoms()):
            pos = original_conf.GetAtomPosition(i)
            original_coords.append([pos.x, pos.y, pos.z])
        original_coords = np.array(original_coords)
        logger.info(f"Original ligand atoms: {original_mol.GetNumAtoms()}")
        logger.info(f"Original ligand center: {np.mean(original_coords, axis=0)}")
        logger.info(f"Original ligand extent: {np.max(original_coords, axis=0) - np.min(original_coords, axis=0)}")
    
    # 3. Get reference ligand
    logger.info("=== STEP 2: REFERENCE LIGAND COORDINATES ===")
    ref_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{system_id}/{system_id}_ligand.sdf"
    ref_mol = Chem.SDMolSupplier(ref_ligand_path)[0]
    
    if ref_mol:
        ref_conf = ref_mol.GetConformer()
        ref_coords = []
        for i in range(ref_mol.GetNumAtoms()):
            pos = ref_conf.GetAtomPosition(i)
            ref_coords.append([pos.x, pos.y, pos.z])
        ref_coords = np.array(ref_coords)
        logger.info(f"Reference ligand atoms: {ref_mol.GetNumAtoms()}")
        logger.info(f"Reference ligand center: {np.mean(ref_coords, axis=0)}")
        logger.info(f"Reference ligand extent: {np.max(ref_coords, axis=0) - np.min(ref_coords, axis=0)}")
    
    # 4. Calculate RMSD before any alignment
    logger.info("=== STEP 3: RMSD BEFORE ALIGNMENT ===")
    try:
        raw_rmsd = rdMolAlign.AlignMol(original_mol, ref_mol)
        logger.info(f"Raw RMSD (before any alignment): {raw_rmsd:.3f} Å")
    except Exception as e:
        logger.error(f"Raw RMSD calculation failed: {e}")
    
    # 5. Test protein alignment
    logger.info("=== STEP 4: PROTEIN ALIGNMENT ===")
    ref_protein_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{system_id}/{system_id}_protein.pdb"
    
    if os.path.exists(ref_protein_path):
        rotation_matrix, translation_vector = approach._align_protein_structures(pred_path, ref_protein_path)
        
        if rotation_matrix is not None:
            logger.info(f"Protein alignment successful")
            logger.info(f"Rotation matrix shape: {rotation_matrix.shape}")
            logger.info(f"Translation vector: {translation_vector}")
            logger.info(f"Rotation matrix determinant: {np.linalg.det(rotation_matrix):.6f}")
        else:
            logger.error("Protein alignment failed")
            return
    else:
        logger.error(f"Reference protein not found: {ref_protein_path}")
        return
    
    # 6. Apply transformation manually and check step by step
    logger.info("=== STEP 5: MANUAL COORDINATE TRANSFORMATION ===")
    
    # Apply transformation to original coordinates
    transformed_coords = []
    for coord in original_coords:
        # Apply rotation and translation
        new_coord = np.dot(coord, rotation_matrix.T) + translation_vector
        transformed_coords.append(new_coord)
    transformed_coords = np.array(transformed_coords)
    
    logger.info(f"Transformed ligand center: {np.mean(transformed_coords, axis=0)}")
    logger.info(f"Center displacement: {np.linalg.norm(np.mean(transformed_coords, axis=0) - np.mean(original_coords, axis=0)):.3f} Å")
    
    # 7. Create a new molecule with transformed coordinates
    logger.info("=== STEP 6: CREATING TRANSFORMED MOLECULE ===")
    transformed_mol = Chem.Mol(original_mol)
    transformed_conf = transformed_mol.GetConformer()
    
    for i in range(transformed_mol.GetNumAtoms()):
        transformed_conf.SetAtomPosition(i, transformed_coords[i])
    
    # 8. Calculate RMSD after transformation
    logger.info("=== STEP 7: RMSD AFTER TRANSFORMATION ===")
    try:
        transformed_rmsd = rdMolAlign.AlignMol(transformed_mol, ref_mol)
        logger.info(f"RMSD after transformation: {transformed_rmsd:.3f} Å")
        
        # Also try without further alignment
        manual_rmsd = np.sqrt(np.mean(np.sum((transformed_coords[:len(ref_coords)] - ref_coords[:len(transformed_coords)])**2, axis=1)))
        logger.info(f"Manual RMSD (no further alignment): {manual_rmsd:.3f} Å")
        
    except Exception as e:
        logger.error(f"Transformed RMSD calculation failed: {e}")
    
    # 9. Compare with the approach's own transformation
    logger.info("=== STEP 8: APPROACH'S OWN TRANSFORMATION ===")
    approach_mol = approach._apply_transformation_to_ligand(original_mol, rotation_matrix, translation_vector)
    
    if approach_mol:
        try:
            approach_rmsd = rdMolAlign.AlignMol(approach_mol, ref_mol)
            logger.info(f"Approach's transformation RMSD: {approach_rmsd:.3f} Å")
        except Exception as e:
            logger.error(f"Approach RMSD calculation failed: {e}")
    
    # 10. Check if there's a coordinate system issue
    logger.info("=== STEP 9: COORDINATE SYSTEM ANALYSIS ===")
    logger.info(f"Original coordinates range:")
    logger.info(f"  X: {np.min(original_coords[:, 0]):.2f} to {np.max(original_coords[:, 0]):.2f}")
    logger.info(f"  Y: {np.min(original_coords[:, 1]):.2f} to {np.max(original_coords[:, 1]):.2f}")
    logger.info(f"  Z: {np.min(original_coords[:, 2]):.2f} to {np.max(original_coords[:, 2]):.2f}")
    
    logger.info(f"Reference coordinates range:")
    logger.info(f"  X: {np.min(ref_coords[:, 0]):.2f} to {np.max(ref_coords[:, 0]):.2f}")
    logger.info(f"  Y: {np.min(ref_coords[:, 1]):.2f} to {np.max(ref_coords[:, 1]):.2f}")
    logger.info(f"  Z: {np.min(ref_coords[:, 2]):.2f} to {np.max(ref_coords[:, 2]):.2f}")
    
    logger.info(f"Transformed coordinates range:")
    logger.info(f"  X: {np.min(transformed_coords[:, 0]):.2f} to {np.max(transformed_coords[:, 0]):.2f}")
    logger.info(f"  Y: {np.min(transformed_coords[:, 1]):.2f} to {np.max(transformed_coords[:, 1]):.2f}")
    logger.info(f"  Z: {np.min(transformed_coords[:, 2]):.2f} to {np.max(transformed_coords[:, 2]):.2f}")
    
    # 11. Test the existing aligned file
    logger.info("=== STEP 10: EXISTING ALIGNED FILE ===")
    aligned_sdf_path = pred_path.replace('.pdb', '_ligand_protein_aligned.sdf')
    if os.path.exists(aligned_sdf_path):
        try:
            aligned_mol = Chem.SDMolSupplier(aligned_sdf_path)[0]
            if aligned_mol:
                existing_rmsd = rdMolAlign.AlignMol(aligned_mol, ref_mol)
                logger.info(f"Existing aligned file RMSD: {existing_rmsd:.3f} Å")
                
                # Check its coordinates
                aligned_conf = aligned_mol.GetConformer()
                aligned_coords = []
                for i in range(aligned_mol.GetNumAtoms()):
                    pos = aligned_conf.GetAtomPosition(i)
                    aligned_coords.append([pos.x, pos.y, pos.z])
                aligned_coords = np.array(aligned_coords)
                logger.info(f"Existing aligned center: {np.mean(aligned_coords, axis=0)}")
            else:
                logger.error("Could not read existing aligned molecule")
        except Exception as e:
            logger.error(f"Error reading existing aligned file: {e}")
    else:
        logger.warning(f"Existing aligned file not found: {aligned_sdf_path}")

if __name__ == "__main__":
    debug_coordinate_transformation()

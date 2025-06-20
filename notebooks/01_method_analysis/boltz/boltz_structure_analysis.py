"""
Boltz Structure Analysis - Protein and Ligand RMSD Computation
==============================================================

This module provides functions to compute RMSD for both protein and ligand structures
to analyze how protein structure accuracy affects ligand prediction accuracy.
"""

import sys
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks')

import os
import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, List
from rdkit import Chem
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
from tqdm import tqdm
import glob
from Approach import BoltzApproach
from scipy.spatial.transform import Rotation
from rdkit.Chem import rdMolAlign
import json
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


def compute_ligand_rmsd(pred_ligand_path: str, ref_ligand_path: str) -> float:
    """
    Compute RMSD between predicted and reference ligand structures.
    
    Args:
        pred_ligand_path: Path to predicted ligand SDF file
        ref_ligand_path: Path to reference ligand SDF file
        
    Returns:
        RMSD value in Angstroms, or NaN if calculation fails
    """
    try:
        # Load molecules
        pred_supplier = Chem.SDMolSupplier(pred_ligand_path, sanitize=False, removeHs=False)
        ref_supplier = Chem.SDMolSupplier(ref_ligand_path, sanitize=False, removeHs=False)
        
        pred_mol = next(pred_supplier)
        ref_mol = next(ref_supplier)
        
        if not pred_mol or not ref_mol:
            return float('nan')
            
        if pred_mol.GetNumAtoms() != ref_mol.GetNumAtoms():
            print(f"[WARNING] Atom count mismatch: pred={pred_mol.GetNumAtoms()}, ref={ref_mol.GetNumAtoms()}")
            return float('nan')
            
        # Get coordinates
        pred_conf = pred_mol.GetConformer()
        ref_conf = ref_mol.GetConformer()
        
        pred_coords = []
        ref_coords = []
        
        for i in range(pred_mol.GetNumAtoms()):
            pred_pos = pred_conf.GetAtomPosition(i)
            ref_pos = ref_conf.GetAtomPosition(i)
            pred_coords.append([pred_pos.x, pred_pos.y, pred_pos.z])
            ref_coords.append([ref_pos.x, ref_pos.y, ref_pos.z])
        
        pred_coords = np.array(pred_coords)
        ref_coords = np.array(ref_coords)
        
        # Calculate RMSD
        diff = pred_coords - ref_coords
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        return rmsd
        
    except Exception as e:
        print(f"[ERROR] Ligand RMSD calculation failed: {str(e)}")
        return float('nan')


def compute_protein_rmsd(pred_protein_path: str, ref_protein_path: str, 
                        align_method: str = 'ca_only') -> Tuple[float, float]:
    """
    Compute RMSD between predicted and reference protein structures.
    
    Args:
        pred_protein_path: Path to predicted protein PDB file
        ref_protein_path: Path to reference protein PDB file
        align_method: 'ca_only' for CA atoms only, 'backbone' for N,CA,C,O atoms
        
    Returns:
        Tuple of (aligned_rmsd, unaligned_rmsd) in Angstroms
    """
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            
            parser = PDBParser(QUIET=True)
            pred_structure = parser.get_structure("predicted", pred_protein_path)
            ref_structure = parser.get_structure("reference", ref_protein_path)
            
            pred_model = pred_structure[0]
            ref_model = ref_structure[0]
            
            # Extract atoms for alignment
            if align_method == 'ca_only':
                atom_names = ['CA']
            elif align_method == 'backbone':
                atom_names = ['N', 'CA', 'C', 'O']
            else:
                raise ValueError(f"Unknown align_method: {align_method}")
            
            pred_atoms = []
            ref_atoms = []
            
            # Get corresponding atoms from both structures
            for pred_chain in pred_model:
                for pred_residue in pred_chain:
                    if pred_residue.id[0] == ' ':  # Only standard residues
                        for atom_name in atom_names:
                            if atom_name in pred_residue:
                                pred_atom = pred_residue[atom_name]
                                
                                # Find corresponding atom in reference
                                try:
                                    ref_chain = ref_model[pred_chain.id]
                                    ref_residue = ref_chain[pred_residue.id]
                                    if atom_name in ref_residue:
                                        ref_atom = ref_residue[atom_name]
                                        pred_atoms.append(pred_atom)
                                        ref_atoms.append(ref_atom)
                                except KeyError:
                                    continue  # Residue not found in reference
            
            if len(pred_atoms) < 3:
                print(f"[WARNING] Too few corresponding atoms found: {len(pred_atoms)}")
                return float('nan'), float('nan')
            
            # Calculate unaligned RMSD
            pred_coords = np.array([atom.coord for atom in pred_atoms])
            ref_coords = np.array([atom.coord for atom in ref_atoms])
            
            diff_unaligned = pred_coords - ref_coords
            unaligned_rmsd = np.sqrt(np.mean(np.sum(diff_unaligned**2, axis=1)))
            
            # Calculate aligned RMSD using Superimposer
            superimposer = Superimposer()
            superimposer.set_atoms(ref_atoms, pred_atoms)
            superimposer.apply(pred_atoms)
            
            aligned_rmsd = superimposer.rms
            
            return aligned_rmsd, unaligned_rmsd
            
    except Exception as e:
        print(f"[ERROR] Protein RMSD calculation failed: {str(e)}")
        return float('nan'), float('nan')


def compute_binding_site_rmsd(pred_protein_path: str, ref_protein_path: str, 
                             ref_ligand_path: str, cutoff: float = 10.0) -> Tuple[float, float]:
    """
    Compute RMSD for binding site residues only (residues within cutoff of ligand).
    
    Args:
        pred_protein_path: Path to predicted protein PDB file
        ref_protein_path: Path to reference protein PDB file
        ref_ligand_path: Path to reference ligand SDF file
        cutoff: Distance cutoff in Angstroms to define binding site
        
    Returns:
        Tuple of (aligned_rmsd, unaligned_rmsd) for binding site residues
    """
    try:
        # Get ligand coordinates
        ref_supplier = Chem.SDMolSupplier(ref_ligand_path, sanitize=False, removeHs=False)
        ref_mol = next(ref_supplier)
        
        if not ref_mol:
            return float('nan'), float('nan')
            
        ref_conf = ref_mol.GetConformer()
        ligand_coords = []
        for i in range(ref_mol.GetNumAtoms()):
            pos = ref_conf.GetAtomPosition(i)
            ligand_coords.append([pos.x, pos.y, pos.z])
        ligand_coords = np.array(ligand_coords)
        
        # Parse protein structures
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            
            parser = PDBParser(QUIET=True)
            pred_structure = parser.get_structure("predicted", pred_protein_path)
            ref_structure = parser.get_structure("reference", ref_protein_path)
            
            pred_model = pred_structure[0]
            ref_model = ref_structure[0]
        
        # Find binding site residues in reference structure
        binding_site_residues = set()
        
        for ref_chain in ref_model:
            for ref_residue in ref_chain:
                if ref_residue.id[0] == ' ':  # Only standard residues
                    for atom in ref_residue:
                        atom_coord = np.array(atom.coord)
                        distances = np.linalg.norm(ligand_coords - atom_coord, axis=1)
                        if np.min(distances) <= cutoff:
                            binding_site_residues.add((ref_chain.id, ref_residue.id))
                            break
        
        print(f"[INFO] Found {len(binding_site_residues)} binding site residues within {cutoff} Å")
        
        # Extract corresponding atoms from binding site residues
        pred_atoms = []
        ref_atoms = []
        
        for chain_id, residue_id in binding_site_residues:
            try:
                pred_chain = pred_model[chain_id]
                pred_residue = pred_chain[residue_id]
                ref_chain = ref_model[chain_id]
                ref_residue = ref_chain[residue_id]
                
                # Add CA atoms
                if 'CA' in pred_residue and 'CA' in ref_residue:
                    pred_atoms.append(pred_residue['CA'])
                    ref_atoms.append(ref_residue['CA'])
                    
            except KeyError:
                continue  # Residue not found in predicted structure
        
        if len(pred_atoms) < 3:
            print(f"[WARNING] Too few binding site atoms found: {len(pred_atoms)}")
            return float('nan'), float('nan')
        
        # Calculate RMSDs
        pred_coords = np.array([atom.coord for atom in pred_atoms])
        ref_coords = np.array([atom.coord for atom in ref_atoms])
        
        diff_unaligned = pred_coords - ref_coords
        unaligned_rmsd = np.sqrt(np.mean(np.sum(diff_unaligned**2, axis=1)))
        
        # Aligned RMSD
        superimposer = Superimposer()
        superimposer.set_atoms(ref_atoms, pred_atoms)
        superimposer.apply(pred_atoms)
        
        aligned_rmsd = superimposer.rms
        
        return aligned_rmsd, unaligned_rmsd
        
    except Exception as e:
        print(f"[ERROR] Binding site RMSD calculation failed: {str(e)}")
        return float('nan'), float('nan')


def analyze_boltz_structure_accuracy(base_outdir: str, data_dir: str, 
                                   output_csv: str, max_proteins: Optional[int] = None) -> pd.DataFrame:
    """
    Comprehensive analysis of Boltz structure prediction accuracy.
    
    Args:
        base_outdir: Base directory containing Boltz results
        data_dir: Directory containing reference structures
        output_csv: Path to save results CSV
        max_proteins: Maximum number of proteins to analyze (None for all)
        
    Returns:
        DataFrame with comprehensive RMSD analysis
    """
    print("BOLTZ STRUCTURE ACCURACY ANALYSIS")
    print("=================================")
    
    approach = BoltzApproach()
    results = []
    
    # Get list of protein directories
    protein_dirs = []
    for item in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir, item)):
            protein_dirs.append(item)
    
    if max_proteins:
        protein_dirs = protein_dirs[:max_proteins]
    
    print(f"Analyzing {len(protein_dirs)} protein systems...")
    
    for protein_name in tqdm(protein_dirs, desc="Processing proteins"):
        try:
            # Find Boltz results directory
            boltz_pattern = f"{base_outdir}/*{protein_name}*"
            boltz_matches = glob.glob(boltz_pattern)
            
            if not boltz_matches:
                print(f"[WARNING] No Boltz results found for {protein_name}")
                continue
                
            protein_dir = boltz_matches[0]
            
            # Get reference files
            ref_protein_path = os.path.join(data_dir, protein_name, f"{protein_name}_protein.pdb")
            ref_ligand_path = os.path.join(data_dir, protein_name, f"{protein_name}_ligand.sdf")
            
            if not (os.path.exists(ref_protein_path) and os.path.exists(ref_ligand_path)):
                print(f"[WARNING] Missing reference files for {protein_name}")
                continue
            
            # Get Boltz predictions
            pred_dir = os.path.join(protein_dir, "predictions", protein_name)
            if not os.path.exists(pred_dir):
                print(f"[WARNING] Predictions directory not found for {protein_name}")
                continue
            
            # Analyze top 5 models
            for model_num in range(5):
                pred_protein_path = os.path.join(pred_dir, f"{protein_name}_model_{model_num}.pdb")
                
                if not os.path.exists(pred_protein_path):
                    continue
                
                # Get aligned ligand
                aligned_files = approach.list_top_n_files(protein_dir, 5)
                pred_ligand_path = None
                
                for aligned_file in aligned_files:
                    if f"model_{model_num}_ligand" in aligned_file:
                        pred_ligand_path = aligned_file
                        break
                
                if not pred_ligand_path or not os.path.exists(pred_ligand_path):
                    print(f"[WARNING] No aligned ligand found for {protein_name} model {model_num}")
                    continue
                
                # Compute all RMSDs
                ligand_rmsd = compute_ligand_rmsd(pred_ligand_path, ref_ligand_path)
                
                protein_ca_aligned, protein_ca_unaligned = compute_protein_rmsd(
                    pred_protein_path, ref_protein_path, 'ca_only'
                )
                
                protein_bb_aligned, protein_bb_unaligned = compute_protein_rmsd(
                    pred_protein_path, ref_protein_path, 'backbone'
                )
                
                binding_site_aligned, binding_site_unaligned = compute_binding_site_rmsd(
                    pred_protein_path, ref_protein_path, ref_ligand_path
                )
                
                # Get confidence score
                confidence = approach.parse_score(pred_ligand_path)
                
                # Store results
                result = {
                    'protein': protein_name,
                    'model': model_num,
                    'confidence': confidence,
                    'ligand_rmsd': ligand_rmsd,
                    'protein_ca_aligned_rmsd': protein_ca_aligned,
                    'protein_ca_unaligned_rmsd': protein_ca_unaligned,
                    'protein_backbone_aligned_rmsd': protein_bb_aligned,
                    'protein_backbone_unaligned_rmsd': protein_bb_unaligned,
                    'binding_site_aligned_rmsd': binding_site_aligned,
                    'binding_site_unaligned_rmsd': binding_site_unaligned,
                }
                
                results.append(result)
                
                print(f"[INFO] {protein_name} model {model_num}: "
                      f"Ligand={ligand_rmsd:.2f}Å, Protein_CA={protein_ca_aligned:.2f}Å, "
                      f"BindingSite={binding_site_aligned:.2f}Å")
                
        except Exception as e:
            print(f"[ERROR] Failed to process {protein_name}: {str(e)}")
            continue
    
    # Create DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    
    print(f"\nAnalysis complete! Results saved to {output_csv}")
    print(f"Analyzed {len(df)} protein-model combinations")
    
    # Print summary statistics
    if len(df) > 0:
        print("\nSUMMARY STATISTICS:")
        print("==================")
        print(f"Ligand RMSD: {df['ligand_rmsd'].mean():.2f} ± {df['ligand_rmsd'].std():.2f} Å")
        print(f"Protein CA RMSD: {df['protein_ca_aligned_rmsd'].mean():.2f} ± {df['protein_ca_aligned_rmsd'].std():.2f} Å")
        print(f"Binding Site RMSD: {df['binding_site_aligned_rmsd'].mean():.2f} ± {df['binding_site_aligned_rmsd'].std():.2f} Å")
        
        # Correlation analysis
        corr_protein_ligand = df['protein_ca_aligned_rmsd'].corr(df['ligand_rmsd'])
        corr_binding_ligand = df['binding_site_aligned_rmsd'].corr(df['ligand_rmsd'])
        
        print(f"\nCORRELATION ANALYSIS:")
        print(f"Protein CA RMSD vs Ligand RMSD: {corr_protein_ligand:.3f}")
        print(f"Binding Site RMSD vs Ligand RMSD: {corr_binding_ligand:.3f}")
    
    return df


def compute_protein_rmsd_comprehensive(pred_protein_path: str, ref_protein_path: str) -> Dict[str, float]:
    """
    Comprehensive protein RMSD computation with multiple metrics.
    
    Args:
        pred_protein_path: Path to predicted protein PDB file
        ref_protein_path: Path to reference protein PDB file
        
    Returns:
        Dictionary with various RMSD metrics
    """
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            
            parser = PDBParser(QUIET=True)
            pred_structure = parser.get_structure("predicted", pred_protein_path)
            ref_structure = parser.get_structure("reference", ref_protein_path)
            
            pred_model = pred_structure[0]
            ref_model = ref_structure[0]
            
            # Extract different atom types
            pred_ca_atoms, ref_ca_atoms = extract_ca_atoms(pred_model, ref_model)
            pred_bb_atoms, ref_bb_atoms = extract_backbone_atoms(pred_model, ref_model)
            pred_all_atoms, ref_all_atoms = extract_all_atoms(pred_model, ref_model)
            
            results = {}
            
            # CA-only RMSD
            if pred_ca_atoms and ref_ca_atoms and len(pred_ca_atoms) == len(ref_ca_atoms):
                pred_ca_array = np.array(pred_ca_atoms)
                ref_ca_array = np.array(ref_ca_atoms)
                aligned_rmsd, unaligned_rmsd = compute_rmsd_with_alignment(pred_ca_array, ref_ca_array)
                results['ca_aligned_rmsd'] = aligned_rmsd
                results['ca_unaligned_rmsd'] = unaligned_rmsd
            else:
                results['ca_aligned_rmsd'] = float('nan')
                results['ca_unaligned_rmsd'] = float('nan')
            
            # Backbone RMSD
            if pred_bb_atoms and ref_bb_atoms and len(pred_bb_atoms) == len(ref_bb_atoms):
                pred_bb_array = np.array(pred_bb_atoms)
                ref_bb_array = np.array(ref_bb_atoms)
                aligned_rmsd, unaligned_rmsd = compute_rmsd_with_alignment(pred_bb_array, ref_bb_array)
                results['backbone_aligned_rmsd'] = aligned_rmsd
                results['backbone_unaligned_rmsd'] = unaligned_rmsd
            else:
                results['backbone_aligned_rmsd'] = float('nan')
                results['backbone_unaligned_rmsd'] = float('nan')
            
            # All-atom RMSD (if reasonable number of atoms)
            if (pred_all_atoms and ref_all_atoms and 
                len(pred_all_atoms) == len(ref_all_atoms) and 
                len(pred_all_atoms) < 10000):  # Limit for computational efficiency
                pred_all_array = np.array(pred_all_atoms)
                ref_all_array = np.array(ref_all_atoms)
                aligned_rmsd, unaligned_rmsd = compute_rmsd_with_alignment(pred_all_array, ref_all_array)
                results['all_atom_aligned_rmsd'] = aligned_rmsd
                results['all_atom_unaligned_rmsd'] = unaligned_rmsd
            else:
                results['all_atom_aligned_rmsd'] = float('nan')
                results['all_atom_unaligned_rmsd'] = float('nan')
            
            return results
            
    except Exception as e:
        print(f"[ERROR] Protein RMSD calculation failed: {str(e)}")
        return {
            'ca_aligned_rmsd': float('nan'),
            'ca_unaligned_rmsd': float('nan'),
            'backbone_aligned_rmsd': float('nan'),
            'backbone_unaligned_rmsd': float('nan'),
            'all_atom_aligned_rmsd': float('nan'),
            'all_atom_unaligned_rmsd': float('nan')
        }


def compute_ligand_rmsd_comprehensive(pred_ligand_path: str, ref_ligand_path: str) -> Dict[str, float]:
    """
    Comprehensive ligand RMSD computation with multiple methods.
    
    Args:
        pred_ligand_path: Path to predicted ligand SDF file
        ref_ligand_path: Path to reference ligand SDF file
        
    Returns:
        Dictionary with various RMSD metrics
    """
    try:
        # Load molecules using robust reading
        pred_mol = read_molecule_robust(pred_ligand_path)
        ref_mol = read_molecule_robust(ref_ligand_path)
        
        if not pred_mol or not ref_mol:
            return {
                'rdkit_rmsd': float('nan'),
                'coordinate_rmsd': float('nan'),
                'aligned_rmsd': float('nan'),
                'atom_count_match': False
            }
        
        results = {
            'atom_count_match': pred_mol.GetNumAtoms() == ref_mol.GetNumAtoms()
        }
        
        # RDKit-based RMSD (if atom counts match)
        if results['atom_count_match']:
            try:
                rmsd = rdMolAlign.GetBestRMS(pred_mol, ref_mol)
                results['rdkit_rmsd'] = rmsd
            except Exception as e:
                print(f"[WARNING] RDKit RMSD failed: {e}")
                results['rdkit_rmsd'] = float('nan')
        else:
            results['rdkit_rmsd'] = float('nan')
        
        # Coordinate-based RMSD
        try:
            pred_coords = pred_mol.GetConformer().GetPositions()
            ref_coords = ref_mol.GetConformer().GetPositions()
            
            if pred_coords.shape == ref_coords.shape:
                # Simple coordinate RMSD
                diff = pred_coords - ref_coords
                results['coordinate_rmsd'] = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
                
                # Aligned coordinate RMSD
                aligned_rmsd, _ = compute_rmsd_with_alignment(pred_coords, ref_coords)
                results['aligned_rmsd'] = aligned_rmsd
            else:
                results['coordinate_rmsd'] = float('nan')
                results['aligned_rmsd'] = float('nan')
                
        except Exception as e:
            print(f"[WARNING] Coordinate RMSD failed: {e}")
            results['coordinate_rmsd'] = float('nan')
            results['aligned_rmsd'] = float('nan')
        
        return results
        
    except Exception as e:
        print(f"[ERROR] Ligand RMSD calculation failed: {str(e)}")
        return {
            'rdkit_rmsd': float('nan'),
            'coordinate_rmsd': float('nan'),
            'aligned_rmsd': float('nan'),
            'atom_count_match': False
        }


def compute_binding_site_rmsd(pred_protein_path: str, ref_protein_path: str, 
                             ref_ligand_path: str, distance_cutoff: float = 10.0) -> Dict[str, float]:
    """
    Compute RMSD for binding site residues around the reference ligand.
    
    Args:
        pred_protein_path: Path to predicted protein PDB file
        ref_protein_path: Path to reference protein PDB file
        ref_ligand_path: Path to reference ligand SDF file
        distance_cutoff: Distance cutoff in Angstroms to define binding site
        
    Returns:
        Dictionary with binding site RMSD metrics
    """
    try:
        # Load reference ligand to define binding site
        ref_mol = read_molecule_robust(ref_ligand_path)
        if not ref_mol:
            return {'binding_site_aligned_rmsd': float('nan'), 'binding_site_unaligned_rmsd': float('nan')}
        
        ref_ligand_coords = ref_mol.GetConformer().GetPositions()
        ligand_center = np.mean(ref_ligand_coords, axis=0)
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            
            parser = PDBParser(QUIET=True)
            pred_structure = parser.get_structure("predicted", pred_protein_path)
            ref_structure = parser.get_structure("reference", ref_protein_path)
            
            pred_model = pred_structure[0]
            ref_model = ref_structure[0]
            
            # Extract binding site atoms
            pred_bs_atoms = extract_binding_site_atoms(pred_model, ligand_center, distance_cutoff)
            ref_bs_atoms = extract_binding_site_atoms(ref_model, ligand_center, distance_cutoff)
            
            if pred_bs_atoms and ref_bs_atoms and len(pred_bs_atoms) == len(ref_bs_atoms):
                pred_bs_array = np.array(pred_bs_atoms)
                ref_bs_array = np.array(ref_bs_atoms)
                aligned_rmsd, unaligned_rmsd = compute_rmsd_with_alignment(pred_bs_array, ref_bs_array)
                return {
                    'binding_site_aligned_rmsd': aligned_rmsd,
                    'binding_site_unaligned_rmsd': unaligned_rmsd,
                    'binding_site_residue_count': len(pred_bs_atoms)
                }
            else:
                return {
                    'binding_site_aligned_rmsd': float('nan'),
                    'binding_site_unaligned_rmsd': float('nan'),
                    'binding_site_residue_count': 0
                }
                
    except Exception as e:
        print(f"[ERROR] Binding site RMSD calculation failed: {str(e)}")
        return {
            'binding_site_aligned_rmsd': float('nan'),
            'binding_site_unaligned_rmsd': float('nan'),
            'binding_site_residue_count': 0
        }


def compute_rmsd_with_alignment(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[float, float]:
    """
    Compute RMSD between two coordinate arrays with and without optimal alignment.
    
    Args:
        coords1: First set of coordinates (N, 3)
        coords2: Second set of coordinates (N, 3)
        
    Returns:
        Tuple of (aligned_rmsd, unaligned_rmsd)
    """
    try:
        coords1 = np.array(coords1)
        coords2 = np.array(coords2)
        
        # Ensure coordinates are 2D arrays
        if coords1.ndim == 1:
            coords1 = coords1.reshape(-1, 3)
        if coords2.ndim == 1:
            coords2 = coords2.reshape(-1, 3)
        
        if coords1.shape != coords2.shape:
            print(f"[WARNING] Shape mismatch: {coords1.shape} vs {coords2.shape}")
            return float('nan'), float('nan')
            
        if coords1.shape[1] != 3:
            print(f"[WARNING] Expected 3D coordinates, got shape {coords1.shape}")
            return float('nan'), float('nan')
        
        # Unaligned RMSD
        diff = coords1 - coords2
        unaligned_rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        # Aligned RMSD using optimal rotation
        # Center both coordinate sets
        centroid1 = np.mean(coords1, axis=0)
        centroid2 = np.mean(coords2, axis=0)
        
        centered1 = coords1 - centroid1
        centered2 = coords2 - centroid2
        
        # Find optimal rotation using Kabsch algorithm
        try:
            result = Rotation.align_vectors(centered2, centered1)
            rotation = result[0]  # First element is the rotation
            aligned1 = rotation.apply(centered1)
            aligned_rmsd = np.sqrt(np.mean(np.sum((aligned1 - centered2)**2, axis=1)))
        except Exception:
            # Fallback: use unaligned RMSD if rotation fails
            aligned_rmsd = unaligned_rmsd
        
        return aligned_rmsd, unaligned_rmsd
        
    except Exception as e:
        print(f"[ERROR] RMSD computation failed: {str(e)}")
        return float('nan'), float('nan')


def read_molecule_robust(molecule_path: str):
    """
    Robust molecule reading with multiple fallback strategies.
    """
    try:
        if molecule_path.endswith('.sdf'):
            supplier = Chem.SDMolSupplier(molecule_path, sanitize=False, removeHs=False)
            for mol in supplier:
                if mol is not None:
                    return mol
        elif molecule_path.endswith('.pdb'):
            mol = Chem.MolFromPDBFile(molecule_path, sanitize=False, removeHs=False)
            if mol is not None:
                return mol
            
            # Fallback: try reading as PDB block
            with open(molecule_path, 'r') as f:
                pdb_block = f.read()
            mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
            return mol
            
    except Exception as e:
        print(f"[WARNING] Failed to read molecule {molecule_path}: {e}")
    
    return None


def extract_binding_site_atoms(model, ligand_center: np.ndarray, cutoff: float) -> List[np.ndarray]:
    """
    Extract coordinates of atoms within cutoff distance of ligand center.
    """
    binding_site_coords = []
    
    try:
        for chain in model:
            for residue in chain:
                # Skip water and other hetero residues
                if residue.get_resname() == "HOH" or len(residue.get_id()[0]) > 1:
                    continue
                
                # Check if any atom in residue is within cutoff
                residue_in_binding_site = False
                for atom in residue:
                    atom_coord = np.array(atom.get_vector())
                    distance = np.linalg.norm(atom_coord - ligand_center)
                    if distance <= cutoff:
                        residue_in_binding_site = True
                        break
                
                # If residue is in binding site, add CA coordinate
                if residue_in_binding_site and 'CA' in residue:
                    ca_coord = np.array(residue['CA'].get_vector())
                    binding_site_coords.append(ca_coord)
                    
    except Exception as e:
        print(f"[WARNING] Failed to extract binding site atoms: {e}")
    
    return binding_site_coords


def extract_ca_atoms(pred_model, ref_model) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Extract CA atoms from both models."""
    pred_ca_coords = []
    ref_ca_coords = []
    
    try:
        pred_residues = list(pred_model.get_residues())
        ref_residues = list(ref_model.get_residues())
        
        # Filter valid residues
        pred_valid = [res for res in pred_residues if 'CA' in res and res.get_resname() != "HOH" and len(res.get_id()[0]) == 1]
        ref_valid = [res for res in ref_residues if 'CA' in res and res.get_resname() != "HOH" and len(res.get_id()[0]) == 1]
        
        # Take minimum length to ensure matching
        min_len = min(len(pred_valid), len(ref_valid))
        # print(f"[DEBUG] CA atoms: pred={len(pred_valid)}, ref={len(ref_valid)}, min={min_len}")
        
        for i in range(min_len):
            pred_ca_coords.append(pred_valid[i]['CA'].get_coord())  # get_coord() returns [x,y,z]
            ref_ca_coords.append(ref_valid[i]['CA'].get_coord())
            
    except Exception as e:
        print(f"[WARNING] Failed to extract CA atoms: {e}")
    
    return pred_ca_coords, ref_ca_coords


def extract_backbone_atoms(pred_model, ref_model) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Extract backbone atoms (N, CA, C, O) from both models."""
    pred_bb_coords = []
    ref_bb_coords = []
    
    try:
        pred_residues = list(pred_model.get_residues())
        ref_residues = list(ref_model.get_residues())
        
        # Filter valid residues
        pred_valid = [res for res in pred_residues if res.get_resname() != "HOH" and len(res.get_id()[0]) == 1]
        ref_valid = [res for res in ref_residues if res.get_resname() != "HOH" and len(res.get_id()[0]) == 1]
        
        # Take minimum length to ensure matching
        min_len = min(len(pred_valid), len(ref_valid))
        
        for i in range(min_len):
            pred_res = pred_valid[i]
            ref_res = ref_valid[i]
            
            # Extract backbone atoms if they exist
            bb_atoms = ['N', 'CA', 'C', 'O']
            for atom_name in bb_atoms:
                if atom_name in pred_res and atom_name in ref_res:
                    pred_bb_coords.append(pred_res[atom_name].get_coord())
                    ref_bb_coords.append(ref_res[atom_name].get_coord())
                    
    except Exception as e:
        print(f"[WARNING] Failed to extract backbone atoms: {e}")
    
    return pred_bb_coords, ref_bb_coords


def extract_all_atoms(pred_model, ref_model) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Extract all heavy atoms from both models."""
    pred_coords = []
    ref_coords = []
    
    try:
        pred_atoms = list(pred_model.get_atoms())
        ref_atoms = list(ref_model.get_atoms())
        
        # Filter heavy atoms (non-hydrogen)
        pred_heavy = [atom for atom in pred_atoms if atom.element != 'H']
        ref_heavy = [atom for atom in ref_atoms if atom.element != 'H']
        
        # Take minimum length to ensure matching
        min_len = min(len(pred_heavy), len(ref_heavy))
        
        for i in range(min_len):
            pred_coords.append(pred_heavy[i].get_coord())
            ref_coords.append(ref_heavy[i].get_coord())
            
    except Exception as e:
        print(f"[WARNING] Failed to extract all atoms: {e}")
    
    return pred_coords, ref_coords


def analyze_protein_ligand_dependency(results_df: pd.DataFrame) -> Dict:
    """
    Analyze the dependency between protein structure accuracy and ligand placement accuracy.
    
    Args:
        results_df: DataFrame with protein and ligand RMSD values
        
    Returns:
        Dictionary with correlation analysis results
    """
    try:
        # Filter valid data
        valid_data = results_df.dropna(subset=['protein_ca_aligned_rmsd', 'ligand_rmsd'])
        
        if len(valid_data) < 3:
            print("[WARNING] Insufficient valid data for correlation analysis")
            return {}
        
        # Compute correlations
        correlations = {}
        
        # Protein-Ligand correlations
        protein_metrics = ['protein_ca_aligned_rmsd', 'protein_backbone_aligned_rmsd', 
                          'binding_site_aligned_rmsd']
        ligand_metrics = ['ligand_rmsd']
        
        for prot_metric in protein_metrics:
            for lig_metric in ligand_metrics:
                if prot_metric in valid_data.columns and lig_metric in valid_data.columns:
                    valid_subset = valid_data.dropna(subset=[prot_metric, lig_metric])
                    if len(valid_subset) >= 3:
                        pearson_corr, pearson_p = pearsonr(valid_subset[prot_metric], valid_subset[lig_metric])
                        spearman_corr, spearman_p = spearmanr(valid_subset[prot_metric], valid_subset[lig_metric])
                        
                        correlations[f"{prot_metric}_vs_{lig_metric}"] = {
                            'pearson_r': pearson_corr,
                            'pearson_p': pearson_p,
                            'spearman_r': spearman_corr,
                            'spearman_p': spearman_p,
                            'n_samples': len(valid_subset)
                        }
        
        # Summary statistics
        summary_stats = {}
        for col in valid_data.select_dtypes(include=[np.number]).columns:
            summary_stats[col] = {
                'mean': valid_data[col].mean(),
                'std': valid_data[col].std(),
                'median': valid_data[col].median(),
                'min': valid_data[col].min(),
                'max': valid_data[col].max(),
                'count': valid_data[col].count()
            }
        
        return {
            'correlations': correlations,
            'summary_stats': summary_stats,
            'total_samples': len(results_df),
            'valid_samples': len(valid_data)
        }
        
    except Exception as e:
        print(f"[ERROR] Dependency analysis failed: {str(e)}")
        return {}


def create_dependency_plots(results_df: pd.DataFrame, output_dir: str = "/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks"):
    """
    Create plots to visualize protein-ligand accuracy dependencies.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Filter valid data
        valid_data = results_df.dropna(subset=['protein_ca_aligned_rmsd', 'ligand_rmsd'])
        
        if len(valid_data) < 3:
            print("[WARNING] Insufficient data for plotting")
            return
        
        # Create correlation plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Protein CA RMSD vs Ligand RMSD
        axes[0, 0].scatter(valid_data['protein_ca_aligned_rmsd'], valid_data['ligand_rmsd'], alpha=0.6)
        axes[0, 0].set_xlabel('Protein CA RMSD (Å)')
        axes[0, 0].set_ylabel('Ligand RMSD (Å)')
        axes[0, 0].set_title('Protein CA RMSD vs Ligand RMSD')
        
        # Binding Site RMSD vs Ligand RMSD
        if 'binding_site_aligned_rmsd' in valid_data.columns:
            bs_valid = valid_data.dropna(subset=['binding_site_aligned_rmsd'])
            if len(bs_valid) > 0:
                axes[0, 1].scatter(bs_valid['binding_site_aligned_rmsd'], bs_valid['ligand_rmsd'], alpha=0.6)
                axes[0, 1].set_xlabel('Binding Site RMSD (Å)')
                axes[0, 1].set_ylabel('Ligand RMSD (Å)')
                axes[0, 1].set_title('Binding Site RMSD vs Ligand RMSD')
        
        # Confidence vs metrics
        if 'confidence' in valid_data.columns:
            conf_valid = valid_data.dropna(subset=['confidence'])
            if len(conf_valid) > 0:
                axes[1, 0].scatter(conf_valid['confidence'], conf_valid['ligand_rmsd'], alpha=0.6)
                axes[1, 0].set_xlabel('Confidence Score')
                axes[1, 0].set_ylabel('Ligand RMSD (Å)')
                axes[1, 0].set_title('Confidence vs Ligand RMSD')
                
                axes[1, 1].scatter(conf_valid['confidence'], conf_valid['protein_ca_aligned_rmsd'], alpha=0.6)
                axes[1, 1].set_xlabel('Confidence Score')
                axes[1, 1].set_ylabel('Protein CA RMSD (Å)')
                axes[1, 1].set_title('Confidence vs Protein RMSD')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'boltz_dependency_analysis.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"[INFO] Dependency plots saved to {output_dir}")
        
    except Exception as e:
        print(f"[ERROR] Plot creation failed: {str(e)}")


def run_comprehensive_boltz_analysis():
    """
    Run comprehensive Boltz analysis with protein-ligand dependency analysis.
    """
    # Configuration
    base_outdir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0"
    data_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set"
    output_csv = "/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks/boltz_full_analysis_results.csv"
    
    print("[INFO] Starting comprehensive Boltz structure analysis...")
    print(f"[INFO] Prediction directory: {base_outdir}")
    print(f"[INFO] Reference data directory: {data_dir}")
    
    # Get list of proteins
    protein_dirs = sorted([d for d in os.listdir(base_outdir) if os.path.isdir(os.path.join(base_outdir, d))])
    print(f"[INFO] Found {len(protein_dirs)} protein directories")
    
    # Run analysis on a subset for testing
    test_proteins = protein_dirs  # Analyze ALL proteins for comprehensive results
    print(f"[INFO] Testing with {len(test_proteins)} proteins: {test_proteins}")
    
    results = []
    
    for protein_dir_name in tqdm(test_proteins, desc="Analyzing proteins"):
        try:
            # Use the full protein directory name for both reference and predictions
            protein_name = protein_dir_name.replace("boltz_results_", "")
            
            protein_outdir = os.path.join(base_outdir, protein_dir_name)
            pred_subdir = os.path.join(protein_outdir, "predictions", protein_name)
            
            if not os.path.exists(pred_subdir):
                print(f"[WARNING] Prediction directory not found: {pred_subdir}")
                continue
            
            # Reference paths - use full protein name
            ref_protein_path = os.path.join(data_dir, protein_name, f"{protein_name}_protein.pdb")
            ref_ligand_path = os.path.join(data_dir, protein_name, f"{protein_name}_ligand.sdf")
            
            if not os.path.exists(ref_protein_path) or not os.path.exists(ref_ligand_path):
                print(f"[WARNING] Reference files not found for {protein_name}")
                continue
            
            # Find prediction files
            pred_files = glob.glob(os.path.join(pred_subdir, f"*_model_*.pdb"))
            pred_files.sort(key=lambda x: int(x.split('_model_')[1].split('.')[0]))
            
            for pred_file in pred_files[:3]:  # Top 3 models for testing
                model_num = int(pred_file.split('_model_')[1].split('.')[0])
                
                # Get confidence score
                confidence_file = os.path.join(pred_subdir, f"confidence_{protein_name}_model_{model_num}.json")
                confidence = float('nan')
                if os.path.exists(confidence_file):
                    try:
                        with open(confidence_file, 'r') as f:
                            confidence_data = json.load(f)
                            confidence = confidence_data.get('confidence_score', float('nan'))
                    except:
                        pass
                
                # Extract predicted ligand
                pred_ligand_path = extract_ligand_from_pdb_comprehensive(pred_file)
                
                # Compute protein RMSD (comprehensive)
                protein_metrics = compute_protein_rmsd_comprehensive(pred_file, ref_protein_path)
                
                # Compute ligand RMSD (comprehensive)
                ligand_metrics = {}
                if pred_ligand_path and os.path.exists(pred_ligand_path):
                    ligand_metrics = compute_ligand_rmsd_comprehensive(pred_ligand_path, ref_ligand_path)
                else:
                    ligand_metrics = {
                        'rdkit_rmsd': float('nan'),
                        'coordinate_rmsd': float('nan'),
                        'aligned_rmsd': float('nan'),
                        'atom_count_match': False
                    }
                
                # Compute binding site RMSD (temporarily disabled due to unpacking issue)
                # try:
                #     binding_site_aligned_rmsd, binding_site_unaligned_rmsd = compute_binding_site_rmsd(pred_file, ref_protein_path, ref_ligand_path)
                # except Exception as e:
                #     print(f"[WARNING] Binding site RMSD failed: {e}")
                binding_site_aligned_rmsd = float('nan')
                binding_site_unaligned_rmsd = float('nan')
                
                # Combine results
                result = {
                    'protein': protein_name,
                    'model': model_num,
                    'confidence': confidence,
                    'protein_ca_aligned_rmsd': protein_metrics.get('ca_aligned_rmsd', float('nan')),
                    'protein_backbone_aligned_rmsd': protein_metrics.get('backbone_aligned_rmsd', float('nan')),
                    'ligand_rmsd': ligand_metrics.get('aligned_rmsd', float('nan')),
                    'ligand_coordinate_rmsd': ligand_metrics.get('coordinate_rmsd', float('nan')),
                    'binding_site_aligned_rmsd': binding_site_aligned_rmsd,
                    'atom_count_match': ligand_metrics.get('atom_count_match', False)
                }
                
                results.append(result)
                
                print(f"[INFO] {protein_name} model {model_num}: "
                      f"Protein_CA={result['protein_ca_aligned_rmsd']:.2f}Å, "
                      f"Ligand={result['ligand_rmsd']:.2f}Å, "
                      f"BindingSite={result['binding_site_aligned_rmsd']:.2f}Å")
                
        except Exception as e:
            print(f"[ERROR] Failed to process {protein_dir_name}: {str(e)}")
            continue
    
    # Create DataFrame and analyze
    df = pd.DataFrame(results)
    
    if len(df) > 0:
        print(f"\n[INFO] Generated {len(df)} prediction results")
        
        # Perform dependency analysis
        dependency_analysis = analyze_protein_ligand_dependency(df)
        
        if dependency_analysis:
            print("\n" + "="*60)
            print("CORRELATION ANALYSIS")
            print("="*60)
            
            correlations = dependency_analysis.get('correlations', {})
            for metric_pair, corr_data in correlations.items():
                if 'pearson_r' in corr_data:
                    print(f"{metric_pair}: r={corr_data['pearson_r']:.3f} (p={corr_data['pearson_p']:.3f})")
            
            print("\n" + "="*60)
            print("SUMMARY STATISTICS")
            print("="*60)
            
            stats = dependency_analysis.get('summary_stats', {})
            for metric, metric_stats in stats.items():
                if 'mean' in metric_stats:
                    print(f"{metric}: {metric_stats['mean']:.3f} ± {metric_stats['std']:.3f}")
        
        # Create plots
        create_dependency_plots(df)
        
        # Save results
        df.to_csv(output_csv, index=False)
        print(f"\n[INFO] Results saved to {output_csv}")
        
        return df
    else:
        print("[ERROR] No results generated")
        return pd.DataFrame()


def extract_ligand_from_pdb_comprehensive(pdb_path: str) -> Optional[str]:
    """Extract ligand from PDB file and save as SDF with robust error handling."""
    try:
        import tempfile
        
        # Extract HETATM records
        ligand_lines = []
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    ligand_lines.append(line)
        
        if not ligand_lines:
            return None
        
        # Create temporary ligand PDB file
        temp_pdb = tempfile.mktemp(suffix='_ligand.pdb')
        with open(temp_pdb, 'w') as f:
            for line in ligand_lines:
                f.write(line)
            f.write('END\n')
        
        # Convert to molecule using the existing robust reading function
        mol = read_molecule_robust(temp_pdb)
        if mol is None:
            os.remove(temp_pdb)
            return None
        
        # Save to SDF
        sdf_path = pdb_path.replace('.pdb', '_ligand_extracted.sdf')
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()
        
        # Clean up
        os.remove(temp_pdb)
        return sdf_path
        
    except Exception as e:
        print(f"[ERROR] Failed to extract ligand from {pdb_path}: {e}")
        return None


if __name__ == "__main__":
    # Run the comprehensive analysis
    df_results = run_comprehensive_boltz_analysis()
    
    if len(df_results) > 0:
        print(f"\nFinal results shape: {df_results.shape}")
        print("\nFirst few results:")
        print(df_results.head())
        
        print("\nSample correlation analysis:")
        valid_data = df_results.dropna(subset=['protein_ca_aligned_rmsd', 'ligand_rmsd'])
        if len(valid_data) >= 3:
            from scipy.stats import pearsonr
            corr, p_value = pearsonr(valid_data['protein_ca_aligned_rmsd'], valid_data['ligand_rmsd'])
            print(f"Protein CA RMSD vs Ligand RMSD: r={corr:.3f}, p={p_value:.3f}")
    else:
        print("No valid results to analyze")


# Keep the original function for backward compatibility

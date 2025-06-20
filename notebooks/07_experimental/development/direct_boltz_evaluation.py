#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
try:
    from rdkit.Chem import rdMolDescriptors, rdMolAlign
except ImportError:
    # Import from alternative locations if needed
    import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
    import rdkit.Chem.rdMolAlign as rdMolAlign

import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_direct_evaluation():
    """Test Boltz evaluation by using existing processed files directly"""
    logger.info("Starting direct Boltz evaluation...")
    
    approach = BoltzApproach()
    
    # Test systems with their full system names
    test_systems = {
        '1afb__1__1.A__1.D_1.F': '1afb__1__1.A__1.D_1.F',
        '1b5d__1__1.A_1.B__1.D': '1b5d__1__1.A_1.B__1.D',
        '1bcj__1__1.B__1.I_1.K': '1bcj__1__1.B__1.I_1.K',
        '1ci0__1__1.A_1.B__1.D': '1ci0__1__1.A_1.B__1.D'
    }
    
    results = []
    
    for system_id, reference_id in test_systems.items():
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing system: {system_id}")
        logger.info(f"{'='*50}")
        
        try:
            # Get the processed ligand SDF file that already exists
            pred_sdf = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}/predictions/{system_id}/{system_id}_model_0_ligand_protein_aligned.sdf"
            ref_sdf = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{reference_id}/{reference_id}_ligand.sdf"
            
            logger.info(f"Prediction SDF: {pred_sdf}")
            logger.info(f"Reference SDF: {ref_sdf}")
            logger.info(f"Pred exists: {os.path.exists(pred_sdf)}")
            logger.info(f"Ref exists: {os.path.exists(ref_sdf)}")
            
            if os.path.exists(pred_sdf) and os.path.exists(ref_sdf):
                try:
                    # Read molecules directly
                    pred_mol = Chem.SDMolSupplier(pred_sdf)[0]
                    ref_mol = Chem.SDMolSupplier(ref_sdf)[0]
                    
                    if pred_mol and ref_mol:
                        logger.info(f"Successfully read both molecules")
                        logger.info(f"Pred atoms: {pred_mol.GetNumAtoms()}, Ref atoms: {ref_mol.GetNumAtoms()}")
                        
                        # Get molecular formulas
                        pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                        ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                        formula_match = pred_formula == ref_formula
                        
                        logger.info(f"Pred formula: {pred_formula}, Ref formula: {ref_formula}")
                        logger.info(f"Formula match: {formula_match}")
                        
                        # Calculate RMSD
                        try:
                            complex_rmsd = rdMolAlign.AlignMol(pred_mol, ref_mol)
                            logger.info(f"Complex RMSD: {complex_rmsd:.3f} Å")
                            
                            # Calculate heavy-atom RMSD for comparison
                            heavy_rmsd = calculate_heavy_atom_rmsd(pred_mol, ref_mol)
                            logger.info(f"Heavy-atom RMSD: {heavy_rmsd:.3f} Å")
                            
                            # Extract the PDB ID from the system name
                            pdb_id = system_id.split('__')[0]
                            
                            results.append({
                                'system_id': system_id,
                                'pdb_id': pdb_id,
                                'complex_rmsd': complex_rmsd,
                                'heavy_atom_rmsd': heavy_rmsd,
                                'pred_formula': pred_formula,
                                'ref_formula': ref_formula,
                                'formula_match': formula_match,
                                'pred_atoms': pred_mol.GetNumAtoms(),
                                'ref_atoms': ref_mol.GetNumAtoms()
                            })
                            
                        except Exception as e:
                            logger.error(f"RMSD calculation failed: {str(e)}")
                    else:
                        logger.error(f"Failed to read molecules: pred={pred_mol is not None}, ref={ref_mol is not None}")
                        
                except Exception as e:
                    logger.error(f"Error reading SDF files: {str(e)}")
            else:
                logger.warning(f"Files not found for {system_id}")
                
        except Exception as e:
            logger.error(f"Error processing {system_id}: {str(e)}")
    
    # Summary statistics
    logger.info(f"\n{'='*50}")
    logger.info("DIRECT EVALUATION SUMMARY")
    logger.info(f"{'='*50}")
    logger.info(f"Systems successfully processed: {len(results)}")
    
    if results:
        df = pd.DataFrame(results)
        
        logger.info(f"\nResults summary:")
        logger.info(f"Mean complex RMSD: {df['complex_rmsd'].mean():.3f} ± {df['complex_rmsd'].std():.3f} Å")
        logger.info(f"Mean heavy-atom RMSD: {df['heavy_atom_rmsd'].mean():.3f} ± {df['heavy_atom_rmsd'].std():.3f} Å")
        logger.info(f"Systems with complex RMSD < 2 Å: {(df['complex_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Systems with heavy-atom RMSD < 2 Å: {(df['heavy_atom_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Formula matches: {df['formula_match'].sum()}/{len(df)}")
        
        # Success rate analysis
        success_threshold = 2.0
        complex_success_rate = (df['complex_rmsd'] < success_threshold).mean() * 100
        heavy_success_rate = (df['heavy_atom_rmsd'] < success_threshold).mean() * 100
        
        logger.info(f"\nSuccess rates (RMSD < {success_threshold} Å):")
        logger.info(f"Complex RMSD success rate: {complex_success_rate:.1f}%")
        logger.info(f"Heavy-atom RMSD success rate: {heavy_success_rate:.1f}%")
        
        # Detailed analysis
        logger.info(f"\nDetailed results:")
        for _, row in df.iterrows():
            logger.info(f"{row['pdb_id']} ({row['system_id']}): "
                       f"Complex RMSD={row['complex_rmsd']:.3f}Å, "
                       f"Heavy RMSD={row['heavy_atom_rmsd']:.3f}Å, "
                       f"Formula match={row['formula_match']}, "
                       f"Atoms: {row['pred_atoms']}/{row['ref_atoms']}")
        
        # Save results
        output_path = '/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks/direct_boltz_evaluation_results.csv'
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        return df
    else:
        logger.warning("No results obtained")
        return None

def calculate_heavy_atom_rmsd(mol1, mol2):
    """Calculate RMSD using only heavy atoms"""
    try:
        from scipy.spatial.distance import cdist
        
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        # Get heavy atom indices
        heavy_atoms1 = [i for i in range(mol1.GetNumAtoms()) if mol1.GetAtomWithIdx(i).GetSymbol() != 'H']
        heavy_atoms2 = [i for i in range(mol2.GetNumAtoms()) if mol2.GetAtomWithIdx(i).GetSymbol() != 'H']
        
        if len(heavy_atoms1) != len(heavy_atoms2):
            logger.warning(f"Different heavy atom counts: {len(heavy_atoms1)} vs {len(heavy_atoms2)}")
            # Use minimum for comparison
            min_atoms = min(len(heavy_atoms1), len(heavy_atoms2))
            heavy_atoms1 = heavy_atoms1[:min_atoms]
            heavy_atoms2 = heavy_atoms2[:min_atoms]
        
        if len(heavy_atoms1) == 0:
            return float('inf')
        
        # Get coordinates
        coords1 = np.array([conf1.GetAtomPosition(i) for i in heavy_atoms1])
        coords2 = np.array([conf2.GetAtomPosition(i) for i in heavy_atoms2])
        
        # Simple distance-based assignment for now
        distances = cdist(coords1, coords2)
        used_j = set()
        total_dist_sq = 0
        assignments = 0
        
        for i in range(len(coords1)):
            best_j = None
            best_dist = float('inf')
            for j in range(len(coords2)):
                if j not in used_j and distances[i][j] < best_dist:
                    best_dist = distances[i][j]
                    best_j = j
                    
            if best_j is not None:
                used_j.add(best_j)
                total_dist_sq += best_dist**2
                assignments += 1
                
        if assignments == 0:
            return float('inf')
            
        return np.sqrt(total_dist_sq / assignments)
        
    except Exception as e:
        logger.error(f"Heavy atom RMSD calculation failed: {str(e)}")
        return float('inf')

if __name__ == "__main__":
    test_direct_evaluation()

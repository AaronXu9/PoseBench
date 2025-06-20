#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
try:
    from rdkit.Chem import rdMolDescriptors, rdMolAlign
except ImportError:
    import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
    import rdkit.Chem.rdMolAlign as rdMolAlign

import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

try:
    from notebooks.alignment import calculate_complex_rmsd
except ImportError:
    # Fallback implementation if alignment module is not available
    def calculate_complex_rmsd(pred_path, ref_path):
        """Simplified RMSD calculation"""
        try:
            pred_mol = Chem.SDMolSupplier(pred_path)[0]
            ref_mol = Chem.SDMolSupplier(ref_path)[0]
            if pred_mol and ref_mol:
                return rdMolAlign.AlignMol(pred_mol, ref_mol)
            return float('inf')
        except:
            return float('inf')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_heavy_atom_rmsd(mol1, mol2):
    """Calculate RMSD using only heavy atoms"""
    if not mol1 or not mol2:
        return float('inf')
        
    try:
        # Get heavy atom coordinates
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        heavy_atoms1 = [i for i in range(mol1.GetNumAtoms()) if mol1.GetAtomWithIdx(i).GetSymbol() != 'H']
        heavy_atoms2 = [i for i in range(mol2.GetNumAtoms()) if mol2.GetAtomWithIdx(i).GetSymbol() != 'H']
        
        # Use minimum number of atoms
        min_atoms = min(len(heavy_atoms1), len(heavy_atoms2))
        if min_atoms == 0:
            return float('inf')
            
        # Get coordinates of first min_atoms heavy atoms
        coords1 = np.array([conf1.GetAtomPosition(heavy_atoms1[i]) for i in range(min_atoms)])
        coords2 = np.array([conf2.GetAtomPosition(heavy_atoms2[i]) for i in range(min_atoms)])
        
        # Simple alignment - just calculate RMSD without optimization
        return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
        
    except Exception as e:
        logger.error(f"Heavy atom RMSD calculation failed: {str(e)}")
        return float('inf')

def test_existing_boltz_results():
    """Test using existing processed Boltz results"""
    logger.info("Testing existing Boltz processed results...")
    
    # Systems with existing processed files
    test_systems = ['1afb__1__1.A__1.D_1.F', '1b5d__1__1.A_1.B__1.D', '1bcj__1__1.B__1.I_1.K', '1ci0__1__1.A_1.B__1.D']
    
    results = []
    for system_id in test_systems:
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing system: {system_id}")
        logger.info(f"{'='*50}")
        
        try:
            # Path to Boltz results
            base_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}/predictions/{system_id}"
            
            # Try different processed ligand files
            ligand_files = [
                f"{system_id}_model_0_ligand.sdf",
                f"{system_id}_model_0_ligand_protein_aligned.sdf"
            ]
            
            processed_ligand_path = None
            for ligand_file in ligand_files:
                test_path = os.path.join(base_path, ligand_file)
                if os.path.exists(test_path):
                    processed_ligand_path = test_path
                    break
                    
            if not processed_ligand_path:
                logger.warning(f"No processed ligand found for {system_id}")
                continue
                
            logger.info(f"Using processed ligand: {processed_ligand_path}")
            
            # Get reference ligand
            pdb_id = system_id.split('__')[0]  # Extract PDB ID
            ref_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{pdb_id}/{pdb_id}_ligand.sdf"
            
            if not os.path.exists(ref_ligand_path):
                logger.warning(f"Reference ligand not found: {ref_ligand_path}")
                continue
                
            # Load molecules
            try:
                pred_mol = Chem.SDMolSupplier(processed_ligand_path)[0]
                ref_mol = Chem.SDMolSupplier(ref_ligand_path)[0]
                
                if not pred_mol or not ref_mol:
                    logger.error(f"Could not load molecules for {system_id}")
                    continue
                    
                # Calculate RMSD
                try:
                    complex_rmsd = calculate_complex_rmsd(processed_ligand_path, ref_ligand_path)
                except:
                    complex_rmsd = float('inf')
                    
                heavy_rmsd = calculate_heavy_atom_rmsd(pred_mol, ref_mol)
                
                # Get molecular formulas
                pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                
                logger.info(f"Complex RMSD: {complex_rmsd:.3f} Å")
                logger.info(f"Heavy-atom RMSD: {heavy_rmsd:.3f} Å")
                logger.info(f"Predicted formula: {pred_formula}")
                logger.info(f"Reference formula: {ref_formula}")
                
                results.append({
                    'system_id': system_id,
                    'pdb_id': pdb_id,
                    'complex_rmsd': complex_rmsd,
                    'heavy_atom_rmsd': heavy_rmsd,
                    'pred_formula': pred_formula,
                    'ref_formula': ref_formula,
                    'formula_match': pred_formula == ref_formula,
                    'ligand_file_used': os.path.basename(processed_ligand_path)
                })
                
            except Exception as e:
                logger.error(f"Error processing molecules for {system_id}: {str(e)}")
                
        except Exception as e:
            logger.error(f"Error processing {system_id}: {str(e)}")
            
    # Summary
    logger.info(f"\n{'='*50}")
    logger.info("BOLTZ RESULTS SUMMARY")
    logger.info(f"{'='*50}")
    logger.info(f"Systems processed: {len(results)}")
    
    if results:
        df = pd.DataFrame(results)
        
        logger.info(f"\nResults summary:")
        logger.info(f"Mean complex RMSD: {df['complex_rmsd'].mean():.3f} ± {df['complex_rmsd'].std():.3f} Å")
        logger.info(f"Mean heavy-atom RMSD: {df['heavy_atom_rmsd'].mean():.3f} ± {df['heavy_atom_rmsd'].std():.3f} Å")
        logger.info(f"Systems with complex RMSD < 2 Å: {(df['complex_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Systems with heavy-atom RMSD < 2 Å: {(df['heavy_atom_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Formula matches: {df['formula_match'].sum()}/{len(df)}")
        
        # Detailed analysis
        logger.info(f"\nDetailed analysis:")
        for _, row in df.iterrows():
            logger.info(f"{row['system_id']}: Complex RMSD={row['complex_rmsd']:.3f}Å, "
                       f"Heavy RMSD={row['heavy_atom_rmsd']:.3f}Å, "
                       f"Formula match={row['formula_match']}, File={row['ligand_file_used']}")
        
        # Save results
        output_path = '/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks/boltz_direct_results.csv'
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        return df
    else:
        logger.warning("No results obtained")
        return None

if __name__ == "__main__":
    test_existing_boltz_results()

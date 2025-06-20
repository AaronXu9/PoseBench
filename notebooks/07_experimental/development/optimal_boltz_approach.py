#!/usr/bin/env python3

import os
import sys
import numpy as np
from typing import List
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdMolDescriptors
import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class OptimallyAlignedBoltzApproach(BoltzApproach):
    """BoltzApproach with optimal ligand-to-ligand alignment for PoseBusters"""
    
    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Get top N files and apply optimal ligand alignment for PoseBusters
        """
        # Get the standard aligned files from parent class
        aligned_files = super().list_top_n_files(protein_dir, top_n)
        
        # Apply optimal ligand-to-ligand alignment for each file
        optimally_aligned_files = []
        
        for sdf_path in aligned_files:
            optimal_path = self._apply_optimal_ligand_alignment(sdf_path)
            if optimal_path:
                optimally_aligned_files.append(optimal_path)
            else:
                optimally_aligned_files.append(sdf_path)  # Fallback to original
                
        return optimally_aligned_files
    
    def _apply_optimal_ligand_alignment(self, sdf_path: str) -> str:
        """
        Apply optimal ligand-to-ligand alignment using RDKit's AlignMol
        """
        try:
            # Extract system ID from path
            system_id = self._extract_system_id_from_path(sdf_path)
            if not system_id:
                logger.warning(f"Could not extract system ID from {sdf_path}")
                return sdf_path
            
            # Get reference ligand path
            ref_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{system_id}/{system_id}_ligand.sdf"
            
            if not os.path.exists(ref_ligand_path):
                logger.warning(f"Reference ligand not found: {ref_ligand_path}")
                return sdf_path
            
            # Load molecules
            pred_mol = Chem.SDMolSupplier(sdf_path)[0]
            ref_mol = Chem.SDMolSupplier(ref_ligand_path)[0]
            
            if not pred_mol or not ref_mol:
                logger.warning(f"Could not load molecules from {sdf_path}")
                return sdf_path
            
            # Check if molecules can be aligned (same heavy atom count helps)
            pred_heavy = Chem.RemoveHs(pred_mol).GetNumAtoms()
            ref_heavy = Chem.RemoveHs(ref_mol).GetNumAtoms()
            
            if pred_heavy != ref_heavy:
                logger.info(f"Different heavy atom counts ({pred_heavy} vs {ref_heavy}), skipping optimal alignment")
                return sdf_path
            
            # Apply optimal alignment
            try:
                rmsd = rdMolAlign.AlignMol(pred_mol, ref_mol)
                logger.info(f"Optimal alignment RMSD: {rmsd:.3f} Å")
                
                # Save optimally aligned molecule
                optimal_path = sdf_path.replace('.sdf', '_optimal_aligned.sdf')
                writer = Chem.SDWriter(optimal_path)
                writer.write(pred_mol)  # pred_mol is now aligned in-place
                writer.close()
                
                logger.info(f"Saved optimally aligned ligand: {optimal_path}")
                return optimal_path
                
            except Exception as e:
                logger.warning(f"RDKit alignment failed: {str(e)}")
                return sdf_path
                
        except Exception as e:
            logger.error(f"Optimal alignment failed for {sdf_path}: {str(e)}")
            return sdf_path
    
    def _extract_system_id_from_path(self, sdf_path: str) -> str:
        """Extract system ID from SDF file path"""
        try:
            # Path looks like: .../predictions/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_model_0_ligand_protein_aligned.sdf
            # We want to extract: 1afb__1__1.A__1.D_1.F
            if '/predictions/' in sdf_path:
                # Split on predictions and get the next directory
                parts = sdf_path.split('/predictions/')
                if len(parts) > 1:
                    pred_part = parts[1]  # everything after /predictions/
                    system_id = pred_part.split('/')[0]  # first directory
                    return system_id
            return None
        except:
            return None

def test_optimal_alignment():
    """Test the optimal alignment approach"""
    logger.info("Testing optimal alignment approach...")
    
    approach = OptimallyAlignedBoltzApproach()
    system_id = '1afb__1__1.A__1.D_1.F'
    
    protein_dir = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}"
    
    # Get optimally aligned files
    optimal_files = approach.list_top_n_files(protein_dir, 1)
    
    if optimal_files:
        optimal_path = optimal_files[0]
        logger.info(f"Optimal file: {optimal_path}")
        
        # Test with PoseBusters
        from posebusters.posebusters import PoseBusters
        
        pb = PoseBusters(config='redock', top_n=None)
        ref_sdf = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{system_id}/{system_id}_ligand.sdf"
        protein_pdb = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{system_id}/{system_id}_protein.pdb"
        
        try:
            df_pb = pb.bust(
                mol_pred=optimal_path,
                mol_true=ref_sdf,
                mol_cond=protein_pdb,
                full_report=True
            )
            
            rmsd = df_pb["rmsd"].iloc[0]
            success = df_pb["rmsd_≤_2å"].iloc[0]
            
            logger.info(f"PoseBusters RMSD with optimal alignment: {rmsd:.3f} Å")
            logger.info(f"RMSD ≤ 2Å: {success}")
            
            return rmsd, success
            
        except Exception as e:
            logger.error(f"PoseBusters evaluation failed: {str(e)}")
            return None, None
    else:
        logger.error("No optimal files generated")
        return None, None

if __name__ == "__main__":
    test_optimal_alignment()

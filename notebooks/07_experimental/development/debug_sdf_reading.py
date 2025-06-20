#!/usr/bin/env python3

import os
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def debug_sdf_reading():
    """Debug SDF file reading issues"""
    
    # Test a specific SDF file that was generated
    sdf_path = '/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F/predictions/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_model_0_ligand_protein_aligned.sdf'
    
    logger.info(f"Testing SDF file: {sdf_path}")
    logger.info(f"File exists: {os.path.exists(sdf_path)}")
    
    if os.path.exists(sdf_path):
        # Test different ways of reading
        logger.info("Testing different reading methods...")
        
        # Method 1: SDMolSupplier with removeHs=False
        try:
            supplier1 = Chem.SDMolSupplier(sdf_path, removeHs=False)
            mol1 = supplier1[0]
            logger.info(f"Method 1 (removeHs=False): {mol1 is not None}")
            if mol1:
                logger.info(f"  Atoms: {mol1.GetNumAtoms()}")
                logger.info(f"  Formula: {rdMolDescriptors.CalcMolFormula(mol1)}")
        except Exception as e:
            logger.error(f"Method 1 failed: {str(e)}")
        
        # Method 2: SDMolSupplier with removeHs=True
        try:
            supplier2 = Chem.SDMolSupplier(sdf_path, removeHs=True)
            mol2 = supplier2[0]
            logger.info(f"Method 2 (removeHs=True): {mol2 is not None}")
            if mol2:
                logger.info(f"  Atoms: {mol2.GetNumAtoms()}")
                logger.info(f"  Formula: {rdMolDescriptors.CalcMolFormula(mol2)}")
        except Exception as e:
            logger.error(f"Method 2 failed: {str(e)}")
        
        # Method 3: SDMolSupplier with defaults
        try:
            supplier3 = Chem.SDMolSupplier(sdf_path)
            mol3 = supplier3[0]
            logger.info(f"Method 3 (default): {mol3 is not None}")
            if mol3:
                logger.info(f"  Atoms: {mol3.GetNumAtoms()}")
                logger.info(f"  Formula: {rdMolDescriptors.CalcMolFormula(mol3)}")
        except Exception as e:
            logger.error(f"Method 3 failed: {str(e)}")
        
        # Method 4: MolFromMolFile
        try:
            mol4 = Chem.MolFromMolFile(sdf_path, removeHs=False)
            logger.info(f"Method 4 (MolFromMolFile): {mol4 is not None}")
            if mol4:
                logger.info(f"  Atoms: {mol4.GetNumAtoms()}")
                logger.info(f"  Formula: {rdMolDescriptors.CalcMolFormula(mol4)}")
        except Exception as e:
            logger.error(f"Method 4 failed: {str(e)}")
    
    # Test the BoltzApproach list_top_n_files method
    logger.info("\nTesting BoltzApproach.list_top_n_files...")
    
    approach = BoltzApproach()
    
    # Test the system that had the SDF file
    system_dir = '/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F'
    logger.info(f"Testing system directory: {system_dir}")
    
    try:
        files = approach.list_top_n_files(system_dir, 1)
        logger.info(f"Found {len(files)} files: {files}")
        
        if files:
            pred_file = files[0]
            logger.info(f"Prediction file: {pred_file}")
            logger.info(f"File exists: {os.path.exists(pred_file)}")
            
            # Test reading this file
            if os.path.exists(pred_file):
                mol = Chem.MolFromPDBFile(pred_file, removeHs=False)
                logger.info(f"Can read PDB: {mol is not None}")
                if mol:
                    logger.info(f"  PDB Atoms: {mol.GetNumAtoms()}")
                    logger.info(f"  PDB Formula: {rdMolDescriptors.CalcMolFormula(mol)}")
    except Exception as e:
        logger.error(f"BoltzApproach test failed: {str(e)}")

if __name__ == "__main__":
    debug_sdf_reading()

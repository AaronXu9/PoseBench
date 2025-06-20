#!/usr/bin/env python3

import os
import sys
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_sdf_reading_writing():
    """Test SDF file reading/writing to identify the issue"""
    
    # Path to a sample Boltz prediction
    boltz_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0"
    
    # Find a test file
    for root, dirs, files in os.walk(boltz_dir):
        for file in files:
            if file.endswith("_ligand_0.sdf"):
                test_file = os.path.join(root, file)
                logger.info(f"Testing with file: {test_file}")
                
                # Test reading the original file
                logger.info("1. Reading original SDF file...")
                try:
                    supplier = Chem.SDMolSupplier(test_file)
                    mol = supplier[0]
                    if mol:
                        logger.info(f"Successfully read molecule: {Chem.MolToSmiles(mol)}")
                        logger.info(f"Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
                        logger.info(f"Num atoms: {mol.GetNumAtoms()}")
                    else:
                        logger.error("Failed to read molecule from original file")
                        continue
                except Exception as e:
                    logger.error(f"Error reading original file: {e}")
                    continue
                
                # Test writing and re-reading
                logger.info("2. Writing and re-reading SDF file...")
                temp_file = "/tmp/test_output.sdf"
                try:
                    writer = Chem.SDWriter(temp_file)
                    writer.write(mol)
                    writer.close()
                    
                    # Try to read it back
                    supplier2 = Chem.SDMolSupplier(temp_file)
                    mol2 = supplier2[0]
                    if mol2:
                        logger.info("Successfully wrote and read back molecule")
                        logger.info(f"Re-read SMILES: {Chem.MolToSmiles(mol2)}")
                    else:
                        logger.error("Failed to read back written molecule")
                except Exception as e:
                    logger.error(f"Error in write/read test: {e}")
                
                # Test with hydrogen manipulation
                logger.info("3. Testing with hydrogen manipulation...")
                try:
                    mol_no_h = Chem.RemoveHs(mol)
                    mol_with_h = Chem.AddHs(mol_no_h)
                    
                    temp_file2 = "/tmp/test_output_h.sdf"
                    writer = Chem.SDWriter(temp_file2)
                    writer.write(mol_with_h)
                    writer.close()
                    
                    supplier3 = Chem.SDMolSupplier(temp_file2)
                    mol3 = supplier3[0]
                    if mol3:
                        logger.info("Successfully manipulated hydrogens and wrote/read")
                        logger.info(f"H-manipulated SMILES: {Chem.MolToSmiles(mol3)}")
                    else:
                        logger.error("Failed to read hydrogen-manipulated molecule")
                        
                except Exception as e:
                    logger.error(f"Error in hydrogen manipulation test: {e}")
                
                return  # Just test the first file
    
    logger.error("No suitable test files found")

if __name__ == "__main__":
    test_sdf_reading_writing()

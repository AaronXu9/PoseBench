#!/usr/bin/env python3

import os
from rdkit import Chem

def simple_analysis():
    """Simple analysis of reference vs predicted ligand"""
    
    print("Starting simple ligand analysis...")
    
    # Paths
    protein_name = "1afb__1__1.A__1.D_1.F"
    ref_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
    
    print(f"Checking reference ligand: {ref_ligand_path}")
    
    if os.path.exists(ref_ligand_path):
        print("File exists, trying to read...")
        try:
            suppl = Chem.SDMolSupplier(ref_ligand_path, removeHs=False)
            mol = None
            for m in suppl:
                if m is not None:
                    mol = m
                    break
            
            if mol:
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                print(f"Reference formula: {formula}")
                print(f"Reference atoms: {mol.GetNumAtoms()}")
            else:
                print("Could not read molecule from SDF")
        except Exception as e:
            print(f"Error reading SDF: {e}")
    else:
        print("Reference file does not exist")

if __name__ == "__main__":
    simple_analysis()

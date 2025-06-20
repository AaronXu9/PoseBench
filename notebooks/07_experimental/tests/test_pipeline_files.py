#!/usr/bin/env python3

import os
import sys
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach

def test_pipeline_files():
    """Test what files the pipeline is actually getting"""
    
    approach = BoltzApproach()
    system_id = '1afb__1__1.A__1.D_1.F'
    
    # This mimics what the pipeline does
    protein_dir = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}"
    
    print(f"Testing pipeline file discovery for: {system_id}")
    print(f"Protein directory: {protein_dir}")
    print()
    
    # Get top 5 files (same as pipeline)
    sdf_paths = approach.list_top_n_files(protein_dir, 5)
    
    print(f"Files returned by list_top_n_files ({len(sdf_paths)}):")
    for i, path in enumerate(sdf_paths, 1):
        print(f"{i}. {path}")
        print(f"   Exists: {os.path.exists(path)}")
        if os.path.exists(path):
            # Check file size to see if it's reasonable
            size = os.path.getsize(path)
            print(f"   Size: {size} bytes")
            
            # Try to read it and check atom count
            try:
                from rdkit import Chem
                mol = Chem.SDMolSupplier(path)[0]
                if mol:
                    print(f"   Atoms: {mol.GetNumAtoms()}")
                else:
                    print(f"   Could not read molecule")
            except Exception as e:
                print(f"   Read error: {e}")
        print()

if __name__ == "__main__":
    test_pipeline_files()

#!/usr/bin/env python3

import os
import sys
from Approach import BoltzApproach

def test_boltz_alignment():
    """Test BoltzApproach alignment functionality on a single structure"""
    
    print("Starting Boltz alignment test...")
    
    # Initialize the approach
    try:
        boltz = BoltzApproach()
        print(f"✓ Successfully initialized {boltz.get_name()} approach")
    except Exception as e:
        print(f"✗ Failed to initialize BoltzApproach: {e}")
        return
    
    # Test on a specific protein structure
    protein_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F"
    
    print(f"Checking test directory: {protein_dir}")
    if not os.path.exists(protein_dir):
        print(f"✗ Test directory not found: {protein_dir}")
        # Let's try to find any boltz results directory
        base_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0"
        if os.path.exists(base_dir):
            print(f"Base directory exists, looking for any boltz_results_* directory...")
            dirs = [d for d in os.listdir(base_dir) if d.startswith('boltz_results_')]
            if dirs:
                protein_dir = os.path.join(base_dir, dirs[0])
                print(f"Using: {protein_dir}")
            else:
                print(f"No boltz_results_* directories found in {base_dir}")
                return
        else:
            print(f"Base directory does not exist: {base_dir}")
            return
    else:
        print(f"✓ Test directory found")
    
    print(f"Testing on: {protein_dir}")
    
    # Test list_top_n_files function with alignment
    try:
        print("Calling list_top_n_files...")
        sdf_files = boltz.list_top_n_files(protein_dir, top_n=2)
        print(f"✓ Found {len(sdf_files)} aligned SDF files:")
        for i, sdf_path in enumerate(sdf_files):
            print(f"  {i+1}. {sdf_path}")
            if os.path.exists(sdf_path):
                print(f"     File exists: ✓")
                # Test score parsing
                score = boltz.parse_score(sdf_path)
                print(f"     Score: {score}")
            else:
                print(f"     File exists: ✗")
    except Exception as e:
        print(f"✗ Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_boltz_alignment()

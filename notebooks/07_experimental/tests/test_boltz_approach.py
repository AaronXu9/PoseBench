#!/usr/bin/env python3

import sys
import os
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks')

from Approach import BoltzApproach

def test_boltz_approach():
    """Test the BoltzApproach class"""
    approach = BoltzApproach()
    
    # Test basic functionality
    print(f"Method name: {approach.get_name()}")
    
    # Test with a real directory
    test_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F"
    
    if os.path.exists(test_dir):
        print(f"Testing with directory: {test_dir}")
        
        # List top files
        files = approach.list_top_n_files(test_dir, 5)
        print(f"Found {len(files)} SDF files:")
        for i, file in enumerate(files):
            print(f"  {i+1}: {file}")
            
            # Test score parsing
            if os.path.exists(file):
                score = approach.parse_score(file)
                print(f"     Score: {score}")
    else:
        print(f"Test directory not found: {test_dir}")

if __name__ == "__main__":
    test_boltz_approach()

#!/usr/bin/env python3

from posebusters.posebusters import PoseBusters
from rdkit import Chem
import os
import pandas as pd
import re 
from typing import List
from tqdm import tqdm
import glob

from Approach import BoltzApproach

def run_posebusters_test(
    approach,
    base_outdir: str,
    data_dir: str,
    top_n: int = 2,
    max_proteins: int = 3,
) -> pd.DataFrame:
    """
    Test PoseBusters on a limited number of proteins to validate alignment
    """
    pb = PoseBusters(config="redock", top_n=None)
    method_name = approach.get_name()
    all_rows = []
    
    protein_count = 0
    protein_dirs = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d))]
    
    for protein_name in tqdm(protein_dirs[:max_proteins]):
        protein_count += 1
        print(f"\n=== Processing {protein_count}/{max_proteins}: {protein_name} ===")
        
        try: 
            protein_dir = glob.glob(f"{base_outdir}/*{protein_name}*")[0]
        except: 
            print(f"Could not find {protein_name} in {base_outdir}")
            continue
            
        if not os.path.isdir(protein_dir):
            continue

        # Retrieve up to top-N .sdf file paths
        print(f"Getting top {top_n} files...")
        sdf_paths = approach.list_top_n_files(protein_dir, top_n)
        if not sdf_paths:
            print(f"[{method_name}] No top-{top_n} SDF files found for {protein_name}")
            continue

        # References
        true_ligand = os.path.join(data_dir, protein_name, f"{protein_name}_ligand.sdf")
        protein_pdb = os.path.join(data_dir, protein_name, f"{protein_name}_protein.pdb")
        if not (os.path.isfile(true_ligand) and os.path.isfile(protein_pdb)):
            print(f"[{method_name}] Missing reference for {protein_name}")
            continue

        rank_counter = 1
        for sdf_path in sdf_paths:
            try:
                print(f"  Processing rank {rank_counter}: {os.path.basename(sdf_path)}")
                
                # Check if alignment worked (file should contain 'aligned' in name)
                is_aligned = 'aligned' in sdf_path
                print(f"    Aligned: {'✓' if is_aligned else '✗'}")
                
                # Parse score
                numeric_score = approach.parse_score(sdf_path)
                print(f"    Score: {numeric_score}")
                
                # Run PoseBusters with proper parameters
                print(f"    Running PoseBusters...")
                df_pb = pb.bust(
                    mol_pred=sdf_path,
                    mol_true=true_ligand,
                    mol_cond=protein_pdb,
                    full_report=True
                )
                
                # Extract key metrics
                rmsd = df_pb.get('rmsd', [float('nan')])[0] if 'rmsd' in df_pb.columns else float('nan')
                rmsd_lt2 = df_pb.get('rmsd_≤_2å', [False])[0] if 'rmsd_≤_2å' in df_pb.columns else False
                
                print(f"    RMSD: {rmsd:.3f}")
                print(f"    RMSD ≤ 2Å: {'✓' if rmsd_lt2 else '✗'}")
                
                # Add metadata
                df_pb["score"] = numeric_score
                df_pb["method"] = method_name
                df_pb["protein"] = protein_name
                df_pb["rank"] = rank_counter
                df_pb["aligned"] = is_aligned
                rank_counter += 1

                all_rows.append(df_pb)
                
            except Exception as e:
                print(f"[{method_name}] [ERROR] PoseBusters failed for {protein_name} rank {rank_counter}: {e}")

    if not all_rows:
        return pd.DataFrame()
    return pd.concat(all_rows, ignore_index=True)

def main():
    base_outdir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks"
    exp_name = "plinder_set_0"
    data_dir = "/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set"
    
    print("=== Testing Boltz Alignment in PoseBusters Pipeline ===")
    
    # Test with BoltzApproach
    approach = BoltzApproach()
    method_name = approach.get_name()
    base_dir = f"{base_outdir}/boltz/inference/{exp_name}"
    
    print(f"Method: {method_name}")
    print(f"Base directory: {base_dir}")
    print(f"Data directory: {data_dir}")
    
    df_method = run_posebusters_test(
        approach,
        base_outdir=base_dir,
        data_dir=data_dir,
        top_n=2,
        max_proteins=3
    )
    
    if not df_method.empty:
        print(f"\n=== Results Summary ===")
        print(f"Total results: {len(df_method)}")
        print(f"Methods: {df_method['method'].unique()}")
        print(f"Proteins: {df_method['protein'].unique()}")
        print(f"Average RMSD: {df_method['rmsd'].mean():.3f}")
        print(f"RMSD ≤ 2Å success rate: {df_method['rmsd_≤_2å'].mean():.2%}")
        
        # Save results
        output_file = f"boltz_alignment_test_results.csv"
        df_method.to_csv(output_file, index=False)
        print(f"Results saved to: {output_file}")
    else:
        print("No results generated")

if __name__ == "__main__":
    main()

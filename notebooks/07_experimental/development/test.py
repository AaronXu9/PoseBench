import os
import os
import pickle

import pandas as pd
import MDAnalysis as mda
import prolif as plf
from rdkit import DataStructs
from rdkit import Chem
from prolif.molecule import Molecule

import numpy as np

from analysis import ConformationAnalyzer, analyze_all_complexes, compute_reference_fingerprint
from plotting import plot_single_approach_posebusters, plot_overall_bar, plot_cumulative_ecdf

BASE_DIR = '/Users/aoxu/projects/DrugDiscovery/PoseBench'
methods = ["icm", "surfdock", "gnina", "diffdock_pocket_only"]

def load_df():
    BASE_DIR = '/Users/aoxu/projects/DrugDiscovery/PoseBench'

    df = pd.read_csv(f"{BASE_DIR}/notebooks/posebusters_results_filtered_with_descriptors_predicted_ligands.csv")

    for col in ["protein_pdb", "predicted_ligand", "true_ligand"]:
        df[col] = df[col].str.replace("/Users/aoxu/projects/DrugDiscovery/PoseBench", BASE_DIR)

    return df

def compute_interactions(df, method):
    for method in methods:
        method_df = df[df["method"] == method]

        # Group by protein and select the row with minimal RMSD for each protein.
        best_indices = method_df.groupby("protein")["rmsd"].idxmin()
        best_df = method_df.loc[best_indices].reset_index(drop=True)
        selected_dfs[method] = best_df

        print(f"Selected best conformations for method: {method}")
        # display(best_df)

        for split in ["failure", "success"]:
            input_df = best_df[best_df["rmsd_≤_2å"]] if split == "success" else best_df[~best_df["rmsd_≤_2å"]]
            # display(input_df['rmsd_≤_2å'].count())
            all_results = analyze_all_complexes(input_df, [method])

            # Print or save the results. For example, print overall comparisons for each complex:
            for complex_id, complex_data in all_results.items():
                print(f"\nComplex {complex_id}:")
                for method, result_dict in complex_data.items():
                    print(f"  Method: {method}")
                    print(result_dict["overall"])

            with open(os.path.join(BASE_DIR, f"PoseBench/notebooks/{method}_{split}.pkl"), "wb") as f:
                pickle.dump(all_results, f)
            
            break

def load_interactions(method):
    with open(os.path.join(BASE_DIR, f"PoseBench/notebooks/{method}_success.pkl"), "rb") as f:
        success_results = pickle.load(f)
    with open(os.path.join(BASE_DIR, f"PoseBench/notebooks/{method}_failure.pkl"), "rb") as f:
        failure_results = pickle.load(f)
    return success_results, failure_results

if __name__ == "__main__":
    df = load_df()
    selected_dfs = {}




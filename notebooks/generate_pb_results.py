from posebusters.posebusters import PoseBusters
from rdkit import Chem
import os
import pandas as pd
import re 
from typing import List
from tqdm import tqdm
import glob

from rdkit import Chem
import os
import pandas as pd
import re 
from typing import List
from Approach import DockingApproach, DiffDockApproach, ICMApproach, ChaiApproach, VinaApproach, GninaApproach, SurfDockApproach, DiffDockPocketApproach, BoltzApproach

def run_posebusters(
    approach: DockingApproach,
    base_outdir: str,
    data_dir: str,
    top_n: int = 5,
    docking: bool = False,
) -> pd.DataFrame:
    """
    For each protein subdir in base_outdir:
      - list up to top_n SDF files (method-specific naming)
      - run PoseBusters on each
      - parse a numeric score or confidence if available
      - collect results in a DataFrame
    """
    pb = PoseBusters(config="redock", top_n=None)
    method_name = approach.get_name()
    all_rows = []

    for protein_name in tqdm(os.listdir(data_dir)):
        try: 
            if "plinder" in data_dir:
                # protein_dir = glob.glob(f"{base_outdir}/*{protein_name.replace('.','_')}*")[0] ## for ICM specific
                protein_dir = glob.glob(f"{base_outdir}/*{protein_name}*")[0]
            else:
                protein_dir = glob.glob(f"{base_outdir}/*{protein_name}*")[0]
        # protein_dir = os.path.join(base_outdir, protein_name)
        except: 
            print(f"Could not find {protein_name} in {base_outdir}")
            continue
        if not os.path.isdir(protein_dir):
            continue

        # Retrieve up to top-N .sdf file paths
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
                df_pb = pb.bust(
                    mol_pred=sdf_path,
                    mol_true=true_ligand,
                    mol_cond=protein_pdb,
                    full_report=True
                )
                # parse numeric score or confidence if available
                numeric_score = approach.parse_score(sdf_path)

                df_pb["score"] = numeric_score  # or "confidence_score" or "docking_score"
                df_pb["method"] = method_name
                df_pb["protein"] = protein_name
                df_pb["rank"] = rank_counter
                rank_counter += 1

                all_rows.append(df_pb)
            except Exception as e:
                print(f"[{method_name}] [ERROR] PoseBusters failed for {protein_name}: {e}")

    if not all_rows:
        return pd.DataFrame()
    return pd.concat(all_rows, ignore_index=True)

# The folder containing subdirectories like 5S8I_2LY, 5SD5_HWI, etc.
base_outdir = "/home/aoxu/projects/PoseBench/forks"

exp_name = "plinder_set_0"

# Initialize your approaches
approaches = [
    # ICMApproach(),
    # DiffDockApproach(),
    # SurfDockApproach(),
    # DiffDockPocketApproach(), 
    # GninaApproach(), 
    # ChaiApproach(),
    # VinaApproach()
    BoltzApproach(),
    # ... add more if needed
]

BASE_DIRS = {
    "icm": f"{base_outdir}/ICM/inference/plinder_set_0",
    "diffdock": f"{base_outdir}/DiffDock/inference/diffdock_plinder_benchmark_output_0",
    "diffdock_pocket_only": f"{base_outdir}/DiffDock/inference/diffdock_pocket_only_plinder_benchmark_output_1",
    # "chai-1": f"{base_outdir}/chai-lab/inference/chai-lab_{exp_name}",
    # "vina": f"{base_outdir}/Vina/inference/GNINA_plinder_output_0", 
    "gnina":f"{base_outdir}/GNINA/inference/GNINA_plinder_output_0",
    "surfdock":f"/home/aoxu/projects/docking_result/plinder_benchmark/SurfDock_docking_result",
    "boltz": f"{base_outdir}/boltz/inference/{exp_name}",
    # ...
}

# The folder containing the real (crystal) ligand and protein PDB:
#  PoseBench/data/posebusters_benchmark_set/<protein>/
DATA_DIR = "/home/aoxu/projects/PoseBench/data/plinder_set"


df_all = []

for approach in approaches:
    method_name = approach.get_name()
    base_outdir = BASE_DIRS[method_name]
    df_method = run_posebusters(
        approach,
        base_outdir=base_outdir,
        data_dir=DATA_DIR,
        top_n=5
    )
    df_all.append(df_method)
    df_method.to_csv(f"./notebooks/{method_name}_{exp_name}_results.csv", index=False)
df_combined = pd.concat(df_all, ignore_index=True)
print(df_combined.shape)
df_combined.head()

# Now you have a single DataFrame with columns:
#   - 'protein'
#   - 'method' (icm / diffdock / chai-1)
#   - 'rank' (1..5)
#   - 'score' (NaN or a float)
#   - 'rmsd', 'rmsd_≤_2å', etc. from PoseBusters
# You can groupby or pivot as you wish.
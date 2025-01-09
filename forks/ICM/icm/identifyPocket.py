import os
import multiprocessing
import gc
import time 
import glob
import tqdm
from dock import ICMBatchDocking

# Copy the current environment, then override or remove any ICM references
env = os.environ.copy()

# If you had an ICMHOME pointing to BrowserPro, remove or update it:
env.pop("ICMHOME", None)  # remove if it exists
env["ICMHOME"] = "/Applications/MolsoftICM64.app/Contents/Resources/icm"

# Configure these paths according to your setup
BASE_DIR = "/Users/aoxu/projects/DrugDiscovery/PoseBench/data/posebusters_benchmark_set"
SUBSET_DIR = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/DiffDock/inference/diffdock_posebusters_benchmark_output_orig_structure_1/"
ICM_EXECUTABLE = "/Applications/MolsoftICM64.app/Contents/MacOS/icm64"
ICM_DOCKING_DIR = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/ICM"
DOCKING_MAPS = "manual"
dock_scan_path="/Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan",

def prepare_docking():
    ICM_SCRIPT_TEMPLATE = "icm_docking_scripts_template/icm_pocket_identification_template.icm"
    docking_processor = ICMBatchDocking(
        base_dir=BASE_DIR,
        icm_script_template=ICM_SCRIPT_TEMPLATE,
        icm_executable=ICM_EXECUTABLE,
        subset_dir=SUBSET_DIR,
        icm_docking_dir=ICM_DOCKING_DIR,
        docking_maps=DOCKING_MAPS
    )

    STAGE = "pocket_identification"
    docking_processor.run_docking_preparation(num_conf=100, stage=STAGE)

def dock():
    # Create docking processor instance
    docking_processor = ICMBatchDocking(
        base_dir=BASE_DIR,
        icm_executable=ICM_EXECUTABLE,
        icm_dockscan_path=dock_scan_path,
        icm_script_template=None,
        subset_dir=SUBSET_DIR,
        icm_docking_dir=ICM_DOCKING_DIR,
        docking_maps=DOCKING_MAPS
    )
    
    # Run docking with specified number of conformations
    docking_processor.run_docking(num_conf=10, thorough=10)

if __name__ == "__main__":
    # prepare_docking()
    dock()
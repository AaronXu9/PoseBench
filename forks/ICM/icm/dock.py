import os
import subprocess
import shutil
import platform

class ICMBatchDocking:
    def __init__(self, data_dir, icm_script_template, icm_executable, icm_docking_dir, subset_dir=None, docking_maps="manual", icm_dockscan_path=None, repeat_index=0):
        """
        Initialize the batch docking processor
        
        Args:
            data_dir (str): Path to the base directory containing protein subdirectories
            icm_script_template (str): Path to the ICM script template filecm_executable (str): Path to ICM executable
        """
        self.data_dir = data_dir
        self.project_name = os.path.basename(os.path.dirname(data_dir))
        self.icm_script_template = icm_script_template
        self.icm_executable = icm_executable
        self.icm_dockscan_path = icm_dockscan_path
        self.subset_dir = subset_dir
        self.icm_docking_dir = icm_docking_dir
        self.docking_maps = docking_maps
        self.repeat_index = repeat_index
    
        # Copy current environment
        self.env = os.environ.copy()

        # Conditionally set ICMHOME only on macOS
        if platform.system() == "Darwin":
            self.env.pop("ICMHOME", None)
            self.env["ICMHOME"] = "/Applications/MolsoftICM64.app/Contents/Resources/icm"

            # For convenience, store the path to _dockScan
            # On Linux you might have /home/gpcr/icm/bin/icm64 or something similar
            # On macOS: "/Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan"
            self.icm_dockscan_path = "/Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan"

    def is_valid_protein_directory(self, dirname):
        """
        Validate if directory name follows the required format:
        - Must have exactly two parts separated by '_'
        - First part must contain a number
        
        Args:
            dirname (str): Directory name to validate
            
        Returns:
            bool: True if valid, False otherwise
        """
        # Check if directory name contains exactly one underscore
        parts = dirname.split('_')
        # if len(parts) != 2:
        #     return False
            
        # Check if first part contains at least one number
        if not any(char.isdigit() for char in parts[0]):
            return False
            
        return True
        
    def get_protein_ligand_pairs(self):
        """Get all protein-ligand pairs from the protein directory"""
        pairs = []
        for dirname in os.listdir(self.data_dir if self.subset_dir is None else self.subset_dir):   
            dir_path = os.path.join(self.data_dir, dirname)
            
            # Skip if not a directory or doesn't match required format
            if not os.path.isdir(dir_path) or not self.is_valid_protein_directory(dirname):
                print(f"Skipping '{dirname}': Not a valid protein directory format")
                continue
            
            protein_name = dirname  # Get protein name (e.g., 8BPL from 8BPL_RFO)
            ligand_name = dirname  # Full name (e.g., 8BPL_RFO)
            
            # Check if required files exist
            pdb_path = os.path.join(dir_path, f"{dirname}_protein.pdb")
            sdf_path = os.path.join(dir_path, f"{dirname}_ligand.sdf")
            
            if not os.path.exists(pdb_path):
                print(f"Skipping '{dirname}': Missing PDB file {pdb_path}")
                continue
                
            if not os.path.exists(sdf_path):
                print(f"Skipping '{dirname}': Missing SDF file {sdf_path}")
                continue
            
            pairs.append({
                'protein_name': protein_name,
                'ligand_name': ligand_name,
                'pdb_path': pdb_path,
                'sdf_path': sdf_path,
                'dir_path': dir_path
            })
            
        print(f"Found {len(pairs)} valid protein-ligand pairs")
        return sorted(pairs, key=lambda x: x['protein_name'])

    def create_modified_script(self, data_dir, template_path, script_out_path, icb_out_dir, protein_name):
        """
        Create a modified ICM script with the correct protein name.
        
        For proteins in the "plinder_set" dataset, if protein_name contains dots ('.'),
        the file loading commands (openFile lines) will use the original protein_name,
        while the rest of the script uses a 'safe' version (with dots removed) for naming.
        For other datasets, protein_name is used throughout.
        """
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # First, replace the dataset placeholder in the entire template.
        template_content = template_content.replace("DatasetNameHolder", self.project_name)
        template_content = template_content.replace("ICBOutDirHolder", icb_out_dir)
        template_content = template_content.replace("DataDirHolder", data_dir)
        # Check if the current dataset is "plinder_set"
        if self.project_name == "plinder_set":
            # Create a safe version of the protein name (remove '.')
            safe_protein_name = protein_name.replace('.', '_')
            modified_lines = []
            for line in template_content.splitlines():
                # If this is a file-loading line, use the original protein_name.
                if line.lstrip().startswith("openFile"):
                    modified_lines.append(line.replace("ProteinNameHolder", protein_name))
                else:
                    # For all other lines, use the safe version.
                    modified_lines.append(line.replace("ProteinNameHolder", safe_protein_name))
            modified_content = "\n".join(modified_lines)
        else:
            # For other datasets, use the original protein_name everywhere.
            modified_content = template_content.replace("ProteinNameHolder", protein_name)
        
        with open(script_out_path, 'w') as f:
            f.write(modified_content)


    def _run_icm_dockscan(
        self,
        protein_name: str,
        protein_project_dir: str,
        outdir_path: str,
        ligand_path: str,
        confs: int = 3,
        thorough: float = 3.0,
        effort: int = 3
    ):
        """
        Internal helper to run `_dockScan` for one protein-ligand pair.

        Args:
            protein_project_dir (str): The directory "protein_<protein_name>" containing the map files.
            confs (int): Number of conformations, e.g. 3 or 10
            thorough (float): Thoroughness level for docking
            effort (int): The "effort" parameter
        """
        # Example command:  icm64 _dockScan protein_5S8I_2LY -a -S confs=3 effort=3
        # We'll pass those arguments carefully
        icm_ligand_name = f"v{protein_name}_ligand"

        cmd = [
            self.icm_executable,       # e.g. "/Applications/MolsoftICM64.app/Contents/MacOS/icm64"
            self.icm_dockscan_path,    # e.g. "/Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan"
            "-s",                      # run silently
            "-a", f"thorough={thorough}.",
            f"name={icm_ligand_name}",
            f"input={ligand_path}",
            "-S", f"confs={confs}",
            f"protein_{protein_name}",
            f"output={os.path.join(outdir_path, protein_name, f'answers_{protein_name}.sdf')}",
        ]
        
        print("Running command (dockScan):", " ".join(cmd))
        
        # We can either change directory or specify `cwd`.
        # We'll do: cwd=protein_project_dir so that dockScan sees local map files
        subprocess.run(cmd, check=True, env=self.env, cwd=protein_project_dir)


    def run_docking(self, num_conf=10, thorough=3.0, effort=3):
        """Run docking for all protein-ligand pairs."""
        pairs = self.get_protein_ligand_pairs()
        if not pairs:
            print("No valid protein-ligand pairs found.")
            return
        sufix = 'p' if self.project_name == 'plinder_set' else ''
        for pair in pairs:
            protein_name = pair["protein_name"]  # e.g. "5S8I_2LY"
            
            # The project folder named "protein_<protein_name>"
            project_dir = os.path.join(
                self.icm_docking_dir,
                f"ICM_{self.docking_maps}_docking_maps",
                f"{sufix}_{protein_name}"
            )
            
            out_dir = os.path.join(
                self.icm_docking_dir, 
                "inference",
                f"{self.project_name}_docking_results",
            )
            os.makedirs(out_dir, exist_ok=True)

            # Ensure we have the project_dir
            if not os.path.isdir(project_dir):
                print(f"Project directory {project_dir} not found. Did you generate the maps first?")
                continue
            
            if self.project_name == 'posebuster_benchmark'  and protein_name < "7LCU_XTA":
                continue
            # Example command to run docking in that folder
            try:
                self._run_icm_dockscan(
                    protein_name=protein_name,
                    protein_project_dir=project_dir,
                    outdir_path=out_dir,
                    ligand_path=pair["sdf_path"],
                    confs=num_conf,
                    thorough=thorough,
                    effort=effort
                )
                print(f"Docking completed successfully for {protein_name}.")
            
            except subprocess.CalledProcessError as e:
                print(f"Error during docking for {protein_name}: {str(e)}")


    def run_docking_preparation(self, num_conf=100, stage="preprocessing"):
        """Run docking for all protein-ligand pairs"""
        pairs = self.get_protein_ligand_pairs()
        
        if not pairs:
            print("No valid protein-ligand pairs found. Please check your directory structure and file names.")
            return
            
        for pair in pairs:
            print(f"\nProcessing {pair['ligand_name']}...")
            
            # Create working directory for this pair
            # work_dir = os.path.join(pair['dir_path'], 'docking')
            work_dir = os.path.join(os.path.curdir, "icm_docking_scripts", stage)
            os.makedirs(work_dir, exist_ok=True)
            
            # Create modixfied script for this protein
            script_path = os.path.join(work_dir, f"{pair['ligand_name']}_script.icm")
            
            self.create_modified_script(data_dir, self.icm_script_template, script_path, self.config['icb_out_dir'], pair['protein_name'])
            
            # Copy necessary files to working directory
            # shutil.copy2(pair['pdb_path'], work_dir)
            # shutil.copy2(pair['sdf_path'], work_dir)
            
            # create the subdirectories for docking maps
            os.makedirs(os.path.join(self.icm_docking_dir, f"ICM_{self.docking_maps}_docking_maps", 
                                     f"{self.project_name}", f"protein_{pair['protein_name']}"), exist_ok=True)

            # Run ICM command
            try:
                # Run the modified ICM script
                cmd = [f"{self.icm_executable}",
                    #    "-p",
                    #    "/Applications/MolsoftICM64.app/Contents/Resources/icm",
                        "-g",
                        f"{script_path}"]
                print(f"Running command: {cmd}")
                # subprocess.run(cmd, shell=True, check=True)
                subprocess.run(cmd, check=True, env=self.env)
                
                print(f"Successfully completed docking for {pair['ligand_name']}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error processing {pair['ligand_name']}: {str(e)}")
                continue
        
# Example usage
if __name__ == "__main__":
    # Configure these paths according to your setup
    data_dir = "/home/aoxu/projects/PoseBench/data/posebusters_benchmark_set"
    SUBSET_DIR = "/home/aoxu/projects/PoseBench/forks/DiffDock/inference/diffdock_posebusters_benchmark_output_orig_structure_1/"
    ICM_SCRIPT_TEMPLATE = "icm_docking_scripts_template/protein_template_auto_pocket.icm"
    ICM_EXECUTABLE = "/home/aoxu/icm-3.9-4/icm64"
    ICM_DOCKING_DIR = "/home/aoxu/projects/PoseBench/forks/ICM"
    DOCKING_MAPS = "manual"
    dock_scan_path="/home/aoxu/icm-3.9-4/_dockScan",
    
    # Create docking processor instance
    docking_processor = ICMBatchDocking(
        data_dir=data_dir,
        icm_script_template=ICM_SCRIPT_TEMPLATE,
        icm_executable=ICM_EXECUTABLE,
        subset_dir=SUBSET_DIR,
        icm_docking_dir=ICM_DOCKING_DIR,
        docking_maps=DOCKING_MAPS
    )
    
    # Run docking with specified number of conformations
    docking_processor.run_docking(num_conf=10, thorough=10)
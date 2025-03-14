from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List
import os
import json
import yaml
import platform
import subprocess
import shutil
from pathlib import Path

import rootutils
rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)


@dataclass
class ICMDockingConfig:
    """Configuration class for ICM Batch Docking"""
    # Required parameters
    dataset: str
    data_dir: str
    icm_script_template: str
    icm_executable: str
    icm_docking_dir: str
    icm_map_dir: str

    # Optional parameters with defaults
    subset_dir: Optional[str] = None
    docking_maps: str = "manual"
    icm_dockscan_path: Optional[str] = None
    repeat_index: int = 0
    icb_out_dir: str = ""
    skip_existing: bool = False
    
    # Docking parameters
    docking_params: Dict[str, Any] = field(default_factory=lambda: {
        "num_conf": 10,
        "thorough": 3.0,
        "effort": 3
    })
    
    # Platform-specific configurations
    platform_configs: Dict[str, Dict[str, str]] = field(default_factory=lambda: {
        "Darwin": {
            "icm_executable": "/Applications/MolsoftICM64.app/Contents/Resources/icm",
            "icm_dockscan_path": "/Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan"
        },
        "Linux": {
            "icm_executable": "/home/aoxu/icm-3.9-4/icm64",  # Will use system default
            "icm_dockscan_path": "/home/aoxu/icm-3.9-4/_dockScan"  # Will be set based on executable path
        }
    })
    
    logging: Dict[str, Any] = field(default_factory=lambda: {
        "level": "INFO",
        "file": None,
        "console": True
    })

    def __post_init__(self):
        """Initialize derived values after the dataclass is instantiated"""
        # Set project name based on data_dir
        self.icb_out_dir = self.icb_out_dir + "_" + str(self.repeat_index)
        # Apply platform-specific configuration
        system = platform.system()
        if system in self.platform_configs:
            platform_config = self.platform_configs[system]
            
            # Set ICMHOME if specified for this platform
            if platform_config["icm_executable"]:
                self.icm_executable = platform_config["icm_executable"]
            
            # Set dockscan path if not explicitly provided
            if not self.icm_dockscan_path and platform_config["icm_dockscan_path"]:
                self.icm_dockscan_path = platform_config["icm_dockscan_path"]
    
    @classmethod
    def from_file(cls, config_path: str) -> 'ICMDockingConfig':
        """Load configuration from a YAML or JSON file"""
        with open(config_path, 'r') as f:
            if config_path.endswith('.yaml') or config_path.endswith('.yml'):
                config_data = yaml.safe_load(f)
            else:
                config_data = json.load(f)
        
        return cls(**config_data)
    
    def to_file(self, config_path: str) -> None:
        """Save configuration to a YAML or JSON file"""
        # Convert dataclass to dictionary
        config_dict = {k: v for k, v in self.__dict__.items()
                      if not k.startswith('_') and k != 'platform_configs'}
        
        with open(config_path, 'w') as f:
            if config_path.endswith('.yaml') or config_path.endswith('.yml'):
                yaml.dump(config_dict, f, default_flow_style=False)
            else:
                json.dump(config_dict, f, indent=2)


class ICMBatchDocking:
    def __init__(self, config: ICMDockingConfig):
        """
        Initialize the batch docking processor using a configuration object
        
        Args:
            config (ICMDockingConfig): Configuration object with all necessary parameters
        """
        self.config = config
        self.stage = None
        # Copy current environment
        self.env = os.environ.copy()
        # Initialize logger
        log_settings = config.logging
        log_file = log_settings.get('file')
        log_level = log_settings.get('level', 'INFO')
        
        # # Import and set up logger if the module is available
        # try:
        #     from icm_logger import ICMLogger
        #     self.logger = ICMLogger(
        #         log_file=log_file,
        #         console_level=log_level,
        #         file_level='DEBUG'
        #     )
        # except ImportError:
        #     # Fallback to print statements if custom logger is not available
        #     self.logger = None
        #     print(f"Initializing ICMBatchDocking with project: {config.dataset}")

        # Apply platform-specific environment settings
        if platform.system() == "Darwin" and config.icm_executable:
            self.env.pop("ICMHOME", None)
            self.env["ICMHOME"] = config.icm_executable
    
    def is_valid_protein_directory(self, dirname):
        """
        Validate if directory name follows the required format
        
        Args:
            dirname (str): Directory name to validate
            
        Returns:
            bool: True if valid, False otherwise
        """
        # Check if first part contains at least one number
        parts = dirname.split('_')
        if not any(char.isdigit() for char in parts[0]):
            return False
            
        return True
        
    def get_protein_ligand_pairs(self):
        """Get all protein-ligand pairs from the protein directory"""
        pairs = []
        data_source = self.config.subset_dir if self.config.subset_dir else self.config.data_dir
        
        for dirname in os.listdir(data_source):   
            dir_path = os.path.join(self.config.data_dir, dirname)
            
            # Skip if not a directory or doesn't match required format
            if not os.path.isdir(dir_path) or not self.is_valid_protein_directory(dirname):
                print(f"Skipping '{dirname}': Not a valid protein directory format")
                continue
            
            protein_name = dirname
            ligand_name = dirname
            
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
        """Create a modified ICM script with the correct protein name."""
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # First, replace the dataset placeholder in the entire template.
        template_content = template_content.replace("DatasetNameHolder", self.config.dataset)
        template_content = template_content.replace("ICBOutDirHolder", icb_out_dir)
        template_content = template_content.replace("DataDirHolder", data_dir)
        template_content = template_content.replace("MapsOutDirHolder", self.config.icm_map_dir)
        # Check if the current dataset is "plinder_set"
        if self.config.dataset == "plinder_set":
            # Create a safe version of the protein name (remove '.')
            safe_protein_name = protein_name.replace('.', '_')
            modified_lines = []
            for line in template_content.splitlines():
                # If this is a file-loading line, use the original protein_name.
                if line.lstrip().startswith("openFile") and self.stage == 'preprocessing':
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
        **docking_params
    ):
        """
        Internal helper to run _dockScan for one protein-ligand pair.
        
        Args:
            protein_name (str): Name of the protein
            protein_project_dir (str): The directory containing the map files
            outdir_path (str): Output directory path
            ligand_path (str): Path to the ligand SDF file
            **docking_params: Additional docking parameters (confs, thorough, effort)
        """
        # Get docking parameters with defaults
        confs = docking_params.get('confs', self.config.docking_params['num_conf'])
        thorough = docking_params.get('thorough', self.config.docking_params['thorough'])
        
        icm_ligand_name = f"v{protein_name}_ligand"
        output_file = os.path.join(outdir_path, protein_name, f'answers_{protein_name}.sdf')
        
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        cmd = [
            self.config.icm_executable,
            self.config.icm_dockscan_path,
            "-s",                      # run silently
            "-a", f"thorough={thorough}.",
            f"name={icm_ligand_name}",
            f"input={ligand_path}",
            "-S", f"confs={confs}",
            f"protein_{protein_name}",
            f"output={output_file}",
        ]
        
        print("Running command (dockScan):", " ".join(cmd))
        subprocess.run(cmd, check=True, env=self.env, cwd=protein_project_dir)

    def run_docking(self, **docking_params):
        """Run docking for all protein-ligand pairs with optional parameter overrides."""
        # Merge config params with provided overrides
        params = {**self.config.docking_params, **docking_params}
        
        pairs = self.get_protein_ligand_pairs()
        if not pairs:
            print("No valid protein-ligand pairs found.")
            return

        for pair in pairs:
            protein_name = pair["protein_name"]
            
            # Skip specific proteins for posebuster_benchmark if needed
            if self.config.dataset == 'posebuster_benchmark' and protein_name < "7LCU_XTA":
                continue
                
            # The project folder named "protein_<protein_name>"
            project_dir = os.path.join(
                self.config.icm_docking_dir,
                f"ICM_{self.config.docking_maps}_docking_maps",
                f"protein_{protein_name}"
            )
            
            out_dir = os.path.join(
                self.config.icm_docking_dir, 
                "inference",
                f"{self.config.dataset}_docking_results",
            )
            os.makedirs(out_dir, exist_ok=True)

            # Ensure we have the project_dir
            if not os.path.isdir(project_dir):
                print(f"Project directory {project_dir} not found. Did you generate the maps first?")
                continue
            
            try:
                self._run_icm_dockscan(
                    protein_name=protein_name,
                    protein_project_dir=project_dir,
                    outdir_path=out_dir,
                    ligand_path=pair["sdf_path"],
                    **params
                )
                print(f"Docking completed successfully for {protein_name}.")
            
            except subprocess.CalledProcessError as e:
                print(f"Error during docking for {protein_name}: {str(e)}")

    def run_docking_preparation(self, stage="preprocessing"):
        """Run docking preparation for all protein-ligand pairs"""
        self.stage = stage
        pairs = self.get_protein_ligand_pairs()
        
        if not pairs:
            print("No valid protein-ligand pairs found. Please check your directory structure and file names.")
            return
        
        # Create working directory for this pair
        work_dir = os.path.join(os.path.curdir, "forks/ICM",
                                "icm_docking_scripts", stage)
        os.makedirs(work_dir, exist_ok=True)
        
        for pair in pairs:
            print(f"\nProcessing {pair['ligand_name']}...")
            

            os.makedirs(work_dir, exist_ok=True)
            
            # Create modified script for this protein
            script_path = os.path.join(work_dir, f"{pair['ligand_name']}_script.icm")
            
            self.create_modified_script(
                self.config.data_dir, 
                self.config.icm_script_template, 
                script_path, 
                self.config.icb_out_dir, 
                pair['protein_name']
            )
            
            # Create the subdirectories for docking maps
            os.makedirs(
                os.path.join(
                    self.config.icm_docking_dir, 
                    f"ICM_{self.config.docking_maps}_docking_maps", 
                    f"{self.config.dataset}", 
                    f"protein_{pair['protein_name']}"
                ), 
                exist_ok=True
            )

            # Run ICM command
            try:
                cmd = [
                    self.config.icm_executable,
                    "-g",
                    script_path
                ]
                print(f"Running command: {cmd}")
                subprocess.run(cmd, check=True, env=self.env)
                
                print(f"Successfully completed docking for {pair['ligand_name']}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error processing {pair['ligand_name']}: {str(e)}")


# Example usage
if __name__ == "__main__":
    # Load configuration from file

    # config = ICMDockingConfig.from_file('cogligand_config/model/icm_inference.yaml')
    config = ICMDockingConfig.from_file('cogligand_config/data/icm_identify_pocket.yaml')
    
    # # Or create configuration programmatically
    # config = ICMDockingConfig(
    #     data_dir="/home/aoxu/projects/PoseBench/data/posebusters_benchmark_set",
    #     icm_script_template="icm_docking_scripts_template/protein_template_auto_pocket.icm",
    #     icm_executable="/home/aoxu/icm-3.9-4/icm64",
    #     icm_docking_dir="/home/aoxu/projects/PoseBench/forks/ICM",
    #     subset_dir="/home/aoxu/projects/PoseBench/forks/DiffDock/inference/diffdock_posebusters_benchmark_output_orig_structure_1/",
    #     docking_maps="manual",
    #     icm_dockscan_path="/home/aoxu/icm-3.9-4/_dockScan",
    # )
    
    # Save the configuration for future use
    # config.to_file('icm_docking_config.yaml')
    
    # Create docking processor with configuration
    docking_processor = ICMBatchDocking(config)
    
    # Run docking with default parameters from config
    docking_processor.run_docking_preparation(stage="identify_pocket")
    
    # Or override specific parameters
    # docking_processor.run_docking(num_conf=20, thorough=10.0)
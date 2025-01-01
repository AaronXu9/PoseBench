import os
import subprocess
import logging
from enum import Enum
from dataclasses import dataclass
from typing import Dict, List, Optional

class ProcessingStage(Enum):
    PREPROCESS = "preprocess"
    POCKET_FINDER = "pocket_finder"
    DOCK = "dock"

@dataclass
class StageConfig:
    required_files: List[str]
    output_files: List[str]
    script_name: Optional[str] = None

class ICMWorkflowManager:
    def __init__(self, 
                 base_dir: str, 
                 icm_executable: str,
                 dock_scan_path: str,
                 output_dir: str):
        """
        Initialize the workflow manager
        
        Args:
            base_dir (str): Base directory containing protein directories
            icm_executable (str): Path to ICM executable (icm64)
            dock_scan_path (str): Path to _dockScan script
            output_dir (str): Directory for docking outputs
        """
        self.base_dir = base_dir
        self.icm_executable = icm_executable
        self.dock_scan_path = dock_scan_path
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('icm_workflow.log'),
                logging.StreamHandler()
            ]
        )
        
        # Define stage configurations
        self.stage_configs = {
            ProcessingStage.PREPROCESS: StageConfig(
                script_name="preprocess.icm",
                required_files=["protein.pdb", "ligand.sdf"],
                output_files=["processed_protein.pdb"]
            ),
            ProcessingStage.POCKET_FINDER: StageConfig(
                script_name="pocket_finder.icm",
                required_files=["processed_protein.pdb"],
                output_files=["pocket.icm"]
            ),
            ProcessingStage.DOCK: StageConfig(
                required_files=["processed_protein.pdb", "ligand.sdf"],
                output_files=["answers.sdf"]
            )
        }

    def get_protein_name(self, dirname: str) -> str:
        """Extract protein name from directory name"""
        return dirname.split('_')[0]

    def run_docking(self, protein_dir: str, protein_name: str) -> bool:
        """
        Run docking using _dockScan command
        
        Args:
            protein_dir (str): Path to protein directory
            protein_name (str): Name of the protein (e.g., 5S8I)
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Construct paths
            ligand_sdf = os.path.join(protein_dir, f"{protein_name}_ligand.sdf")
            processed_protein = os.path.join(protein_dir, "processed_protein.pdb")
            output_sdf = os.path.join(self.output_dir, f"answers_{protein_name}.sdf")
            
            # Verify input files exist
            if not os.path.exists(ligand_sdf):
                logging.error(f"Missing ligand file: {ligand_sdf}")
                return False
            if not os.path.exists(processed_protein):
                logging.error(f"Missing processed protein file: {processed_protein}")
                return False
            
            # Construct docking command
            cmd = [
                self.icm_executable,
                self.dock_scan_path,
                "-s",
                "-a",
                f"thorough=3.",
                f"name=v{protein_name}_ligand",
                f"input={ligand_sdf}",
                "-S",
                "confs=5",
                processed_protein,
                f"outdir={self.output_dir}",
                f"output=answers_{protein_name}.sdf"
            ]
            
            # Log command
            logging.info(f"Running docking command: {' '.join(cmd)}")
            
            # Run docking
            process = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            # Verify output file exists
            if os.path.exists(output_sdf):
                logging.info(f"Successfully created docking output: {output_sdf}")
                return True
            else:
                logging.error(f"Docking completed but output file not found: {output_sdf}")
                return False
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Docking failed for {protein_name}: {e.stderr}")
            return False
        except Exception as e:
            logging.error(f"Unexpected error during docking for {protein_name}: {str(e)}")
            return False

    def process_protein(self, protein_dir: str) -> bool:
        """Process a single protein through all stages"""
        protein_name = os.path.basename(protein_dir).split('_')[0]
        logging.info(f"\nProcessing protein: {protein_name}")
        
        # Run preprocessing if needed
        if not os.path.exists(os.path.join(protein_dir, "processed_protein.pdb")):
            if not self.run_script(protein_dir, ProcessingStage.PREPROCESS):
                return False
        
        # Run docking
        return self.run_docking(protein_dir, protein_name)

    def process_all_proteins(self) -> Dict[str, bool]:
        """Process all proteins in the base directory"""
        results = {}
        
        for dirname in os.listdir(self.base_dir):
            protein_dir = os.path.join(self.base_dir, dirname)
            
            # Skip if not a directory or doesn't match naming format
            if not os.path.isdir(protein_dir) or '_' not in dirname:
                continue
            
            success = self.process_protein(protein_dir)
            results[dirname] = success
            
        return results

# Example usage
if __name__ == "__main__":
    # Configure paths
    BASE_DIR = "/path/to/protein/directory"
    ICM_EXECUTABLE = "/home/aoxu/icm-3.9-4/icm64"
    DOCK_SCAN_PATH = "/home/aoxu/icm-3.9-4/_dockScan"
    OUTPUT_DIR = "/path/to/output/directory"
    
    # Create workflow manager
    workflow = ICMWorkflowManager(
        base_dir=BASE_DIR,
        icm_executable=ICM_EXECUTABLE,
        dock_scan_path=DOCK_SCAN_PATH,
        output_dir=OUTPUT_DIR
    )
    
    # Process all proteins
    results = workflow.process_all_proteins()
    
    # Print summary
    total = len(results)
    successful = sum(1 for success in results.values() if success)
    print(f"\nProcessing complete: {successful}/{total} successful")
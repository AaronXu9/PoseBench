import os
import re
import json
from abc import ABC, abstractmethod
from typing import List
import pandas as pd
from rdkit import Chem
from posebusters.posebusters import PoseBusters

class DockingApproach(ABC):
    @abstractmethod
    def get_name(self) -> str:
        """A short identifier for this method (e.g. 'icm', 'diffdock', 'chai')."""
        pass

    @abstractmethod
    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Return up to top_n SDF file paths for that protein_dir, 
        in descending or ascending 'rank' order, as appropriate.
        """
        pass

    def parse_score(self, sdf_path: str) -> float:
        """
        If this approach has a numeric score to parse, override this method.
        If there's no numeric score, return None or float('nan').
        """
        return float('nan')  # default: no score
    
class ICMApproach(DockingApproach):
    def get_name(self) -> str:
        return "icm"
    
    def top_n_conformation(self, protein_dir: str):
        protein_name = os.path.basename(protein_dir)
        input_sdf = os.path.join(protein_dir, f"answers_{protein_name}.sdf")
        """select the top n conformations from the input sdf file"""
        try: 
            suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)

            mols = []
            for mol in suppl:
                if mol is not None:
                    score = float(mol.GetProp("Score"))  # Adjust property name as needed
                    mols.append((mol, score))

            # Sort ascending (lowest score first)
            mols.sort(key=lambda x: x[1])

            # Write rank1.sdf
            if mols:
                w1 = Chem.SDWriter(os.path.join(protein_dir, "rank1.sdf"))
                w1.write(mols[0][0])
                w1.close()

            # Write rank2.sdf to rank5.sdf
            for i in range(1, 5):
                if i < len(mols):
                    writer = Chem.SDWriter(os.path.join(protein_dir, f"rank{i+1}.sdf"))
                    writer.write(mols[i][0])
                    writer.close()
        except:
            print(f"Could not find {input_sdf}")

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Assume ICM wrote 'rank1.sdf', 'rank2.sdf', etc.
        We'll list all files matching 'rankX.sdf', sort by X, return top_n.
        """
        self.top_n_conformation(protein_dir)
        all_files = os.listdir(protein_dir)
        sdf_files = [f for f in all_files if f.startswith("rank") and f.endswith(".sdf")]
        
        def extract_rank(fname: str) -> int:
            # "rank(\d+).sdf"
            match = re.match(r"rank(\d+)\.sdf", fname)
            if match:
                return int(match.group(1))
            return 999999
        
        sdf_files.sort(key=extract_rank)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]

    def parse_score(self, sdf_path: str) -> float:
        """
        ICM stores the docking score in the SDF property 'Score'.
        We'll read the single pose from the file and extract that property.
        """
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        for mol in suppl:
            if mol is not None and "Score" in mol.GetPropNames():
                return float(mol.GetProp("Score"))
        return float('nan')

class DiffDockApproach(DockingApproach):
    def get_name(self) -> str:
        return "diffdock"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        DiffDock outputs might be named 'rank1.sdf' or 
        'rank2_confidence-0.16.sdf', etc.
        We'll parse 'rank(\d+)' to get the rank, sort, and return top_n.
        """
        all_files = os.listdir(protein_dir)
        sdf_files = [
            f for f in all_files 
            if f.startswith("rank") and f.endswith(".sdf")
        ]

        def extract_rank(fname: str) -> int:
            match = re.search(r"rank(\d+)", fname)
            if match:
                return int(match.group(1))
            return 999999
        
        sdf_files.sort(key=extract_rank)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]

    def parse_score(self, sdf_path: str) -> float:
        """
        DiffDock doesn't store 'Score' in the SDF properties.
        Instead, there's a confidence value in the filename
        like 'rank2_confidence-0.16.sdf'.
        We'll parse that out. If not present, return NaN.
        """
        fname = os.path.basename(sdf_path)
        # Look for something like _confidence-0.16
        match = re.search(r"_confidence-([\d\.]+)", fname)
        if match:
            conf_str = match.group(1)
            # Remove trailing '.' if present, e.g. "1.40." => "1.40"
            conf_str = conf_str.rstrip('.')
            try:
                return float(conf_str)
            except ValueError:
                print(f"[WARNING] Could not parse confidence '{conf_str}' from {fname}, returning NaN.")
                return float('nan')
        return float('nan')
    
class DiffDockPocketApproach(DockingApproach):
    def get_name(self) -> str:
        return "diffdock_pocket_only"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        DiffDock doesn't store 'Score' in the SDF properties.
        Instead, there's a confidence value in the filename
        like 'rank2_confidence-0.16.sdf'.
        We'll parse that out. If not present, return NaN.
        """
        all_files = os.listdir(protein_dir)
        sdf_files = [
            f for f in all_files 
            if f.startswith("rank") and f.endswith(".sdf")
        ]

        def extract_rank(fname: str) -> int:
            match = re.search(r"rank(\d+)", fname)
            if match:
                return int(match.group(1))
            return 999999
        
        sdf_files.sort(key=extract_rank)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]

    def parse_score(self, sdf_path: str) -> float:
        """
        DiffDock doesn't store 'Score' in the SDF properties.
        Instead, there's a confidence value in the filename
        like 'rank2_confidence-0.16.sdf'.
        We'll parse that out. If not present, return NaN.
        """
        fname = os.path.basename(sdf_path)
        # Look for something like _confidence-0.16
        match = re.search(r"_confidence-([\d\.]+)", fname)
        if match:
            conf_str = match.group(1)
            # Remove trailing '.' if present, e.g. "1.40." => "1.40"
            conf_str = conf_str.rstrip('.')
            try:
                return float(conf_str)
            except ValueError:
                print(f"[WARNING] Could not parse confidence '{conf_str}' from {fname}, returning NaN.")
                return float('nan')
        return float('nan')
    
class ChaiApproach(DockingApproach):
    def get_name(self) -> str:
        return "chai-1"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Example filenames: pred.model_idx_0_ligand_aligned.sdf,
                          pred.model_idx_1_ligand_aligned.sdf, ...
        We'll parse the model_idx to define an order. 
        """
        all_files = os.listdir(protein_dir)
        sdf_files = [
            f for f in all_files
            if f.startswith("pred.model_idx_") and f.endswith("_ligand_aligned.sdf")
        ]

        def extract_model_idx(fname: str) -> int:
            m = re.search(r"model_idx_(\d+)", fname)
            if m:
                return int(m.group(1))
            return 999999

        sdf_files.sort(key=extract_model_idx)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]

    # No score or confidence stored
    def parse_score(self, sdf_path: str) -> float:
        return float('nan')
    
class VinaApproach(DockingApproach):
    def get_name(self) -> str:
        return "vina"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Vina outputs are named like '5SAK_ZRY_pose3_score-7.97.sdf'
        We'll parse 'pose(\d+)' to get the rank, sort, and return top_n.
        """
        all_files = os.listdir(protein_dir)
        sdf_files = [f for f in all_files if '_pose' in f and f.endswith('.sdf')]

        def extract_pose_num(fname: str) -> int:
            match = re.search(r'pose(\d+)', fname)
            if match:
                return int(match.group(1))
            return 999999
        
        sdf_files.sort(key=extract_pose_num)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]


    def parse_score(self, sdf_path: str) -> float:
        """
        Score is in filename like '5SAK_ZRY_pose3_score-7.97.sdf'
        """
        fname = os.path.basename(sdf_path)
        match = re.search(r'score-([\-\d\.]+)', fname)
        if match:
            score_str = match.group(1)
            score_str = score_str.rstrip('.')
            try:
                return float(score_str)
            except ValueError:
                print(f"[WARNING] Could not parse score '{score_str}' from {fname}")
                return float('nan')
        return float('nan')
    
class GninaApproach(DockingApproach):
    def get_name(self) -> str:
        return "gnina"
    
    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        GNINA outputs are named like 'rank1_score-7.97.sdf'
        We'll parse 'rank(\d+)' to get the rank, sort, and return top_n.
        """
        all_files = os.listdir(protein_dir)
        sdf_files = [f for f in all_files if f.startswith('rank') and f.endswith('.sdf')]
        
        def extract_rank_num(fname: str) -> int:
            match = re.search(r'rank(\d+)', fname)
            if match:
                return int(match.group(1))
            return 999999
        
        sdf_files.sort(key=extract_rank_num)
        return [os.path.join(protein_dir, f) for f in sdf_files[:top_n]]
    
    def parse_score(self, sdf_path: str) -> float:
        """
        Score is in filename like 'rank1_score-7.97.sdf'
        """
        fname = os.path.basename(sdf_path)
        match = re.search(r'score([\-\d\.]+)', fname)
        if match:
            score_str = match.group(1)
            score_str = score_str.rstrip('.')
            try:
                return float(score_str)
            except ValueError:
                print(f"[WARNING] Could not parse score '{score_str}' from {fname}")
                return float('nan')
        return float('nan')
    
class SurfDockApproach(DockingApproach):
    def get_name(self) -> str:
        return "surfdock"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Rename the directory from something like:
          5SAK_ZRY_protein_8A_5SAK_ZRY_ligand
        to:
          5SAK_ZRY
        Then locate any SDF files and return the top_n entries (sorted if needed).
        """
        # Example directory rename (adjust as necessary)
        new_dir = re.sub(r"_protein_.*_ligand$", "", protein_dir)
        if os.path.exists(protein_dir) and not os.path.exists(new_dir):
            os.rename(protein_dir, new_dir)
        
        all_files = os.listdir(new_dir)
        sdf_files = [f for f in all_files if f.endswith(".sdf")]

        # Sort by rank if needed, for example:
        # parse something like ..._rank_11_... 
        def extract_rank(fname: str) -> int:
            match = re.search(r"rank_(\d+)", fname)
            if match:
                return int(match.group(1))
            return 999999

        sdf_files.sort(key=extract_rank)
        return [os.path.join(new_dir, f) for f in sdf_files[:top_n]]

    def parse_score(self, sdf_path: str) -> float:
        """
        Extract RMSD and confidence from filenames like:
          5SAK_5SAK_ZRY_A_404_5SAK_ZRY_ligand.sdf_file_inner_idx_0_sample_idx_19_rank_11_rmsd_0.38968_confidence_280.8584.sdf
        Parse both RMSD and confidence, return confidence as the 'score'.
        """
        fname = os.path.basename(sdf_path)
        
        # Parse RMSD
        match_rmsd = re.search(r"rmsd_([\-\d\.]+)", fname)
        if match_rmsd:
            try:
                rmsd = float(match_rmsd.group(1))
            except ValueError:
                print(f"[WARNING] Could not parse RMSD from {fname}")
                rmsd = float('nan')
        else:
            rmsd = float('nan')
        
        # Parse confidence
        match_conf = re.search(r"confidence_([\-\d\.]+)", fname)
        if match_conf:
            try:
                confidence = float(match_conf.group(1).rstrip('.'))
            except ValueError:
                print(f"[WARNING] Could not parse confidence from {fname}")
                confidence = float('nan')
        else:
            confidence = float('nan')
        
        print(f"Parsed from {fname}: RMSD={rmsd}, confidence={confidence}")
        # Return confidence as the 'score' for consistency
        return confidence

class BoltzApproach(DockingApproach):
    def get_name(self) -> str:
        return "boltz"

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Boltz outputs are in predictions/<protein_name>/<protein_name>_model_X.pdb
        We'll convert PDB files to SDF format and return paths to top_n files.
        """
        # Navigate to the predictions subdirectory
        protein_name = os.path.basename(protein_dir).replace("boltz_results_", "")
        pred_dir = os.path.join(protein_dir, "predictions", protein_name)
        
        if not os.path.exists(pred_dir):
            print(f"[WARNING] Predictions directory not found: {pred_dir}")
            return []
        
        # Find all PDB model files
        all_files = os.listdir(pred_dir)
        pdb_files = [f for f in all_files if f.endswith('.pdb') and '_model_' in f]
        
        def extract_model_num(fname: str) -> int:
            match = re.search(r'_model_(\d+)\.pdb', fname)
            if match:
                return int(match.group(1))
            return 999999
        
        # Sort by model number
        pdb_files.sort(key=extract_model_num)
        
        # Convert PDB files to SDF and return paths
        sdf_paths = []
        for pdb_file in pdb_files[:top_n]:
            pdb_path = os.path.join(pred_dir, pdb_file)
            sdf_path = self._convert_pdb_to_sdf(pdb_path)
            if sdf_path:
                sdf_paths.append(sdf_path)
        
        return sdf_paths

    def _convert_pdb_to_sdf(self, pdb_path: str) -> str:
        """
        Extract ligand coordinates from PDB file (HETATM records) and convert to SDF.
        """
        try:
            # Read PDB file and extract HETATM records for ligands
            ligand_lines = []
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('HETATM') and 'LIG' in line:
                        ligand_lines.append(line)
            
            if not ligand_lines:
                print(f"[WARNING] No ligand HETATM records found in {pdb_path}")
                return None
            
            # Create a temporary PDB file with only ligand atoms
            temp_pdb = pdb_path.replace('.pdb', '_ligand_temp.pdb')
            with open(temp_pdb, 'w') as f:
                f.writelines(ligand_lines)
                f.write('END\n')
            
            # Convert to SDF using RDKit
            sdf_path = pdb_path.replace('.pdb', '_ligand.sdf')
            
            # Try to read the ligand PDB with RDKit
            mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False)
            if mol is not None:
                writer = Chem.SDWriter(sdf_path)
                writer.write(mol)
                writer.close()
                
                # Clean up temporary file
                os.remove(temp_pdb)
                return sdf_path
            else:
                print(f"[WARNING] Could not parse ligand from {pdb_path}")
                os.remove(temp_pdb)
                return None
                
        except Exception as e:
            print(f"[ERROR] Failed to convert {pdb_path} to SDF: {str(e)}")
            return None

    def parse_score(self, sdf_path: str) -> float:
        """
        Parse confidence score from the corresponding JSON file.
        Also try to get affinity value as an alternative score.
        """
        try:
            # Get the base path and model number from SDF path
            base_path = sdf_path.replace('_ligand.sdf', '')
            pred_dir = os.path.dirname(sdf_path)
            
            # Extract model number
            match = re.search(r'_model_(\d+)_ligand\.sdf', sdf_path)
            if not match:
                return float('nan')
            
            model_num = match.group(1)
            
            # Look for confidence JSON file
            json_pattern = f"confidence_*_model_{model_num}.json"
            json_files = [f for f in os.listdir(pred_dir) if re.match(json_pattern.replace('*', '.*'), f)]
            
            if json_files:
                json_path = os.path.join(pred_dir, json_files[0])
                with open(json_path, 'r') as f:
                    data = json.load(f)
                    # Return confidence score as primary metric
                    return data.get('confidence_score', float('nan'))
            
            # If no confidence file, try affinity file
            affinity_files = [f for f in os.listdir(pred_dir) if f.startswith('affinity_') and f.endswith('.json')]
            if affinity_files:
                affinity_path = os.path.join(pred_dir, affinity_files[0])
                with open(affinity_path, 'r') as f:
                    data = json.load(f)
                    # Return affinity prediction value
                    return data.get('affinity_pred_value', float('nan'))
            
            return float('nan')
            
        except Exception as e:
            print(f"[ERROR] Failed to parse score from {sdf_path}: {str(e)}")
            return float('nan')
import os
import re
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

    def list_top_n_files(self, protein_dir: str, top_n: int) -> List[str]:
        """
        Assume ICM wrote 'rank1.sdf', 'rank2.sdf', etc.
        We'll list all files matching 'rankX.sdf', sort by X, return top_n.
        """
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


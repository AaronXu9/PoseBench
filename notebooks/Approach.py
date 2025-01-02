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
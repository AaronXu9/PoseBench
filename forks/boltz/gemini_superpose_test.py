import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS
import traceback
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolTransforms # <-- Import the new, required module
import traceback
import numpy as np

def calculate_ligand_rmsd_definitive(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates the ligand RMSD after correctly superimposing the protein backbones.

    This version uses the standard, documented API call for rdMolAlign.CalcRMS,
    which is correct for a clean conda-forge RDKit installation.
    """
    try:
        # Steps 1-5: These are all correct and remain unchanged.
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_atoms = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }
        
        protein_atom_map = []
        for pred_atom in pred_complex.GetAtoms():
            if is_ca(pred_atom):
                pred_info = pred_atom.GetPDBResidueInfo()
                pred_key = (pred_info.GetChainId().strip(), pred_info.GetResidueNumber())
                ref_idx = ref_ca_atoms.get(pred_key)
                if ref_idx is not None:
                    protein_atom_map.append((pred_atom.GetIdx(), ref_idx))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms.")
        
        rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        pred_lig_atom_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        mcs_result = rdFMCS.FindMCS([ref_lig, pred_lig], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs_result.numAtoms == 0:
            raise RuntimeError("No common substructure found between reference and predicted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = tuple(ref_lig.GetSubstructMatch(mcs_mol))
        pred_match = tuple(pred_lig.GetSubstructMatch(mcs_mol))

        if not ref_match or not pred_match:
            raise RuntimeError("Failed to match MCS to ligands.")

        ligand_atom_map = list(zip(pred_match, ref_match))
        
        #
        # --- THE CORRECT CALL FOR A STANDARD RDKIT INSTALLATION ---
        #
        ligand_rmsd = rdMolAlign.CalcRMS(
            prbMol=pred_lig,
            refMol=ref_lig,
            atomMap=ligand_atom_map
        )

        return {
            'ligand_rmsd': ligand_rmsd,
            'success': True,
            'error': None
        }

    except Exception as e:
        return {
            'ligand_rmsd': None,
            'success': False,
            'error': f"{type(e).__name__}: {e}",
            'traceback': traceback.format_exc()
        }

def calculate_ligand_rmsd_final_attempt(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates the ligand RMSD after correctly superimposing the protein backbones.

    Final Attempt: This version uses a hyper-literal interpretation of the C++ error
    signature for rdMolAlign.CalcRMS, providing the first four arguments positionally
    and using the 'map' keyword.
    """
    try:
        # Steps 1-5, which are known to be working correctly, remain unchanged.
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_atoms = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }
        
        protein_atom_map = []
        for pred_atom in pred_complex.GetAtoms():
            if is_ca(pred_atom):
                pred_info = pred_atom.GetPDBResidueInfo()
                pred_key = (pred_info.GetChainId().strip(), pred_info.GetResidueNumber())
                ref_idx = ref_ca_atoms.get(pred_key)
                if ref_idx is not None:
                    protein_atom_map.append((pred_atom.GetIdx(), ref_idx))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms.")
        
        rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        pred_lig_atom_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        mcs_result = rdFMCS.FindMCS([ref_lig, pred_lig], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs_result.numAtoms == 0:
            raise RuntimeError("No common substructure found between reference and predicted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = tuple(ref_lig.GetSubstructMatch(mcs_mol))
        pred_match = tuple(pred_lig.GetSubstructMatch(mcs_mol))

        if not ref_match or not pred_match:
            raise RuntimeError("Failed to match MCS to ligands.")

        ligand_atom_map = list(zip(pred_match, ref_match))
        
        #
        # --- THE CRITICAL LINE ---
        # Call CalcRMS providing the arguments exactly as the C++ signature dictates:
        # (molecule, molecule, int, int, map=...)
        #
        ligand_rmsd = rdMolAlign.CalcRMS(
            pred_lig,        # Positional Arg 1: prbMol
            ref_lig,         # Positional Arg 2: refMol
            -1,              # Positional Arg 3: prbId
            -1,              # Positional Arg 4: refId
            map=ligand_atom_map # Keyword Arg 5: map
        )

        return {
            'ligand_rmsd': ligand_rmsd,
            'success': True,
            'error': None
        }

    except Exception as e:
        return {
            'ligand_rmsd': None,
            'success': False,
            'error': f"{type(e).__name__}: {e}",
            'traceback': traceback.format_exc()
        }

def calculate_ligand_rmsd_manual(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates ligand RMSD after protein superposition.

    This definitive version bypasses the faulty RDKit RMSD functions and
    calculates the RMSD manually using NumPy to guarantee a result.
    """
    try:
        # Steps 1-5: Loading, Protein Superposition, Ligand Extraction, and MCS.
        # This logic is sound and remains unchanged.
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_atoms = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }
        
        protein_atom_map = []
        for pred_atom in pred_complex.GetAtoms():
            if is_ca(pred_atom):
                pred_info = pred_atom.GetPDBResidueInfo()
                pred_key = (pred_info.GetChainId().strip(), pred_info.GetResidueNumber())
                ref_idx = ref_ca_atoms.get(pred_key)
                if ref_idx is not None:
                    protein_atom_map.append((pred_atom.GetIdx(), ref_idx))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms.")
        
        rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        pred_lig_atom_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        mcs_result = rdFMCS.FindMCS([ref_lig, pred_lig], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs_result.numAtoms == 0:
            raise RuntimeError("No common substructure found between reference and predicted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = tuple(ref_lig.GetSubstructMatch(mcs_mol))
        pred_match = tuple(pred_lig.GetSubstructMatch(mcs_mol))

        if not ref_match or not pred_match:
            raise RuntimeError("Failed to match MCS to ligands.")

        # This map is correct: [(pred_idx_1, ref_idx_1), (pred_idx_2, ref_idx_2), ...]
        ligand_atom_map = list(zip(pred_match, ref_match))
        
        #
        # --- THE MANUAL RMSD CALCULATION ---
        #
        
        # 6. Get coordinates for the mapped atoms from each ligand
        ref_conf = ref_lig.GetConformer()
        pred_conf = pred_lig.GetConformer()

        squared_deviations = []
        for pred_idx, ref_idx in ligand_atom_map:
            pred_pos = np.array(pred_conf.GetAtomPosition(pred_idx))
            ref_pos = np.array(ref_conf.GetAtomPosition(ref_idx))
            
            # Calculate the squared distance between the two atoms
            dist_sq = np.sum((pred_pos - ref_pos) ** 2)
            squared_deviations.append(dist_sq)

        # 7. Calculate the final RMSD value
        if not squared_deviations:
            raise ValueError("No atom pairs to calculate RMSD.")
            
        mean_squared_error = np.mean(squared_deviations)
        ligand_rmsd = np.sqrt(mean_squared_error)

        return {
            'ligand_rmsd': ligand_rmsd,
            'n_ligand_atoms_matched': len(ligand_atom_map),
            'success': True,
            'error': None
        }

    except Exception as e:
        return {
            'ligand_rmsd': None,
            'n_ligand_atoms_matched': None,
            'success': False,
            'error': f"{type(e).__name__}: {e}",
            'traceback': traceback.format_exc()
        }

def calculate_ligand_rmsd_symmetric(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates a symmetry-corrected ligand RMSD after protein superposition.

    This version uses a specific call to rdMolAlign.CalcRMS that is tailored
    to the user's environment and leverages RDKit's internal, symmetry-aware
    atom matching.
    """
    try:
        # Steps 1-4: Loading, Protein Superposition, and Ligand Extraction.
        # This logic is sound and remains unchanged.
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_atoms = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }
        
        protein_atom_map = []
        for pred_atom in pred_complex.GetAtoms():
            if is_ca(pred_atom):
                pred_info = pred_atom.GetPDBResidueInfo()
                pred_key = (pred_info.GetChainId().strip(), pred_info.GetResidueNumber())
                ref_idx = ref_ca_atoms.get(pred_key)
                if ref_idx is not None:
                    protein_atom_map.append((pred_atom.GetIdx(), ref_idx))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms.")
        
        # This superimposes the complex and transforms the coordinates of pred_complex
        rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        # Extract the ligand from the *now aligned* complex
        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        pred_lig_atom_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        #
        # --- THE SYMMETRY-AWARE RMSD CALCULATION ---
        #
        
        # 5. Calculate the RMSD using RDKit's internal atom-matcher.
        # This call provides the positional integer arguments required by your
        # environment and allows the function to handle atom mapping and symmetry.
        # This function does NOT perform a new alignment.
        ligand_rmsd = rdMolAlign.CalcRMS(
            pred_lig,        # Positional Arg 1: prbMol
            ref_lig,         # Positional Arg 2: refMol
            -1,              # Positional Arg 3: prbId
            -1               # Positional Arg 4: refId
        )

        return {
            'ligand_rmsd': ligand_rmsd,
            'success': True,
            'error': None
        }

    except Exception as e:
        return {
            'ligand_rmsd': None,
            'success': False,
            'error': f"{type(e).__name__}: {e}",
            'traceback': traceback.format_exc()
        }
    # --- Final & Correct Script: MCS-Templated RMSD ---

def calculate_all_rmsd_metrics(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates and returns all relevant RMSD metrics to give a complete
    picture of model performance.
    """
    try:
        # 1. Load structures and create the unified reference complex
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])
        n_prot_atoms = ref_prot.GetNumAtoms()
        ref_complex = Chem.CombineMols(ref_prot, ref_lig)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)

        # 2. Perform and verify the protein alignment
        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_map = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_complex.GetAtoms() if is_ca(a) and a.GetIdx() < n_prot_atoms
        }
        
        protein_atom_map = []
        for atom in pred_complex.GetAtoms():
            if is_ca(atom):
                info = atom.GetPDBResidueInfo()
                key = (info.GetChainId().strip(), info.GetResidueNumber())
                if key in ref_ca_map:
                    protein_atom_map.append((atom.GetIdx(), ref_ca_map[key]))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms.")
        
        protein_rmsd = AllChem.AlignMol(prbMol=pred_complex, refMol=ref_complex, atomMap=protein_atom_map)

        # 3. Extract and sanitize ligands from the ALIGNED complexes
        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        
        pred_lig_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues]
        if not pred_lig_indices: raise RuntimeError("No ligand atoms found in predicted complex.")
        pred_lig_mol = Chem.PathToSubmol(pred_complex, pred_lig_indices)
        Chem.SanitizeMol(pred_lig_mol)

        ref_lig_indices = list(range(n_prot_atoms, ref_complex.GetNumAtoms()))
        ref_lig_mol = Chem.PathToSubmol(ref_complex, ref_lig_indices)

        # 4. Calculate Pose RMSD (manual, non-symmetric)
        mcs_result = rdFMCS.FindMCS([ref_lig_mol, pred_lig_mol], timeout=10)
        if mcs_result.numAtoms == 0: raise RuntimeError("MCS failed between extracted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = ref_lig_mol.GetSubstructMatch(mcs_mol)
        pred_match = pred_lig_mol.GetSubstructMatch(mcs_mol)

        ref_coords = np.array([ref_lig_mol.GetConformer().GetAtomPosition(i) for i in ref_match])
        pred_coords = np.array([pred_lig_mol.GetConformer().GetAtomPosition(i) for i in pred_match])
        
        diff = ref_coords - pred_coords
        pose_rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))

        # 5. Calculate Minimum RMSD (symmetry-aware, ignores pose)
        # This performs a new alignment of only the ligands.
        minimum_rmsd = AllChem.GetBestRMS(pred_lig_mol, ref_lig_mol)

        return {
            'pose_rmsd': pose_rmsd,
            'minimum_rmsd': minimum_rmsd,
            'protein_alignment_rmsd': protein_rmsd,
            'n_atoms_matched': mcs_result.numAtoms,
            'success': True
        }

    except Exception as e:
        return { 'success': False, 'error': str(e), 'traceback': traceback.format_exc() }

def calculate_ligand_rmsd_final_answer(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates the true ligand RMSD by first locating the ligand in the
    predicted complex via a substructure search, making the process robust
    to malformed PDB files from ML models.
    """
    try:
        # 1. Load the reference ligand (our search query) and the reference protein.
        ref_lig_template = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        
        # 2. Load the predicted complex and sanitize it to enable substructure search.
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        Chem.SanitizeMol(pred_complex)

        # --- THE ULTIMATE FIX: Find the ligand via substructure search ---
        # This is guaranteed to find the correct atoms if they exist.
        pred_lig_indices = pred_complex.GetSubstructMatch(ref_lig_template)

        if not pred_lig_indices:
            raise RuntimeError("CRITICAL FAILURE: The reference ligand's chemical structure could not be found within the predicted complex PDB file.")

        # 4. Now that we know the ligand atoms, we can align the proteins.
        # The protein is everything EXCEPT the identified ligand atoms.
        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        # Create a lookup map for the reference protein C-alphas
        ref_ca_map = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }

        # Find matching C-alphas in the prediction, ensuring they aren't part of the ligand
        protein_atom_map = []
        for atom in pred_complex.GetAtoms():
            if atom.GetIdx() in pred_lig_indices:
                continue # Skip ligand atoms
            
            if is_ca(atom):
                info = atom.GetPDBResidueInfo()
                key = (info.GetChainId().strip(), info.GetResidueNumber())
                if key in ref_ca_map:
                    protein_atom_map.append((atom.GetIdx(), ref_ca_map[key]))

        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms for alignment.")

        # 5. Align the predicted complex to the reference protein
        # We align to ref_prot, not a combined complex, for simplicity.
        protein_rmsd = AllChem.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)
        print(f"  [DIAGNOSTIC] Protein C-alpha alignment RMSD: {protein_rmsd:.3f} Å")


        # 6. Get coordinates from the aligned prediction and the original reference ligand.
        pred_conf = pred_complex.GetConformer()
        ref_conf = ref_lig_template.GetConformer()
        
        # The substructure match gives us indices for pred_complex. The reference indices are just 0, 1, 2...
        pred_coords = np.array([pred_conf.GetAtomPosition(i) for i in pred_lig_indices])
        ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(ref_lig_template.GetNumAtoms())])
        
        # 7. Calculate the final manual RMSD.
        diff = ref_coords - pred_coords
        ligand_rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
        
        return {
            'ligand_rmsd': ligand_rmsd,
            'protein_alignment_rmsd': protein_rmsd,
            'n_atoms_matched': len(pred_lig_indices),
            'success': True
        }

    except Exception as e:
        return { 'success': False, 'error': str(e), 'traceback': traceback.format_exc() }


def calculate_ligand_rmsd_production(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    This definitive script calculates the true ligand pose RMSD. It is robust
    to "dirty" PDB files containing extraneous HETATM records by performing a
    substructure search to guarantee chemical identity before calculation.
    """
    try:
        # 1. Load reference structures.
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        
        # 2. Load the predicted complex.
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)

        if not all([ref_lig, ref_prot, pred_complex]):
            raise ValueError("Failed to load one or more input files.")

        # --- THE ULTIMATE LIGAND IDENTIFICATION METHOD ---
        # a. First, get all HETATM records as a starting point.
        all_het_indices = [a.GetIdx() for a in pred_complex.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetIsHeteroAtom()]
        if not all_het_indices:
            raise RuntimeError("No HETATM records found in the predicted PDB file.")
            
        # b. Create a temporary molecule from ALL HETATMs and sanitize it.
        all_hets_mol = Chem.PathToSubmol(pred_complex, all_het_indices)
        Chem.SanitizeMol(all_hets_mol)

        # c. Perform a substructure search to find the TRUE ligand within the HETATMs.
        true_lig_match_in_hets = all_hets_mol.GetSubstructMatch(ref_lig)
        if not true_lig_match_in_hets:
            raise ValueError("Chemical Mismatch! The true ligand structure could not be found among the HETATM records of the predicted PDB.")
        
        # d. Map the indices from the temporary molecule back to the full complex.
        pred_lig_indices = [all_het_indices[i] for i in true_lig_match_in_hets]
        
        # 4. Now that we have the TRUE ligand atoms, align the proteins.
        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        ref_ca_map = {
            (a.GetPDBResidueInfo().GetChainId().strip(), a.GetPDBResidueInfo().GetResidueNumber()): a.GetIdx()
            for a in ref_prot.GetAtoms() if is_ca(a)
        }
        protein_atom_map = []
        for atom in pred_complex.GetAtoms():
            if atom.GetIdx() not in pred_lig_indices and is_ca(atom): # Must not be part of the identified ligand
                info = atom.GetPDBResidueInfo()
                key = (info.GetChainId().strip(), info.GetResidueNumber())
                if key in ref_ca_map:
                    protein_atom_map.append((atom.GetIdx(), ref_ca_map[key]))
        
        if len(protein_atom_map) < 10:
            raise ValueError(f"Found only {len(protein_atom_map)} common C-alpha atoms for alignment.")
        
        protein_rmsd = AllChem.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)
        print(f"  [DIAGNOSTIC] Protein C-alpha alignment RMSD: {protein_rmsd:.3f} Å")


        # 5. Get coordinates from the aligned prediction and the original reference ligand.
        pred_conf = pred_complex.GetConformer()
        ref_conf = ref_lig.GetConformer()
        
        # We know the atom order in `pred_lig_indices` now corresponds to the atom order in `ref_lig`
        # because GetSubstructMatch preserves it.
        pred_coords = np.array([pred_conf.GetAtomPosition(i) for i in pred_lig_indices])
        ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(ref_lig.GetNumAtoms())])
        
        # 6. Calculate the final manual RMSD.
        diff = ref_coords - pred_coords
        ligand_rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
        
        return {
            'ligand_pose_rmsd': ligand_rmsd,
            'protein_ca_rmsd': protein_rmsd,
            'n_atoms_matched': len(pred_lig_indices),
            'success': True
        }

    except Exception as e:
        return { 'success': False, 'error': str(e), 'traceback': traceback.format_exc() }
    
def compute_for_al():
    BASE_DIR = '/home/aoxu/projects/PoseBench/data/plinder_set/'
    data_dir = '/home/aoxu/projects/PoseBench/forks/boltz/inference/plinder_set_0/'
    rmsds = []
    for protein_name in os.listdir(data_dir):
        try:
            ref_prot_path = os.path.join(BASE_DIR, protein_name, f"{protein_name}_protein.pdb")
            pred_complex_path = os.path.join(data_dir, protein_name, f"predictions/{protein_name}/{protein_name}_model_0.pdb")
            ref_lig_path = os.path.join(BASE_DIR, protein_name, f"{protein_name}_ligand.sdf")

            result = calculate_ligand_rmsd_final_robust(ref_prot_path, pred_complex_path, ref_lig_path)
            if result['success']:
                print(f"{protein_name}: Ligand RMSD = {result['ligand_rmsd']:.3f} Å")
            else:
                print(f"{protein_name}: Error - {result['error']}")
            rmsds.append(result)
        except Exception as e:
            print(f"Error processing {protein_name}: {e}")

    print("\n--- Summary of RMSDs ---")
    for result in rmsds:
        if result.get('success') and result.get('ligand_rmsd') is not None and result['ligand_rmsd'] < 2:
            print(f"{result}")
    print(f"\nNumber of systems with ligand RMSD < 2 Å: {sum(1 for r in rmsds if r.get('success') and r.get('ligand_rmsd') is not None and r['ligand_rmsd'] < 2)}")

# --- Example Usage ---
# To test this, you would need to have actual files.
# Let's assume you have files named:
# 'reference_protein.pdb', 'predicted_complex.pdb', 'reference_ligand.sdf'

# Create dummy files for demonstration if they don't exist
# NOTE: This part is just to make the script runnable. Replace with your actual file paths.
try:
    # This is a placeholder. You should use your actual files.
    # We'll skip creating files and just show how to call the function.
    print("Please replace the placeholder file paths with your actual files to run a real test.")

    
    result = calculate_ligand_rmsd_production(
        ref_prot_path='/home/aoxu/projects/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_protein.pdb',
        pred_complex_path='/home/aoxu/projects/PoseBench/forks/boltz/inference//plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F/predictions/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_model_0.pdb',
        ref_lig_path='/home/aoxu/projects/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_ligand.sdf'
    )
    
    if result['success']:
         print(f"\n--- Final Result ---")
         print(f"Ligand RMSD: {result['ligand_rmsd']:.3f} Å")
    else:
         print(f"\n--- Calculation Failed ---")
         print(f"Error: {result['error']}")

except Exception as e:
    print(f"Could not run example: {e}")

# compute_for_al()
    
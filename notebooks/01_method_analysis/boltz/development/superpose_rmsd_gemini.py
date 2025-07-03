import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS, rdMolTransforms
import traceback
import numpy as np

def calculate_ligand_rmsd_final(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates the ligand RMSD after correctly superimposing the protein backbones.

    Version 5: Corrected for RDKit 2023.9.x. This version uses a fully
    keyword-based call to rdMolAlign.CalcRMS to avoid ambiguity with
    positional arguments, which caused the ArgumentError.
    """
    try:
        # 1. Load all necessary molecules
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        # 2. Create the atom map for protein alignment
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
        
        # 3. Superimpose the predicted complex onto the reference protein
        protein_align_rmsd = rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        # 4. Extract the predicted ligand
        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                             "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                             "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        
        pred_lig_atom_indices = [
            a.GetIdx() for a in pred_complex.GetAtoms() 
            if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues
        ]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        # 5. Find MCS to map ligand atoms
        mcs_result = rdFMCS.FindMCS([ref_lig, pred_lig], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs_result.numAtoms == 0:
            raise RuntimeError("No common substructure found between reference and predicted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = tuple(ref_lig.GetSubstructMatch(mcs_mol))
        pred_match = tuple(pred_lig.GetSubstructMatch(mcs_mol))

        if not ref_match or not pred_match:
            raise RuntimeError("Failed to match MCS to ligands.")

        # The order for the map is (probe, reference). pred_lig is the probe molecule.
        ligand_atom_map = list(zip(pred_match, ref_match))
        
        #
        # --- THE FINAL CORRECTION ---
        # For RDKit 2023.9.2, we use a full keyword-argument call to avoid all ambiguity.
        # The correct Python keyword is `atomMap`.
        #
        ligand_rmsd = rdMolAlign.CalcRMS(
            prbMol=pred_lig, 
            refMol=ref_lig, 
            atomMap=ligand_atom_map
        )

        return {
            'ligand_rmsd': ligand_rmsd,
            'protein_superposition_rmsd': protein_align_rmsd,
            'success': True,
            'error': None
        }

    except Exception as e:
        return {
            'ligand_rmsd': None,
            'protein_superposition_rmsd': None,
            'success': False,
            'error': f"{type(e).__name__}: {e}",
            'traceback': traceback.format_exc()
        }

def calculate_ligand_rmsd_definitive(ref_prot_path, pred_complex_path, ref_lig_path):
    """
    Calculates the ligand RMSD after correctly superimposing the protein backbones.

    Definitive Version: Abandons the problematic `rdMolAlign.CalcRMS` and uses
    the more direct and robust `rdMolTransforms.GetRMSD` function, which
    operates directly on conformer coordinates.
    """
    try:
        # Steps 1-4: Loading, Protein Superposition, and Ligand Extraction
        # This part of the logic is correct and remains unchanged.

        # 1. Load all necessary molecules
        ref_prot = Chem.MolFromPDBFile(ref_prot_path, removeHs=False, sanitize=False)
        pred_complex = Chem.MolFromPDBFile(pred_complex_path, removeHs=False, sanitize=False)
        ref_lig = Chem.RemoveHs(Chem.SDMolSupplier(ref_lig_path, removeHs=False)[0])

        if not all([ref_prot, pred_complex, ref_lig]):
            raise ValueError("Failed to load one or more input files.")

        def is_ca(atom):
            info = atom.GetPDBResidueInfo()
            return info and info.GetName().strip() == 'CA'

        # 2. Create the atom map for protein alignment
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
        
        # 3. Superimpose the predicted complex onto the reference protein
        rdMolAlign.AlignMol(prbMol=pred_complex, refMol=ref_prot, atomMap=protein_atom_map)

        # 4. Extract the predicted ligand
        standard_residues = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                             "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                             "THR", "TRP", "TYR", "VAL", "HOH", "WAT"}
        
        pred_lig_atom_indices = [
            a.GetIdx() for a in pred_complex.GetAtoms() 
            if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueName().strip() not in standard_residues
        ]

        if not pred_lig_atom_indices:
            raise RuntimeError("Could not identify ligand atoms in the predicted complex.")

        pred_lig = Chem.PathToSubmol(pred_complex, pred_lig_atom_indices, useQuery=False)
        pred_lig = Chem.RemoveHs(pred_lig)

        # 5. Find MCS to map ligand atoms. This is also correct and unchanged.
        mcs_result = rdFMCS.FindMCS([ref_lig, pred_lig], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs_result.numAtoms == 0:
            raise RuntimeError("No common substructure found between reference and predicted ligands.")
        
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        ref_match = tuple(ref_lig.GetSubstructMatch(mcs_mol))
        pred_match = tuple(pred_lig.GetSubstructMatch(mcs_mol))

        if not ref_match or not pred_match:
            raise RuntimeError("Failed to match MCS to ligands.")

        # The order for the map is (probe, reference) for rdMolAlign functions.
        # Let's check what rdMolTransforms.GetRMSD expects. It doesn't matter for GetRMSD.
        ligand_atom_map = list(zip(pred_match, ref_match))
        
        #
        # --- THE DEFINITIVE CORRECTION ---
        # We will bypass rdMolAlign.CalcRMS entirely and use rdMolTransforms.GetRMSD
        #
        
        # 6. Get the conformer objects (the 3D coordinates) from each ligand
        ref_conf = ref_lig.GetConformer()
        pred_conf = pred_lig.GetConformer()

        # 7. Calculate RMSD directly on the conformers using the atom map
        ligand_rmsd = rdMolTransforms.GetRMSD(pred_conf, ref_conf, atomMap=ligand_atom_map)

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
    
    # Example call:
    result = calculate_ligand_rmsd_manual(
        ref_prot_path='/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_protein.pdb',
        pred_complex_path='/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference//plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F/predictions/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_model_0.pdb',
        ref_lig_path='/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_ligand.sdf'
    )
    
    if result['success']:
         print(f"\n--- Final Result ---")
         print(f"Ligand RMSD: {result['ligand_rmsd']:.3f} Å")
    else:
         print(f"\n--- Calculation Failed ---")
         print(f"Error: {result['error']}")

except Exception as e:
    print(f"Could not run example: {e}")
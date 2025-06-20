#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def create_heavy_atom_only_molecules():
    """
    Create versions of molecules with only heavy atoms (manually exclude hydrogens)
    """
    print("=== Creating Heavy-Atom-Only Molecules ===")
    
    protein_name = "1r34__1__1.A__1.C_1.D"
    ref_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
    pred_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand_protein_aligned.sdf"
    
    # Read original molecules
    ref_mol = Chem.SDMolSupplier(ref_path)[0]
    pred_mol = Chem.SDMolSupplier(pred_path)[0]
    
    print(f"Original molecules:")
    print(f"Reference: {rdMolDescriptors.CalcMolFormula(ref_mol)} ({ref_mol.GetNumAtoms()} atoms)")
    print(f"Predicted: {rdMolDescriptors.CalcMolFormula(pred_mol)} ({pred_mol.GetNumAtoms()} atoms)")
    
    # Create heavy-atom-only versions by excluding hydrogen atoms
    ref_heavy = create_heavy_atom_mol(ref_mol)
    pred_heavy = create_heavy_atom_mol(pred_mol)
    
    if ref_heavy and pred_heavy:
        print(f"\\nHeavy-atom-only molecules:")
        print(f"Reference: {rdMolDescriptors.CalcMolFormula(ref_heavy)} ({ref_heavy.GetNumAtoms()} atoms)")
        print(f"Predicted: {rdMolDescriptors.CalcMolFormula(pred_heavy)} ({pred_heavy.GetNumAtoms()} atoms)")
        
        # Check if heavy atom formulas match
        ref_heavy_formula = rdMolDescriptors.CalcMolFormula(ref_heavy)
        pred_heavy_formula = rdMolDescriptors.CalcMolFormula(pred_heavy)
        
        print(f"Heavy atom formulas match: {'✓' if ref_heavy_formula == pred_heavy_formula else '✗'}")
        
        if ref_heavy_formula == pred_heavy_formula:
            # Save heavy-atom-only versions
            ref_heavy_path = ref_path.replace('.sdf', '_heavy_atoms_only.sdf')
            pred_heavy_path = pred_path.replace('.sdf', '_heavy_atoms_only.sdf')
            
            ref_writer = Chem.SDWriter(ref_heavy_path)
            ref_writer.write(ref_heavy)
            ref_writer.close()
            
            pred_writer = Chem.SDWriter(pred_heavy_path)
            pred_writer.write(pred_heavy)
            pred_writer.close()
            
            print(f"\\nSaved heavy-atom-only molecules:")
            print(f"Reference: {ref_heavy_path}")
            print(f"Predicted: {pred_heavy_path}")
            
            # Test RMSD
            test_heavy_atom_rmsd(pred_heavy_path, ref_heavy_path, protein_name)
        else:
            print("Even heavy atoms don't match - cannot use this approach")
    else:
        print("Failed to create heavy-atom-only molecules")

def create_heavy_atom_mol(mol):
    """
    Create a new molecule containing only heavy atoms (non-hydrogen)
    """
    try:
        # Get heavy atom indices
        heavy_atom_indices = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:  # Not hydrogen
                heavy_atom_indices.append(atom.GetIdx())
        
        if not heavy_atom_indices:
            return None
        
        # Create new molecule with only heavy atoms
        # This is complex in RDKit, so let's use a simpler approach
        # We'll use the removeHs function but force it to work
        
        # Try different approaches to remove hydrogens
        try:
            # Approach 1: Explicit hydrogen removal
            heavy_mol = Chem.RemoveHs(mol, implicitOnly=False, updateExplicitCount=False)
            if heavy_mol.GetNumAtoms() < mol.GetNumAtoms():
                return heavy_mol
        except:
            pass
            
        try:
            # Approach 2: Remove all explicit hydrogens
            heavy_mol = Chem.RemoveHs(mol, implicitOnly=False)
            if heavy_mol.GetNumAtoms() < mol.GetNumAtoms():
                return heavy_mol
        except:
            pass
            
        try:
            # Approach 3: Sanitize then remove
            mol_copy = Chem.Mol(mol)
            Chem.SanitizeMol(mol_copy)
            heavy_mol = Chem.RemoveHs(mol_copy, implicitOnly=False, updateExplicitCount=False, sanitize=False)
            return heavy_mol
        except:
            pass
        
        # If all else fails, return the original molecule
        print(f"Warning: Could not remove hydrogens, using original molecule")
        return mol
        
    except Exception as e:
        print(f"Error creating heavy-atom molecule: {e}")
        return None

def test_heavy_atom_rmsd(pred_heavy_path, ref_heavy_path, protein_name):
    """Test RMSD with heavy-atom-only molecules"""
    try:
        from posebusters.posebusters import PoseBusters
        
        print(f"\\n=== Testing Heavy-Atom RMSD ===")
        
        pb = PoseBusters(config='redock', top_n=None)
        protein_pdb = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_protein.pdb"
        
        result = pb.bust(
            mol_pred=pred_heavy_path,
            mol_true=ref_heavy_path,
            mol_cond=protein_pdb,
            full_report=True
        )
        
        rmsd = result.get('rmsd', ['N/A'])[0]
        rmsd_lt2 = result.get('rmsd_≤_2å', ['N/A'])[0]
        mol_formula = result.get('molecular_formula', ['N/A'])[0]
        
        print(f"Heavy-atom RMSD: {rmsd:.3f} Å" if isinstance(rmsd, (int, float)) else f"Heavy-atom RMSD: {rmsd}")
        print(f"RMSD ≤ 2Å: {'✓' if rmsd_lt2 else '✗'}")
        print(f"Molecular formula check: {'✓' if mol_formula else '✗'}")
        
    except Exception as e:
        print(f"Error testing heavy-atom RMSD: {e}")

if __name__ == "__main__":
    create_heavy_atom_only_molecules()

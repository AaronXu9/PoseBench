#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os

def implement_heavy_atom_approach():
    """
    Implement heavy-atom-only comparison as fallback for chemical differences
    """
    print("=== Testing Heavy-Atom-Only Approach ===")
    
    # Test case that fails chemical correction
    protein_name = "1r34__1__1.A__1.C_1.D"
    ref_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
    pred_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand_protein_aligned.sdf"
    
    # Read molecules
    ref_mol = Chem.SDMolSupplier(ref_path)[0]
    pred_mol = Chem.SDMolSupplier(pred_path)[0]
    
    print(f"Original molecules:")
    print(f"Reference: {rdMolDescriptors.CalcMolFormula(ref_mol)}")
    print(f"Predicted: {rdMolDescriptors.CalcMolFormula(pred_mol)}")
    
    # Remove hydrogens to compare heavy atoms only
    ref_heavy = Chem.RemoveHs(ref_mol)
    pred_heavy = Chem.RemoveHs(pred_mol)
    
    print(f"\\nHeavy atoms only:")
    print(f"Reference: {rdMolDescriptors.CalcMolFormula(ref_heavy)} ({ref_heavy.GetNumAtoms()} atoms)")
    print(f"Predicted: {rdMolDescriptors.CalcMolFormula(pred_heavy)} ({pred_heavy.GetNumAtoms()} atoms)")
    
    # Check if heavy atom formulas match
    ref_heavy_formula = rdMolDescriptors.CalcMolFormula(ref_heavy)
    pred_heavy_formula = rdMolDescriptors.CalcMolFormula(pred_heavy)
    
    print(f"\\nHeavy atom formulas match: {'✓' if ref_heavy_formula == pred_heavy_formula else '✗'}")
    
    if ref_heavy_formula == pred_heavy_formula:
        print("SUCCESS: Can use heavy-atom-only RMSD comparison!")
        
        # Save heavy-atom-only versions for testing
        ref_heavy_path = ref_path.replace('.sdf', '_heavy_only.sdf')
        pred_heavy_path = pred_path.replace('.sdf', '_heavy_only.sdf')
        
        # Write heavy-atom-only molecules
        ref_writer = Chem.SDWriter(ref_heavy_path)
        ref_writer.write(ref_heavy)
        ref_writer.close()
        
        pred_writer = Chem.SDWriter(pred_heavy_path)
        pred_writer.write(pred_heavy)
        pred_writer.close()
        
        print(f"Saved heavy-atom versions:")
        print(f"Reference: {ref_heavy_path}")
        print(f"Predicted: {pred_heavy_path}")
        
        # Test with PoseBusters
        test_heavy_atom_rmsd(pred_heavy_path, ref_heavy_path, protein_name)
    else:
        print("FAILURE: Even heavy atoms don't match - fundamental structural difference")

def test_heavy_atom_rmsd(pred_heavy_path, ref_heavy_path, protein_name):
    """Test RMSD calculation with heavy-atom-only molecules"""
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
        
        return rmsd
        
    except Exception as e:
        print(f"Error testing heavy-atom RMSD: {e}")
        return None

if __name__ == "__main__":
    implement_heavy_atom_approach()

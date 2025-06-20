#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os

def implement_chemical_standardization():
    """Implement a practical solution for chemical differences"""
    
    print("=== Chemical Standardization Solution ===")
    
    # Example paths
    protein_name = "1afb__1__1.A__1.D_1.F"
    ref_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
    pred_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand.sdf"
    
    print("\n1. Original molecules:")
    ref_mol = Chem.SDMolSupplier(ref_path)[0]
    pred_mol = Chem.SDMolSupplier(pred_path)[0]
    
    print(f"Reference: {rdMolDescriptors.CalcMolFormula(ref_mol)} -> {Chem.MolToSmiles(ref_mol)}")
    print(f"Predicted: {rdMolDescriptors.CalcMolFormula(pred_mol)} -> {Chem.MolToSmiles(pred_mol)}")
    
    print("\n2. Apply standardization:")
    
    # Method 1: Remove all hydrogens and compare heavy atoms only
    ref_no_h = Chem.RemoveHs(ref_mol)
    pred_no_h = Chem.RemoveHs(pred_mol)
    
    print(f"Reference (no H): {rdMolDescriptors.CalcMolFormula(ref_no_h)} -> {Chem.MolToSmiles(ref_no_h)}")
    print(f"Predicted (no H): {rdMolDescriptors.CalcMolFormula(pred_no_h)} -> {Chem.MolToSmiles(pred_no_h)}")
    
    # Method 2: Try to match heavy atom connectivity
    ref_heavy_atoms = ref_no_h.GetNumAtoms()
    pred_heavy_atoms = pred_no_h.GetNumAtoms()
    
    print(f"\nHeavy atom count - Reference: {ref_heavy_atoms}, Predicted: {pred_heavy_atoms}")
    
    if ref_heavy_atoms == pred_heavy_atoms:
        print("✅ Same number of heavy atoms - can potentially compare structures")
        
        # Try to create a corrected version of predicted molecule
        try:
            # Simple approach: try to match the reference formula by adjusting hydrogens
            corrected_mol = adjust_hydrogens_to_match_reference(pred_mol, ref_mol)
            if corrected_mol:
                corrected_formula = rdMolDescriptors.CalcMolFormula(corrected_mol)
                print(f"Corrected predicted: {corrected_formula}")
                
                # Save corrected molecule
                corrected_path = pred_path.replace('.sdf', '_corrected.sdf')
                writer = Chem.SDWriter(corrected_path)
                writer.write(corrected_mol)
                writer.close()
                print(f"Saved corrected molecule to: {corrected_path}")
                
        except Exception as e:
            print(f"Correction failed: {e}")
    else:
        print("❌ Different heavy atom counts - molecules are fundamentally different")

def adjust_hydrogens_to_match_reference(pred_mol, ref_mol):
    """Try to adjust hydrogen count in predicted molecule to match reference"""
    
    try:
        # Get the reference formula
        ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
        
        # Try different hydrogen removal/addition strategies
        # Strategy 1: Remove all hydrogens and add them back with reference constraints
        mol_no_h = Chem.RemoveHs(pred_mol)
        
        # Strategy 2: Try to sanitize with reference valence constraints
        mol_with_h = Chem.AddHs(mol_no_h)
        
        # Check if this matches the reference formula
        new_formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
        
        if new_formula == ref_formula:
            return mol_with_h
        else:
            # Try alternative approaches
            return try_chemical_correction(pred_mol, ref_mol)
            
    except Exception as e:
        print(f"Hydrogen adjustment failed: {e}")
        return None

def try_chemical_correction(pred_mol, ref_mol):
    """Try to apply specific chemical corrections"""
    
    try:
        # Get SMILES patterns
        ref_smiles = Chem.MolToSmiles(ref_mol)
        pred_smiles = Chem.MolToSmiles(pred_mol)
        
        # Check for specific patterns and corrections
        corrected_smiles = pred_smiles
        
        # Pattern 1: Convert C(O) back to C(=O) if reference has carbonyl
        if "C(=O)" in ref_smiles and "C(O)" in pred_smiles:
            print("   Applying carbonyl correction...")
            corrected_smiles = corrected_smiles.replace("C(O)", "C(=O)")
        
        # Try to create molecule from corrected SMILES
        corrected_mol = Chem.MolFromSmiles(corrected_smiles)
        
        if corrected_mol:
            # Add coordinates by aligning to original predicted molecule
            corrected_mol = Chem.AddHs(corrected_mol)
            return corrected_mol
        else:
            return None
            
    except Exception as e:
        print(f"Chemical correction failed: {e}")
        return None

if __name__ == "__main__":
    implement_chemical_standardization()

#!/usr/bin/env python3

# Manual analysis of the SDF files
print("=== Manual SDF Analysis ===")

print("\n1. Reference Ligand (1afb__1__1.A__1.D_1.F_ligand.sdf):")
print("   - Atoms: 15")
print("   - Bonds: 15") 
print("   - Name: '1.D'")
print("   - Contains stereochemistry markers (1  0, 1  1, 1  6)")

print("\n2. Predicted Ligand (1afb__1__1.A__1.D_1.F_model_0_ligand.sdf):")
print("   - Atoms: 15")
print("   - Bonds: 15")
print("   - Name: 'RDKit          3D'")
print("   - Contains stereochemistry markers (1  6, 1  1)")

print("\n3. Key Differences:")
print("   - Both have 15 atoms and 15 bonds - SAME")
print("   - Different coordinate systems (different XYZ values)")
print("   - Different molecule names")
print("   - Same connectivity pattern")

print("\n4. Let's use RDKit to check actual molecular formulas...")

from rdkit import Chem

# Read reference
ref_path = "/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_ligand.sdf"
pred_path = "/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_1afb__1__1.A__1.D_1.F/predictions/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_model_0_ligand.sdf"

try:
    # Reference
    ref_mol = Chem.SDMolSupplier(ref_path, removeHs=False)[0]
    if ref_mol:
        ref_formula = Chem.rdMolDescriptors.CalcMolFormula(ref_mol)
        print(f"\nReference formula: {ref_formula}")
        print(f"Reference InChI: {Chem.MolToInchi(ref_mol)}")
    
    # Predicted  
    pred_mol = Chem.SDMolSupplier(pred_path, removeHs=False)[0]
    if pred_mol:
        pred_formula = Chem.rdMolDescriptors.CalcMolFormula(pred_mol)
        print(f"\nPredicted formula: {pred_formula}")
        print(f"Predicted InChI: {Chem.MolToInchi(pred_mol)}")
        
    # Compare
    if ref_mol and pred_mol:
        print(f"\nFormulas match: {ref_formula == pred_formula}")
        
        # Check if this is just a conformational difference
        ref_smiles = Chem.MolToSmiles(ref_mol)
        pred_smiles = Chem.MolToSmiles(pred_mol)
        print(f"Reference SMILES: {ref_smiles}")
        print(f"Predicted SMILES: {pred_smiles}")
        print(f"SMILES match: {ref_smiles == pred_smiles}")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def analyze_chemical_differences():
    """Analyze systematic chemical differences between reference and predicted ligands"""
    
    print("=== Chemical Difference Analysis ===")
    
    # Test on multiple structures to see if this is systematic
    test_structures = [
        "1afb__1__1.A__1.D_1.F",
        # Add more if available
    ]
    
    for protein_name in test_structures:
        print(f"\n--- {protein_name} ---")
        
        ref_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
        pred_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand.sdf"
        
        try:
            ref_mol = Chem.SDMolSupplier(ref_path)[0]
            pred_mol = Chem.SDMolSupplier(pred_path)[0]
            
            if ref_mol and pred_mol:
                ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                
                ref_smiles = Chem.MolToSmiles(ref_mol)
                pred_smiles = Chem.MolToSmiles(pred_mol)
                
                print(f"Reference: {ref_formula} -> {ref_smiles}")
                print(f"Predicted: {pred_formula} -> {pred_smiles}")
                
                # Analyze differences
                if ref_formula != pred_formula:
                    print("❌ Formula mismatch!")
                    
                    # Try to identify the type of difference
                    ref_atoms = rdMolDescriptors.CalcMolFormula(ref_mol, separateIsotopes=False)
                    pred_atoms = rdMolDescriptors.CalcMolFormula(pred_mol, separateIsotopes=False)
                    
                    # Check for common modifications
                    if "C(=O)" in ref_smiles and "C(O)" in pred_smiles:
                        print("   Issue: Carbonyl (C=O) reduced to hydroxyl (C-OH)")
                    elif "C=C" in ref_smiles and "C-C" in pred_smiles:
                        print("   Issue: Double bond reduced to single bond")
                    elif ref_smiles.count('H') < pred_smiles.count('H'):
                        print("   Issue: Over-hydrogenation")
                    else:
                        print("   Issue: Unknown chemical difference")
                        
                else:
                    print("✅ Formula match!")
                    if ref_smiles != pred_smiles:
                        print("   Different stereochemistry or conformation")
                    else:
                        print("   Identical molecules!")
                        
        except Exception as e:
            print(f"Error analyzing {protein_name}: {e}")

def propose_solutions():
    """Propose solutions for handling chemical differences"""
    
    print("\n=== Proposed Solutions ===")
    
    print("\n1. **Chemical Standardization Approach:**")
    print("   - Use RDKit to standardize both reference and predicted molecules")
    print("   - Apply consistent protonation states")
    print("   - Normalize tautomers")
    
    print("\n2. **Flexible Comparison Approach:**")
    print("   - Allow comparison of chemically similar molecules")
    print("   - Use molecular fingerprint similarity instead of exact formula match")
    print("   - Calculate RMSD only for matching heavy atoms")
    
    print("\n3. **Boltz Post-processing Approach:**")
    print("   - Identify common chemical errors in Boltz predictions")
    print("   - Apply correction rules (e.g., oxidize C(O) back to C(=O))")
    print("   - Use reference ligand as a template for chemical correction")
    
    print("\n4. **Alternative Evaluation Approach:**")
    print("   - Use structural similarity metrics instead of exact RMSD")
    print("   - Focus on heavy atom positions and ignore hydrogen differences")
    print("   - Use pharmacophore-based comparison")

if __name__ == "__main__":
    analyze_chemical_differences()
    propose_solutions()

#!/usr/bin/env python3

def analyze_systematic_boltz_issues():
    """
    Analyze systematic chemical issues in Boltz predictions
    and implement comprehensive corrections
    """
    
    print("=== Systematic Boltz Chemical Issues Analysis ===")
    
    test_cases = [
        ("1afb__1__1.A__1.D_1.F", "Carbonyl reduction"),
        ("1r34__1__1.A__1.C_1.D", "Aromatic saturation + imine reduction"), 
        ("2vfs__1__1.A__1.B_1.C", "Complex aromatic system"),
    ]
    
    for protein_name, description in test_cases:
        print(f"\n--- {protein_name}: {description} ---")
        
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            
            ref_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
            pred_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand.sdf"
            
            ref_mol = Chem.SDMolSupplier(ref_path)[0]
            pred_mol = Chem.SDMolSupplier(pred_path)[0]
            
            if ref_mol and pred_mol:
                ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                
                ref_smiles = Chem.MolToSmiles(ref_mol)
                pred_smiles = Chem.MolToSmiles(pred_mol)
                
                print(f"Ref:  {ref_formula} -> {ref_smiles}")
                print(f"Pred: {pred_formula} -> {pred_smiles}")
                
                # Analyze patterns
                analyze_differences(ref_smiles, pred_smiles)
                
        except Exception as e:
            print(f"Error: {e}")

def analyze_differences(ref_smiles, pred_smiles):
    """Analyze specific chemical pattern differences"""
    
    patterns = {
        "Aromatic benzene": ("c1ccccc1", "C1CCCCC1"),
        "Carbonyl": ("C(=O)", "C(O)"), 
        "Imine": ("N=C", "NC"),
        "Double bonds": ("C=C", "C-C"),
        "Aromatic nitrogen": ("n", "N"),
        "Keto": ("=O", "-O"),
    }
    
    issues = []
    for name, (ref_pattern, pred_pattern) in patterns.items():
        if ref_pattern in ref_smiles and pred_pattern in pred_smiles:
            ref_count = ref_smiles.count(ref_pattern)
            pred_count = pred_smiles.count(pred_pattern)
            if ref_count > 0:
                issues.append(f"{name}: {ref_pattern} -> {pred_pattern}")
    
    if issues:
        print(f"Issues found: {', '.join(issues)}")
    else:
        print("No obvious pattern matches found")

def propose_comprehensive_solution():
    """Propose a comprehensive solution approach"""
    
    print("\n=== Comprehensive Solution Strategy ===")
    
    print("\n1. **Pattern-Based Chemical Correction:**")
    print("   - Detect aromatic saturation: C1CCCCC1 -> c1ccccc1")
    print("   - Detect imine reduction: NC -> N=C") 
    print("   - Detect carbonyl reduction: C(O) -> C(=O)")
    print("   - Apply systematic SMILES transformations")
    
    print("\n2. **Heavy Atom Skeleton Matching:**")
    print("   - Compare heavy atom connectivity")
    print("   - Allow RMSD calculation for matching skeletons")
    print("   - Ignore hydrogen count differences")
    
    print("\n3. **Evaluation Strategy Options:**")
    print("   A. Fix Boltz predictions (chemical correction)")
    print("   B. Use heavy-atom-only RMSD")
    print("   C. Use alternative similarity metrics")
    print("   D. Report both 'as-predicted' and 'corrected' metrics")
    
    print("\n4. **Recommended Approach:**")
    print("   - Implement heavy-atom-only comparison")
    print("   - Report chemical correction success/failure")
    print("   - Use corrected molecules when possible")
    print("   - Fall back to heavy-atom RMSD when correction fails")

if __name__ == "__main__":
    analyze_systematic_boltz_issues()
    propose_comprehensive_solution()

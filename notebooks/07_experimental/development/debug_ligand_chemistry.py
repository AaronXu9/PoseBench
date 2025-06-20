#!/usr/bin/env python3

import os
from rdkit import Chem

def analyze_ligand_chemistry():
    """Analyze the chemical differences between reference and predicted ligands"""
    
    # Paths
    protein_name = "1afb__1__1.A__1.D_1.F"
    ref_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/{protein_name}/{protein_name}_ligand.sdf"
    pred_ligand_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0_ligand.sdf"
    pred_pdb_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{protein_name}/predictions/{protein_name}/{protein_name}_model_0.pdb"
    
    print("=== Analyzing Chemical Differences ===")
    
    # 1. Read reference ligand
    print("\n1. Reference Ligand Analysis:")
    if os.path.exists(ref_ligand_path):
        ref_suppl = Chem.SDMolSupplier(ref_ligand_path, removeHs=False)
        ref_mol = None
        for mol in ref_suppl:
            if mol is not None:
                ref_mol = mol
                break
        
        if ref_mol:
            print(f"   File: {ref_ligand_path}")
            print(f"   Formula: {Chem.rdMolDescriptors.CalcMolFormula(ref_mol)}")
            print(f"   Num atoms: {ref_mol.GetNumAtoms()}")
            print(f"   Num bonds: {ref_mol.GetNumBonds()}")
            print(f"   Num heavy atoms: {ref_mol.GetNumHeavyAtoms()}")
            print(f"   InChI Key: {Chem.inchi.MolToInchiKey(ref_mol)}")
            
            # Check aromatic bonds
            aromatic_bonds = sum(1 for bond in ref_mol.GetBonds() if bond.GetIsAromatic())
            print(f"   Aromatic bonds: {aromatic_bonds}")
        else:
            print("   ERROR: Could not read reference molecule")
    else:
        print(f"   ERROR: Reference file not found: {ref_ligand_path}")
    
    # 2. Read predicted ligand (if exists)
    print("\n2. Predicted Ligand Analysis (SDF):")
    if os.path.exists(pred_ligand_path):
        pred_suppl = Chem.SDMolSupplier(pred_ligand_path, removeHs=False)
        pred_mol = None
        for mol in pred_suppl:
            if mol is not None:
                pred_mol = mol
                break
        
        if pred_mol:
            print(f"   File: {pred_ligand_path}")
            print(f"   Formula: {Chem.rdMolDescriptors.CalcMolFormula(pred_mol)}")
            print(f"   Num atoms: {pred_mol.GetNumAtoms()}")
            print(f"   Num bonds: {pred_mol.GetNumBonds()}")
            print(f"   Num heavy atoms: {pred_mol.GetNumHeavyAtoms()}")
            print(f"   InChI Key: {Chem.inchi.MolToInchiKey(pred_mol)}")
            
            # Check aromatic bonds
            aromatic_bonds = sum(1 for bond in pred_mol.GetBonds() if bond.GetIsAromatic())
            print(f"   Aromatic bonds: {aromatic_bonds}")
        else:
            print("   ERROR: Could not read predicted molecule")
    else:
        print(f"   INFO: Predicted SDF not found (will be created)")
    
    # 3. Analyze the raw PDB ligand data
    print("\n3. Raw PDB Ligand Analysis:")
    if os.path.exists(pred_pdb_path):
        print(f"   File: {pred_pdb_path}")
        
        # Extract ligand lines
        ligand_lines = []
        atom_count = 0
        elements = {}
        
        with open(pred_pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') and 'LIG' in line:
                    ligand_lines.append(line)
                    atom_count += 1
                    
                    # Parse element from PDB line
                    element = line[76:78].strip() or line[12:14].strip()
                    element = element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()
                    elements[element] = elements.get(element, 0) + 1
        
        print(f"   Total ligand atoms in PDB: {atom_count}")
        print(f"   Element composition: {elements}")
        
        # Calculate expected formula from PDB
        if elements:
            formula_parts = []
            for elem in sorted(elements.keys()):
                if elements[elem] > 1:
                    formula_parts.append(f"{elem}{elements[elem]}")
                else:
                    formula_parts.append(elem)
            pdb_formula = "".join(formula_parts)
            print(f"   Expected formula from PDB: {pdb_formula}")
        
        # Show first few ligand lines
        print(f"   First few ligand lines:")
        for i, line in enumerate(ligand_lines[:5]):
            print(f"     {line.strip()}")
        
    else:
        print(f"   ERROR: PDB file not found: {pred_pdb_path}")
    
    # 4. Test different conversion approaches
    print("\n4. Testing Conversion Approaches:")
    if os.path.exists(pred_pdb_path):
        
        # Method 1: Current approach (removeHs=False)
        print("\n   Method 1: Current approach (removeHs=False)")
        test_conversion_method(pred_pdb_path, "method1", removeHs=False)
        
        # Method 2: Remove hydrogens
        print("\n   Method 2: Remove hydrogens (removeHs=True)")
        test_conversion_method(pred_pdb_path, "method2", removeHs=True)
        
        # Method 3: Sanitize molecule
        print("\n   Method 3: With sanitization")
        test_conversion_method(pred_pdb_path, "method3", removeHs=False, sanitize=True)

def test_conversion_method(pdb_path, method_name, removeHs=False, sanitize=False):
    """Test different approaches to convert PDB ligand to SDF"""
    try:
        # Extract ligand lines
        ligand_lines = []
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') and 'LIG' in line:
                    ligand_lines.append(line)
        
        # Create temporary PDB
        temp_pdb = f"/tmp/test_ligand_{method_name}.pdb"
        with open(temp_pdb, 'w') as f:
            f.writelines(ligand_lines)
            f.write('END\n')
        
        # Convert using RDKit
        mol = Chem.MolFromPDBFile(temp_pdb, removeHs=removeHs)
        if mol is not None:
            if sanitize:
                try:
                    Chem.SanitizeMol(mol)
                except:
                    print(f"     Sanitization failed")
                    
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            print(f"     Formula: {formula}")
            print(f"     Num atoms: {mol.GetNumAtoms()}")
            print(f"     Num heavy atoms: {mol.GetNumHeavyAtoms()}")
            
            # Check if we can write to SDF
            test_sdf = f"/tmp/test_ligand_{method_name}.sdf"
            writer = Chem.SDWriter(test_sdf)
            writer.write(mol)
            writer.close()
            print(f"     SDF written successfully")
            
            # Clean up
            os.remove(temp_pdb)
            os.remove(test_sdf)
        else:
            print(f"     ERROR: Could not parse molecule")
            os.remove(temp_pdb)
            
    except Exception as e:
        print(f"     ERROR: {str(e)}")

if __name__ == "__main__":
    analyze_ligand_chemistry()

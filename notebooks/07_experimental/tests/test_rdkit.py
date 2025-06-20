#!/usr/bin/env python3

print("Testing basic functionality...")

try:
    from rdkit import Chem
    print("RDKit imported")
    
    ref_path = "/Users/aoxu/projects/DrugDiscovery/PoseBench/data/plinder_set/1afb__1__1.A__1.D_1.F/1afb__1__1.A__1.D_1.F_ligand.sdf"
    print(f"Trying to read: {ref_path}")
    
    suppl = Chem.SDMolSupplier(ref_path)
    print("SDMolSupplier created")
    
    mol = next(suppl)
    print("Molecule extracted")
    
    if mol:
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        print(f"Formula: {formula}")
    else:
        print("Molecule is None")
        
except Exception as e:
    print(f"Exception: {e}")
    import traceback
    traceback.print_exc()

print("Script finished")

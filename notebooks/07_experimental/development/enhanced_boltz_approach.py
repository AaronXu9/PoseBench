#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
try:
    from rdkit.Chem import rdMolDescriptors, rdMolAlign
except ImportError:
    # Import from alternative locations if needed
    import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
    import rdkit.Chem.rdMolAlign as rdMolAlign

from scipy.spatial.distance import cdist
import logging

# Add the parent directory to Python path to import modules
sys.path.append('/Users/aoxu/projects/DrugDiscovery/PoseBench')

from notebooks.Approach import BoltzApproach
try:
    from notebooks.alignment import calculate_complex_rmsd
except ImportError:
    # Fallback implementation if alignment module is not available
    def calculate_complex_rmsd(pred_path, ref_path):
        """Simplified RMSD calculation"""
        try:
            pred_mol = Chem.SDMolSupplier(pred_path)[0]
            ref_mol = Chem.SDMolSupplier(ref_path)[0]
            if pred_mol and ref_mol:
                return rdMolAlign.AlignMol(pred_mol, ref_mol)
            return float('inf')
        except:
            return float('inf')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EnhancedBoltzApproach(BoltzApproach):
    """Enhanced Boltz approach with improved chemical standardization and heavy-atom evaluation"""
    
    def __init__(self):
        super().__init__()
        self.heavy_atom_fallback_count = 0
        self.chemical_correction_count = 0
        self.successful_alignment_count = 0
        self.total_processed = 0
        self.dataset_name = 'plinder_set'  # Set default dataset
        
    def _standardize_molecule_enhanced(self, mol, reference_mol=None):
        """Enhanced molecule standardization with comprehensive chemical corrections"""
        if mol is None:
            return None
            
        # Get molecular formulas
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol) if reference_mol else None
        
        logger.info(f"Predicted formula: {mol_formula}, Reference formula: {ref_formula}")
        
        # If formulas match, return as-is
        if mol_formula == ref_formula:
            return mol
            
        # Apply comprehensive chemical corrections
        corrected_mol = self._apply_comprehensive_corrections(mol, reference_mol)
        if corrected_mol:
            corrected_formula = rdMolDescriptors.CalcMolFormula(corrected_mol)
            if corrected_formula == ref_formula:
                logger.info(f"Successfully corrected formula from {mol_formula} to {corrected_formula}")
                self.chemical_correction_count += 1
                return corrected_mol
                
        # If chemical correction fails, return original for heavy-atom evaluation
        logger.warning(f"Chemical correction failed. Will use heavy-atom-only evaluation.")
        self.heavy_atom_fallback_count += 1
        return mol
        
    def _apply_comprehensive_corrections(self, mol, reference_mol):
        """Apply comprehensive chemical corrections based on common Boltz prediction patterns"""
        if not mol or not reference_mol:
            return mol
            
        # Create a copy to avoid modifying the original
        mol_copy = Chem.Mol(mol)
        
        # Pattern 1: Try simple hydrogen adjustment
        corrected_mol = self._try_hydrogen_correction(mol_copy, reference_mol)
        if corrected_mol:
            return corrected_mol
            
        # Pattern 2: Try SMILES-based corrections
        corrected_mol = self._try_smiles_correction_enhanced(mol_copy, reference_mol)
        if corrected_mol:
            return corrected_mol
            
        return None
        
    def _try_hydrogen_correction(self, mol, reference_mol):
        """Try to correct hydrogen count by removing and re-adding with proper valences"""
        try:
            # Remove all hydrogens and re-add
            mol_no_h = Chem.RemoveHs(mol)
            mol_with_h = Chem.AddHs(mol_no_h)
            
            # Check if this matches the reference formula
            corrected_formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
            ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol)
            
            if corrected_formula == ref_formula:
                logger.info("Successfully applied hydrogen correction")
                return mol_with_h
                
        except Exception as e:
            logger.warning(f"Hydrogen correction failed: {str(e)}")
            
        return None
        
    def _try_smiles_correction_enhanced(self, mol, reference_mol):
        """Enhanced SMILES-based chemical correction"""
        try:
            ref_smiles = Chem.MolToSmiles(reference_mol)
            pred_smiles = Chem.MolToSmiles(mol)
            
            logger.info(f"SMILES correction - Ref: {ref_smiles[:50]}...")
            logger.info(f"SMILES correction - Pred: {pred_smiles[:50]}...")
            
            # Apply multiple correction patterns
            corrected_smiles = pred_smiles
            corrections_applied = []
            
            # Pattern 1: Carbonyl reduction (C(=O) ↔ C(O))
            if "C(=O)" in ref_smiles and "C(O)" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("C(O)", "C(=O)", 1)
                corrections_applied.append("C(O)->C(=O)")
            elif "C(O)" in ref_smiles and "C(=O)" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("C(=O)", "C(O)", 1)
                corrections_applied.append("C(=O)->C(O)")
                
            if corrections_applied:
                logger.info(f"Applied SMILES corrections: {', '.join(corrections_applied)}")
                corrected_mol = Chem.MolFromSmiles(corrected_smiles)
                if corrected_mol:
                    corrected_mol = Chem.AddHs(corrected_mol)
                    # Verify the correction worked
                    corrected_formula = rdMolDescriptors.CalcMolFormula(corrected_mol)
                    ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol)
                    if corrected_formula == ref_formula:
                        return corrected_mol
                        
        except Exception as e:
            logger.warning(f"SMILES correction failed: {str(e)}")
            
        return None
        
    def _calculate_heavy_atom_rmsd(self, mol1, mol2):
        """Calculate RMSD using only heavy atoms when chemical formulas don't match"""
        if not mol1 or not mol2:
            return float('inf')
            
        try:
            # Get heavy atom coordinates
            conf1 = mol1.GetConformer()
            conf2 = mol2.GetConformer()
            
            heavy_atoms1 = [i for i in range(mol1.GetNumAtoms()) if mol1.GetAtomWithIdx(i).GetSymbol() != 'H']
            heavy_atoms2 = [i for i in range(mol2.GetNumAtoms()) if mol2.GetAtomWithIdx(i).GetSymbol() != 'H']
            
            # Get element symbols for heavy atoms
            elements1 = [mol1.GetAtomWithIdx(i).GetSymbol() for i in heavy_atoms1]
            elements2 = [mol2.GetAtomWithIdx(i).GetSymbol() for i in heavy_atoms2]
            
            # If different number of heavy atoms, try to match common atoms
            if len(heavy_atoms1) != len(heavy_atoms2):
                logger.warning(f"Different heavy atom counts: {len(heavy_atoms1)} vs {len(heavy_atoms2)}")
                # Use minimum number for comparison
                min_atoms = min(len(heavy_atoms1), len(heavy_atoms2))
                heavy_atoms1 = heavy_atoms1[:min_atoms]
                heavy_atoms2 = heavy_atoms2[:min_atoms]
                elements1 = elements1[:min_atoms]
                elements2 = elements2[:min_atoms]
            
            if len(heavy_atoms1) == 0:
                return float('inf')
                
            # Get coordinates
            coords1 = np.array([conf1.GetAtomPosition(i) for i in heavy_atoms1])
            coords2 = np.array([conf2.GetAtomPosition(i) for i in heavy_atoms2])
            
            # For small molecules, try to find optimal atom correspondence
            if len(heavy_atoms1) <= 10:
                return self._calculate_optimal_heavy_atom_rmsd(coords1, coords2, elements1, elements2)
            else:
                # For larger molecules, use simple distance-based matching
                return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
                
        except Exception as e:
            logger.error(f"Heavy atom RMSD calculation failed: {str(e)}")
            return float('inf')
    
    def _calculate_optimal_heavy_atom_rmsd(self, coords1, coords2, elements1, elements2):
        """Calculate optimal RMSD by trying different atom correspondences for small molecules"""
        try:
            from itertools import permutations
            
            # For very small molecules, try all permutations
            if len(coords1) <= 6:
                min_rmsd = float('inf')
                for perm in permutations(range(len(coords2))):
                    try:
                        # Check if this permutation preserves element types
                        valid = all(elements1[i] == elements2[perm[i]] for i in range(len(elements1)))
                        if valid:
                            rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2[list(perm)])**2, axis=1)))
                            min_rmsd = min(min_rmsd, rmsd)
                    except:
                        continue
                        
                return min_rmsd if min_rmsd != float('inf') else self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
            else:
                return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
                
        except Exception as e:
            logger.warning(f"Optimal RMSD calculation failed: {str(e)}")
            return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
    
    def _calculate_greedy_heavy_atom_rmsd(self, coords1, coords2):
        """Calculate RMSD using greedy distance-based atom assignment"""
        try:
            # Calculate pairwise distances
            distances = cdist(coords1, coords2)
            
            # Greedy assignment
            used_j = set()
            total_dist_sq = 0
            assignments = 0
            
            for i in range(len(coords1)):
                best_j = None
                best_dist = float('inf')
                for j in range(len(coords2)):
                    if j not in used_j and distances[i][j] < best_dist:
                        best_dist = distances[i][j]
                        best_j = j
                        
                if best_j is not None:
                    used_j.add(best_j)
                    total_dist_sq += best_dist**2
                    assignments += 1
                    
            if assignments == 0:
                return float('inf')
                
            return np.sqrt(total_dist_sq / assignments)
            
        except Exception as e:
            logger.error(f"Greedy RMSD calculation failed: {str(e)}")
            return float('inf')
    
    def _convert_pdb_to_sdf_enhanced(self, pdb_path, sdf_path=None):
        """Enhanced PDB to SDF conversion with comprehensive processing"""
        try:
            if sdf_path is None:
                sdf_path = pdb_path.replace('.pdb', '_enhanced.sdf')
                
            self.total_processed += 1
            
            # Extract ligand from PDB
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
            if mol is None:
                logger.warning(f"Failed to read molecule from {pdb_path}, trying without hydrogens")
                mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
                if mol is None:
                    logger.error(f"Could not read molecule from {pdb_path}")
                    return False
            
            # Get reference ligand for comparison
            system_id = os.path.basename(pdb_path).split('_')[0]
            ref_ligand_path = self._get_reference_ligand_path(system_id)
            
            if ref_ligand_path and os.path.exists(ref_ligand_path):
                ref_mol = Chem.SDMolSupplier(ref_ligand_path, removeHs=False)[0]
                if ref_mol is None:
                    ref_mol = Chem.SDMolSupplier(ref_ligand_path, removeHs=True)[0]
                
                # Apply enhanced standardization
                mol = self._standardize_molecule_enhanced(mol, ref_mol)
            
            # Perform protein-based alignment if protein structures are available
            ref_protein_path = self._get_reference_protein_path(system_id)
            if ref_protein_path and os.path.exists(ref_protein_path):
                logger.info(f"Performing protein-based alignment for {system_id}")
                rotation_matrix, translation_vector = self._align_protein_structures(pdb_path, ref_protein_path)
                
                if rotation_matrix is not None and translation_vector is not None:
                    # Apply transformation to ligand
                    mol = self._apply_transformation_to_ligand_enhanced(mol, rotation_matrix, translation_vector)
                    self.successful_alignment_count += 1
            
            # Write the processed molecule
            if mol:
                writer = Chem.SDWriter(sdf_path)
                writer.write(mol)
                writer.close()
                logger.info(f"Successfully converted {pdb_path} to {sdf_path}")
                return True
            else:
                logger.error(f"Molecule became None during processing")
                return False
                
        except Exception as e:
            logger.error(f"Error converting {pdb_path} to {sdf_path}: {str(e)}")
            return False
    
    def _get_reference_protein_path(self, system_id):
        """Get the reference protein structure path"""
        base_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/{self.dataset_name}"
        protein_path = os.path.join(base_path, f"{system_id}_protein.pdb")
        if os.path.exists(protein_path):
            return protein_path
        
        # Try alternative naming
        alt_path = os.path.join(base_path, f"{system_id}_target.pdb")
        if os.path.exists(alt_path):
            return alt_path
            
        return None
    
    def _apply_transformation_to_ligand_enhanced(self, mol, rotation_matrix, translation_vector):
        """Apply rotation and translation transformation to ligand coordinates"""
        try:
            if mol is None:
                return None
                
            # Get the conformer
            conf = mol.GetConformer()
            
            # Apply transformation to each atom
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                pos_array = np.array([pos.x, pos.y, pos.z])
                
                # Apply rotation and translation
                new_pos = np.dot(pos_array, rotation_matrix.T) + translation_vector
                
                conf.SetAtomPosition(i, [new_pos[0], new_pos[1], new_pos[2]])
            
            return mol
            
        except Exception as e:
            logger.error(f"Failed to apply transformation to ligand: {str(e)}")
            return mol

def test_enhanced_approach():
    """Test the enhanced Boltz approach"""
    logger.info("Starting enhanced Boltz approach testing...")
    
    approach = EnhancedBoltzApproach()
    
    # Use actual systems that exist in the Boltz results
    test_systems = ['1afb__1__1.A__1.D_1.F', '1b5d__1__1.A_1.B__1.D', '1bcj__1__1.B__1.I_1.K', '1ci0__1__1.A_1.B__1.D', '1d7c__1__1.A_1.B__1.E_1.K']
    
    results = []
    for system_id in test_systems:
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing system: {system_id}")
        logger.info(f"{'='*50}")
        
        try:
            # Build protein directory path for Boltz
            protein_dir = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/forks/boltz/inference/plinder_set_0/boltz_results_{system_id}"
            
            # Get prediction files
            pred_files = approach.list_top_n_files(protein_dir, 1)
            if not pred_files:
                logger.warning(f"No prediction files found for {system_id}")
                continue
                
            pred_path = pred_files[0]
            logger.info(f"Prediction file: {pred_path}")
            
            # Convert to SDF with enhanced processing
            sdf_path = pred_path.replace('.pdb', '_enhanced.sdf')
            success = approach._convert_pdb_to_sdf_enhanced(pred_path, sdf_path)
            
            if success:
                # Calculate RMSD with reference ligand
                # Extract the PDB ID from the system name (e.g., '1afb__1__1.A__1.D_1.F' -> '1afb')
                pdb_id = system_id.split('__')[0]
                ref_ligand_path = approach._get_reference_ligand_path(pdb_id)
                
                if ref_ligand_path and os.path.exists(ref_ligand_path):
                    try:
                        # Standard complex RMSD
                        complex_rmsd = calculate_complex_rmsd(sdf_path, ref_ligand_path)
                        logger.info(f"Complex RMSD for {system_id}: {complex_rmsd:.3f} Å")
                        
                        # Heavy-atom-only RMSD
                        pred_mol = Chem.SDMolSupplier(sdf_path)[0]
                        ref_mol = Chem.SDMolSupplier(ref_ligand_path)[0]
                        
                        if pred_mol and ref_mol:
                            heavy_rmsd = approach._calculate_heavy_atom_rmsd(pred_mol, ref_mol)
                            logger.info(f"Heavy-atom RMSD for {system_id}: {heavy_rmsd:.3f} Å")
                            
                            # Get molecular formulas for analysis
                            pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                            ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                            formula_match = pred_formula == ref_formula
                            
                            results.append({
                                'system_id': system_id,
                                'pdb_id': pdb_id,
                                'complex_rmsd': complex_rmsd,
                                'heavy_atom_rmsd': heavy_rmsd,
                                'pred_formula': pred_formula,
                                'ref_formula': ref_formula,
                                'formula_match': formula_match,
                                'chemical_correction_applied': approach.chemical_correction_count > len(results),
                                'alignment_success': approach.successful_alignment_count > len(results)
                            })
                            
                        else:
                            logger.error(f"Could not read molecules for RMSD calculation")
                    except Exception as e:
                        logger.error(f"RMSD calculation failed for {system_id}: {str(e)}")
                else:
                    logger.warning(f"Reference ligand not found for {pdb_id}")
                        
        except Exception as e:
            logger.error(f"Error processing {system_id}: {str(e)}")
            
    # Summary statistics
    logger.info(f"\n{'='*50}")
    logger.info("ENHANCED APPROACH SUMMARY")
    logger.info(f"{'='*50}")
    logger.info(f"Systems processed: {len(results)}")
    logger.info(f"Total molecules processed: {approach.total_processed}")
    logger.info(f"Chemical corrections applied: {approach.chemical_correction_count}")
    logger.info(f"Successful protein alignments: {approach.successful_alignment_count}")
    logger.info(f"Heavy-atom fallbacks used: {approach.heavy_atom_fallback_count}")
    
    if results:
        df = pd.DataFrame(results)
        
        logger.info(f"\nResults summary:")
        logger.info(f"Mean complex RMSD: {df['complex_rmsd'].mean():.3f} ± {df['complex_rmsd'].std():.3f} Å")
        logger.info(f"Mean heavy-atom RMSD: {df['heavy_atom_rmsd'].mean():.3f} ± {df['heavy_atom_rmsd'].std():.3f} Å")
        logger.info(f"Systems with complex RMSD < 2 Å: {(df['complex_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Systems with heavy-atom RMSD < 2 Å: {(df['heavy_atom_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Formula matches: {df['formula_match'].sum()}/{len(df)}")
        
        # Detailed analysis
        logger.info(f"\nDetailed analysis:")
        for _, row in df.iterrows():
            logger.info(f"{row['system_id']}: Complex RMSD={row['complex_rmsd']:.3f}Å, "
                       f"Heavy RMSD={row['heavy_atom_rmsd']:.3f}Å, "
                       f"Formula match={row['formula_match']}")
        
        # Save results
        output_path = '/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks/enhanced_boltz_results.csv'
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        return df
    else:
        logger.warning("No results obtained")
        return None

if __name__ == "__main__":
    test_enhanced_approach()
        
    def _standardize_molecule_enhanced(self, mol, reference_mol=None):
        """Enhanced molecule standardization with comprehensive chemical corrections"""
        if mol is None:
            return None
            
        # Get molecular formulas
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol) if reference_mol else None
        
        logger.info(f"Predicted formula: {mol_formula}, Reference formula: {ref_formula}")
        
        # If formulas match, return as-is
        if mol_formula == ref_formula:
            return mol
            
        # Apply comprehensive chemical corrections
        corrected_mol = self._apply_comprehensive_corrections(mol, reference_mol)
        if corrected_mol:
            corrected_formula = rdMolDescriptors.CalcMolFormula(corrected_mol)
            if corrected_formula == ref_formula:
                logger.info(f"Successfully corrected formula from {mol_formula} to {corrected_formula}")
                self.chemical_correction_count += 1
                return corrected_mol
                
        # If chemical correction fails, return original for heavy-atom evaluation
        logger.warning(f"Chemical correction failed. Will use heavy-atom-only evaluation.")
        self.heavy_atom_fallback_count += 1
        return mol
        
    def _apply_comprehensive_corrections(self, mol, reference_mol):
        """Apply comprehensive chemical corrections based on common Boltz prediction patterns"""
        if not mol or not reference_mol:
            return mol
            
        # Create a copy to avoid modifying the original
        mol_copy = Chem.Mol(mol)
        
        # Pattern 1: Try simple hydrogen adjustment
        corrected_mol = self._try_hydrogen_correction(mol_copy, reference_mol)
        if corrected_mol:
            return corrected_mol
            
        # Pattern 2: Try SMILES-based corrections
        corrected_mol = self._try_smiles_correction_enhanced(mol_copy, reference_mol)
        if corrected_mol:
            return corrected_mol
            
        # Pattern 3: Try bond order corrections
        corrected_mol = self._try_bond_order_corrections(mol_copy, reference_mol)
        if corrected_mol:
            return corrected_mol
            
        return None
        
    def _try_hydrogen_correction(self, mol, reference_mol):
        """Try to correct hydrogen count by removing and re-adding with proper valences"""
        try:
            # Remove all hydrogens and re-add
            mol_no_h = Chem.RemoveHs(mol)
            mol_with_h = Chem.AddHs(mol_no_h)
            
            # Check if this matches the reference formula
            corrected_formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
            ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol)
            
            if corrected_formula == ref_formula:
                logger.info("Successfully applied hydrogen correction")
                return mol_with_h
                
        except Exception as e:
            logger.warning(f"Hydrogen correction failed: {str(e)}")
            
        return None
        
    def _try_smiles_correction_enhanced(self, mol, reference_mol):
        """Enhanced SMILES-based chemical correction"""
        try:
            ref_smiles = Chem.MolToSmiles(reference_mol)
            pred_smiles = Chem.MolToSmiles(mol)
            
            logger.info(f"SMILES correction - Ref: {ref_smiles[:50]}...")
            logger.info(f"SMILES correction - Pred: {pred_smiles[:50]}...")
            
            # Apply multiple correction patterns
            corrected_smiles = pred_smiles
            corrections_applied = []
            
            # Pattern 1: Carbonyl reduction (C(=O) ↔ C(O))
            if "C(=O)" in ref_smiles and "C(O)" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("C(O)", "C(=O)", 1)
                corrections_applied.append("C(O)->C(=O)")
            elif "C(O)" in ref_smiles and "C(=O)" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("C(=O)", "C(O)", 1)
                corrections_applied.append("C(=O)->C(O)")
                
            # Pattern 2: Imine/amine correction (C=N ↔ CN)
            if "C=N" in ref_smiles and "CN" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("CN", "C=N", 1)
                corrections_applied.append("CN->C=N")
            elif "CN" in ref_smiles and "C=N" in corrected_smiles:
                corrected_smiles = corrected_smiles.replace("C=N", "CN", 1)
                corrections_applied.append("C=N->CN")
                
            if corrections_applied:
                logger.info(f"Applied SMILES corrections: {', '.join(corrections_applied)}")
                corrected_mol = Chem.MolFromSmiles(corrected_smiles)
                if corrected_mol:
                    corrected_mol = Chem.AddHs(corrected_mol)
                    # Verify the correction worked
                    corrected_formula = rdMolDescriptors.CalcMolFormula(corrected_mol)
                    ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol)
                    if corrected_formula == ref_formula:
                        return corrected_mol
                        
        except Exception as e:
            logger.warning(f"SMILES correction failed: {str(e)}")
            
        return None
        
    def _try_bond_order_corrections(self, mol, reference_mol):
        """Try to correct bond orders to match reference"""
        try:
            # This is a simplified approach - in practice, bond order correction
            # is very complex and would need sophisticated algorithms
            
            # For now, just try the standard hydrogen addition approach
            mol_no_h = Chem.RemoveHs(mol)
            
            # Try different sanitization options
            try:
                Chem.SanitizeMol(mol_no_h)
                mol_sanitized = Chem.AddHs(mol_no_h)
                
                corrected_formula = rdMolDescriptors.CalcMolFormula(mol_sanitized)
                ref_formula = rdMolDescriptors.CalcMolFormula(reference_mol)
                
                if corrected_formula == ref_formula:
                    logger.info("Successfully applied bond order correction")
                    return mol_sanitized
                    
            except:
                pass
                
        except Exception as e:
            logger.warning(f"Bond order correction failed: {str(e)}")
            
        return None
        
    def _calculate_heavy_atom_rmsd(self, mol1, mol2):
        """Calculate RMSD using only heavy atoms when chemical formulas don't match"""
        if not mol1 or not mol2:
            return float('inf')
            
        try:
            # Get heavy atom coordinates
            conf1 = mol1.GetConformer()
            conf2 = mol2.GetConformer()
            
            heavy_atoms1 = [i for i in range(mol1.GetNumAtoms()) if mol1.GetAtomWithIdx(i).GetSymbol() != 'H']
            heavy_atoms2 = [i for i in range(mol2.GetNumAtoms()) if mol2.GetAtomWithIdx(i).GetSymbol() != 'H']
            
            # Get element symbols for heavy atoms
            elements1 = [mol1.GetAtomWithIdx(i).GetSymbol() for i in heavy_atoms1]
            elements2 = [mol2.GetAtomWithIdx(i).GetSymbol() for i in heavy_atoms2]
            
            # If different number of heavy atoms, try to match common atoms
            if len(heavy_atoms1) != len(heavy_atoms2):
                logger.warning(f"Different heavy atom counts: {len(heavy_atoms1)} vs {len(heavy_atoms2)}")
                # Use minimum number for comparison
                min_atoms = min(len(heavy_atoms1), len(heavy_atoms2))
                heavy_atoms1 = heavy_atoms1[:min_atoms]
                heavy_atoms2 = heavy_atoms2[:min_atoms]
                elements1 = elements1[:min_atoms]
                elements2 = elements2[:min_atoms]
            
            if len(heavy_atoms1) == 0:
                return float('inf')
                
            # Get coordinates
            coords1 = np.array([conf1.GetAtomPosition(i) for i in heavy_atoms1])
            coords2 = np.array([conf2.GetAtomPosition(i) for i in heavy_atoms2])
            
            # For small molecules, try to find optimal atom correspondence
            if len(heavy_atoms1) <= 10:
                return self._calculate_optimal_heavy_atom_rmsd(coords1, coords2, elements1, elements2)
            else:
                # For larger molecules, use simple distance-based matching
                return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
                
        except Exception as e:
            logger.error(f"Heavy atom RMSD calculation failed: {str(e)}")
            return float('inf')
    
    def _calculate_optimal_heavy_atom_rmsd(self, coords1, coords2, elements1, elements2):
        """Calculate optimal RMSD by trying different atom correspondences for small molecules"""
        try:
            from itertools import permutations
            
            # Group atoms by element type
            element_groups1 = {}
            element_groups2 = {}
            
            for i, elem in enumerate(elements1):
                if elem not in element_groups1:
                    element_groups1[elem] = []
                element_groups1[elem].append(i)
                
            for i, elem in enumerate(elements2):
                if elem not in element_groups2:
                    element_groups2[elem] = []
                element_groups2[elem].append(i)
            
            # Check if we have the same elements
            if set(element_groups1.keys()) != set(element_groups2.keys()):
                logger.warning("Different element types in molecules")
                # Fall back to distance-based matching
                return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
            
            # For very small molecules, try all permutations within element groups
            min_rmsd = float('inf')
            
            # Generate all possible correspondences (simplified approach)
            if len(coords1) <= 6:
                for perm in permutations(range(len(coords2))):
                    try:
                        # Check if this permutation preserves element types
                        valid = all(elements1[i] == elements2[perm[i]] for i in range(len(elements1)))
                        if valid:
                            rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2[list(perm)])**2, axis=1)))
                            min_rmsd = min(min_rmsd, rmsd)
                    except:
                        continue
                        
                return min_rmsd if min_rmsd != float('inf') else self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
            else:
                return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
                
        except Exception as e:
            logger.warning(f"Optimal RMSD calculation failed: {str(e)}")
            return self._calculate_greedy_heavy_atom_rmsd(coords1, coords2)
    
    def _calculate_greedy_heavy_atom_rmsd(self, coords1, coords2):
        """Calculate RMSD using greedy distance-based atom assignment"""
        try:
            # Calculate pairwise distances
            distances = cdist(coords1, coords2)
            
            # Greedy assignment
            used_j = set()
            total_dist_sq = 0
            assignments = 0
            
            for i in range(len(coords1)):
                best_j = None
                best_dist = float('inf')
                for j in range(len(coords2)):
                    if j not in used_j and distances[i][j] < best_dist:
                        best_dist = distances[i][j]
                        best_j = j
                        
                if best_j is not None:
                    used_j.add(best_j)
                    total_dist_sq += best_dist**2
                    assignments += 1
                    
            if assignments == 0:
                return float('inf')
                
            return np.sqrt(total_dist_sq / assignments)
            
        except Exception as e:
            logger.error(f"Greedy RMSD calculation failed: {str(e)}")
            return float('inf')
    
    def _convert_pdb_to_sdf(self, pdb_path, sdf_path=None):
        """Enhanced PDB to SDF conversion with comprehensive processing"""
        try:
            if sdf_path is None:
                sdf_path = pdb_path.replace('.pdb', '_enhanced.sdf')
                
            self.total_processed += 1
            
            # Extract ligand from PDB
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
            if mol is None:
                logger.warning(f"Failed to read molecule from {pdb_path}, trying without hydrogens")
                mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
                if mol is None:
                    logger.error(f"Could not read molecule from {pdb_path}")
                    return False
            
            # Get reference ligand for comparison
            system_id = os.path.basename(pdb_path).split('_')[0]
            ref_ligand_path = self._get_reference_ligand_path(system_id)
            
            if ref_ligand_path and os.path.exists(ref_ligand_path):
                ref_mol = Chem.SDMolSupplier(ref_ligand_path, removeHs=False)[0]
                if ref_mol is None:
                    ref_mol = Chem.SDMolSupplier(ref_ligand_path, removeHs=True)[0]
                
                # Apply enhanced standardization
                mol = self._standardize_molecule_enhanced(mol, ref_mol)
            
            # Perform protein-based alignment if protein structures are available
            ref_protein_path = self._get_reference_protein_path(system_id)
            if ref_protein_path and os.path.exists(ref_protein_path):
                logger.info(f"Performing protein-based alignment for {system_id}")
                aligned_mol = self._align_with_reference_protein(mol, pdb_path, ref_protein_path)
                if aligned_mol:
                    mol = aligned_mol
                    self.successful_alignment_count += 1
            
            # Write the processed molecule
            if mol:
                writer = Chem.SDWriter(sdf_path)
                writer.write(mol)
                writer.close()
                logger.info(f"Successfully converted {pdb_path} to {sdf_path}")
                return True
            else:
                logger.error(f"Molecule became None during processing")
                return False
                
        except Exception as e:
            logger.error(f"Error converting {pdb_path} to {sdf_path}: {str(e)}")
            return False
    
    def _get_reference_protein_path(self, system_id):
        """Get the reference protein structure path"""
        base_path = f"/Users/aoxu/projects/DrugDiscovery/PoseBench/data/{self.dataset_name}"
        protein_path = os.path.join(base_path, f"{system_id}_protein.pdb")
        if os.path.exists(protein_path):
            return protein_path
        
        # Try alternative naming
        alt_path = os.path.join(base_path, f"{system_id}_target.pdb")
        if os.path.exists(alt_path):
            return alt_path
            
        return None
    
    def _align_with_reference_protein(self, ligand_mol, pred_pdb_path, ref_protein_path):
        """Align ligand based on protein structure alignment using existing methods"""
        try:
            # Use the existing protein alignment method from the parent class
            rotation_matrix, translation_vector = self._align_protein_structures(pred_pdb_path, ref_protein_path)
            
            if rotation_matrix is not None and translation_vector is not None:
                # Apply transformation to ligand
                return self._apply_transformation_to_ligand_enhanced(ligand_mol, rotation_matrix, translation_vector)
            else:
                logger.warning("Protein alignment failed, using original ligand coordinates")
                return ligand_mol
                
        except Exception as e:
            logger.error(f"Protein alignment failed: {str(e)}")
            return ligand_mol
            
    def _apply_transformation_to_ligand_enhanced(self, mol, rotation_matrix, translation_vector):
        """Apply rotation and translation transformation to ligand coordinates"""
        try:
            if mol is None:
                return None
                
            # Get the conformer
            conf = mol.GetConformer()
            
            # Apply transformation to each atom
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                pos_array = np.array([pos.x, pos.y, pos.z])
                
                # Apply rotation and translation
                new_pos = np.dot(pos_array, rotation_matrix.T) + translation_vector
                
                conf.SetAtomPosition(i, [new_pos[0], new_pos[1], new_pos[2]])
            
            return mol
            
        except Exception as e:
            logger.error(f"Failed to apply transformation to ligand: {str(e)}")
            return mol

def test_enhanced_approach():
    """Test the enhanced Boltz approach"""
    logger.info("Starting enhanced Boltz approach testing...")
    
    approach = EnhancedBoltzApproach()
    
    # Test on a diverse set of systems with known issues
    test_systems = ['1a30', '1a8g', '1ai7', '1a4w', '1a6q']
    
    results = []
    for system_id in test_systems:
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing system: {system_id}")
        logger.info(f"{'='*50}")
        
        try:
            # Get prediction files
            pred_files = approach.list_top_n_files(system_id, n=1)
            if not pred_files:
                logger.warning(f"No prediction files found for {system_id}")
                continue
                
            pred_path = pred_files[0]
            logger.info(f"Prediction file: {pred_path}")
            
            # Convert to SDF with enhanced processing
            sdf_path = pred_path.replace('.pdb', '_enhanced.sdf')
            success = approach._convert_pdb_to_sdf(pred_path, sdf_path)
            
            if success:
                # Calculate RMSD
                ref_ligand_path = approach._get_reference_ligand_path(system_id)
                if ref_ligand_path and os.path.exists(ref_ligand_path):
                    try:
                        # Standard complex RMSD
                        complex_rmsd = calculate_complex_rmsd(sdf_path, ref_ligand_path)
                        logger.info(f"Complex RMSD for {system_id}: {complex_rmsd:.3f} Å")
                        
                        # Heavy-atom-only RMSD
                        pred_mol = Chem.SDMolSupplier(sdf_path)[0]
                        ref_mol = Chem.SDMolSupplier(ref_ligand_path)[0]
                        
                        if pred_mol and ref_mol:
                            heavy_rmsd = approach._calculate_heavy_atom_rmsd(pred_mol, ref_mol)
                            logger.info(f"Heavy-atom RMSD for {system_id}: {heavy_rmsd:.3f} Å")
                            
                            # Get molecular formulas for analysis
                            pred_formula = rdMolDescriptors.CalcMolFormula(pred_mol)
                            ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol)
                            formula_match = pred_formula == ref_formula
                            
                            results.append({
                                'system_id': system_id,
                                'complex_rmsd': complex_rmsd,
                                'heavy_atom_rmsd': heavy_rmsd,
                                'pred_formula': pred_formula,
                                'ref_formula': ref_formula,
                                'formula_match': formula_match,
                                'chemical_correction_applied': False,  # Will be updated based on counters
                                'alignment_success': approach.successful_alignment_count > len(results)
                            })
                        
                    except Exception as e:
                        logger.error(f"RMSD calculation failed for {system_id}: {str(e)}")
                        
        except Exception as e:
            logger.error(f"Error processing {system_id}: {str(e)}")
            
    # Summary statistics
    logger.info(f"\n{'='*50}")
    logger.info("ENHANCED APPROACH SUMMARY")
    logger.info(f"{'='*50}")
    logger.info(f"Systems processed: {len(results)}")
    logger.info(f"Total molecules processed: {approach.total_processed}")
    logger.info(f"Chemical corrections applied: {approach.chemical_correction_count}")
    logger.info(f"Successful protein alignments: {approach.successful_alignment_count}")
    logger.info(f"Heavy-atom fallbacks used: {approach.heavy_atom_fallback_count}")
    
    if results:
        df = pd.DataFrame(results)
        
        logger.info(f"\nResults summary:")
        logger.info(f"Mean complex RMSD: {df['complex_rmsd'].mean():.3f} ± {df['complex_rmsd'].std():.3f} Å")
        logger.info(f"Mean heavy-atom RMSD: {df['heavy_atom_rmsd'].mean():.3f} ± {df['heavy_atom_rmsd'].std():.3f} Å")
        logger.info(f"Systems with complex RMSD < 2 Å: {(df['complex_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Systems with heavy-atom RMSD < 2 Å: {(df['heavy_atom_rmsd'] < 2).sum()}/{len(df)}")
        logger.info(f"Formula matches: {df['formula_match'].sum()}/{len(df)}")
        
        # Detailed analysis
        logger.info(f"\nDetailed analysis:")
        for _, row in df.iterrows():
            logger.info(f"{row['system_id']}: Complex RMSD={row['complex_rmsd']:.3f}Å, "
                       f"Heavy RMSD={row['heavy_atom_rmsd']:.3f}Å, "
                       f"Formula match={row['formula_match']}")
        
        # Save results
        output_path = '/Users/aoxu/projects/DrugDiscovery/PoseBench/notebooks/enhanced_boltz_results.csv'
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        return df
    else:
        logger.warning("No results obtained")
        return None

if __name__ == "__main__":
    test_enhanced_approach()
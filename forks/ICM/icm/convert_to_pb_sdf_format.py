from rdkit import Chem
import os

def convert_one_sdf_format(input_file, output_file):
    """
    Convert a complex SDF file to a simpler format by extracting just the first molecule
    and its basic structure information using RDKit.
    
    Args:
        input_file (str): Path to the input SDF file
        output_file (str): Path to save the converted SDF file
    """
    try:
        # Read the input SDF file
        supplier = Chem.SDMolSupplier(input_file)
        
        # Get the first molecule
        mol = next(supplier)
        
        if mol is None:
            print("Error: No valid molecule found in the input file")
            return
            
        # Create an SDWriter for the output file
        writer = Chem.SDWriter(output_file)
        
        # Write the molecule to the new file
        writer.write(mol)
        
        # Close the writer
        writer.close()
        
        print(f"Successfully converted {input_file} to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file {input_file}")
    except Exception as e:
        print(f"Error occurred during conversion: {str(e)}")

def convert_sdf_format(base_dir):
    for dirname in os.listdir(base_dir):
        dir_path = os.path.join(base_dir, dirname)
        sdf_path = os.path.join(dir_path, f"v{dirname}_ligand_docked.sdf")
        output_path = os.path.join(dir_path, f"{dirname}.sdf")
        
        if not os.path.exists(sdf_path):
            print(f"Skipping '{dirname}': Missing SDF file {sdf_path}")
            continue

        convert_one_sdf_format(sdf_path, output_path)
    return

def main():
    # Example usage
    # base_dir = "7TOM_5AD"  # Change this to your directory path
    # input_file = os.path.join(base_dir, "v7TOM_5AD_ligand_docked.sdf")
    # output_file = os.path.join(base_dir, "7TOM_5AD.sdf")
    
    # # Create directory if it doesn't exist
    # os.makedirs(base_dir, exist_ok=True)
    
    # # Convert the file
    # convert_sdf_format(input_file, output_file)
    base_dir = "/home/aoxu/projects/PoseBench/forks/ICM/inference/icm_auto_posebusters_benchmark_outputs_1"
    convert_sdf_format(base_dir)

if __name__ == "__main__":
    main()
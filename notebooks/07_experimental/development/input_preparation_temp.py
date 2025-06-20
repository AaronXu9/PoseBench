import csv

input_file = '/home/aoxu/projects/PoseBench/forks/DiffDock/inference/diffdock_posebusters_benchmark_inputs_orig.csv'
output_file = '/home/aoxu/projects/PoseBench/forks/DiffDock/inference/diffdock_pocket_only_posebusters_benchmark_inputs_orig.csv'

with open(input_file, mode='r') as infile, open(output_file, mode='w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    # Write the header
    header = next(reader)
    writer.writerow(header)
    
    for row in reader:
        # Modify the protein_path
        complex_name = row[0]
        new_protein_path = f"data/posebusters_benchmark_set/posebusters_benchmark_set_orig_bs_cropped/{complex_name}_protein_bs_cropped.pdb"
        row[1] = new_protein_path
        
        # Write the modified row to the new file
        writer.writerow(row)
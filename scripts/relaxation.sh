#!/bin/bash

cd /home/aoxu/projects/PoseBench
# Redirect all output to a log file
exec > >(tee -i logs/relaxation_script_output.log)
exec 2>&1

# Diffdock relax and create the bust_results by analyze

# dynamicbind relax and create the bust_results by analyze
python3 posebench/models/inference_relaxation.py method=dynamicbind dataset=posebusters_benchmark remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=dynamicbind dataset=posebusters_benchmark repeat_index=4
python3 posebench/models/inference_relaxation.py method=dynamicbind dataset=posebusters_benchmark remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=5
python3 posebench/analysis/inference_analysis.py method=dynamicbind dataset=posebusters_benchmark repeat_index=5

# fabind relax and create the bust_results by analyze
python3 posebench/models/inference_relaxation.py method=fabind dataset=posebusters_benchmark remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=5
python3 posebench/analysis/inference_analysis.py method=fabind dataset=posebusters_benchmark repeat_index=5

# neuralplexer relax and create the bust_results by analyze
python3 posebench/models/inference_relaxation.py method=neuralplexer dataset=posebusters_benchmark num_processes=1 remove_initial_protein_hydrogens=true assign_partial_charges_manually=true cache_files=false repeat_index=4
python3 posebench/analysis/complex_alignment.py method=neuralplexer dataset=posebusters_benchmark repeat_index=4
python3 posebench/analysis/inference_analysis.py method=neuralplexer dataset=posebusters_benchmark repeat_index=4

python3 posebench/models/inference_relaxation.py method=neuralplexer dataset=posebusters_benchmark num_processes=1 remove_initial_protein_hydrogens=true assign_partial_charges_manually=true cache_files=false repeat_index=5
python3 posebench/analysis/complex_alignment.py method=neuralplexer dataset=posebusters_benchmark repeat_index=5
python3 posebench/analysis/inference_analysis.py method=neuralplexer dataset=posebusters_benchmark repeat_index=5

# autodock vina relax and create the bust_results by analyze
mkdir -p forks/Vina/inference/vina_diffdock_posebusters_benchmark_outputs_4 && cp -r data/test_cases/posebusters_benchmark/vina_diffdock_posebusters_benchmark_outputs_4/* forks/Vina/inference/vina_diffdock_posebusters_benchmark_outputs_4
python3 posebench/models/inference_relaxation.py method=vina vina_binding_site_method=diffdock dataset=posebusters_benchmark remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=vina vina_binding_site_method=diffdock dataset=posebusters_benchmark repeat_index=4

python3 posebench/models/vina_inference.py dataset=posebusters_benchmark method=diffdock repeat_index=5 # NOTE: DiffDock-L's binding pockets are recommended as the default Vina input
mkdir -p forks/Vina/inference/vina_diffdock_posebusters_benchmark_outputs_5 && cp -r data/test_cases/posebusters_benchmark/vina_diffdock_posebusters_benchmark_outputs_5/* forks/Vina/inference/vina_diffdock_posebusters_benchmark_outputs_5
python3 posebench/models/inference_relaxation.py method=vina vina_binding_site_method=diffdock dataset=posebusters_benchmark remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=5
python3 posebench/analysis/inference_analysis.py method=vina vina_binding_site_method=diffdock dataset=posebusters_benchmark repeat_index=5

# tulip relax and create the bust_results by analyze
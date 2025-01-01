#!/bin/bash

cd /home/aoxu/projects/PoseBench
# Redirect all output to a log file
exec > >(tee -i logs/relaxation_script_output.log)
exec 2>&1

# Diffdock relax and create the bust_results by analyze
# python3 posebench/models/diffdock_inference.py dataset=astex_diverse repeat_index=4
python3 posebench/models/inference_relaxation.py method=diffdock dataset=astex_diverse remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=diffdock dataset=astex_diverse repeat_index=4

# dynamicbind relax and create the bust_results by analyze
python3 posebench/data/dynamicbind_input_preparation.py dataset=astex_diverse
python3 posebench/models/dynamicbind_inference.py dataset=astex_diverse repeat_index=4
python3 posebench/models/inference_relaxation.py method=dynamicbind dataset=astex_diverse remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=dynamicbind dataset=astex_diverse repeat_index=4

# fabind relax and create the bust_results by analyze
python3 posebench/data/fabind_input_preparation.py dataset=astex_diverse
python3 posebench/models/fabind_inference.py dataset=astex_diverse repeat_index=4
python3 posebench/models/inference_relaxation.py method=fabind dataset=astex_diverse remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=fabind dataset=astex_diverse repeat_index=4

# neuralplexer relax and create the bust_results by analyze
python3 posebench/data/neuralplexer_input_preparation.py dataset=astex_diverse
python3 posebench/models/neuralplexer_inference.py dataset=astex_diverse repeat_index=4
python3 posebench/models/inference_relaxation.py method=neuralplexer dataset=astex_diverse num_processes=1 remove_initial_protein_hydrogens=true assign_partial_charges_manually=true cache_files=false repeat_index=4
python3 posebench/analysis/complex_alignment.py method=neuralplexer dataset=astex_diverse repeat_index=4python3 posebench/analysis/inference_analysis.py method=neuralplexer dataset=astex_diverse repeat_index=4

# autodock vina relax and create the bust_results by analyze
cp forks/DiffDock/inference/diffdock_astex_diverse_inputs.csv forks/Vina/inference/vina_astex_diverse_inputs.csv
python3 posebench/models/vina_inference.py dataset=astex_diverse method=diffdock repeat_index=4
mkdir -p forks/Vina/inference/vina_diffdock_astex_diverse_outputs_4 && cp -r data/test_cases/astex_diverse/vina_diffdock_astex_diverse_outputs_4/* forks/Vina/inference/vina_diffdock_astex_diverse_outputs_4
python3 posebench/models/inference_relaxation.py method=vina vina_binding_site_method=diffdock dataset=astex_diverse remove_initial_protein_hydrogens=true assign_partial_charges_manually=true num_processes=1 repeat_index=4
python3 posebench/analysis/inference_analysis.py method=vina vina_binding_site_method=diffdock dataset=astex_diverse repeat_index=4

# tulip relax and create the bust_results by analyze
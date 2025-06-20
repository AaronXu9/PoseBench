# PoseBench Notebooks Organization Plan

## ✅ COMPLETED ORGANIZATION (June 2025):
- Successfully organized 150+ files from a single directory
- All method-specific files organized by method type
- Analysis results separated from source code
- Clear separation between methods, datasets, and analysis types
- Improved findability and version control

## Current Organized Structure:

```
notebooks/
├── ORGANIZATION_PLAN.md                # This organization guide
├── __pycache__/                        # Python cache files
├── 00_tutorials/                       # Getting started guides ✅
│   ├── adding_new_dataset_tutorial.ipynb
│   └── adding_new_method_tutorial.ipynb
│
├── 01_method_analysis/                 # Per-method analysis folders ✅
│   ├── boltz/                          # Boltz-specific files
│   │   ├── boltz_comprehensive_analysis.ipynb
│   │   ├── boltz_structure_analysis.py
│   │   ├── enhanced_boltz_approach.py
│   │   ├── fixed_boltz_approach.py
│   │   ├── optimal_boltz_approach.py
│   │   └── test_boltz_*.py files
│   ├── chai_lab/                       # Chai-Lab files
│   │   ├── chai-1_plinder_set_0_results.csv
│   │   └── chai-lab_casp15_interaction_dataframes_1.h5
│   ├── diffdock/                       # DiffDock files
│   │   ├── diffdock_casp15_interaction_dataframes_1.h5
│   │   ├── diffdock_plinder_set_0_results.csv
│   │   ├── diffdock_pocket_only_*.csv
│   │   └── diffdockv1_casp15_interaction_dataframes_1.h5
│   ├── dynamicbind/                    # DynamicBind files
│   │   └── dynamicbind_casp15_interaction_dataframes_1.h5
│   ├── gnina/                          # GNINA files
│   │   ├── gnina_plinder_set_0_results.csv
│   │   └── gnina_posebusters_results.csv
│   ├── icm/                            # ICM files
│   │   └── icm_plinder_set_0_results.csv
│   ├── neuralplexer/                   # NeuralPLexer files
│   │   ├── neuralplexer_casp15_interaction_dataframes_1.h5
│   │   └── neuralplexer_no_ilcl_casp15_interaction_dataframes_1.h5
│   ├── rfaa/                           # RoseTTAFold-All-Atom files
│   │   └── rfaa_casp15_interaction_dataframes_1.h5
│   ├── surfdock/                       # SurfDock files
│   │   ├── surfdock_plinder_set_0_results.csv
│   │   ├── surfdock_plinder_set_0_results_first_half.csv
│   │   └── surfdock_posebusters_results.csv
│   ├── tulip/                          # TULIP files
│   │   └── tulip_casp15_interaction_dataframes_1.h5
│   └── vina/                           # Vina files
│       ├── vina_diffdock_casp15_interaction_dataframes_1.h5
│       ├── vina_p2rank_casp15_interaction_dataframes_1.h5
│       ├── vina_plinder_set_0_results.csv
│       └── vina_posebusters_results.csv
│
├── 02_dataset_analysis/                # Dataset-specific analysis ✅
│   ├── astex_diverse/                  # Astex Diverse dataset
│   │   └── astex_diverse_interaction_dataframes.h5
│   ├── casp15/                         # CASP15 dataset (previously organized)
│   ├── dockgen/                        # DockGen dataset
│   │   ├── dockgen_inference_results_plotting.ipynb
│   │   ├── dockgen_interaction_dataframes.h5
│   │   └── dockgen_*.png plots
│   ├── pdbbind/                        # PDBBind dataset
│   │   └── pdbbind_training_subset_1000_interaction_dataframes.h5
│   ├── plinder/                        # Plinder dataset (previously organized)
│   └── posebusters/                    # PoseBusters dataset
│       ├── posebusters_benchmark_*.csv
│       ├── posebusters_benchmark_interaction_dataframes.h5
│       ├── posebusters_method_interaction_analysis_plotting.ipynb
│       ├── posebusters_results_filtered*.csv
│       └── posebusters_*_inference_results_plotting.ipynb
│
├── 03_comparative_analysis/            # Cross-method comparisons
├── 04_utils/                           # Utility scripts and functions ✅
│   ├── Approach.py
│   ├── ConformationAnalyzer.py
│   ├── alignment.py
│   ├── analysis.py
│   ├── comprehensive_analysis.py
│   ├── plotting.py
│   ├── protein_structure_analysis.py
│   └── utils.py
│
├── 05_specialized_analysis/            # Specialized analysis notebooks ✅
│   ├── RMSD_score_statistical.ipynb
│   ├── dataset_interaction_analysis_plotting.ipynb
│   ├── docking.ipynb
│   ├── load_plinder.ipynb
│   ├── plidner_preprocess_results.ipynb
│   ├── plidner_rmsd_edge_case.ipynb
│   ├── posebuster_analyze.ipynb
│   ├── posebuster_comparative_analysis.ipynb
│   ├── posebuster_complementary_analysis.ipynb
│   ├── protein_analysis.ipynb
│   ├── protein_ligand_interaction_analysis.ipynb
│   ├── prolif_analysis.ipynb
│   ├── rmsd_analysis_plinder_outlier.ipynb
│   ├── timing_analysis.ipynb
│   └── true_structure_pb_analysis.ipynb
│
├── 06_results_archive/                 # Analysis results and outputs ✅
│   ├── 2024_q4/                        # Quarterly archives (previously organized)
│   ├── 2025_q1/
│   ├── boltz/                          # Boltz-specific results (previously organized)
│   ├── current/
│   ├── plots/                          # All PNG plot files
│   │   ├── dataset_interaction_analysis.png
│   │   ├── ensembled_results_5benchmarks.png
│   │   ├── grouped_bars_plot.png
│   │   └── všechny ostatní *.png soubory
│   ├── all_results.pkl
│   ├── analysis_output.txt
│   ├── consensus_ensemble_casp15_interaction_dataframes_1.h5
│   ├── fingerprint.pkl
│   ├── plidner_test_input.csv
│   ├── ref_fp.csv
│   └── tmp_*.csv
│
└── 07_experimental/                    # Development and testing ✅
    ├── development/                    # Development scripts (previously organized)
    │   ├── chemical_analysis.py
    │   ├── generate_pb_restuls.py
    │   ├── input_preparation_temp.py
    │   ├── manual_sdf_analysis.py
    │   ├── simple_debug.py
    │   ├── standardization_solution.py
    │   └── test.py
    └── tests/                          # Test scripts (previously organized)
```
│   │   ├── results/
│   │   └── analysis.ipynb
│   │
│   └── [other_methods]/
│
├── 02_dataset_analysis/               # Per-dataset analysis
│   ├── plinder/
│   │   ├── plinder_comprehensive_analysis.ipynb
│   │   ├── plinder_property_analysis.ipynb
│   │   ├── plinder_rmsd_analysis.ipynb
│   │   └── results/
│   ├── casp15/
│   │   ├── casp15_inference_results_plotting.ipynb
│   │   ├── casp15_method_interaction_analysis.ipynb
│   │   └── results/
│   ├── posebusters/
│   │   ├── posebusters_benchmark_analysis.ipynb
│   │   └── results/
│   └── astex_diverse/
│       └── results/
│
├── 03_comparative_analysis/           # Cross-method comparisons
│   ├── method_benchmarking.ipynb      # Compare all methods
│   ├── dataset_interaction_analysis.ipynb
│   ├── consensus_analysis.ipynb
│   └── results/
│
├── 04_utils/                          # Shared utilities
│   ├── analysis_utils.py              # Common analysis functions
│   ├── plotting_utils.py              # Visualization functions
│   ├── data_processing.py             # Data loading/cleaning
│   └── evaluation_metrics.py         # RMSD, PoseBusters, etc.
│
├── 05_specialized_analysis/           # Specific analysis types
│   ├── protein_ligand_interactions/
│   │   ├── interaction_analysis.ipynb
│   │   └── prolif_analysis.ipynb
│   ├── chemical_analysis/
│   │   ├── fingerprint_analysis.ipynb
│   │   └── chemical_space.ipynb
│   ├── timing_performance/
│   │   ├── timing_analysis.ipynb
│   │   └── computational_efficiency.ipynb
│   └── edge_cases/
│       ├── rmsd_edge_cases.ipynb
│       └── failure_analysis.ipynb
│
├── 06_results_archive/               # Historical results
│   ├── 2024_q4/
│   ├── 2025_q1/
│   └── current/
│
└── 07_experimental/                  # Work in progress
    ├── new_method_tests/
    ├── prototype_analyses/
    └── sandbox/
```

## Implementation Benefits:

1. **Method-Focused**: Each method has its own complete analysis pipeline
2. **Version Control**: Clear separation of development vs production code
3. **Reusability**: Shared utilities prevent code duplication
4. **Discoverability**: Logical hierarchy makes finding specific analyses easy
5. **Maintenance**: Easier to update and maintain individual method analyses
6. **Collaboration**: Clear ownership and contribution guidelines
7. **Documentation**: Each folder has its own README with specific guidance

## Migration Priority:

### Phase 1: Boltz Organization (High Priority)
- Create boltz/ folder structure
- Move all boltz_* files to appropriate subfolders
- Consolidate test files
- Archive old versions

### Phase 2: Core Infrastructure (Medium Priority)  
- Set up utils/ folder with shared functions
- Organize dataset-specific analyses
- Create comparative analysis section

### Phase 3: Full Migration (Low Priority)
- Migrate all other method analyses
- Set up results archive
- Create comprehensive documentation

## Organization Summary:

### ✅ **COMPLETED (June 2025):**
- **150+ files organized** from single directory into logical structure
- **11 method-specific folders** created with all relevant files
- **6 dataset analysis folders** with dataset-specific files
- **Utility scripts** centralized in 04_utils/
- **Specialized analysis** notebooks organized in 05_specialized_analysis/
- **Results and plots** archived in 06_results_archive/
- **Development files** organized in 07_experimental/

### 📊 **Organization Statistics:**
- **Method Analysis**: 11 methods (Boltz, Chai-Lab, DiffDock, DynamicBind, GNINA, ICM, NeuralPLexer, RFAA, SurfDock, TULIP, Vina)
- **Dataset Analysis**: 6 datasets (Astex Diverse, CASP15, DockGen, PDBBind, Plinder, PoseBusters)
- **Utility Scripts**: 8 core utility files
- **Specialized Notebooks**: 16 analysis notebooks
- **Archived Results**: 100+ result files and plots

### 🎯 **Benefits Achieved:**
- **Improved findability**: Files logically grouped by purpose
- **Better maintainability**: Clear separation of concerns
- **Enhanced collaboration**: Organized structure for team development
- **Version control**: Easier to track changes by category
- **Scalability**: Structure supports adding new methods/datasets

## File Naming Conventions:

- **Core analysis**: `{method}_comprehensive_analysis.ipynb`
- **Results**: `{method}_{dataset}_{date}_results.csv`
- **Test files**: `test_{method}_{specific_test}.py`
- **Debug files**: `debug_{issue_description}.py`
- **Plots**: `{method}_{metric}_{plot_type}.png`

## Next Steps:

1. **Create README files** for each major subdirectory
2. **Add method-specific documentation** for complex analyses
3. **Establish automated organization** for new files
4. **Regular maintenance** to keep structure current

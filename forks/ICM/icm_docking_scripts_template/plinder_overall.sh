
icm/def> openFile "/home/aoxu/projects/PoseBench/forks/ICM/inference/plinder_test_cases/1afb__1__1_A__1_D_1_F/1afb__1__1_A__1_D_1_F.icb" 0 yes no no no " append" 
1 1 icb 


icm/p_1afb__1__1_A__1_D_1_F_rec> S_out = Sarray( s_projectsDir + "icmdock/" ) 
icm/p_1afb__1__1_A__1_D_1_F_rec> cool Mol( as_graph )
icm/p_1afb__1__1_A__1_D_1_F_rec>  s_out = ( Nof(as_graph & a_1afb__1__1_A__1_D_1_F_protein. ) > 0 )?" Error> the pocket is NOT EMPTY, move ligand out first!":"" 
icm/p_1afb__1__1_A__1_D_1_F_rec>  as_graph = Sphere(as_graph a_1afb__1__1_A__1_D_1_F_protein. 5. ) 
icm/p_1afb__1__1_A__1_D_1_F_rec>  print s_out 
 
icm/p_1afb__1__1_A__1_D_1_F_rec> if ( ! Exist( "/home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0/" ) ) make directory "/home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0/" 
icm/p_1afb__1__1_A__1_D_1_F_rec> set directory "/home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0/" 
 Info> current directory is set to /home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0
icm/p_1afb__1__1_A__1_D_1_F_rec> currentDockProj.data[8] = "yes" 
icm/p_1afb__1__1_A__1_D_1_F_rec> if( Type( as_graph ) == "aselection" ) tempsel = as_graph & a_1afb__1__1_A__1_D_1_F_protein. 
icm/p_1afb__1__1_A__1_D_1_F_rec> if( Type( as_graph ) == "unknown" ) tempsel = as_graph & a_1afb__1__1_A__1_D_1_F_protein. 
icm/p_1afb__1__1_A__1_D_1_F_rec> dock2SetupReceptor "p_1afb__1__1_A__1_D_1_F" a_1afb__1__1_A__1_D_1_F_protein. tempsel yes no ? "enhancedBorn" : "none" 
Two following receptor setup steps are: 
1. adjustment of the initial ligand position; 2. adjustment of the box size/position. 
1. If necessary, re-orient the yellow probe. Hold SHIFT for global rotation. 
Press 'ENTER' or click 'Go' to continue.
If necessary, adjust the size/position of the box around the binding site 
(hold LEFT MOUSE BUTTON with the cursor at any corner of the box). 
Press 'ENTER' or click 'Go' to continue.

Receptor object has been created: a_p_1afb__1__1_A__1_D_1_F_rec.
Next:
Setup ligand using any of the following menu options 
        'Docking.Ligand Setup.From Loaded ICM Object' 
        'Docking.Ligand Setup.From Loaded non-ICM Object' 
        'Docking.Ligand Setup.From File (ICM)' 
        'Docking.Ligand Setup.From File (Mol/Mol2)' 
        'Docking.Ligand Setup.From Database' 
Summary appended to the log file /home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0//p_1afb__1__1_A__1_D_1_F_LOG.htm
Started:  03-15-2025 01:19:26
Completed:  03-15-2025 01:19:3
Total time, min: 0.15
icm/p_1afb__1__1_A__1_D_1_F_rec> if( yes & currentDockProj.l_readyReceptor ) dock5CalcMaps "p_1afb__1__1_A__1_D_1_F" 0.5 4.0 no 
Calculation of maps with grid size 0.50 A and terms gh,gc,ge,gb,gs,gl...
 m_gh... done 
 Info> Map m_gh created. GridStep=0.50 Dimensions: 32 40 40, Size=51200

Maps ready: grid size: 0.50 A, terms: gh,gc,ge,gb,gs,gl.
run virtual screening ICM script : _dockScan
To screen the database 
NOTSET
run the following UNIX command: 
    icm /home/aoxu/icm-3.9-4/_dockScan p_1afb__1__1_A__1_D_1_F & 
(make sure that the database
 NOTSET 
is indexed before your start). 
Summary appended to the log file /home/aoxu/projects/PoseBench/forks/ICM/ICM_manual_docking_maps/plinder_set_0//p_1afb__1__1_A__1_D_1_F_LOG.htm
Started:  03-15-2025 01:19:35
Completed:  03-15-2025 01:19:3
Total time, min: 0.02
icm/p_1afb__1__1_A__1_D_1_F_rec> currentDockProj.data[1] = "p_1afb__1__1_A__1_D_1_F" 
icm/p_1afb__1__1_A__1_D_1_F_rec> s_out = Path() 
icm/p_1afb__1__1_A__1_D_1_F_rec>  S_out = Name( model, "DOCKING" ) 
icm/p_1afb__1__1_A__1_D_1_F_rec> dockMakeInputTable currentDockProj.data[1] Name( v1afb__1__1_A__1_D_1_F_ligand table ) 
 Info> 1 column  'docked' added to table 'v1afb__1__1_A__1_D_1_F_ligand'
 Info> 1 column  'IX' added to table 'v1afb__1__1_A__1_D_1_F_ligand'
 Info> 1 header  'toolsPanel' added to table 'v1afb__1__1_A__1_D_1_F_ligand'
 Info> 1 header  's_results' added to table 'v1afb__1__1_A__1_D_1_F_ligand'
 Info> 1 header  's_projName' added to table 'v1afb__1__1_A__1_D_1_F_ligand'
 Info> table  v1afb__1__1_A__1_D_1_F_ligand  (s_out) (4 headers, 4 arrays (i_out) of 1 rows) has been created
  Headers: v1afb__1__1_A__1_D_1_F_ligand.toolsPanel v1afb__1__1_A__1_D_1_F_ligand.s_results v1afb__1__1_A__1_D_1_F_ligand.s_projName v1afb__1__1_A__1_D_1_F_ligand.l_busyDocking
  Arrays : v1afb__1__1_A__1_D_1_F_ligand.mol v1afb__1__1_A__1_D_1_F_ligand.NAME_ v1afb__1__1_A__1_D_1_F_ligand.docked v1afb__1__1_A__1_D_1_F_ligand.IX
 Info> s_projName,s_cursor,s_chemtable,l_embed temp.variables deleted
icm/p_1afb__1__1_A__1_D_1_F_rec> dockInputTable Name( v1afb__1__1_A__1_D_1_F_ligand table ) 
/home/aoxu/icm-3.9-4/icm64 /home/aoxu/icm-3.9-4/_dockScan -s -a thorough=3. name=v1afb__1__1_A__1_D_1_F_ligand input=v1afb__1__1_A__1_D_1_F_ligand_input.sdf -S confs=3 p_1afb__1__1_A__1_D_1_F 
icm/p_1afb__1__1_A__1_D_1_F_rec> v1afb__1__1_A__1_D_1_F_ligand.l_busyDocking = no
icm/p_1afb__1__1_A__1_D_1_F_rec> dockProcessAnswersFromInputTable "v1afb__1__1_A__1_D_1_F_ligand"
Total number of solutions to browse:  3 
Processing file p_1afb__1__1_A__1_D_1_F_v1afb__1__1_A__1_D_1_F_ligand1.ob ..
 
  no
icm/p_1afb__1__1_A__1_D_1_F_rec> 

# Running the model (Computational Model of the Whisker Thalamocortical Pathway)

## Description
Detailed computational model of the whisker thalamocortical pathway of the mouse developed in NetPyNE ([www.netpyne.org](https://www.netpyne.org/)).

The model is described in the following publication:

- Moreira JVS, Borges FS, Atherton Z, Crandall S, Varela C, Dura-Bernal S (2025) TITLE. JOURNAL.

## Cloning the model
- If you follow the steps in `README_installation.md`, you should be able to clone and run the computational model in your local machine (not recommended), or in an HPC environment (recommended).  
- We do not recommend attempting to run the full scale model locally, because it is very computationally expensive. A single run of the model using a integration step `dt = 0.025 ms` for `16 seconds` requires on average 2h35 minutes using 64 cores, and 128Gb of RAM memory.

## Compiling the MOD files
- We provide a custom script to compile the MOD files used to create the ionic currents and synaptic mechanisms (`compile_unified.sh` or `compile_unified_HPC.sh`, for MAC local computer or Linux HPC, respectively), but the user can also perform the compilation manualy by running the command `nrnivmodl` in the terminal from within the `\mod` folder.  
    (1) open a terminal window inside `/mod`  
    (2) run `nrnivmodl`  
    or  
    (2) run `./compile_unified.sh`  
- Make sure the compilation was successful by checking the output of the `nrnivmodl` for any messages indicating a failure to compile any mechanism.
- **Pro tip**: if you are not familiar with compilation of MOD files, AI assistants can be very helpful in this proccess. Paste the output of the compilation into the chat and have it help you debug the possible mistakes in the compilation.

## Running the model to reproduce the data

To run the model and reproduce the figures associated with [this publication](add link), you can run the command:
- `python barreloid_batch.py`  
This will run all the simulations necessary to produce the data used to make the figures. However, this should not be performed in a local computer, since it will run 203 simulations, each requiring 128Gb of RAM and with a runtime of 2h35min using 64 cores.  

Following that, the user should run the script `./run_copy_files_to_folders.sh` ( `./run_copy_files_to_folders.sh ../data/batch_sims/PATH_TO_DATASET`) followed by the path of the data folder to generate the subfolder `sim_output`. This folder has copies of the lighter datasets saved during the simulation, and it is used to gather the data from the HPC and generate the figures locally.

- e.g.: `./run_copy_files_to_folders.sh ../data/batch_sims/barreloid_batch_fig02_sleep_wake/`
- output: the following subfolder will be created, with all the datasets `../data/batch_sims/barreloid_batch_fig02_sleep_wake/sim_output`

To skip this process of having to run the simulations and post-process the data, a simplified dataset with the data used to make the figures can be found in the folder `/sim_data`.  
- e.g.: `../sim_data/paper_000/paper_000_barreloid_batch_fig02_sleep_wake/sim_output/`

The scripts to run the plotting functions can be found in the folder `\sim_analysis`, and all figures can be generated simultaneously by running the script `run_analysis.sh`. This script runs all the plotting scripts found in `\sim_analysis`, and generates the figures for each simulation.

This way, the model can be verified both by (1) re-running the simulations from scrath, or (2) using the datasets provided (see `README_analysis.md`).

## (1) Running the model to make modifications/improvements

This model is an open source platform, and we encourage that other use this to make further explorations and test their scientific questions.

To re-run the model and produce the datasets used in this study, the user can re-run the script `barreloid_batch.py`. 
- The file was simplified in a way that all analysis are executed at once, using 60 cores in an HPC environment.
    - each analysis is submited for running using the commands at the end of the script:
    ```
    # RUN SIMULATIONS FOR FIGURE 02
    # # === Sleep-Wake === #
    batchLabel = 'barreloid_batch_fig02_sleep_wake' ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    # RUN SIMULATIONS FOR FIGURE 05
    # === Awake === #
    batchLabel = 'barreloid_batch_fig05_awake'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    # RUN SIMULATIONS FOR FIGURE 06
    # === Sleep === #
    batchLabel = 'barreloid_batch_fig06_sleep'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    ## RUN SIMULATIONS FOR FIGURE 07
    # === Awake CT modulation === # - Closed-loop CT
    #  -  Uniform  -  #
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_0'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_10'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_25'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_50'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_100'    ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Mixed  -  #
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_0'        ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_10'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_25'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_50'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_100'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Closed  -  #
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_0'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_10'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_25'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_50'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_100'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    
    # === Awake CT modulation === # - UNIFORM CT
    #  -  Uniform  -  #
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_0'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_10'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_25'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_50'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_100'    ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Mixed  -  #
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_0'        ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_10'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_25'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_50'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_100'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Closed  -  #
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_0'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_10'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_25'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_50'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_100'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)    
    ```

- To add a new analysis/simulation condition in the model, the user should modify the `barreloid_batch.py` script and configure a new entry into the `network_state_dict` dictionary.
    - e.g.:  
    >elif batch_label == 'barreloid_batch_fig02_sleep_wake':
    >>network_state_dict={
    >>>'sleep__baseline': {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[6.0,  6.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_upDown2, 'replace_ct_dict': None},
    >>>'wake__baseline_10Hz':{'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None},  
    >>}
- obs.: the names of the keys in the dictionary must either include a `wake__` or `sleep__` flag, so that the code knows which set of parameters to use. These can be further configured to suit the needs of your specific implementation.

- you should also select the connectivity strategy used, as defined in the `conn_schema_flags` (uniform, open-loop, closed-loop, etc). These are predefined and loaded from datasets that store the connectivity template for each network arrangement.
- e.g.: 
    > if batch_label == 'barreloid_batch_fig02_sleep_wake':
    >>conn_schema_flags=[     
    >>>['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],  
    >>]
    
The user can also run the `barreloid_init.py`, which uses default paramters and will output a single run of the model. We do not encourage using this approach, since the parameter tuning was performed using the batch file, and all updated parameters are stored there. But if that is the choice, then make sure to modify the parameters to reproduce the condition desired.

## (2) Running the datasets to reproduce the model's data
See `README_analysis.md`

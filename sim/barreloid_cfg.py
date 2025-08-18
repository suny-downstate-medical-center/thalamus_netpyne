"""
barreloid_cfg.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""

from netpyne import specs
import pickle
import NetPyNE_BBP
import json
import sys
import os
import numpy as np
import Build_Net as BN

saveFigFlag = '' # declared here, but altered afterwards


#------------------------------------------------------------------------------
#   Main config
#------------------------------------------------------------------------------
cfg       = specs.SimConfig()     # object of class SimConfig to store simulation configuration

cfg.model_checkpoint = 'stable_version__04__2024_10_15__post_SfN'
cfg.model_checkpoint_comments = 'Current version compares 4 types of intrathalamic feedback - uniform, topological (open), closed, and mixed.\nWe also reproduce Waxing-waning by increasing wGABA and slightly depolarizing VPM and TRN neurons\nThis version was presented at SfN 2024'

# (Bazhenov, 2002) - Rescaling of Ca2+-dependent K+ current
cfg.map_targetRMP = { 'VPL_TC': -68, 'Rt_RC': -70}

####################################################################################################################################################################################
# --- Modify currents for all Thalamus
####################################################################################################################################################################################
cfg.rescale_SK_E2   = 1 # default (None - no rescaling) - should not be changed, except when testing this parameter
cfg.rescale_pas     = 1 # default (None - no rescaling) - should not be changed, except when testing this parameter

# new variables
cfg.rescale_iT      = 1
cfg.rescale_ih      = 1
cfg.rescale_iA      = 1
cfg.rescale_iNap    = 1

cfg.add_KLeak       = 0 # (default: 0 - model has no active KLeak)

cfg.add_iT_shift__Thalamus  = 0 # (mV) - shifts iT activation and inactivation curves for all thalamic neurons

####################################################################################################################################################################################
# --- Modify currents by Thalamic pop
####################################################################################################################################################################################
cfg.rescale_SK_E2__VPM  = 1
cfg.rescale_SK_E2__TRN  = 1

cfg.rescale_ih_gmax__VPM= 1
cfg.rescale_ih_gmax__TRN= 1
cfg.add_ih_shift__VPM   = 0

cfg.add_iT_shift__VPM       = 0 # (mV) - shifts iT activation and inactivation curves for VPM neurons
cfg.add_iT_shift__TRN       = 0 # (mV) - shifts iT activation and inactivation curves for TRN neurons

cfg.rescale_pas__VPM    = 1
cfg.modify_pas_e__VPM   = -80                       # default value (-80 mV)
cfg.modify_pas_g__VPM   = 3.4702549429081374e-05    # default value (3.4702549429081374e-05 )
cfg.rescale_pas__TRN    = 1
cfg.modify_pas_e__TRN   = -80                       # default value (-80 mV)
cfg.modify_pas_g__TRN   = 8.617446501142974e-05     # default value (8.617446501142974e-05 )

# Oscillating:  (for a reference rmp of -55 mV)
cfg.add_KLeak__VPM = 0 # (default:0 - added as a posteriori modification) scaling factor for default current of (g = 1.0e-5 (S/cm2)) - for comparison, default g_pas ("g": 3.4702549429081374e-05)
# # Tonic:        (for a reference rmp of -55 mV)
# cfg.add_KLeak__VPM = 2 # (default:0 - added as a posteriori modification) scaling factor for default current of (g = 1.0e-5 (S/cm2)) - for comparison, default g_pas ("g": 3.4702549429081374e-05)

cfg.add_KLeak__TRN = 0 # (default:0 - added as a posteriori modification) scaling factor for default current of (g = 1.0e-5 (S/cm2)) - for comparison, default g_pas ("g": 8.617446501142974e-05)

####################################################################################################################################################################################

# cfg.target_angle = None
# cfg.target_angle = 0
# cfg.target_angle = 45
# cfg.target_angle = 90
cfg.target_angle = 135
# cfg.target_angle = 180 
# cfg.target_angle = 225
# cfg.target_angle = 270
# cfg.target_angle = 315

cfg.sim_tag = 'bWgt_9241__post_SfN_tests__'+str(cfg.target_angle)+''

#------------------------------------------------------------------------------
#   Simulation options
#------------------------------------------------------------------------------

cfg.duration        = 16000 # (ms)
cfg.dt              = 0.025                 # Internal integration timestep to use
cfg.printRunTime    = 0.1 # (s)
cfg.printPopAvgRates = True

cfg.hParams         = {'v_init': -60, 'celsius': 34} # changing to -60 mV to remove initial bursts in "in vivo"-like sim
cfg.verbose         = False                 # Show detailed messages

# --- Must be (connRandomSecFromList=True) and (distributeSynsUniformly=False), so that conns are picked from a secList, and the list changes everytime
cfg.connRandomSecFromList   = True  # default, so that 'inputs__distal', 'inputs__intermediate', 'inputs__proximal' are randomized
cfg.distributeSynsUniformly = False

cfg.use_fast_imem = True  # use CVode fast_imem to record membrane voltage via i_membrane_

cfg.base_random_seed = 100000

cfg.rand123GlobalIndex = cfg.base_random_seed

cfg.seeds = {'conn': cfg.base_random_seed, 
             'stim': cfg.base_random_seed, 
             'loc':  cfg.base_random_seed, 
             'cell': cfg.base_random_seed}

cfg.distributedSaving = False

cfg.includeParamsLabel = True

#------------------------------------------------------------------------------
# Recording
#------------------------------------------------------------------------------
cfg.center_point        = 500           # (X-Z center of the network)

cfg.recordStep = 0.1              # Step size in ms to save data (eg. V traces, LFP, etc)
# cfg.recordStep = cfg.dt         # Step size in ms to save data (eg. V traces, LFP, etc)

cfg.recordTraces = {    'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'},
                        # 'V_ptr': {'var': 'ptr'},
                        }

cfg.recordCellMembraneI=False
if cfg.recordCellMembraneI: cfg.recordTraces.update({'I_soma':       {'sec':'soma_0',  'loc':0.5, 'var':'i_membrane_'}})

cfg.recordCurrents=False
rec_curr = [('SK_E2','ik'),('TC_HH','ina'),('TC_HH','ik'),('TC_Nap_Et2','ina'),
            ('TC_iA','ik'),('TC_iL','ica'),('TC_iT_Des98','ica'),('TC_ih_Bud97','ih')]
if cfg.recordCurrents:
    for curr in rec_curr: cfg.recordTraces.update({'i__soma_0__'+curr[0]+'__'+curr[1]:{'sec':'soma_0','loc':0.5,'mech':curr[0],'var':curr[1]},})

'''
get set of currents:
mech_names=[]
for cell_name in netParams.cellParams.keys():
    for sec in netParams.cellParams[cell_name]['secs'].keys():
        mech_names.extend(list(netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys()))
mech_names_set=list(set(mech_names))
out:    
    ['TC_Nap_Et2', 'TC_cad', 'SK_E2', 'TC_Kleak', 'TC_iA', 'TC_iL', 'TC_iT_Des98', 'TC_ih_Bud97', 'TC_HH', 'pas']
'''


cfg.plotLFP=True

# Removing LFP recording because it is breaking the code when more morphologies are added
cfg.recordLFP = [
    [cfg.center_point,2675,cfg.center_point],
    [cfg.center_point,3650,cfg.center_point]
]
cfg.saveLFPPops = ['VPM__pop','TRN__pop']

#------------------------------------------------------------------------------
# Path configuration
#------------------------------------------------------------------------------
cfg.convertCellMorphologies = False
cfg.loadCellModel           = True
cfg.convertCellModel        = False
cfg.saveCellModel           = True
cfg.plotSynLocation         = True

# cfg.base_dir = '/Users/joao'
cfg.base_dir=os.path.expanduser("~")

cfg.BBP_rootFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'
cfg.morphologyFolder_asc            = cfg.BBP_rootFolder+'/sonata/morphologies/morphologies_asc'

cfg.NetPyNE_rootFolder              = cfg.base_dir+'/Research/Models/BBP/thalamus_netpyne'
cfg.NetPyNE_JSON_cells              = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies'
cfg.NetPyNE_templateCells           = cfg.NetPyNE_rootFolder+'/mod'
cfg.NetPyNE_exportedCells           = cfg.NetPyNE_rootFolder+'/cells/morphologies_swc'
cfg.NetPyNE_L6A_JSON_cells          = cfg.NetPyNE_rootFolder+'/cells/S1_BBP_cells/'
cfg.NetPyNE_network_template        = cfg.NetPyNE_rootFolder+'/conn/barreloid_network_template/network_template.json'

cfg.prepare_HPC_mode = False # prepares the project to run in HPC environments
cfg.run_HPC_mode     = True # configure the script to run in HPC environments
if cfg.prepare_HPC_mode or cfg.run_HPC_mode: 
    # cfg.NetPyNE_JSON_cells+='_barreloid'
    cfg.NetPyNE_JSON_cells=cfg.NetPyNE_rootFolder+'/cells/netpyne_validated_thalamus_morphologies'

cfg.loadCircuitProperties           = True
cfg.saveCircuitProperties           = True
cfg.stored_circuit_path             = cfg.NetPyNE_rootFolder+'/conn/bbp_circuit_propeties/circuit_dict.pkl'

cfg.BBP_conn_properties             = cfg.NetPyNE_rootFolder+'/conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'

cfg.dump_cell_properties            = cfg.NetPyNE_rootFolder+'/network_template/barreloid_network_cell_properties.json' # file that stores the network template used in <updateNetwork> and <ResampleParameters.processSimObj>

# --- Path to BBP data files
cfg.Fig_4_1_dataFolder              = cfg.NetPyNE_rootFolder+'/conn/bbp_input_data/Fig4_wakefulness_spontaneous'
cfg.Fig_4_1_spikes                  = cfg.Fig_4_1_dataFolder+'/out.h5'
# --- Load Artificial CT Spike Times
cfg.ct_spk_files_dict={
                    'uniform':      '../stims/Artificial/CT/v3_spiking_dict__nCells_818__sampling_25_uS__freq_5_Hz__downTimes__None__.json',
                }

####################################
cfg.replace_spk_file = None # (Default for not replacing spikes)
# # # Replace with 360 angular tuned CT cells
# cfg.replace_spk_file = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes4_deflectionAngle|180_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'

####################################

cfg.ct_downTimes = None

# --- Alternative CT spikes
cfg.ct_spk_file             = '../stims/Artificial/CT/spiking_dict__nCells_1000__sampling_25_uS__freq_10_Hz__downTimes__1000|1500_2000|2500_3000|3500__.json'

#------------------------------------------------------------------------------
# Simulation Configuration
#------------------------------------------------------------------------------

cfg.modType = 'Prob_original'

# cfg.randomizeNetwork=True
cfg.randomizeNetwork= False

cfg.simplify_gids   = True # simplifying the number of morphologies added

calibrated_currents_files=[ '../cells/calibrated_currents/calibrated_currents.json',
                            '../cells/calibrated_currents/calibrated_currents2.json',
                            '../cells/calibrated_currents/calibrated_currents3.json',
                            '../cells/calibrated_currents/calibrated_currents4.json',
                            ]
map_ThresholdCurrent={}
for calibrated_currents_file in calibrated_currents_files:
    try:
        with open(calibrated_currents_file, 'r') as file: map_ThresholdCurrent_temp = json.load(file)
        for new_cell in map_ThresholdCurrent_temp.keys():
            if new_cell not in map_ThresholdCurrent.keys():map_ThresholdCurrent.update({new_cell:map_ThresholdCurrent_temp[new_cell]})
    except:
        print('Failed to load file ',calibrated_currents_file)      

# # --- Simplifying mophologies to test the effects of connectivity
# # # TC: 36219
# # TC: 36963
# # RT: 31864 (cNAD) - non adapting
# select_morphologies=[36963,31864]
# map_ThresholdCurrent_ = {key:map_ThresholdCurrent[key] for key in map_ThresholdCurrent for select_morphology in select_morphologies if str(select_morphology) in key}
# map_ThresholdCurrent=map_ThresholdCurrent_

''' 
Adapting TRN cells - burst and then accomodate firing
31738 31191 32846 33325
'''
# # --- Simplifying mophologies to test the effects of connectivity
# # # TC: 36219
# # TC: 36963
# # RT: 31738 (cAD) -  adapting
# select_morphologies=[36963,31738]
# map_ThresholdCurrent_ = {key:map_ThresholdCurrent[key] for key in map_ThresholdCurrent for select_morphology in select_morphologies if str(select_morphology) in key}
# map_ThresholdCurrent=map_ThresholdCurrent_

'''
cNAD with less afterhyperpolarization
31924, 
29476, 
29934, 30838, 30891, 30905, 31400, 31528, 32267
'''
# # --- Simplifying mophologies to test the effects of connectivity
# # # TC: 36219
# # TC: 36963
# # RT: 31924 (cNAD) -  non-adapting w/ less afterhyperpolarization
# select_morphologies=[36963,31924]
# map_ThresholdCurrent_ = {key:map_ThresholdCurrent[key] for key in map_ThresholdCurrent for select_morphology in select_morphologies if str(select_morphology) in key}
# map_ThresholdCurrent=map_ThresholdCurrent_

# select_morphologies=[36963,31924] # original morphs - work
# select_morphologies=[34349,28925] # original morphs - work
# select_morphologies=[36219,30252] # original morphs - work
# select_morphologies=[37520,31315] # original morphs - work
# select_morphologies=[41912,32049] # original morphs - work

# testing adding more morphologies
select_morphologies=[
                        # VPM
                        36963, 
                        34349, 
                        36219, 37520, 41912,     
                        
                        # TRN
                        31924, 
                        28925, 
                        30252, 31315, 32049,
                     ]     
map_ThresholdCurrent_ = {key:map_ThresholdCurrent[key] for key in map_ThresholdCurrent for select_morphology in select_morphologies if str(select_morphology) in key}
map_ThresholdCurrent=map_ThresholdCurrent_

if cfg.randomizeNetwork:
    NetPyNE_BBP.Prompt.headerMsg('Randomizing network gids')
    cfg.VPM_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='VPL_TC', cfg_file=cfg.sonataConfigFile, gids = cfg.load_vpm_gids)
    cfg.TRN_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='Rt_RC',  cfg_file=cfg.sonataConfigFile, gids = cfg.load_trn_gids)
else:
    NetPyNE_BBP.Prompt.headerMsg('Loading preset network gids')
    if cfg.simplify_gids:
        print('\t>>\tAdding varied network morphologies - best 5 VPM and 5 TRN')

        cfg.VPM_gids = [int(cell_thresh.split('__')[1]) for cell_thresh in map_ThresholdCurrent.keys() if 'VPL_TC' in cell_thresh]
        cfg.TRN_gids = [int(cell_thresh.split('__')[1]) for cell_thresh in map_ThresholdCurrent.keys() if 'Rt_RC'  in cell_thresh]

        # cfg.VPM_gids=[34349, 36043, 36219, 37520, 41912]
        # cfg.TRN_gids=[28925, 
        #             #   29892, # traces look weird
        #               30252, 31315, 32049]

        # # # VPM (318) GIDs picked from: [40352, 34066, 37812, 38737, 35050]
        # # cfg.VPM_gids=[35050, 37812, 38737, 37812, 34066, 38737, 38737, 40352, 37812, 34066, 34066, 40352, 35050, 40352, 34066, 34066, 40352, 35050, 37812, 38737, 34066, 34066, 40352, 38737, 37812, 34066, 37812, 37812, 34066, 38737, 35050, 34066, 34066, 34066, 40352, 38737, 37812, 40352, 40352, 38737, 34066, 38737, 37812, 35050, 35050, 37812, 35050, 34066, 40352, 38737, 35050, 35050, 34066, 40352, 37812, 37812, 34066, 35050, 40352, 38737, 35050, 35050, 37812, 34066, 37812, 38737, 40352, 38737, 35050, 38737, 37812, 34066, 37812, 34066, 40352, 35050, 38737, 34066, 35050, 40352, 35050, 34066, 38737, 38737, 40352, 40352, 38737, 34066, 40352, 35050, 34066, 35050, 38737, 35050, 40352, 34066, 40352, 37812, 34066, 35050, 38737, 40352, 34066, 38737, 38737, 35050, 40352, 40352, 34066, 37812, 40352, 34066, 38737, 38737, 34066, 40352, 34066, 34066, 38737, 35050, 37812, 37812, 38737, 35050, 38737, 40352, 37812, 34066, 37812, 34066, 34066, 35050, 34066, 40352, 34066, 40352, 40352, 37812, 37812, 34066, 40352, 34066, 38737, 38737, 40352, 34066, 40352, 34066, 40352, 34066, 38737, 40352, 40352, 34066, 40352, 34066, 37812, 37812, 38737, 40352, 34066, 35050, 35050, 35050, 37812, 38737, 38737, 38737, 35050, 40352, 38737, 37812, 34066, 40352, 35050, 37812, 37812, 34066, 34066, 34066, 40352, 38737, 40352, 40352, 37812, 40352, 40352, 37812, 35050, 40352, 34066, 34066, 35050, 38737, 37812, 34066, 37812, 35050, 35050, 37812, 38737, 34066, 38737, 35050, 40352, 40352, 37812, 38737, 35050, 35050, 38737, 40352, 37812, 38737, 37812, 40352, 37812, 35050, 38737, 40352, 34066, 37812, 38737, 38737, 35050, 35050, 37812, 40352, 35050, 35050, 40352, 40352, 35050, 35050, 35050, 34066, 34066, 37812, 38737, 40352, 40352, 34066, 35050, 40352, 35050, 37812, 35050, 38737, 38737, 35050, 38737, 38737, 37812, 40352, 37812, 35050, 35050, 35050, 37812, 40352, 40352, 37812, 35050, 35050, 34066, 38737, 35050, 38737, 37812, 34066, 37812, 37812, 38737, 37812, 34066, 37812, 38737, 34066, 38737, 34066, 38737, 34066, 35050, 40352, 35050, 40352, 37812, 40352, 40352, 37812, 34066, 34066, 34066, 35050, 38737, 34066, 38737, 37812, 38737, 37812, 34066, 38737, 35050, 40352, 34066, 40352, 37812, 38737, 40352, 35050, 40352, 34066, 37812, 37812, 38737, 34066, 38737, 37812]
        # # # TRN (101) GIDs picked from: [30404, 32411, 31136, 32952, 29223]
        # # cfg.TRN_gids=[32411, 31136, 32411, 31136, 30404, 29223, 32952, 31136, 32411, 32411, 29223, 32952, 30404, 32952, 31136, 29223, 29223, 31136, 29223, 29223, 29223, 32411, 32952, 32952, 31136, 32952, 32952, 32411, 32952, 32952, 29223, 29223, 30404, 29223, 31136, 30404, 32411, 32952, 32952, 30404, 32411, 31136, 32411, 32952, 32952, 32952, 30404, 30404, 29223, 32952, 32411, 29223, 32411, 29223, 30404, 29223, 29223, 32411, 32411, 32411, 32411, 29223, 30404, 30404, 32952, 32952, 29223, 29223, 30404, 32411, 32411, 30404, 30404, 30404, 31136, 30404, 32411, 32411, 32952, 32411, 32411, 29223, 30404, 30404, 30404, 29223, 32411, 31136, 29223, 32952, 29223, 32952, 32952, 32411, 32952, 32411, 32411, 32411, 31136, 31136, 29223]
    else:
        cfg.VPM_gids=[40352, 34066, 36777, 40166, 37812, 38737, 40296, 39900, 35050, 36666, 42272, 38785, 40826, 42023, 33933, 34466, 38284, 34844, 40465, 39867, 41216, 42326, 39167, 33887, 33559, 36658, 36388, 38608, 36299, 39299, 39251, 41563, 39633, 37330, 34195, 40932, 41549, 37891, 36698, 37851, 37428, 36855, 37478, 40842, 35711, 36579, 40249, 37405, 37301, 42351, 34692, 36310, 35985, 41386, 36520, 41474, 42158, 35215, 42321, 40463, 38254, 41946, 40128, 40039, 40816, 35193, 41190, 40600, 40850, 40988, 39574, 34574, 35455, 36670, 35877, 39030, 40488, 33700, 40106, 37496, 37444, 37426, 34747, 37101, 40674, 35132, 38297, 37758, 41972, 36705, 37135, 38593, 36018, 39293, 34664, 37342, 41990, 41412, 39180, 34903, 40977, 34286, 36937, 37325, 40079, 40269, 34673, 33834, 41329, 35825, 37624, 41490, 37614, 42178, 33677, 36258, 36395, 35841, 34269, 37281, 34671, 39485, 35312, 34693, 37944, 36934, 40523, 38097, 39204, 42120, 41049, 33592, 34605, 38808, 40746, 37217, 34309, 36336, 33772, 37477, 37649, 35219, 39933, 35387, 42283, 41852, 40970, 33988, 33758, 40845, 33689, 40199, 36333, 41191, 34732, 42423, 33737, 40843, 38058, 38726, 34346, 38513, 39488, 40697, 33600, 37144, 38973, 38386, 38718, 40777, 39544, 40444, 34787, 39026, 34423, 37580, 40068, 33920, 36134, 38708, 40036, 38719, 39155, 33905, 37264, 37151, 40453, 41337, 37789, 33562, 39932, 39922, 35064, 38248, 35647, 37979, 40020, 34817, 35795, 36598, 33985, 35998, 42276, 39691, 41416, 36638, 36443, 41706, 41254, 41863, 36222, 36365, 34493, 40836, 33862, 35899, 40611, 37020, 34949, 34376, 34504, 40345, 37422, 35177, 39764, 37443, 36304, 34710, 39791, 34897, 36422, 38384, 38540, 40865, 41027, 36096, 36713, 39218, 39964, 34763, 38743, 39782, 41319, 40184, 42079, 34306, 41217, 41135, 34665, 39008, 34956, 39088, 37048, 41379, 39607, 35124, 41322, 40979, 37712, 36709, 34047, 42106, 36603, 37607, 38261, 35941, 35797, 35232, 34896, 39099, 39425, 39383, 36036, 35978, 39064, 37567, 39386, 40720, 38734, 36073, 42100, 33765, 41042, 37268, 38324, 42065, 37046, 40622, 35461, 38062, 38757, 36794, 38525, 34648, 37965, 36058, 39154, 38321, 39955, 36947, 34892, 36176, 41843, 41270, 35459, 40379, 36575, 41413, 40712, 41249, 34417, 34447, 39220, 37459, 34960, 38394, 34801, 39512]
        cfg.TRN_gids=[31364, 30106, 29095, 30404, 32411, 32689, 31136, 32952, 29223, 29226, 31158, 30716, 30874, 29694, 32244, 28976, 31159, 29987, 29060, 28963, 31917, 33413, 29656, 28828, 30883, 33262, 28657, 31312, 30553, 32956, 33496, 31127, 30686, 32101, 28629, 29205, 29588, 32390, 31050, 31497, 32746, 31493, 32515, 30601, 30161, 32655, 29546, 28867, 31388, 30041, 30876, 33378, 33405, 29303, 32878, 29152, 28751, 30860, 30265, 32556, 33156, 32386, 30801, 29146, 29527, 32576, 29159, 32174, 29217, 29736, 31557, 32851, 29401, 31530, 31839, 30846, 30879, 31170, 32483, 31756, 29935, 29172, 30493, 29094, 33437, 30013, 30440, 28812, 30307, 29266, 31156, 33485, 30093, 30987, 32245, 31512, 29027, 31093, 32087, 29125, 32921]
    
    if cfg.prepare_HPC_mode: 
        NetPyNE_BBP.Prompt.headerMsg('UPDATING CELL MORPHOLOGY FILES')
        cfg.loadCellModel = False
        cfg.convertCellModel= True
    elif cfg.run_HPC_mode:   cfg.loadCellModel = True
    else:
        cfg.loadCellModel           = True
        cfg.convertCellMorphologies = False
cfg.select_thal_gids=cfg.VPM_gids+cfg.TRN_gids
cfg.gids_lists = [cfg.VPM_gids,cfg.TRN_gids]

cfg.printCellProperties=True

#------------------------------------------------------------------------------
#   Cells
#------------------------------------------------------------------------------
cfg.re_rescale  = 1
cfg.cao_secs    = 1.2
cfg.rescaleUSE  = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file

cfg.modifyNMDAratio_L6A_VPM = 2.0
cfg.modifyNMDAratio_L6A_TRN = 1.0

cfg.addUSE_MLe      = 0.0 # (default: 0.0) adds a fixed value to all u_syn values in each conn
cfg.rescaleFac_MLe  = 1.0 # (default: 1.0) multiplies the Facilitation constant of all synapses in each conn
cfg.rescaleDep_MLe  = 1.0 # (default: 1.0) multiplies the Depression   constant of all synapses in each conn

cfg.max_gaps_perCell = None

#------------------------------------------------------------------------------
#   Network
#------------------------------------------------------------------------------
# - General pop configs
cfg.pop_shape           = 'cylinder'    # (Options: 'cylinder' or 'cube')
cfg.singleCellPops      = False         # (Options: bool or list of pops)

cfg.scale_density       = 1.0           # (Default: 1.0)

# --- Adds pops to test network step-wise
cfg.includeBS           = True          # (adds Brainstem pop  (MLe))
cfg.includeTH           = True          # (adds Thalamus  pops (VPM and TRN))
cfg.includeL6           = True          # (adds Cortical  pops (CT, CC and IN))

# --- Adds pops to test network step-wise
# Either all False or a single True
cfg.simulateThalonly    = False         # (removes all pops that are not from TH    to tune connectivity)
cfg.simulateThalMle     = False         # (removes all pops that are not from TH/BS to tune connectivity)
cfg.simulateThalMle_vCT = True          # (removes all pops that are not from TH/BS and Virtual CT to tune connectivity)
cfg.simulateL6only      = False         # (removes all pops that are not from L6    to tune connectivity)
cfg.simulateMleonly     = False         # (removes all pops that are not from TH/BS to tune connectivity)
cfg.simulateMleVPM      = False         # (removes all pops that are not from TH/BS to tune connectivity)
cfg.simulateTRNgap      = False         # (removes all pops that are not from TH/BS to tune connectivity)

# - Configuring TH
cfg.TRN_adacent_pops    = ['TRN_ring']
cfg.addTRNadjacentpops  = 'biophysical'
cfg.numTRNadjacentpops  = 4 # 4x the number of cells in the central pop

cfg.align_thal_y=True

cfg.align_virtual_cells=True
# cfg.align_virtual_cells=False

# - Configuring L6
cfg.addL6Adetailed      = False

cfg.L6ACTpops           = ['L6A_activated','L6A_suppressed','L6A_sparse','L6A_silent']
cfg.addL6Asubpops       = True
cfg.L6Asubpops          = ['L6CC_TPC_L1_cAD','L6CC_UTPC_cAD','L6CC_BPC_cAD','L6CC_IPC_cAD','L6IN_LBC_bAC','L6IN_LBC_bNA','L6IN_LBC_cNA','L6IN_MC_bAC','L6IN_MC_bNA','L6IN_MC_cAC']
cfg.simplifyL6A         = False

# - Configuring BS
cfg.align_cells         = False         # (Options: (True: aligns X-positions in crescent order based on angle/depth) / (False: random X-positions))

# --- Connectivity
# - General conn configs
cfg.removeConns         = False # (Options: True / False / 'elec' / 'chem')
cfg.removeESyns         = False
cfg.addL6Ainterconns    = True  # (MERGE with Fernando code)

##########################################################################################################################################################################
# --- Stimulation
##########################################################################################################################################################################
# Values obtained running the script in <calibrate_threshold_currents.py>
# Values are stored in <../cells/calibrated_currents/calibrated_currents.json>
# These are the currents that force the cells into a non-spiking regime (~ -80mV)
# Don't change, otherwise the currents will need to be recalibrated - change cfg.addRiCurrentTargetVoltage instead
# cfg.addHoldingCurrent           = False
cfg.addHoldingCurrent           = True
cfg.addHoldingCurrentPops       = ['VPM__pop','TRN__pop','TRN_ring__pop']
cfg.addHoldingCurrentPercent    = { 'VPM__pop':     100,         
                                    'TRN__pop':     100,       
                                    'TRN_ring__pop':100} # --- (base value is already negative) in percent (e.g 0%, 100%, 150%, ...)

##########################################################################################################################################################################
cfg.addThresholdCurrent         = True
cfg.addThresholdCurrentPops     = ['VPM__pop','TRN__pop','TRN_ring__pop']
# cfg.addThresholdCurrentPercent  = { 'VPM__pop':     400,          
#                                     'TRN__pop':     250,       
#                                     'TRN_ring__pop':250} # --- in percent (e.g 0%, 100%, 150%, ...)
# # cfg.addCellTypeCompensation     = True # Multiplier added to threshold current for each individual cell type

# cfg.map_ThresholdCurrent = {
#                                 'VPL_TC__36043':    {'rmp':[-81.22027408108586, -78.00540339412105, -74.79119872686145, -71.78805400254281, -55.60552391344034  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 560.0]},
#                                 'VPL_TC__34349':    {'rmp':[-81.18859298106919, -77.94358425440531, -74.68609848153007, -71.6622466872409,  -55.36875966501423  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 560.0]},
#                                 'VPL_TC__36219':    {'rmp':[-81.14892418724747, -78.13753987922506, -75.10106409159293, -72.24187020896431, -57.28931991373245  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 600.0]},
#                                 'VPL_TC__37520':    {'rmp':[-79.95278649296449, -75.81918972548455, -71.69588877070348, -67.93279085530456, -55.02821154374243  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 360.0]},
#                                 'VPL_TC__41912':    {'rmp':[-81.06140268393905, -78.05523668723227, -75.03976973952162, -72.20960110213443, -55.384782023891454 ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 600.0]},
#                                 'Rt_RC__29892':     {'rmp':[-79.66866936639987, -74.17962590606305, -69.08657546920023, -64.57382133341058, -60.05802883480706  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 250.0]},
#                                 'Rt_RC__31315':     {'rmp':[-79.88288702059448, -74.95201887540686, -70.131654545683,   -65.4893384346661,  -60.80437117022895  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 250.0]},
#                                 'Rt_RC__30252':     {'rmp':[-80.37284490706574, -75.94629294993345, -71.58608705491505, -67.34981340174933, -63.121416632696736 ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 250.0]},
#                                 'Rt_RC__32049':     {'rmp':[-79.88157898482328, -74.95208140903134, -70.1322659268712,  -65.48894140461844, -60.80111160896434  ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 250.0]},
#                                 'Rt_RC__28925':     {'rmp':[-80.38704639959978, -75.96684470322047, -71.60901794240488, -67.36898967780695, -63.137300292653684 ],    'i_thresh':[50.0, 100.0, 150.0, 200.0, 250.0]},
#                                 }

cfg.map_ThresholdCurrent = map_ThresholdCurrent

cfg.map_targetRMP_byCellTemplate = {}
for cellTemplate in cfg.map_ThresholdCurrent.keys():
    if   'VPL_TC' in cellTemplate: cfg.map_targetRMP_byCellTemplate.update({cellTemplate:cfg.map_targetRMP['VPL_TC']})
    elif 'Rt_RC'  in cellTemplate: cfg.map_targetRMP_byCellTemplate.update({cellTemplate:cfg.map_targetRMP['Rt_RC']})

cfg.addThresholdCurrentPops_byCellTemplate={cellTemplate:NetPyNE_BBP.Utils.cubic_extrapolation(cfg.map_ThresholdCurrent[cellTemplate]['rmp'], cfg.map_ThresholdCurrent[cellTemplate]['i_thresh'], cfg.map_targetRMP_byCellTemplate[cellTemplate]) for cellTemplate in cfg.map_ThresholdCurrent.keys()}
'''
# Code to confirm the RMPs after sim run (calibrated for t = 1000 ms)
for cell_ind,cell in enumerate(sim.net.cells): print(sim.net.cells[cell_ind].gid, '\t',sim.net.cells[cell_ind].tags['label'],'\t', sim.simData['V_soma']['cell_'+str(sim.net.cells[cell_ind].gid)][-1])
'''

# Adds the parameters for the vector.play() method to change the RMP during simulation runtime
# p.s.: to avoid changing the RMP, cfg.shift_t can be defined as < cfg.shift_t=[] >, so that the code in CurrentStim.py will ignore these modifications
cfg.addModifyRMP_vecPlay=True
if cfg.addModifyRMP_vecPlay:
    
    # # --- Avoids modifying the stims
    # cfg.shift_t                 = []

    # --- Modifies the stims
    cfg.shift_t                 = [1,       4000,]
    cfg.target_VPL_TC_shiftRMP  = [-62,     -60, ]
    cfg.target_Rt_RC_shiftRMP   = [-67,     -61, ]

    cfg.modifyRMP_vecPlay_dict={'t':cfg.shift_t}    
    cfg.modifyRMP_vecPlay_dict.update({key:[] for key in cfg.map_ThresholdCurrent.keys()})
    for i in range(len(cfg.shift_t)):
        # store_new_targetI={}
        for cellTemplate in cfg.map_ThresholdCurrent.keys():
            if   'VPL_TC' in cellTemplate: new_targetRMP = cfg.target_VPL_TC_shiftRMP[i]
            elif 'Rt_RC'  in cellTemplate: new_targetRMP = cfg.target_Rt_RC_shiftRMP[i]
            cfg.modifyRMP_vecPlay_dict[cellTemplate].append(NetPyNE_BBP.Utils.cubic_extrapolation(cfg.map_ThresholdCurrent[cellTemplate]['rmp'], cfg.map_ThresholdCurrent[cellTemplate]['i_thresh'], new_targetRMP))

#
cfg.addModifyConnWeight_vecPlay=True
cfg.modifyConnWeight_vecPlay_dict ={
                                    'conn|TRN|VPM|chem':                        [(1,1.75),(4000,0.15),(7000,1.75)],
                                    'conn|TRN@biophysical|TRN@biophysical|chem':[(1,1.75),(4000,0.15),(7000,1.75)], 
                                    }

#------------------------------------------------------------------------------
#   MLe inputs
#------------------------------------------------------------------------------
from GenerateStimDataset import SampleData

cfg.deflection_events = [[2000,2150],[3000,3150],[4000,4250],[5000,5250],[6000,6350],[7000,7350]]

cfg.samplingInterv = 25     # us
# cfg.samplingInterv = 1000   # us

if cfg.target_angle is not None:    
    cfg.deflection_times = cfg.deflection_events+[[10000,10020]]
    cfg.stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in cfg.deflection_times)
else:                               
    cfg.deflection_times = [[10000,10020]]
    cfg.stims_string = '_stims|None'

cfg.stim_source_code_folder    = '/stims/Minnery_2003/fig1'
cfg.deflection_model_folder    = cfg.NetPyNE_rootFolder+cfg.stim_source_code_folder+'/deflection_model/'   # Folder to save the output deflection model and sampled raster
# --- Load deflection model and sampled raster (stored in JSON for readability - small files)
cfg.deflection_dict_path = SampleData.LoadJSON(file_path=   cfg.deflection_model_folder+'deflectionModel_'+str(cfg.target_angle)+'|deg'+'_samplingBins|'+str(cfg.samplingInterv)+'us'+cfg.stims_string+'.json')
cfg.deflection_dataset_path =                cfg.deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|'+str(cfg.target_angle)+'_samplingBins|'+str(cfg.samplingInterv)+'us'+cfg.stims_string+'_simplifiedDataset'+'.pkl'

### 10 Hz ML input that simulates driver input during sleep
# cfg.deflection_dataset_path =                '../stims/Artificial/CT/spiking_dict__nCells_1000__sampling_25_uS__freq_10_Hz__downTimes__None__.json'

cfg.extend_spike_times=True

#------------------------------------------------------------------------------
#   VPM ring inputs
#------------------------------------------------------------------------------
cfg.includeBaselineDriver=True
cfg.baseline_driver_path =                cfg.deflection_model_folder+'spike_dicts/baselineVPMSpikes_deflectionAngle|None_samplingBins|1000us_stims|None_simplifiedDataset.pkl'
#------------------------------------------------------------------------------
#   Stimulation
#------------------------------------------------------------------------------

'''
Currents - from source paper:
One step (200% of firing threshold) on top of a hyperpolarizing current was therefore used to constrain the bursting response, 
while three depolarizing steps on top of a depolarizing holding current (to reach -64 mV) were used to constrain the tonic firing responses, as explained above.
'''

# Tunned currents - not perfec because of membrane sag (specially in TC), but close enough (+-2mV)
cfg.addNoiseIClamp=True 
cfg.NoiseIClampParams={
                        # custom
                        'VPM__pop':             {'amp':0},
                        'TRN__pop':             {'amp':0},
                        'TRN_ring__pop':        {'amp':0},                    
                    }

#------------------------------------------------------------------------------
#   Connectivity
#------------------------------------------------------------------------------

# # --- Connectivity to update Gap conns
# cfg.conn_MLe_VPM    = True     # works
# cfg.conn_TRN_VPM    = True     # works
# cfg.conn_TRN_TRN    = True     # works
# cfg.conn_VPM_TRN    = True     # works

# cfg.conn_L6A_VPM    = True     # works
# cfg.conn_L6A_TRN    = True     # works
# cfg.conn_VPM_L6A    = True     # works

# cfg.conn_TRNe_TRNe_connList = False
# cfg.addCTvirtual    = True

# # --- Removing Connectivity
# cfg.conn_MLe_VPM    = False     # works
# cfg.conn_TRN_VPM    = False     # works
# cfg.conn_TRN_TRN    = False     # works
# cfg.conn_VPM_TRN    = False     # works

# cfg.conn_L6A_VPM    = False     # works
# cfg.conn_L6A_TRN    = False     # works
# cfg.conn_VPM_L6A    = False     # works

# cfg.conn_TRNe_TRNe_connList = False
# cfg.addCTvirtual    = True

# --- Enabling Connectivity
cfg.conn_MLe_VPM    = True     # works
cfg.conn_TRN_VPM    = True     # works
cfg.conn_TRN_TRN    = True     # works
cfg.conn_VPM_TRN    = True     # works

cfg.conn_L6A_VPM    = True     # works
cfg.conn_L6A_TRN    = True     # works
cfg.conn_VPM_L6A    = True     # works

cfg.conn_TRNe_TRNe_connList = True
cfg.addCTvirtual    = True

#------------------------------------------------------------------------------

cfg.delayCTvirtual  = 0
cfg.delayBS         = 0

cfg.reset_cell_rotation=False
extra_rotation = -90    # extra rotation to compensate for the plotting code being rotated 90 deg on the z-axis (depth axis - into the screen)
flip_trn       = 180    # extra 180 deg rotation so that concave part of TRN cells faces the VPM, which is more accurate with thalamic anatomy
cfg.rotation_angle = {  
                        'Rt_RC__28925':     np.deg2rad(-65  + extra_rotation + flip_trn),
                        'Rt_RC__30252':     np.deg2rad(-65  + extra_rotation + flip_trn),
                        'Rt_RC__31315':     np.deg2rad(-85  + extra_rotation + flip_trn),
                        'Rt_RC__31924':     np.deg2rad( 45  + extra_rotation + flip_trn),
                        'Rt_RC__32049':     np.deg2rad(-85  + extra_rotation + flip_trn),
                        
                        'VPL_TC__34349':    np.deg2rad(40   + extra_rotation),
                        # 'VPL_TC__36043':    np.deg2rad(0),
                        'VPL_TC__36219':    np.deg2rad(-15  + extra_rotation),
                        'VPL_TC__36963':    np.deg2rad(45   + extra_rotation),
                        'VPL_TC__37520':    np.deg2rad(40   + extra_rotation),
                        'VPL_TC__41912':    np.deg2rad(0    + extra_rotation),
                        }

# #------------------------------------------------------------------------------

# --- Uniform connectivity

# # Uniform FF / Uniform FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_uniform'; 

# # (new) Uniform FF / Open FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_open'; 

# # (new) Uniform FF / Closed FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_closed'; 

# # (new) Uniform FF / Mixed FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_mixed'; 

# # (new) Uniform FF / Remixed FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_remixed'; 

# --- Topological connectivity

# # (new) Closed FF / Uniform FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_uniform'; 

# # Closed FF / Open FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_open'; 
# # FF_flag = '__VPM_ff_topolog'; FB_flag = '__TRN_fb_topolog'; 

# # Closed FF / Closed FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_closed'; 

# # Closed FF / Mixed FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_mixed'; 

# # (new) Closed FF / Remixed FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_remixed'; 

# # --- Topological connectivity - Closed + Open

# # Mixed FF / Uniform FB
# FF_flag = '__VPM_ff_mixed'; FB_flag = '__TRN_fb_uniform'; 

# # Mixed FF / Open FB
# FF_flag = '__VPM_ff_mixed'; FB_flag = '__TRN_fb_open'; 

# # Mixed FF / Closed FB
# FF_flag = '__VPM_ff_mixed'; FB_flag = '__TRN_fb_closed'; 

# # Mixed FF / Mixed FB
# FF_flag = '__VPM_ff_mixed'; FB_flag = '__TRN_fb_mixed'; 

# # (new) Mixed FF / Remixed FB
# FF_flag = '__VPM_ff_mixed'; FB_flag = '__TRN_fb_remixed'; 

# # --- Topological connectivity - Closed + Uniform

# # (new) Remixed FF / Uniform FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_uniform'; 

# # (new) Remixed FF / Open FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_open'; 

# # (new) Remixed FF / Closed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_closed'; 

# # (new) Remixed FF / Mixed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_mixed'; 

# # (new) Remixed FF / Remixed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_remixed'; 

# #------------------------------------------------------------------------------

# # Uniform FF / Uniform FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_uniform'; ct_FB_flag = '__L6A_fb_uniform'     # uniform

# # Closed FF / Closed FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_closed'; ct_FB_flag = '__L6A_fb_closed'      # closed

# # (new) Remixed FF / Remixed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_remixed'; ct_FB_flag = '__L6A_fb_remixed'     # remixed

# # (new_new) Remixed FF / Remixed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_openRemixed'; ct_FB_flag = '__L6A_fb_remixed'     # remixed

# # (new_new) Closed FF / Open FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_open'; ct_FB_flag = '__L6A_fb_closed'     # open

# # ------- test add uniform CT ------------------------------------------------------------------------------
add_tr=False

# # Uniform FF / Uniform FB
# FF_flag = '__VPM_ff_uniform'; FB_flag = '__TRN_fb_uniform'; ct_FB_flag = '__L6A_fb_uniform'     # uniform

# # Closed FF / Closed FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_closed'; ct_FB_flag = '__L6A_fb_uniform'      # closed

# (new) Remixed FF / Remixed FB
FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_remixed'; ct_FB_flag = '__L6A_fb_uniform'; add_tr=True     # remixed

# # (new_new) Remixed FF / Remixed FB
# FF_flag = '__VPM_ff_remixed'; FB_flag = '__TRN_fb_openRemixed'; ct_FB_flag = '__L6A_fb_uniform'; add_tr=True     # remixed

# # (new_new) Closed FF / Open FB
# FF_flag = '__VPM_ff_closed'; FB_flag = '__TRN_fb_open'; ct_FB_flag = '__L6A_fb_uniform'     # open

# #------------------------------------------------------------------------------

# # Percentage of connections that are topological/uniform
# # cfg.topology_ratio = [0.9,0.1]
# # cfg.topology_ratio = [0.8,0.2]
# # cfg.topology_ratio = [0.7,0.3]
# # cfg.topology_ratio = [0.75,0.25]
# # cfg.topology_ratio = [0.6,0.4]
# cfg.topology_ratio = [0.5,0.5]

# # cfg.topology_ratio = [0.4,0.6]
# # cfg.topology_ratio = [0.3,0.7]
# # cfg.topology_ratio = [0.25,0.75]
# # cfg.topology_ratio = [0.2,0.8]
# # cfg.topology_ratio = [0.1,0.9]

# cfg.topology_ratio = [0.75,0.25]
cfg.topology_ratio = [0.5,0.5]
# cfg.topology_ratio = [0.25,0.75]

if (ct_FB_flag == '__L6A_fb_remixed') or (ct_FB_flag == '__L6A_fb_mixed') or add_tr: ct_FB_flag+='_tr_'+str(int(cfg.topology_ratio[0]*100))+'_'+str(int(cfg.topology_ratio[1]*100))
# if (ct_FB_flag == '__L6A_fb_remixed') or (ct_FB_flag == '__L6A_fb_mixed'): ct_FB_flag+='_tr_'+str(int(cfg.topology_ratio[0]*100))+'_'+str(int(cfg.topology_ratio[1]*100))

# _addVersion=''
_addVersion='_v04'

cfg.dump_cell_properties = cfg.dump_cell_properties.split('.json')[0]+FF_flag+FB_flag+ct_FB_flag+_addVersion+'.json'

#------------------------------------------------------------------------------

def get_synsPerConn(conn_data, pre_pop, post_pop, conn_type):
    try:    
        syns_per_conn = round(conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1
    return syns_per_conn

# --- Connectivity properties
# --- Loads the conn data calculated from the projections in the BBP somatosensory thalamus model
cfg.conn_data = BN.NetworkConversion.convertConnProperties(filePath=cfg.BBP_conn_properties)

#------------------------------------------------------------------------------

if cfg.conn_TRN_TRN == False: 
    NetPyNE_BBP.Prompt.headerMsg('TRN -> TRN CONNECTIONS REMOVED FOR TESTING - BASED ON (Hou, 2016)')
else:
    NetPyNE_BBP.Prompt.headerMsg('TRN -> TRN GABA CONNECTIONS ENABLED')

#------------------------------------------------------------------------------
conversion_factor = 0.001 # must be added to all synapses, except gap junctions

scale_weights = 1     # 
scale_EE = conversion_factor * scale_weights * 1.0
scale_IE = conversion_factor * scale_weights * 1.0 # previous: 1.25 good
scale_EI = conversion_factor * scale_weights * 1.0 # previous: 1 - sim 125 (5 too much sim 127)
scale_II = conversion_factor * scale_weights * 1.0 # previous: 1 - sim 124

weights_dict={
                'MLe_VPM':      conversion_factor*9,
 
                'TRNe_TRNe':    10e-6,
                
                'TRN_VPM':      conversion_factor*3,
                'VPM_TRN':      conversion_factor*4, # bad - reduce and continue

                'TRN_TRN':      conversion_factor*4, # more TRN synchronization
                
                'L6A_VPM':      conversion_factor*4, # temporary, for plotting conns
                'L6A_TRN':      conversion_factor*4, # temporary, for plotting conns

                'VPM_L6A':      conversion_factor*1, # temporary, for plotting conns
                }

cfg.BaselineDriver_weights_VPM_TRN_ring = conversion_factor*2/6

#------------------------------------------------------------------------------
# --- MLe => VPM
''' (Takeuchi, 2017) - Takeuchi, Y., Osaki, H., Yagasaki, Y., Katayama, Y., & Miyata, M. (2017). Afferent Fiber Remodeling in the Somatosensory Thalamus of Mice as a Neural Basis of Somatotopic Reorganization in the Brain and Ectopic Mechanical Hypersensitivity after Peripheral Sensory Nerve Injury. eNeuro, 4(2). https://doi.org/10.1523/ENEURO.0345-16.2017
    number of boutons   / mle fiber: 40.5 +- 25.0 
    number of terminals / mle fiber: 25.4 +- 15.7 
    number of branches  / mle fiber: 22.3 +- 14.1 
'''
# # # cfg.MLe_VPM__syns_per_conn  = 40 # target: 40.5 +- 25.0 boutons per MLe cell (presynaptic) into VPM cells (1 to 3 postsynaptic cells)
# # cfg.MLe_VPM__syns_per_conn  = 55 # target: 40.5 +- 25.0 boutons per MLe cell (presynaptic) into VPM cells (1 to 3 postsynaptic cells)
# cfg.MLe_VPM__syns_per_conn  = 40 # target: 40.5 +- 25.0 boutons per MLe cell (presynaptic) into VPM cells (1 to 3 postsynaptic cells)

# test - 2024_08_26 - convergence = 143.35260115606937
cfg.MLe_VPM__syns_per_conn  = 60 # target: 40.5 +- 25.0 boutons per MLe cell (presynaptic) into VPM cells (1 to 3 postsynaptic cells)
cfg.MLe_VPM__y_tresh        = 0.0065

# # test - 2024_07_04
# cfg.MLe_VPM__syns_per_conn  = 20 # target: 40.5 +- 25.0 boutons per MLe cell (presynaptic) into VPM cells (1 to 3 postsynaptic cells)
# cfg.MLe_VPM__y_tresh        = 0.02

cfg.MLe_VPM__pre_pop        = 'MLe'
cfg.MLe_VPM__post_pop       = 'VPM'
cfg.MLe_VPM__pre_pop_axis   = 'y'
cfg.MLe_VPM__post_pop_axis  = 'y'

cfg.MLe_VPM__syn_mech       = 'syn|MLe|VPM|exc'
cfg.MLe_VPM__conn_type      = 'chem'
cfg.MLe_VPM__weight         = weights_dict['MLe_VPM']
cfg.MLe_VPM__target_secs    = 'inputs__proximal'

#------------------------------------------------------------------------------
# --- Updated conn - including L6A subtypes
# --- L6A => VPM - distanceBasedProbability_1D_exponential

# --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.L6A_VPM__conn_params_dict={
                                # 'conn_prob':2.0827527202658684, 'decay_factor':10, 'baseline_prob': 0.05,
                                'conn_prob':2.0827527202658684, 'decay_factor':8, 'baseline_prob': 0.01,
                                }

cfg.L6A_VPM__pre_pop        = 'all_L6A_subpops'
cfg.L6A_VPM__post_pop       = 'VPM'
cfg.L6A_VPM__pre_pop_axis   = 'theta'
cfg.L6A_VPM__post_pop_axis  = 'y'
cfg.L6A_VPM__conn_type      = 'chem'
cfg.L6A_VPM__syn_mech       = 'syn|L6A|VPM|exc'
cfg.L6A_VPM__target_secs    = 'inputs__distal'
cfg.L6A_VPM__weight         = weights_dict['L6A_VPM']
cfg.L6A_VPM__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop='L6A',post_pop=cfg.L6A_VPM__post_pop,conn_type=cfg.L6A_VPM__conn_type)

#------------------------------------------------------------------------------
# --- Updated conn - including L6A subtypes
# --- L6A => TRN - distanceBasedProbability_1D_exponential

# --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.L6A_TRN__conn_params_dict={
                                # 'conn_prob':0.8010587385637955, 'decay_factor':10, 'baseline_prob': 0.05,
                                'conn_prob':0.8010587385637955, 'decay_factor':8, 'baseline_prob': 0.01,
                                }

cfg.L6A_TRN__pre_pop        = 'all_L6A_subpops'
cfg.L6A_TRN__post_pop       = 'TRN'
cfg.L6A_TRN__pre_pop_axis   = 'theta'
cfg.L6A_TRN__post_pop_axis  = 'y'   # (2023_11_04) - testing adjustment of projections
cfg.L6A_TRN__conn_type      = 'chem'
cfg.L6A_TRN__syn_mech       = 'syn|L6A|TRN|exc'
cfg.L6A_TRN__target_secs    = 'inputs__distal'
cfg.L6A_TRN__weight         = weights_dict['L6A_TRN']
cfg.L6A_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop='L6A',post_pop=cfg.L6A_TRN__post_pop,conn_type=cfg.L6A_TRN__conn_type)

#------------------------------------------------------------------------------
# --- TRN => VPM

cfg.TRN_VPM__ref_pop        = 'MLe'
cfg.TRN_VPM__pre_pop        = 'TRN'
cfg.TRN_VPM__post_pop       = 'VPM'
cfg.TRN_VPM__conn_type      = 'chem'
cfg.TRN_VPM__syn_mech       = 'syn|TRN|VPM|inh'
cfg.TRN_VPM__target_secs    = 'inputs__intermediate'
cfg.TRN_VPM__weight         = weights_dict['TRN_VPM']
cfg.TRN_VPM__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRN_VPM__pre_pop,post_pop=cfg.TRN_VPM__post_pop,conn_type=cfg.TRN_VPM__conn_type)

cfg.TRN_VPM__conn_method = 'probability'
cfg.TRN_VPM__pre_pop_axis   = 'y'
cfg.TRN_VPM__post_pop_axis  = 'y'

# (Off-center) --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.TRN_VPM__conn_params_dict={
                                'conn_prob':1.0, 'decay_factor':15, 'baseline_prob': 0.55,
                                }
cfg.TRNring_VPM__conn_params_dict={
                                'conn_prob':1.0, 'decay_factor':15, 'baseline_prob': 0.55,
                                }

# (On-center) --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.TRN_VPM__conn_params_dict_closed={
                                        'conn_prob':1.0, 'decay_factor':3.75, 'baseline_prob': 0.01,
                                        }
cfg.TRNring_VPM__conn_params_dict_closed={
                                        'conn_prob':1.0, 'decay_factor':3.75, 'baseline_prob': 0.01,
                                        }

#------------------------------------------------------------------------------
# --- TRN => TRN

# --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.TRN_TRN__conn_params_dict={
                                # 'conn_prob':0.4, 'decay_factor':0.022, 'baseline_prob': 0.01,    # bWgt_4059
                                'conn_prob':1.0, 'decay_factor':0.022, 'baseline_prob': 0.01,
                                }

cfg.TRN_TRN__pre_pop        = 'TRN'
cfg.TRN_TRN__post_pop       = 'TRN'
cfg.TRN_TRN__conn_type      = 'chem'
cfg.TRN_TRN__syn_mech       = 'syn|TRN|TRN|inh'
cfg.TRN_TRN__target_secs    = 'inputs__distal'
cfg.TRN_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRN_TRN__pre_pop,post_pop=cfg.TRN_TRN__post_pop,conn_type=cfg.TRN_TRN__conn_type)
cfg.TRN_TRN__conn_method    = 'probability'
cfg.TRN_TRN__weight         = weights_dict['TRN_TRN']
cfg.TRN_TRN__conn_rule      = cfg.conn_data[cfg.TRN_TRN__post_pop][cfg.TRN_TRN__conn_type][cfg.TRN_TRN__pre_pop]['convergence'][0]
print('\t>>\tTRN TRN chemical: ', cfg.TRN_TRN__conn_rule)

#------------------------------------------------------------------------------
# test
if cfg.conn_TRNe_TRNe_connList == False: print('\n\n\n\n\n\n\n ---- WARNING: Gap conns are inactivated!!!!! --- \n\n\n\n\n\n\n')
cfg.individual_gap_conns = True;cfg.simplify_gap_conns=False

cfg.TRNe_TRNe__weight       = weights_dict['TRNe_TRNe']

# # Deprecated
# # --- TRNe => TRNe
# from CalculateProbability import CalculateProbability as CalProb
# '''
# # --- Lee 2014
# gap_number = [2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 12, 12, 16, 17, 18, 21, 21]
# import numpy as np
# np.mean(gap_number) # = 8.625+-5.072905971925756
# '''
# cfg.TRNe_TRNe__pre_pop      = 'TRN'
# cfg.TRNe_TRNe__post_pop     = 'TRN'
# cfg.TRNe_TRNe__conn_type    = 'elec'
# cfg.TRNe_TRNe__syn_mech     = 'gap_nr' # Gap junction mechanism from https://github.com/BlueBrain/neuron_reduce/blob/08bea3520c0f535cdba27ef0c3e4d8f970d08604/tests/TestsFiles/Test_4_LBC_amsalem/mod/halfgap.mod#L4 
# # cfg.TRNe_TRNe__syn_mech     = 'esyn' # Gap junction mechanism from NetPyNE - still breaking with MPI
# cfg.TRNe_TRNe__target_secs  = 'inputs__distal'
# # cfg.TRNe_TRNe__weight       = 0.001
# cfg.TRNe_TRNe__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRNe_TRNe__pre_pop,post_pop=cfg.TRNe_TRNe__post_pop,conn_type=cfg.TRNe_TRNe__conn_type)

# coupled_neurons             = 8.625 # --- 2 to 20 coupled neurons (Iavarone, 2023 - Fig 3)
# # coupling_scaling_factor     = 2.364365671641791 # --- Scaling fator to reach 8.625 coupled neurons - Conn result: 8.422867513611616 -Obtained by running the sim object through CalProb.verifyConns(sim,synsPerConn) to inspect the actual number of conns created
# coupling_scaling_factor     = 2.364365671641791*2.4827586206896552 # --- Scaling fator to reach 8.625 coupled neurons - Conn result: 8.422867513611616 -Obtained by running the sim object through CalProb.verifyConns(sim,synsPerConn) to inspect the actual number of conns created
# max_dist                    = 80
# cell_density_dict           = BN.BuildNetwork.getCellDensity()
# diam,height                 = BN.BuildNetwork.getNetworkSize(cfg.center_point)
# cfg.TRNe_TRNe__conn_method  = 'probability'
# conn_prob = CalProb.calculate_probability_of_connection(  diameter        = diam['TRN_ring'][1]   - diam['TRN_ring'][0],
#                                                                             height          = height[cfg.TRNe_TRNe__pre_pop][1] - height[cfg.TRNe_TRNe__pre_pop][0],
#                                                                             avg_connections = coupled_neurons,
#                                                                             max_distance    = max_dist,
#                                                                             numCells        = 551,
#                                                                             seed            = cfg.base_random_seed,
#                                                                             )
# scaled_conn_prob = conn_prob * coupling_scaling_factor
# cfg.TRNe_TRNe__conn_rule = '%f*(dist_3D<%f)*(dist_y<(%f/5))'% (scaled_conn_prob,max_dist,max_dist)
# # cfg.TRNe_TRNe__conn_rule = '%f*(dist_3D<%f)'% (scaled_conn_prob,max_dist)

# # cfg.TRNe_TRNe__conn_method  = 'convergence'
# # cfg.TRNe_TRNe__conn_rule    = 8.625 # --- 2 to 20 coupled neurons (Iavarone, 2023 - Fig 3)

# # rescale_gap_convergenge     = 0.1
# # cfg.TRNe_TRNe__conn_rule    = rescale_gap_convergenge * cfg.conn_data[cfg.TRNe_TRNe__post_pop][cfg.TRNe_TRNe__conn_type][cfg.TRNe_TRNe__pre_pop]['convergence'][0]
# print('\t>>\tTRN TRN electrical: ', cfg.TRNe_TRNe__conn_rule)
# print('\t\t----------- >   Using GAP junction ', cfg.TRNe_TRNe__syn_mech)

#------------------------------------------------------------------------------
# --- VPM => TRN - distanceBasedProbability_1D_exponential

# --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.VPM_TRN__conn_params_dict={
                                # 'conn_prob':0.2, 'decay_factor':10, 'baseline_prob': 0.05,
                                # 'conn_prob':1.0, 'decay_factor':20, 'baseline_prob': 0.01,    # bWgt_4059
                                'conn_prob':1.0, 'decay_factor':7.5, 'baseline_prob': 0.01,
                                }

cfg.VPM_TRN__pre_pop        = 'VPM'
cfg.VPM_TRN__post_pop       = 'TRN'
cfg.VPM_TRN__pre_pop_axis   = 'y'
cfg.VPM_TRN__post_pop_axis  = 'y'   # (2023_11_04) - testing adjustment of projections
cfg.VPM_TRN__conn_type      = 'chem'
cfg.VPM_TRN__syn_mech       = 'syn|VPM|TRN|exc'
cfg.VPM_TRN__target_secs    = 'inputs__proximal'
cfg.VPM_TRN__weight         = weights_dict['VPM_TRN']
cfg.VPM_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.VPM_TRN__pre_pop,post_pop=cfg.VPM_TRN__post_pop,conn_type=cfg.VPM_TRN__conn_type)

#------------------------------------------------------------------------------
# --- Updated conn - including L6A subtypes
# --- VPM => L6A - distanceBasedProbability_1D_exponential

# --- Creating a dictionary of (probability and decay factor) combinations that result in a similar number of conns
cfg.VPM_L6A__conn_params_dict={
                                # 'conn_prob':0.2, 'decay_factor':10, 'baseline_prob': 0.05,
                                'conn_prob':0.8, 'decay_factor':100, 'baseline_prob': 0.01,
                                }

cfg.VPM_L6A__pre_pop        = 'VPM'
cfg.VPM_L6A__post_pop       = 'L6A_activated'
cfg.VPM_L6A__pre_pop_axis   = 'y'
cfg.VPM_L6A__post_pop_axis  = 'theta'
cfg.VPM_L6A__conn_type      = 'chem'
cfg.VPM_L6A__syn_mech       = 'syn|VPM|L6A|exc'
cfg.VPM_L6A__weight         = weights_dict['VPM_L6A']
cfg.VPM_L6A__TC_syns_per_conn = 9        # from S1-netpyne model (Borges, 2022)

CT_cells            = 110+75 # (activated+suppressed)
CT_bouton_per_cell  = 131      # from (Meyer, 2010) 
TC_cells            = 318
TC_CT_conns = (CT_cells*CT_bouton_per_cell)/(TC_cells*cfg.VPM_L6A__TC_syns_per_conn)
print('\t>>\tTC cell target convergence onto CT neurons: ',round(TC_CT_conns))

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.filename            = 'barr_net_'          # Set file output name
cfg.simLabel            = cfg.sim_tag
cfg.savePickle          = True             # Save params, network and sim output to pickle file
cfg.saveJson            = False              # Save params, network and sim output to JSON file
cfg.saveDataInclude     = ['simConfig','netParams','simData','net']
cfg.saveCellConns       = True
if cfg.saveCellConns == False: NetPyNE_BBP.Prompt.headerMsg('THIS IS A DEVELOPMENT/DEBUG SESSION: CELL CONNS ARE NOT BEING SAVED - for final runs set cfg.saveCellConns = True')
folderName              = 'barreloid_network'

cfg.saveFolder          = '../data/barreloid_sims/'+folderName

#------------------------------------------------------------------------------
# --- Plotting
#------------------------------------------------------------------------------
cfg.saveFigPath = cfg.NetPyNE_rootFolder+'/figs/barreloid_figs'

popColors_dict={}
for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
            'MLe__pop', 
            'L6A_activated__pop',
            'CTvirtual_uniform__pop',
            'CTvirtual_activated__pop',
            ]:
    if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
    elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
    elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
    elif 'CTvirtual_' in    pop:    popColors_dict.update({pop:'r'})
    elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
    elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
    elif 'MLe' in           pop:    popColors_dict.update({pop:'k'})
    else:                           popColors_dict.update({pop:'gray'})

cfg.analysis['plotRaster']  = { 
                                'figSize':[45, 25], 
                                'include':[ 'MLe__pop','VPM__pop', 'TRN__pop', 
                                            'CTvirtual_uniform__pop',
                                            'CTvirtual_activated__pop',
                                           ],
                                'timeRange':[1500,cfg.duration],
                                'orderInverse': True,
                                'orderBy':'y',
                                'spikeHist':True,
                                'popColors':popColors_dict,
                                'saveFig': True,
                                'dpi':300,
                                'tightLayout':False,
                                } # Plot a raster


if cfg.simulateThalonly or cfg.simulateThalMle or cfg.simulateL6only or cfg.simulateTRNgap:
    pops = [
            'VPM__pop', 
            'TRN__pop', 
            'TRN_ring__pop', 
            'SEPop', 
            ]
    record_pops = [(pop,list(np.arange(0,40))) for pop in pops]
    overlay=True
else:
    pops = [
            'VPM__pop', 
            'TRN__pop', 
            'TRN_ring__pop', 
            'L6A_activated__pop', 'L6A_suppressed__pop', 'L6A_sparse__pop', 'L6A_silent__pop',
            'L6CC_TPC_L1_cAD__pop', 'L6CC_UTPC_cAD__pop', 'L6CC_BPC_cAD__pop', 'L6CC_IPC_cAD__pop',
            'L6IN_LBC_bAC__pop', 'L6IN_LBC_bNA__pop', 'L6IN_LBC_cNA__pop', 'L6IN_MC_bAC__pop', 'L6IN_MC_bNA__pop', 'L6IN_MC_cAC__pop',
            'L6IN__pop',
            'TRN_1__pop', 
            'TRN_5__pop', 
            ]
    record_pops = [(pop,[0]) for pop in pops]
    overlay=True

cfg.analysis['plotTraces']  = {
                                'include': ['all'], 
                                'overlay': overlay,
                                'oneFigPer': 'trace',
                                'timeRange':[1500,cfg.duration],
                                'figSize':[40,15],
                                'saveFig': True,
                                }


conn_pre_pops=[ 
                'MLe__pop', 
                'VPM__pop', 
                'TRN__pop', 
                'TRN_ring__pop', 
                'CTvirtual_uniform__pop', 
                'BaselineDriver__pop'
                ]

cfg.analysis['plotConn']    = {'figSize':[10, 10], 
                               'includePre':conn_pre_pops,
                               'includePost':['VPM__pop', 'TRN__pop', 'TRN_ring__pop'],
                               'saveFig':True,
                               } # plot connectivity matrix

# cfg.record_windowed_spectrogram = True # saves the the windowed spectrogram - heavy file - only used in barreloid_batch_fig06_sleep
cfg.record_windowed_spectrogram = False # 

cfg.analysis['plotLFP'] = {'includeAxon':   False, 
                        'figSize':       [40,50], 
                        'timeRange':     [1500,cfg.duration], 
                        'plots':         ['timeSeries', 'spectrogram'],
                        'maxFreq':       75,
                        'saveFig':       True,
                        }  # optional: 'pop': 'E4'

cfg.analysis['plotRatePSD']={
                        'include':          ['all'], 
                        'timeRange':        [1500,cfg.duration], 
                        'binSize':          1, 
                        'minFreq':          1, 
                        'maxFreq':          400, 
                        'transformMethod':  'morlet', 
                        'stepFreq':         1, 
                        'NFFT':             256, 
                        'noverlap':         128, 
                        'smooth':           0, 
                        'norm':             False, 
                        'overlay':          True, 
                        'popColors':        popColors_dict, 
                        'yLogScale':        True, 
                        'ylim':             None, 
                        'figSize':          (20, 20), 
                        'fontSize':         12, 
                        'lineWidth':        1.5, 
                        'saveFig':          True, 
                        }

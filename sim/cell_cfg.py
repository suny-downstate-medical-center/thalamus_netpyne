"""
cfg.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""

from netpyne import specs
import pickle
import NetPyNE_BBP
import numpy as np
import os
saveFigFlag = '' # declared here, but altered afterwards
sim_name = 'cell_validation_000020_'

#------------------------------------------------------------------------------
# Simulation options
#------------------------------------------------------------------------------
cfg                 = specs.SimConfig()     # object of class SimConfig to store simulation configuration
cfg.skipTime        = 5000                  # pre-run simulation for this amount of time 
cfg.duration        = 5.25*1e3 + cfg.skipTime              # Duration of the simulation, in ms
cfg.dt              = 0.025                 # Internal integration timestep to use
cfg.hParams         = {'v_init': -75, 'celsius': 34}
cfg.verbose         = False                 # Show detailed messages
cfg.printRunTime    = 0.1
cfg.connRandomSecFromList = False   # testing post SfN (2023-11-22)

################################################################################################################################################
# Adding patches from issues opened on S1 repo (
# https://github.com/BlueBrain/CoreNeuron/issues/843
# https://github.com/suny-downstate-medical-center/S1_netpyne/commit/513c6db577aa574be6d54c1e576cd28b43f4fcdf
# )
cfg.coreneuron = False
cfg.random123 = True

# testing change random seed of random123 to the same in BlueCfg file
# source: https://github.com/suny-downstate-medical-center/netpyne/blob/3b720f7617cd490d0b3b8923e4a23df450d7d0a5/netpyne/sim/run.py#L75 

cfg.base_random_seed = 100000

cfg.rand123GlobalIndex = cfg.base_random_seed

cfg.seeds = {'conn': cfg.base_random_seed, 
             'stim': cfg.base_random_seed, 
             'loc':  cfg.base_random_seed, 
             'cell': cfg.base_random_seed}

cfg.saveCellSecs        = True # needs to be true if you want to access the membrane current (sim.net.cells[0].secs['soma'](0.5).i_membrane)
cfg.use_fast_imem       = True # RECORD MEMBRANE CURRENT
################################################################################################################################################

cfg.recordTraces = {    'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'},
                        # 'V_ptr': {'var': 'ptr'},
                        }

cfg.recordCurrents=True
rec_curr = [('SK_E2','ik'),('TC_HH','ina'),('TC_HH','ik'),('TC_Nap_Et2','ina'),
            ('TC_iA','ik'),('TC_iL','ica'),('TC_iT_Des98','ica'),('TC_ih_Bud97','ih')]
if cfg.recordCurrents:
    for curr in rec_curr: cfg.recordTraces.update({'i__soma_0__'+curr[0]+'__'+curr[1]:{'sec':'soma_0','loc':0.5,'mech':curr[0],'var':curr[1]},})


# cfg.recordStep = 1              # Step size in ms to save data (eg. V traces, LFP, etc)
cfg.recordStep = cfg.dt         # Step size in ms to save data (eg. V traces, LFP, etc)

# cfg.recordStim = False  
# cfg.recordTime = False

#------------------------------------------------------------------------------
# Load BBP circuit
#------------------------------------------------------------------------------

# --- Path to BBP model files

# cfg.base_dir='/Users/joao'
cfg.base_dir=os.path.expanduser("~")

cfg.BBP_rootFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'
cfg.morphologyFolder_asc            = cfg.BBP_rootFolder+'/sonata/morphologies/morphologies_asc'
cfg.virtual_input_spikes            = cfg.BBP_rootFolder+'/sonata/simulation_spike_files'
cfg.ct_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ct_noise.dat'
cfg.ml_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ml_noise.dat'

# --- Path to BBP data files
cfg.BBP_dataFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_paper_data'
cfg.Fig_4_1_dataFolder              = cfg.BBP_dataFolder+'/Fig4_wakefulness_spontaneous'
cfg.Fig_4_1_spikes                  = cfg.Fig_4_1_dataFolder+'/out.h5'
cfg.Fig_4_1_traces                  = cfg.Fig_4_1_dataFolder+'/soma.bbp.h5'

# --- Path to NetPyNE files
cfg.NetPyNE_rootFolder              = cfg.base_dir+'/Research/Models/BBP/thalamus_netpyne'
cfg.NetPyNE_sim                     = cfg.NetPyNE_rootFolder+'/sim'
cfg.NetPyNE_templateCells           = cfg.NetPyNE_rootFolder+'/mod'
cfg.NetPyNE_exportedCells           = cfg.NetPyNE_rootFolder+'/cells/morphologies_swc'
cfg.NetPyNE_JSON_cells              = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies'
cfg.NetPyNE_conn                    = cfg.NetPyNE_rootFolder+'/conn'
cfg.NetPyNE_data                    = cfg.NetPyNE_rootFolder+'/data'
cfg.NetPyNE_node_pathway_edges      = cfg.NetPyNE_conn+'/node_pathway_edges'
cfg.NetPyNE_input_noise_savePath    = cfg.NetPyNE_rootFolder+'/conn/external_inputs'

cfg.NetPyNE_JSON_cells_validation           = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies_validation'
cfg.NetPyNE_node_pathway_edges_validation   = cfg.NetPyNE_conn+'/node_pathway_edges_validation'


#------------------------------------------------------------------------------
# Sim Properties
#------------------------------------------------------------------------------

cfg.addNoiseIClamp=True 
cfg.scale_threshold_current = 0

#------------------------------------------------------------------------------
# Sim Properties
#------------------------------------------------------------------------------
cfg.select_microcircuit = None
# --- Convert cell models
cfg.convertCellMorphologies=False

from neuron import h
rand = h.Random(cfg.base_random_seed)

# cfg.target = 'VPL_TC'
cfg.target = 'Rt_RC'
if  cfg.target == 'VPL_TC':
    cfg.select_thal_gids = [
                            # 35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847, # sim 01
                            # 33798, 34368, 36219, 39232, 34389, 34242, 38102, 35703, 38487, 41067, 37463, 38468, 36711, 34932, 38346, 34503, 36248, 41454, 36721, 33741, # sim 02
                            # 40602, 34274, 41534, 33640, 36881, 34859, 36169, 38276, 37409, 34707, 38440, 41237, 38052, 36302, 33602, 41247, 38036, 39429, 38474, 35824, # sim 03
                            # 38651, 37968, 40213, 42177, 40168, 40215, 41723, 36655, 38134, 41695, 42422, 42460, 36521, 38775, 35220, 35162, 34349, 36440, 35739, 34954, # sim 04
                            # 37256, 41168, 39751, 38748, 33967, 35343, 40876, 39755, 36185, 41399, 39299, 38971, 37093, 37917, 37599, 34471, 39745, 39477, 42073, 36043, # sim 05
                            # 41388, 38169, 34773, 34401, 41379, 37475, 38090, 40659, 37782, 38709, 42405, 41353, 41307, 40641, 37685, 39390, 39239, 35684, 34363, 37548, # sim 06
                            # 34976, 35398, 34977, 34209, 37751, 39276, 38218, 41138, 37435, 37966, 42345, 35864, 34506, 40105, 38470, 34418, 37141, 39362, 33676, 36674, # sim 07
                            # 36748, 36059, 35158, 40735, 35483, 42198, 34433, 41390, 39229, 40044, 37740, 40122, 36364, 35113, 38793, 40560, 36857, 37553, 41271, 39981, # sim 08 
                            # 41439, 38171, 39183, 41890, 37925, 37824, 38002, 35649, 41579, 38806, 37520, 40430, 33822, 39202, 37863, 41253, 33571, 35332, 35748, 39340, # sim 09
                            33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806, # sim 10
                            ]
elif cfg.target == 'Rt_RC':  
    cfg.select_thal_gids = [
                            # 30550, 30541, 29582, 29174, 29908, 33109, 31521, 32579, 32893, 32954, 32696, 32933, 33334, 31927, 30299, 29934, 30694, 31191, 31989, 32369, # sim 11
                            # 30242, 30823, 29379, 31241, 31793, 31492, 32974, 30653, 29993, 30022, 29770, 32501, 29195, 29892, 30730, 30655, 32740, 32640, 28671, 28831, # sim 12
                            # 28660, 29828, 31704, 28988, 29183, 29690, 31254, 30838, 31637, 30922, 30182, 33200, 28663, 31412, 31625, 31778, 29791, 31120, 30543, 29184, # sim 13
                            # 28612, 30652, 32453, 32047, 29522, 32049, 29342, 31907, 30072, 32729, 29735, 32221, 30986, 33224, 31309, 30551, 31296, 29803, 29007, 30947, # sim 14
                            # 28805, 30849, 33463, 29657, 30946, 32631, 31840, 30892, 31646, 31738, 31315, 29086, 29040, 28852, 29608, 30025, 31528, 32662, 32781, 31170, # sim 15
                            # 32479, 33190, 31420, 28785, 30084, 31972, 30225, 30872, 30506, 32036, 33089, 33362, 32299, 32620, 29371, 32292, 32978, 32313, 32267, 30174, # sim 16
                            # 33014, 30007, 31239, 28733, 32470, 31044, 28694, 29087, 29476, 29687, 30990, 29126, 31800, 28834, 31881, 28925, 30252, 29621, 29094, 29304, # sim 17
                            # 31400, 29526, 31674, 32147, 31113, 29861, 32413, 29052, 30152, 29731, 29205, 31864, 31393, 33031, 30772, 28731, 30090, 33325, 30891, 29863, # sim 18
                            # 30403, 31638, 32406, 33043, 30905, 32926, 30014, 30813, 30854, 29679, 29049, 31751, 31816, 29689, 32540, 29846, 28833, 32411, 32730, 29805, # sim 19
        
                            32846, 29328, 30216, 32641, 29663, 30936, 32371, 29722, 31923, 30609, 32591, 30670, 31012, 31181, 33204, 31924, 32040, 28873, 33230, 31602, # sim 20
        
                            ]
    # cfg.select_thal_gids = [31571]
else:                           
    cfg.target = 'VPL_TC'
    cfg.select_thal_gids = [41337]
    # cfg.target == 'Rt_RC'  
    # cfg.select_thal_gids = [31571]


''''
VPL_TC: 
    good:       35703   34932   34368   
    ok:         36219   
    not bad:    36248   34503   33798   33741
'''

cfg.select_thal_gids.sort()

# --- structure to store the values of the holding currents that will be played with a h.vector in init.py
circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new( cfg_file=cfg.sonataConfigFile)
cfg.store_holdingCurrent={}
cfg.store_thresholdCurrent={}
for thal_gid in cfg.select_thal_gids:
    cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
    cfg.store_holdingCurrent.update(    {str(thal_gid): cell_properties['@dynamics:holding_current']})
    cfg.store_thresholdCurrent.update(  {str(thal_gid): cell_properties['@dynamics:threshold_current']})

# --- Cell Model
cfg.loadCellModel   = False;    cfg.convertCellModel= True;    cfg.saveCellModel   = True; 
cfg.plotSynLocation = False;    cfg.plotCellShape = False

# --- Connectivity Type
cfg.connType        = 'original'    # adds original edge connectivity to the target cell model

# --- Connectivity
cfg.th_updateConns      = False     # (True): generates the edge connectivity to the target cell model | (False): Loads from file | Defaults to True if file doesn't exist
cfg.th_topological      = False      # Redistributes the synapses onto the dendrites according to Driver/IN/RT/CT hierarchy
cfg.th_useTableVals     = False     # Replace edges values for the values in Table S2 from the paper
cfg.saveIndividualNoise = True      # Saves individual input noise files for each postsynaptic cell for later plotting
cfg.th_connOverrides    = True      # Alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model

cfg.removeChemical      = False
cfg.removeElectrical    = True      # wont work with TRN inputs as vecstims

cfg.th_singleTcPop      = False
cfg.modType             = 'Prob_original'
cfg.th_stimType         = 'VecStim'


if cfg.target == 'VPL_TC':
    cfg.th_inter_sources    = ['Rt_RC','VPL_IN']
    cfg.th_select_pathways  = [
                                'MedialLemniscus_projections',
                                'CorticoThalamic_projections', 
                                'thalamus_neurons|Rt_RC',
                                'thalamus_neurons|VPL_IN'
                                ]
elif cfg.target == 'Rt_RC':
    cfg.th_inter_sources    = ['VPL_TC','Rt_RC']
    cfg.th_select_pathways  = [
                                'CorticoThalamic_projections', 
                                'thalamus_neurons|Rt_RC',
                                # 'thalamus_neurons|VPL_IN'
                                'thalamus_neurons|VPL_TC'
                                ]

# --- Extracellular Calcium concentration

# cfg.it_shift    = None # mV --- test: 2023_11_24 (Iavarone, 2019)
# cfg.it_actshift = -9 # mV --- test: 2023_11_24 (Iavarone, 2019)
cfg.it_shift    = None # mV --- test: 2023_11_24 (Iavarone, 2019)
cfg.it_actshift = None # mV --- test: 2023_11_24 (Iavarone, 2019)

cfg.re_rescale = 1
# cfg.re_rescale = 1.55

cfg.cao_secs            = None
# cfg.rescaleUSE          = None # From BlueConfig file

# cfg.cao_secs            = 1.2
# cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file

# cfg.cao_secs            = 1.2
cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file


# --- Change conn weights
base_weight=0.001
cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1*base_weight,
                            'CorticoThalamic_projections':  1*base_weight,
                            'thalamus_neurons|Rt_RC':       1*base_weight,
                            'thalamus_neurons|VPL_IN':      1*base_weight,
                            'thalamus_neurons|VPL_TC':      1*base_weight,
                            }

# --- Change biophysics
cfg.removeExtraCurrents=False

cfg.changeCurrents  = []
cfg.rescaleCurrents = []
# cfg.rescaleCurrents = [('TC_iA','gk_max',0.5)]
# cfg.rescaleCurrents = [('TC_iA','gk_max',0.25)]

# cfg.rescaleCurrents = [('TC_HH','gna_max',0.017196474865155516*100)]


# --- Stims
# cfg.use_BBP_holdingCurrent      = True
# cfg.use_BBP_thresholdCurrent    = not cfg.use_BBP_holdingCurrent # test: use one or the other

cfg.use_BBP_thresholdCurrent    = True # test: use both currents
cfg.use_BBP_holdingCurrent      = False # test: use both currents

cfg.add_current_stims_noise   = False
noise_mean = 0
noise_var = 0.001
# noise_var = abs(cfg.hParams['v_init']*0.001)
cfg.noise_string = 'uniform('+str(noise_mean)+','+str(noise_var)+')'


cfg.add_current_stims   = False
cfg.current_stim_amp    = 0
# cfg.current_stim_amp    = iclamp_amp
cfg.current_stim_start    = 0
cfg.current_stim_duration = cfg.duration

# cfg.current_stim_amp    = 0.061
# cfg.current_stim_amp    = 0.084
# cfg.current_stim_start    = 0

cfg.th_boostFibers = False

# --- Load Spike Times
cfg.th_spikes_file = cfg.Fig_4_1_spikes

# if cfg.th_topological:          saveFigFlag+='_topological'
# if cfg.th_useTableVals:         saveFigFlag+='_tableVals'
# if cfg.th_connOverrides:        saveFigFlag+='_connOverr'
# if cfg.cao_secs is not None:    saveFigFlag+='_ca|'+str(cfg.cao_secs)
# if cfg.removeExtraCurrents:     saveFigFlag+='_noIH_noIA'
# if cfg.add_current_stims_noise:       saveFigFlag+='_preIClampNoise'
# if cfg.add_current_stims:       saveFigFlag+='_IClamp_'+str(cfg.current_stim_amp)
# if len(cfg.changeCurrents)>0:   saveFigFlag+='_changedCurr'
# if len(cfg.rescaleCurrents)>0:  saveFigFlag+='_rescaledCurr'
# if np.prod(list(cfg.rescale_conn_weight.values()))!=1:  saveFigFlag+='_rescaledWeights_'+'_'.join(map(str, list(cfg.rescale_conn_weight.values())))
# if cfg.th_boostFibers:          saveFigFlag+='_boostFibers'
# if cfg.th_singleTcPop:          saveFigFlag+='_singleTcPop'
# if cfg.select_microcircuit:     saveFigFlag+='_mc|'+str(cfg.select_microcircuit)
# if cfg.use_BBP_holdingCurrent:  saveFigFlag+='_holdingCurr'
# if cfg.use_BBP_thresholdCurrent:saveFigFlag+='_thresholdCurr'
# if cfg.it_shift is not None:    saveFigFlag+='_itShift|'+str(cfg.it_shift)
# if cfg.it_actshift is not None: saveFigFlag+='_itActShift|'+str(cfg.it_actshift)

# for pathway in cfg.th_select_pathways:
#     if 'CorticoThalamic' in pathway: saveFigFlag+='_ct'
#     if 'MedialLemniscus' in pathway: saveFigFlag+='_ml'
#     if 'thalamus_neurons|Rt_RC' in pathway: saveFigFlag+='_rt'
#     if 'thalamus_neurons|VPL_IN' in pathway: saveFigFlag+='_in'
#     if 'thalamus_neurons|VPL_TC' in pathway: saveFigFlag+='_tc'

#------------------------------------------------------------------------------
# File Name
#------------------------------------------------------------------------------

cfg.sim_tag = sim_name

# --- Selecting type of MOD file to be used
if cfg.modType == 'Prob_original':saveFigFlag+='_MOD_ThOriginal' # S1 Probabilistic  implementation of the BBP mod files

if cfg.connType == 'original':
    c_name=''
    # if len(cfg.select_thal_gids)==1: c_name = str(cfg.select_thal_gids[0])
    # else:
    #     for cell in cfg.select_thal_gids:c_name+='|'+str(cell)
    outputName = cfg.sim_tag+cfg.target+'__'+c_name+saveFigFlag

fileName = outputName
cellFigName = outputName+'.png'

if   cfg.connType == 'original':    folderName = 'single_cell_inputs'
else:                               folderName = ''

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.filename            = fileName          # Set file output name
cfg.savePickle          = False             # Save params, network and sim output to pickle file
cfg.saveJson            = True              # Save params, network and sim output to JSON file
cfg.saveDataInclude     = ['simData', 'simConfig', 'netParams', 'net']
cfg.saveCellConns       = True

cfg.saveFolder          = '../data/init_sims/'+folderName

#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------

# cfg.analysis['plotRaster']  = { 'orderInverse': True, 'labels': None, 'markerSize':1,'saveFig':'../data/raster_plot.png','dpi':1000}           
cfg.analysis['plotTraces']  = { 
                                'include': list(range(len(cfg.select_thal_gids))),
                                # 'include': [0,1],
                                # 'include': [(pop, 0) for pop in popList],
                                # 'include': ['all'],
                                'timeRange': [0, cfg.duration],
                                'ylim': [-90, 60],
                                'overlay':True,
                                'saveFig':'../init_figs/'+cellFigName,
                                # 'saveFig':True,
                                }
# cfg.analysis['plotShape']   = { 'saveFig':True}

# cfg.analysis['plot2Dnet'] = {'figSize':(8, 20),'saveFig': 'model_2dnet__.png'}                                                # plot 2D cell positions and connections
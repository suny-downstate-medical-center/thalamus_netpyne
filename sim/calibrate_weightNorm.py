import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 30})

import NetPyNE_BBP
import Build_Net as BN

########################################################################################################################################################################################################
# # snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
# from pydoc import source_synopsis
# import sys
# from    matplotlib  import  pyplot  as plt
# print("Matplotlib backend (default): %s" %plt.get_backend())
# modules = []
# for module in sys.modules:
#     if module.startswith('matplotlib'):
#         modules.append(module)
# for module in modules:
#     sys.modules.pop(module)
# import matplotlib
# matplotlib.use("MacOSX")
# from    matplotlib  import  pyplot  as plt
# print("Matplotlib backend (dynamic): %s" %plt.get_backend())

########################################################################################################################################################################################################
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Loading a list of gids")
        # VPL_TC
        cellNumbers = [
                        # VPL_TC
                        # 36963, # selected VPM

                        # 33602, 
                        # 33640, 33676, 33702, 33741, 33774, 33798, 33806, 33822, 33967,
                        # 34274, 34349, 34368, 34389, 34471, 34503, 34569, 34859, 34954, 35103, 
                        # 35162, 35220, 35343, 35703, 35739, 35824, 35864, 35879, 36043, 36169, 
                        # 36185, 36219, 36248, 36302, 36440, 36521, 36655, 36721, 36881, 
                        # 37055, 37082, 37093, 37256, 37409, 37463, 37475, 37548, 37599, 37685,

                        # Rt_RC
                        # 31864, # selected TRN

                        # # cAD RT cells
                        # 31738, 
                        # # 31191, 32846, 33325,

                        # cNAD RT cells with less afterhyperpolarization
                        31924, 
                        # 29476, 
                        # 29934, 30838, 30891, 30905, 31400, 31528, 32267


                        # 32729, 
                        # 32640, 30947, 30225, 
                        # 31412, 30025, 33334, 28612, 29621,
                        # 30609, 28663, 29052, 29328, 29657, 29195, 32662, 32047, 30670, 33463, 
                        # 30072, 32579
                       ]
        
        # # VPL_TC
        # cellNumbers = [
        #                 33602, 33640, 33676, 33702, 33741, 33774, 33798, 33806, 33822, 33967, 
        #                 34274, 34349, 34368, 34389, 34471, 34503, 34569, 34859, 34954, 35103, 
        #                 35162, 35220, 35343, 35703, 35739, 35824, 35864, 35879, 36043, 36169, 
        #                 36185, 36219, 36248, 36302, 36440, 36521, 36655, 36721, 36881, 36963, 
        #                 37055, 37082, 37093, 37256, 37409, 37463, 37475, 37548, 37599, 37685
        #                ]
        
        # # Rt_RC
        # cellNumbers = [
        #                 32729, 32640, 30947, 30225, 31864, 31412, 30025, 33334, 28612, 29621, 
        #                 30609, 28663, 29052, 29328, 29657, 29195, 32662, 32047, 30670, 33463, 
        #                 30072, 32579
        #                ]

        print(cellNumbers)
        # print ("Script need a cell number between 0 and 42510")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=42510:
        print ("Comparing BBP and Netpyne Traces of:")
        cellNumber = int(sys.argv[1])

        cellNumbers=[cellNumber]

    # --- Dict to store the calibrated currents
    map_weightNorm = {}
    # # --- Dict to store the calibrated currents
    # map_scale_gpas = {}

    # --- Run config
    timeStep        = 0.025
    recordStep      = 0.025
    timesimulation  = 750
    # timesimulation  = 1000
    durationstim    = timesimulation
    delayhold       = 0
    delaystim       = 0

    durationhold = timesimulation

    ampPercent_holding_current  =  100.0 # of @dynamics:holding_current
    # KLeak                       =  [0, 5, 10, 15, 20] # of @dynamics:holding_current
    # scale_gpas                 =  [1, 0.8, 0.6, 0.4, 0.2] # of @dynamics:holding_current
    scale_gpas                 =  [1, 5, 10, 15, 20] # of @dynamics:holding_current

    # ampPercent_holding_current   =  [100.0, 100.0, 100.0, 100.0, 100.0] # of @dynamics:holding_current
    # ampPercent_threshold_current =  [50.0,  100.0, 150.0, 200.0, 250.0] # of @dynamics:threshold_current

    calibrated_currents_file='../cells/calibrated_currents/calibrated_currents.json'
    with open(calibrated_currents_file, 'r') as file: map_ThresholdCurrent = json.load(file)

    calibrated_currents_file2='../cells/calibrated_currents/calibrated_currents2.json'
    with open(calibrated_currents_file2, 'r') as file: map_ThresholdCurrent2 = json.load(file)
    for new_cell in map_ThresholdCurrent2.keys():
        if new_cell not in map_ThresholdCurrent.keys():map_ThresholdCurrent.update({new_cell:map_ThresholdCurrent2[new_cell]})

    calibrated_currents_file3='../cells/calibrated_currents/calibrated_currents3.json'
    with open(calibrated_currents_file3, 'r') as file: map_ThresholdCurrent3 = json.load(file)
    for new_cell in map_ThresholdCurrent3.keys():
        if new_cell not in map_ThresholdCurrent.keys():map_ThresholdCurrent.update({new_cell:map_ThresholdCurrent3[new_cell]})

    # --- Simplifying mophologies
    select_morphologies=cellNumbers
    map_ThresholdCurrent_ = {key:map_ThresholdCurrent[key] for key in map_ThresholdCurrent for select_morphology in select_morphologies if str(select_morphology) in key}
    map_ThresholdCurrent=map_ThresholdCurrent_

    ####################################################################################################################################

    from netpyne import sim
    from netpyne import specs
    import pickle

    # cfg = specs.SimConfig()     
    from cfg import cfg
    
    cfg.duration = timesimulation ## Duration of the sim, in ms  
    cfg.dt = timeStep
    cfg.hParams = {'celsius': 34, 'v_init': -70}  
    # cfg.hParams = {'celsius': 34, 'v_init': -65}  
    cfg.verbose = False
    cfg.createNEURONObj = True
    cfg.createPyStruct = True
    # cfg.cvode_active = False
    
    cfg.includeParamsLabel = False
    cfg.printPopAvgRates = True
    cfg.checkErrors = False

    ####################################################################################################################################
    cfg.modType = 'Prob_original'
    ####################################################################################################################################
    cfg.map_targetRMP = { 'VPL_TC': -70, 'Rt_RC': -70}
    cfg.map_ThresholdCurrent = map_ThresholdCurrent
    cfg.map_targetRMP_byCellTemplate = {}
    for cellTemplate in cfg.map_ThresholdCurrent.keys():
        if   'VPL_TC' in cellTemplate: cfg.map_targetRMP_byCellTemplate.update({cellTemplate:cfg.map_targetRMP['VPL_TC']})
        elif 'Rt_RC'  in cellTemplate: cfg.map_targetRMP_byCellTemplate.update({cellTemplate:cfg.map_targetRMP['Rt_RC']})

    cfg.addThresholdCurrentPops_byCellTemplate={cellTemplate:NetPyNE_BBP.Utils.cubic_extrapolation(cfg.map_ThresholdCurrent[cellTemplate]['rmp'], cfg.map_ThresholdCurrent[cellTemplate]['i_thresh'], cfg.map_targetRMP_byCellTemplate[cellTemplate]) for cellTemplate in cfg.map_ThresholdCurrent.keys()}

    #------------------------------------------------------------------------------
    # Analysis and plotting 
    #------------------------------------------------------------------------------
    # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------

    # netParams = specs.NetParams()   # object of class NetParams to store the network parameters
    from netParams import netParams

    circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new(cfg_file=cfg.sonataConfigFile)
    cell_properties_dict={}
    netParams.popParams={}

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # --- Synaptic Mechanisms
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # --- Dummy exc and inh synMechs for testing
    netParams.synMechParams = BN.BuildNetwork.getSynMechParams()

    # --- Selecting type of MOD file to be used
    modAMPANMDA = 'ProbAMPANMDA_EMS_original';  
    modGABA     = 'ProbGABAAB_EMS_original'
    print('\n\t>>\tMOD template\tAMPA: ',   modAMPANMDA, '\tGABA: ', modGABA)

    TableS1_barreloidThalamus   = NetPyNE_BBP.StoreParameters.getTableS1_barreloidThalamus()
    for pre_pop in TableS1_barreloidThalamus.keys():
        if   pre_pop == 'L6A':  mod_file = modAMPANMDA; mech_flag='exc'
        elif pre_pop == 'MLe':  mod_file = modAMPANMDA; mech_flag='exc'
        elif pre_pop == 'VPM':  mod_file = modAMPANMDA; mech_flag='exc'
        else:                   mod_file = modGABA;     mech_flag='inh'
        for post_pop in TableS1_barreloidThalamus[pre_pop].keys():
            syn_values = TableS1_barreloidThalamus[pre_pop][post_pop]
            if cfg.modType == 'Prob':syn_values['paper_reference_values'].update({'n_rrp_vesicles':1}) # only present in the Prob MOD, not in the Det MOD
            edge_dict = NetPyNE_BBP.CreateNetPyNE.modTemplate(syn_values['paper_reference_values'],mod_file)
            netParams.synMechParams.update({'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag:edge_dict})

        # Changing USE to 1, so that syns are always triggered
        for synMechName in netParams.synMechParams.keys():
            if 'Use' in netParams.synMechParams[synMechName].keys():
                print('fixing USE to 1')
                netParams.synMechParams[synMechName]['Use']=1

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # --- Importing cells
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    count_traces=0
    store_trace_names=[]
    for cellNumber_ind, cellNumber in enumerate(cellNumbers):

        print('Processing cell # ', cellNumber_ind, ' \t gid: ', cellNumber, '\t cells left: ', len(cellNumbers)-cellNumber_ind)

        #------------------------------------------------------------------------------
        # Cell parameters
        #------------------------------------------------------------------------------
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,cellNumber)

        cell_properties_dict.update({str(cell_properties.name):cell_properties})

        # map_KLeak.update({str(cell_properties.mtype)+'__'+str(cell_properties.name):{'rmp':[],
        #                                                                                 'kleak':KLeak,
        #                                                                                 }})
        # map_scale_gpas.update({str(cell_properties.mtype)+'__'+str(cell_properties.name):{ 'rmp':[],
        #                                                                                     'scale_gpas':scale_gpas,
        #                                                                                     }})
        
        cell_name = str(cell_properties.mtype)+'__'+str(cell_properties.name)
        map_weightNorm.update({cell_name:{}})

        gid = cell_properties.name
        MorphoName = cell_properties['morphology'] + '.swc'
        hocName = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc'
        cellName = hocName + '_' + str(cell_properties.name)

        cell_label=cell_properties['mtype']+'__'+str(cell_properties.name)
        cell_type = str(cell_properties.name)

        cellRule = netParams.importCellParams(
            label           = cell_label, 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':cell_type},
            # conds           = {'cellType':str(thal_gid)},
            cellArgs        = [cell_properties.name,cfg.NetPyNE_exportedCells,cell_properties['morphology']+'.swc'], 
            importSynMechs  = True, 
            somaAtOrigin    = False, 
            cellInstance    = False,
            )
        
        # --- Fixing up split soma section
        if 'soma_1' in netParams.cellParams[cell_label]['secs'].keys(): 
            print('\t--- Check behavior of model after fixing extra soma_1 section issue')
            netParams.cellParams[cell_label] = NetPyNE_BBP.ConvertMorphology.mergeSomaSections(cell = netParams.cellParams[cell_label])
            single_soma_sec=True
        else: single_soma_sec=False


        #------------------------------------------------------------------------------
        # Select Syn Mechs
        #------------------------------------------------------------------------------

        if   'VPL_TC' in cell_label: 
            # synMechs                = ['syn|MLe|VPM|exc']
            # depolarization_target   = [3.62]
            # synMechs                = ['syn|L6A|VPM|exc']
            # depolarization_target   = [0.071]
            synMechs                = ['syn|TRN|VPM|inh']
            depolarization_target   = [-1.31]
        elif 'Rt_RC'  in cell_label: 
            # synMechs                = ['syn|VPM|TRN|exc']
            # depolarization_target   = [6.79]
            # synMechs                = ['syn|L6A|TRN|exc']
            # depolarization_target   = [2.4*0.071]
            synMechs                = ['syn|TRN|TRN|inh']
            depolarization_target   = [-1.31]

        
        # if   'VPL_TC' in cell_label: 
        #     synMechs                = ['syn|MLe|VPM|exc', 'syn|L6A|VPM|exc', 'syn|TRN|VPM|inh']
        #     depolarization_target   = [             3.62,             0.071,              1.31]
        # elif 'Rt_RC'  in cell_label: 
        #     synMechs                = ['syn|VPM|TRN|exc', 'syn|L6A|TRN|exc', 'syn|TRN|TRN|inh']
        #     depolarization_target   = [             6.79,        2.4*0.071,              1.31]
        #                                                    # Golshani,2001               # TRN->TRN no ref - estimated from TRN->VPM
        
        #------------------------------------------------------------------------------
        # Population parameters
        #------------------------------------------------------------------------------
        for sec in netParams.cellParams[cell_label]['secs'].keys():
        # for sec in ['soma_0']:
        # for sec in ['soma_0','dend_0','dend_1','dend_10','dend_50','dend_77',]:
            map_weightNorm[cell_name].update({sec:{}})

            if 'axon' in sec: continue  # ignores axon section
            for iloc in range(netParams.cellParams[cell_label]['secs'][sec]['geom']['nseg']):
                slot    = (2*iloc)+1
                nslots  = (netParams.cellParams[cell_label]['secs'][sec]['geom']['nseg']*2)
                loc     = slot/(nslots)
                # print(iloc,slot,nslots,loc)

                pop_name = cell_type+'__'+sec+'__'+str(iloc)
                netParams.popParams[pop_name] = {'cellType': cell_type, 'numCells': 1} 

                stim_times      = [  300,     450,    600]
                stim_weights    = [0.001,   0.002,  0.003]
                t_interval      = 100 # ms - slice of the traces that are used to calculate the peaks

                se_stim_delay   = 150

                map_weightNorm[cell_name][sec].update({iloc:{}})
                count_traces+=1
                trace_name=cell_name+'__'+sec+'__'+str(iloc)
                store_trace_names.append(trace_name)

                for stim_ind,stim_time in enumerate(stim_times):
                    
                    # map_weightNorm[cell_name][sec][iloc].update({stim_ind:{}})

                    netParams.stimSourceParams['NStim__'+cell_type+'__'+sec+'__'+str(iloc)+'__'+str(stim_ind)+'__source'] = {'type': 'NetStim', 'interval': 1, 'number': 1, 'start': stim_time, 'noise': 0}
                    netParams.stimTargetParams['NStim__'+cell_type+'__'+sec+'__'+str(iloc)+'__'+str(stim_ind)+'__target'] = {
                            'source':   'NStim__'+cell_type+'__'+sec+'__'+str(iloc)+'__'+str(stim_ind)+'__source',
                            'conds':    {'pop': pop_name},
                            'sec':      sec,
                            'loc':      loc,
                            'weight':   stim_weights[stim_ind],
                            'delay':    0,
                            'synMech':  synMechs[0],
                            }
                    
                    # # won't work

                    # netParams.stimSourceParams['SeStim__'+cell_type+'__'+sec+'__'+str(loc)+'__'+stim_ind+'__source'] = {'type': 'SEClamp', 'interval': 1, 'number': 1, 'start': stim_time, 'noise': 0}
                    # netParams.stimTargetParams['SeStim__'+cell_type+'__'+sec+'__'+str(loc)+'__'+stim_ind+'__target'] = {
                    #         'source':   'SeStim__'+cell_type+'__'+sec+'__'+str(loc)+'__source',
                    #         'sec':      sec,
                    #         'loc':      loc,
                    #         'weight':   stim_weights[stim_ind],
                    #         'delay':    0,
                    #         'conds':    {'pop': pop_name}}
                    

                    # netParams.stimSourceParams['SeStim__'+cell_type+'__'+sec+'__'+str(loc)+'__'+stim_ind+'__source'] = {'type': 'SEClamp', 'dur1': , 'amp1': ,'dur2': , 'amp2': ,'dur3': , 'amp3':,}
                    # netParams.stimTargetParams['SeStim__'+cell_type+'__'+sec+'__'+str(loc)+'__'+stim_ind+'__target'] = {'source': , 'sec':sec, 'loc': loc, 'conds': {'pop': pop_name}}

                    


        #------------------------------------------------------------------------------
        # Current inputs 
        #------------------------------------------------------------------------------

        netParams.stimSourceParams['HClamp__'+cell_type+'__source'] = {'type': 'IClamp', 'delay': 0, 'dur': cfg.duration, 'amp': (ampPercent_holding_current/100)*cell_properties['@dynamics:holding_current']}
        netParams.stimTargetParams['HClamp__'+cell_type+'__target'] = { 'source':   'HClamp__'+cell_type+'__source', 
                                                                        'conds':    {'cellType': cell_type},
                                                                        'sec':      'soma_0', 
                                                                        'loc':      0.5}

        threshold_percent_multiplier = (cfg.addThresholdCurrentPops_byCellTemplate[cell_label])/100
        amp = threshold_percent_multiplier * cell_properties['@dynamics:threshold_current']
        # amp = threshold_percent_multiplier * sim.net.params.store_cell_properties[sim.net.cells[cell_ind].tags['label'][0]]['@dynamics:threshold_current']
        netParams.stimSourceParams['TClamp__'+cell_type+'__source'] = {'type': 'IClamp', 'delay': 0, 'dur': cfg.duration, 'amp': amp}
        netParams.stimTargetParams['TClamp__'+cell_type+'__target'] = { 'source':   'TClamp__'+cell_type+'__source', 
                                                                        'conds':    {'cellType': cell_type},
                                                                        'sec':      'soma_0', 
                                                                        'loc':      0.5}

        stimSourceParams_list = list(netParams.stimSourceParams.keys())
        for key in stimSourceParams_list:
            if ('HClamp__' in key) or ('TClamp__' in key) or ('NStim__' in key):    pass
            else:                                                                   del netParams.stimSourceParams[key] # removing extra keys from default barreloid_netParams
        stimTargetParams_list = list(netParams.stimTargetParams.keys())
        for key in stimTargetParams_list:
            if ('HClamp__' in key) or ('TClamp__' in key) or ('NStim__' in key):    pass
            else:                                                                   del netParams.stimTargetParams[key] # removing extra keys from default barreloid_netParams

        continue

        # netParams.stimSourceParams['TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source'] = {'type': 'IClamp', 'delay': 0, 'dur': cfg.duration, 'amp': (ampPercent_threshold_current[ind]/100)*cell_properties['@dynamics:threshold_current']}
        # netParams.stimTargetParams['TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__target'] = {'source': 'TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source', 'conds': {'pop': str(cell_properties.name)+'_'+str(ind)},'sec': 'soma_0', 'loc': 0.5}


        # if cfg.addHoldingCurrent:
        #     threshold_percent_multiplier = ampPercent_holding_current/100
        #     for pop in cfg.addHoldingCurrentPops:
        #         netParams.stimSourceParams['HoldingCurrent_source__'+pop]   = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': amp}
        #         netParams.stimTargetParams['HoldingCurrent_target__'+pop]   = {'source': 'HoldingCurrent_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}
        # if cfg.addThresholdCurrent:
        #     threshold_percent_multiplier = (cfg.addThresholdCurrentPops_byCellTemplate[sim.net.cells[cell_ind].tags['label'][0]])/100
        #     amp = threshold_percent_multiplier * sim.net.params.store_cell_properties[sim.net.cells[cell_ind].tags['label'][0]]['@dynamics:threshold_current']
        #     for pop in cfg.addThresholdCurrentPops:
        #         netParams.stimSourceParams['ThresholdCurrent_source__'+pop] = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': 0}
        #         netParams.stimTargetParams['ThresholdCurrent_target__'+pop] = {'source': 'ThresholdCurrent_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}
    # sys.exit()
    cfg.connType='' # to avoid running default pops on cfg/netParams file

    #------------------------------------------------------------------------------
    # Recording
    #------------------------------------------------------------------------------

    allpops = list(netParams.popParams.keys())
    # allpops = [str(cell_properties.name)]
    # # allpops = ['L1_1','L1_2','L1_3','L1_4']

    cfg.recordCells=['all']
    # recordPops=[]
    # for ind, curr_percent in enumerate(ampPercent_threshold_current):
    #     recordPops.append(str(cell_properties.name)+'_'+str(ind))
    # cfg.recordCells = recordPops
      # which cells to record from
    # cfg.recordCells = allpops  # which cells to record from
    cfg.recordTraces = {'V_soma_0': {'sec':'soma_0', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
    cfg.recordStim = True
    cfg.recordTime = True
    cfg.recordStep = recordStep            

    cfg.simLabel = 'weightNorm_'+str(cell_properties.name)
    cfg.saveFolder = './validation'
    # cfg.filename =                	## Set file output name
    cfg.savePickle = False         	## Save pkl file
    cfg.saveJson = False           	## Save json file
    cfg.saveDataInclude = ['simConfig', 'netParams'] ## 'simData' , 'simConfig', 'netParams'
    cfg.backupCfgFile = None 		##  
    cfg.gatherOnlySimData = False	##  
    cfg.saveCellSecs = False			##  
    cfg.saveCellConns = False		##  


    # --- Runs all at once
    # sim.createSimulateAnalyze(netParams, cfg)

    # --- Runs stepwise
    sim.initialize(
        simConfig = cfg, 	
        netParams = netParams)  				# create network object and set cfg and net params
    sim.net.createPops()               			# instantiate network populations
    sim.net.createCells()              			# instantiate network cells based on defined populations
    sim.net.connectCells()            			# create connections between cells based on params
    sim.net.addStims() 							# add network stimulation

    # - Add changes to hObj here

    sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
    sim.runSim()                      			# run parallel Neuron simulation  
    sim.gatherData()                  			# gather spiking data and cell info from each node
    sim.analyze()

    t_vec = list(sim.simData['t'])

    timestamps = [(int(t_ind/sim.cfg.dt),int((t_ind+t_interval)/sim.cfg.dt)) for t_ind in stim_times]

    netpyneTraces = []
    netpyneTracesList = []

    # for ind in range(len(scale_gpas)*len(cellNumbers)):
    # for ind in range(len(sim.simData['V_soma_0'])):
    for ind in range(count_traces):
        netpyneTraces.append(np.array(sim.simData['V_soma_0']['cell_'+ str(ind)]))
        netpyneTracesList.append(list(sim.simData['V_soma_0']['cell_'+ str(ind)]))        

    traces_dict={}
    slices_dict={}
    weightNorm_dict={}
    # weightNorm_dict_full={}
    for trace_ind,trace_name in enumerate(store_trace_names):
        
        cellPop,cellNumber,sec,loc = trace_name.split('__')
        if sec not in weightNorm_dict.keys(): weightNorm_dict.update({sec:[]})# initializing dictionary with a scaling factor value of (1)
        # if sec not in weightNorm_dict.keys(): weightNorm_dict_full.update({sec:{}})# initializing dictionary with a scaling factor value of (1)
        voltage_trace = list(sim.simData['V_soma_0']['cell_'+ str(trace_ind)])
        traces_dict.update({trace_name:voltage_trace})

        voltage_trace_slices = [voltage_trace[t_stamp1:t_stamp2] for (t_stamp1,t_stamp2) in timestamps]
        time_trace_slices    = [t_vec[t_stamp1:t_stamp2]         for (t_stamp1,t_stamp2) in timestamps]

        # find peaks
        if   'exc' in synMechs[0]: 
            max_voltages        = [max(voltage_trace_slice) for voltage_trace_slice in voltage_trace_slices]
            max_voltages_time   = [time_trace_slices[v_idx][np.argmax(voltage_trace_slice)] for v_idx,voltage_trace_slice in enumerate(voltage_trace_slices)]
        elif 'inh' in synMechs[0]: 
            max_voltages        = [min(voltage_trace_slice)                          for voltage_trace_slice in voltage_trace_slices]
            max_voltages_time   = [time_trace_slices[v_idx][np.argmin(voltage_trace_slice)] for v_idx,voltage_trace_slice in enumerate(voltage_trace_slices)]
        # voltage at the instant the stim begins
        starting_voltages       = [voltage_trace_slice[0] for voltage_trace_slice in voltage_trace_slices]
        # change in membrane voltage
        PSP_depolarization      = [max_voltage-(starting_voltages[max_voltage_index]) for max_voltage_index,max_voltage in enumerate(max_voltages)]

        '''
        improve the PSP calculation code to take the membrane potential just before the stim, instead of considering it will be at -70mV
        '''
        PSP_depolarization_1_2=PSP_depolarization[1]-PSP_depolarization[0]
        PSP_depolarization_2_3=PSP_depolarization[2]-PSP_depolarization[1]
        PSP_depolarization_1_3=(PSP_depolarization[2]-PSP_depolarization[0])/2

        scaling_factor_1 = depolarization_target[0]/PSP_depolarization_1_2
        scaling_factor_2 = depolarization_target[0]/PSP_depolarization_2_3
        scaling_factor_3 = depolarization_target[0]/PSP_depolarization_1_3

        weightNorm_dict[sec].append(scaling_factor_1)
        # weightNorm_dict_full[trace_name]={  'PSP':{ 
        #                                         'PSP_depolarization_1_2':PSP_depolarization_1_2,
        #                                         'PSP_depolarization_2_3':PSP_depolarization_2_3,
        #                                         'PSP_depolarization_1_3':PSP_depolarization_1_3,
        #                                         },
        #                                     'scaling_factor':{
        #                                         'scaling_factor_1':scaling_factor_1,
        #                                         'scaling_factor_2':scaling_factor_2,
        #                                         'scaling_factor_3':scaling_factor_3,
        #                                         }
        # }

        slice_dict={weight:voltage_trace_slices[weight_ind] for weight_ind,weight in enumerate(stim_weights)}
        slices_dict.update({trace_name:slice_dict})

        colors=['r', 'cyan', 'limegreen']
        plt.figure(figsize=(10,10))
        plt.plot(t_vec,voltage_trace,'-k')
        for voltage_trace_slice_ind, voltage_trace_slice in enumerate(voltage_trace_slices):
            plt.plot(time_trace_slices[voltage_trace_slice_ind],voltage_trace_slices[voltage_trace_slice_ind], '-',c=colors[voltage_trace_slice_ind])
            plt.plot(max_voltages_time,max_voltages, 'o',c='magenta',markersize=10)
        plt.ylim(-75,-65)
        plt.savefig('../figs/figs_weighNorm/'+synMechs[0]+'/test_trace_wnorm_'+trace_name+'_'+synMechs[0]+'.png')
        # plt.savefig('../figs/figs_weighNorm/test_trace_wnorm_'+trace_name+'_'+synMechs[0]+'.png')

    

        # for v_ind, v in enumerate(voltage_trace):
        #     if t_vec[v_ind]
            
            
        #     voltage_trace_slices = 
        # for t in stim_times:

    cellPop     = store_trace_names[0].split('__')[0]
    cellNumber  = store_trace_names[0].split('__')[1]

    sm=synMechs[0].split('|')
    sm_join='_'.join(sm)

    filename = '../cells/calibrated_weightNorm/calibrated_weightNorm__'+cellPop+'__'+cellNumber+'__'+sm_join+'.json'
    with open(filename, "w") as json_file: json.dump(weightNorm_dict, json_file, indent=4)
    
    # sys.exit()
    
    '''
    # for rec_ind,cell_recording in enumerate(cfg.recordCells):
    for trace_ind,trace_name in enumerate(store_trace_names):
        # cell_num,trace_ind = cell_recording.split('_')
        cell_name,sec,iloc= trace_name.split('__')

        # map_scale_gpas[str(cell_properties_dict[str(cell_num)].mtype)+'__'+str(cell_num)]['rmp'].append(netpyneTraces[rec_ind][-1])

    filename = '../cells/calibrated_weightNorm/calibrated_weightNorm.json'
    with open(filename, "w") as json_file: json.dump(map_weightNorm, json_file, indent=4)

    # netpyneTraces = runnetpyne(cellNumbers)






    ####################################################################################################################################






    # for cellInd, cellNumber in enumerate(cellNumbers):
    #     print('Running sim for cell # ', str(cellInd), ' - ', len(cellNumbers)-cellInd, ' cells left')

    #     # --- Import cfg parameters
    #     from cfg import cfg
    #     # --- Load dictionary with thalamic circuit properties
    #     circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new(cfg_file=cfg.sonataConfigFile)
    #     cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,cellNumber)
    #     # sys.exit()

    #     # compareTraces(cell_properties)
    #     # --- Moving Compare Traces outside of function to re-run plots iteractively
    #     netpyneTraces = runnetpyne(cell_properties)
    #     rmp = [trace[-1] for trace in netpyneTraces]
    #     map_ThresholdCurrent.update({str(cell_properties.mtype)+'__'+str(cell_properties.name):{'rmp':rmp,
    #                                                                                             'i_thresh':ampPercent_threshold_current,
    #                                                                                             }})

    #     # # plot both traces overlayed
    #     # timeRange = [0, timesimulation]
    #     # # ~ ylim = [-100, 40]
    #     # figSize = (20,100)
    #     # fig = plt.figure(figsize=figSize)  # Open a new figure
        
    #     # t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
        
    #     # for ind, curr_percent in enumerate(ampPercent_threshold_current):
    #     #     netpyneTrace = netpyneTraces[ind]
    #     #     BBPTrace = BBPTraces[ind]
    #     #     plt.subplot(len(ampPercent_threshold_current), 1, ind+1)
    #     #     plt.ylabel('V (mV)')
    #     #     plt.plot([t[0],t[-1]],[-84,-84],linewidth=2, color='b',linestyle=':')
    #     #     plt.plot([t[0],t[-1]],[-64,-64],linewidth=2, color='r',linestyle=':')
    #     #     plt.plot(t[:len(netpyneTrace)], netpyneTrace, linewidth=5, color='k', label='I = %.1f, NetPyNE' % (ampPercent_threshold_current[ind]))
    #     #     plt.plot(t[:len(BBPTrace)], BBPTrace, linewidth=4, color='grey', label='I = %.1f, NEURON' % (ampPercent_threshold_current[ind]))  # linestyle=':'
    #     #     plt.xlabel('Time (ms)')
    #     #     # plt.xlim(0, timesimulation)
    #     #     plt.xlim(2000, 4000)
    #     #     plt.ylim(-100, 40)
    #     #     plt.title('Percent of threshold current: '+ str(ampPercent_threshold_current[ind]))
    #     #     plt.grid(True)
    #     #     plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
    #     # plt.ion()
    #     # plt.tight_layout()

    #     # MEName = cell_properties['mtype'] + '_' + cell_properties['etype']
    #     # save_fig_path='./validation/nrn-7-8-2_'+str(cell_properties.name)+'_'+MEName+'_comparison_traces_soma_voltage_4steps_'+str(tonic_hold_percent)+'.png'
    #     # plt.savefig(save_fig_path)
    #     # print ('Figure Saved in ', save_fig_path)   

    #     # ########################################################################################################################################################################################################################

    #     # for ind, curr_percent in enumerate(ampPercent_threshold_current):
    #     #     figSize = (20,20)
    #     #     fig = plt.figure(figsize=figSize)  # Open a new figure
    #     #     plt.rcParams.update({'font.size': 40})
    #     #     t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
    #     #     select_inds=[0,1,2,3,8,9,10,15]
    #     #     if ind not in select_inds:continue
    #     #     netpyneTrace = netpyneTraces[ind]
    #     #     BBPTrace = BBPTraces[ind]
    #     #     # plt.subplot(len(ampPercent_threshold_current), 1, ind+1)
    #     #     plt.ylabel('V (mV)')
    #     #     plt.plot([t[0],t[-1]],[-84,-84],linewidth=4, color='b',linestyle=':')
    #     #     plt.plot([t[0],t[-1]],[-64,-64],linewidth=4, color='r',linestyle=':')
    #     #     plt.plot(t[:len(netpyneTrace)], netpyneTrace,   linewidth=10, color='k',        label='I = %.1f, NetPyNE' % (ampPercent_threshold_current[ind]))
    #     #     plt.plot(t[:len(BBPTrace)],     BBPTrace,       linewidth=3, color='fuchsia',  label='I = %.1f, NEURON'  % (ampPercent_threshold_current[ind]))  # linestyle=':'
    #     #     plt.xlabel('Time (ms)')
    #     #     # plt.xlim(0, timesimulation)
    #     #     plt.xlim(2000, 4000)
    #     #     plt.xticks([2000,2500,3000,3500,4000])
    #     #     plt.ylim(-110, 40)
    #     #     plt.title('Percent of threshold current: '+ str(ampPercent_threshold_current[ind]))
    #     #     plt.grid(True)
    #     #     plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
    #     #     plt.ion()
    #     #     plt.tight_layout()

    #     #     MEName = cell_properties['mtype'] + '_' + cell_properties['etype']
    #     #     save_fig_path='./validation/nrn-7-8-2_'+str(cell_properties.name)+'_'+MEName+'_comparison_traces_soma_voltage_4steps_'+str(tonic_hold_percent)+'proposal_'+str(ind)+'.png'
    #     #     plt.savefig(save_fig_path)
    #     #     print ('Figure Saved in ', save_fig_path)


    #     # ########################################################################################################################################################################################################################

    #     # fig = plt.figure(figsize=(30,20))

    #     # for i in range(2):
    #     #     plt.subplot(2,1,i+1)
    #     #     # plt.subplot(1,2,i+1)
    #     #     plt.plot([t[0],t[-1]],[-84,-84],linewidth=2, color='b',linestyle=':')
    #     #     plt.plot([t[0],t[-1]],[-64,-64],linewidth=2, color='r',linestyle=':')
            
            
    #     # plt.subplot(2,1,1,projection='3d')
    #     # # plt.subplot(1,2,1)
    #     # plt.ylabel('V (mV)')
            
    #     # norm = plt.Normalize(0, len(ampPercent_threshold_current))
    #     # colormap=plt.cm.plasma

    #     # for ind_, curr_percent in enumerate(ampPercent_threshold_current):
    #     #     ind=len(ampPercent_threshold_current)-ind_-1
    #     #     netpyneTrace = netpyneTraces[ind]
    #     #     BBPTrace = BBPTraces[ind]
    #     #     color = colormap(norm(ind))

    #     #     plt.subplot(2,1,1,projection='3d')
    #     #     # plt.subplot(1,2,1)
    #     #     plt.xlim(2000, 4000)
    #     #     plt.plot(t[:len(netpyneTrace)], [ind for i in range(len(netpyneTrace))], netpyneTrace, linewidth=2, color=color, label='I = %.1f, NetPyNE' % (ampPercent_threshold_current[ind]))
    #     #     # plt.plot(t[:len(netpyneTrace)], [ind for i in range(len(netpyneTrace))], netpyneTrace, linewidth=2, color='k', label='I = %.1f, NetPyNE' % (ampPercent_threshold_current[ind]))
    #     #     # plt.plot(t[:len(netpyneTrace)], [ind for i in range(len(netpyneTrace))], netpyneTrace, linewidth=5, color='k', label='I = %.1f, NetPyNE' % (ampPercent_threshold_current[ind]),alpha=1-(ind/len(ampPercent_threshold_current)))
            
    #     #     plt.subplot(2,1,2,projection='3d')
    #     #     # plt.subplot(1,2,2)
    #     #     plt.xlim(2000, 4000)
    #     #     plt.plot(t[:len(BBPTrace)], [ind for i in range(len(netpyneTrace))],  BBPTrace, linewidth=4, color=color, label='I = %.1f, NEURON' % (ampPercent_threshold_current[ind]))  # linestyle=':'
        
    #     # plt.savefig('test_bbpVSnet')

    '''
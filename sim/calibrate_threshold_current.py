import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 30})

import NetPyNE_BBP

########################################################################################################################################################################################################
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Loading a list of gids")
        # # VPL_TC
        # cellNumbers = [
        #                 # VPL_TC
        #                 33602, 33640, 33676, 33702, 33741, 33774, 33798, 33806, 33822, 33967,
        #                 34274, 34349, 34368, 34389, 34471, 34503, 34569, 34859, 34954, 35103, 
        #                 35162, 35220, 35343, 35703, 35739, 35824, 35864, 35879, 36043, 36169, 
        #                 36185, 36219, 36248, 36302, 36440, 36521, 36655, 36721, 36881, 36963, 
        #                 37055, 37082, 37093, 37256, 37409, 37463, 37475, 37548, 37599, 37685,

        #                 # Rt_RC
        #                 32729, 32640, 30947, 30225, 31864, 31412, 30025, 33334, 28612, 29621,
        #                 30609, 28663, 29052, 29328, 29657, 29195, 32662, 32047, 30670, 33463, 
        #                 30072, 32579
        #                ]
        cellNumbers = [
                        # VPL_TC

                        # Rt_RC - calibrated_currents2
                        # 31738, 31191, 32846, 33325,

                        # # Rt_RC - calibrated_currents3
                        # 31924, 
                        # 29476, 
                        # # 29934, 30838, 30891, 30905, 31400, 31528, 32267

                        # VPL_TC and Rt_RC - calibrated_currents4 # 2024_10_16
                        37520, 41912,               # VPL_TC
                        28925, 30252, 31315, 32049  # Rt_RC

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
    map_ThresholdCurrent = {}

    # --- Run config
    timeStep        = 0.025
    recordStep      = 0.025
    timesimulation  = 1000
    durationstim    = timesimulation
    delayhold       = 0
    delaystim       = 0

    durationhold = timesimulation

    ampPercent_holding_current   =  [100.0, 100.0, 100.0, 100.0, 100.0] # of @dynamics:holding_current
    ampPercent_threshold_current =  [50.0,  100.0, 150.0, 200.0, 250.0] # of @dynamics:threshold_current


    ####################################################################################################################################

    from netpyne import sim
    from netpyne import specs
    import pickle

    # cfg = specs.SimConfig()     
    from cfg import cfg
    
    cfg.duration = timesimulation ## Duration of the sim, in ms  
    cfg.dt = timeStep
    cfg.hParams = {'celsius': 34, 'v_init': -65}  
    cfg.verbose = False
    cfg.createNEURONObj = True
    cfg.createPyStruct = True
    # cfg.cvode_active = False
    
    cfg.includeParamsLabel = False
    cfg.printPopAvgRates = True
    cfg.checkErrors = False

    #------------------------------------------------------------------------------
    # Analysis and plotting 
    #------------------------------------------------------------------------------
    # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------

    # netParams = specs.NetParams()   # object of class NetParams to store the network parameters
    from netParams import netParams

    cell_properties_dict={}
    netParams.popParams={}
    for cellNumber in cellNumbers:
        #------------------------------------------------------------------------------
        # Cell parameters
        #------------------------------------------------------------------------------
        circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new(cfg_file=cfg.sonataConfigFile)
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,cellNumber)

        cell_properties_dict.update({str(cell_properties.name):cell_properties})

        map_ThresholdCurrent.update({str(cell_properties.mtype)+'__'+str(cell_properties.name):{'rmp':[],
                                                                                        'i_thresh':ampPercent_threshold_current,
                                                                                        }})


        gid = cell_properties.name
        MorphoName = cell_properties['morphology'] + '.swc'
        hocName = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc'
        cellName = hocName + '_' + str(cell_properties.name)

        cell_label=cell_properties['mtype']+'__'+str(cell_properties.name)

        cellRule = netParams.importCellParams(
            label           = cell_label, 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':str(cell_properties.name)},
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



        plotCellShape=True
        if plotCellShape: NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_label],'_validated_cellShape',cell_label=cell_label,savefolder='../validation/validated_morphologies')

        for ind, curr_percent in enumerate(ampPercent_threshold_current):
            #------------------------------------------------------------------------------
            # Population parameters
            #------------------------------------------------------------------------------
            netParams.popParams[str(cell_properties.name)+'_'+str(ind)] = {'cellType': str(cell_properties.name), 'numCells': 1} 

            #------------------------------------------------------------------------------
            # Current inputs 
            #------------------------------------------------------------------------------

            netParams.stimSourceParams['HClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source'] = {'type': 'IClamp', 'delay': 0, 'dur': cfg.duration, 'amp': (ampPercent_holding_current[ind]/100)*cell_properties['@dynamics:holding_current']}
            netParams.stimSourceParams['TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source'] = {'type': 'IClamp', 'delay': 0, 'dur': cfg.duration, 'amp': (ampPercent_threshold_current[ind]/100)*cell_properties['@dynamics:threshold_current']}

            netParams.stimTargetParams['HClamp__'+str(cell_properties.name)+'_'+str(ind)+'__target'] = {'source': 'HClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source', 'conds': {'pop': str(cell_properties.name)+'_'+str(ind)},'sec': 'soma_0', 'loc': 0.5}
            netParams.stimTargetParams['TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__target'] = {'source': 'TClamp__'+str(cell_properties.name)+'_'+str(ind)+'__source', 'conds': {'pop': str(cell_properties.name)+'_'+str(ind)},'sec': 'soma_0', 'loc': 0.5}
    
    cfg.connType='' # to avoid running default pops on cfg/netParams file

    #------------------------------------------------------------------------------
    # Recording
    #------------------------------------------------------------------------------

    allpops = [str(cell_properties.name)]
    # allpops = ['L1_1','L1_2','L1_3','L1_4']

    cfg.recordCells=netParams.popParams.keys()
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

    cfg.simLabel = 'validation_netpyne_neuron_'+str(cell_properties.name)
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

    netpyneTraces = []
    netpyneTracesList = []

    for ind in range(len(ampPercent_threshold_current)*len(cellNumbers)):
        netpyneTraces.append(np.array(sim.simData['V_soma_0']['cell_'+ str(ind)]))
        netpyneTracesList.append(list(sim.simData['V_soma_0']['cell_'+ str(ind)]))        

    for rec_ind,cell_recording in enumerate(cfg.recordCells):
        cell_num,trace_ind = cell_recording.split('_')
        map_ThresholdCurrent[str(cell_properties_dict[str(cell_num)].mtype)+'__'+str(cell_num)]['rmp'].append(netpyneTraces[rec_ind][-1])

    # filename = '../cells/calibrated_currents/calibrated_currents.json'
    # filename = '../cells/calibrated_currents/calibrated_currents2.json'
    # filename = '../cells/calibrated_currents/calibrated_currents3.json'
    filename = '../cells/calibrated_currents/calibrated_currents4.json'
    with open(filename, "w") as json_file: json.dump(map_ThresholdCurrent, json_file, indent=4)

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



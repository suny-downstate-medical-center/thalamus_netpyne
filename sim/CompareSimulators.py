import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 30})
import NetPyNE_BBP
from netpyne import sim

class CompareSimulators():
    def configPaths(cfg):

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
        
        return cfg

    def getCellProperties(cellNumber):
        # --- Import cfg parameters
        
        from netpyne.specs import simConfig as cfg
        
        # from netpyne import specs
        # cfg = specs.simConfig()

        cfg = CompareSimulators.configPaths(cfg)
        # --- Load dictionary with thalamic circuit properties
        circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new(cfg_file=cfg.sonataConfigFile)
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,cellNumber)
        return cell_properties


    def runneuron(cellNumber):
        print(' --- Running NEURON implementation (single cell) --- ')
        from neuron import h
        # --- Import cfg parameters
        
        # from netpyne.specs import simConfig as cfg
        from cell_cfg import cfg
        
        # from netpyne import specs
        # cfg = specs.simConfig()

        cfg = CompareSimulators.configPaths(cfg)
        cell_properties = CompareSimulators.getCellProperties(cellNumber)

        gid = cell_properties.name
        MorphoName      = cell_properties['morphology']+'.swc'
        hocName         = cell_properties['etype']
        MorphoFolder    = cfg.NetPyNE_exportedCells
        hocPath         = cfg.NetPyNE_templateCells+'/'+hocName+'.hoc'
        cellName = hocName + '_' + str(cell_properties.name)

        from validation_cellwrapper import loadCell
        cell = loadCell(hocPath         = hocPath, 
                        MorphoFolder    = MorphoFolder,
                        hocName         = hocName,
                        MorphoName      = MorphoName, 
                        gid             = cell_properties.name,
                        )

        soma_0 = cell.soma[0]

        BBPTraces = []
        BBPTracesList = []

        for protocol_ind in protocols.keys():
            print('Running NEURON sim # ', str(protocol_ind))
            stimulus = h.IClamp(0.5, sec=soma_0)

            stimulus.dur = durationstim # ms
            stimulus.delay = delaystim  # ms    
            
            stimulus.amp = (protocols[protocol_ind]['threshold']/100)*cell_properties['@dynamics:threshold_current']
            # stimulus.amp = (ampPercent_threshold_current[ind]/100)*cell_properties['@dynamics:threshold_current']
            # stimulus.amp = ampstim[ind]

            holding = h.IClamp(0.5, sec=soma_0)
            holding.dur = durationhold # ms
            holding.delay = delayhold  # ms    
            holding.amp = (protocols[protocol_ind]['holding']/100)*cell_properties['@dynamics:holding_current']

            recordings = {}
            recordings['time'] = h.Vector()
            recordings['soma_0(0.5)'] = h.Vector()

            recordings['time'].record(h._ref_t, recordStep)
            recordings['soma_0(0.5)'].record(cell.soma[0](0.5)._ref_v, recordStep)

            h.dt = timeStep
            h.celsius = 34
            h.v_init = -65
            # h.cvode_active(0)
            h.tstop = timesimulation # ms
            h.run()

            time = np.array(recordings['time'])
            soma_voltage = np.array(recordings['soma_0(0.5)'])

            BBPTraces.append(soma_voltage)
            BBPTracesList.append(list(soma_voltage))
        
        return BBPTraces
    
    def runneuron2(cellNumber):
        
        print(' --- Running NEURON implementation (pooled) --- ')

        from neuron import h
        # --- Import cfg parameters
        
        # from netpyne.specs import simConfig as cfg
        from cell_cfg import cfg
        
        # from netpyne import specs
        # cfg = specs.simConfig()

        cfg = CompareSimulators.configPaths(cfg)
        cell_properties = CompareSimulators.getCellProperties(cellNumber)

        gid = cell_properties.name
        MorphoName      = cell_properties['morphology']+'.swc'
        hocName         = cell_properties['etype']
        MorphoFolder    = cfg.NetPyNE_exportedCells
        hocPath         = cfg.NetPyNE_templateCells+'/'+hocName+'.hoc'
        cellName = hocName + '_' + str(cell_properties.name)

        from validation_cellwrapper import loadCell
        cell = loadCell(hocPath         = hocPath, 
                        MorphoFolder    = MorphoFolder,
                        hocName         = hocName,
                        MorphoName      = MorphoName, 
                        gid             = cell_properties.name,
                        )
        # soma_0 = cell.soma[0]
        # cellsList       = [cell for i in range(len(protocols.keys()))]
        cellsList       = [loadCell(hocPath         = hocPath, 
                                    MorphoFolder    = MorphoFolder,
                                    hocName         = hocName,
                                    MorphoName      = MorphoName, 
                                    gid             = cell_properties.name,
                                    ) for i in range(len(protocols.keys()))]
        somaSecsList    = [cellsList[i].soma[0] for i in range(len(protocols.keys()))]
        
        stimulusList    = [h.IClamp(0.5, sec=somaSecsList[i]) for i in range(len(protocols.keys()))]
        holdingStimList = [h.IClamp(0.5, sec=somaSecsList[i]) for i in range(len(protocols.keys()))]

        recordings = {}
        recordings.update({'time':h.Vector()})
        
        v_recs = {'soma_'+str(i)+'(0.5)':h.Vector() for i in range(len(protocols.keys()))}
        recordings.update(v_recs)

        recordings['time'].record(h._ref_t, recordStep)
        for i in range(len(protocols.keys())): recordings['soma_'+str(i)+'(0.5)'].record(cellsList[i].soma[0](0.5)._ref_v, recordStep)
        
        # recordings['soma_0(0.5)'].record(cell.soma[0](0.5)._ref_v, recordStep)

        for i, protocol_ind in enumerate(protocols.keys()):
            stimulusList[i].dur     = durationstim # ms
            stimulusList[i].delay   = delaystim  # ms    
            stimulusList[i].amp     = (protocols[protocol_ind]['threshold']/100)*cell_properties['@dynamics:threshold_current']

            holdingStimList[i].dur  = durationhold # ms
            holdingStimList[i].delay= delayhold  # ms    
            holdingStimList[i].amp  = (protocols[protocol_ind]['holding']/100)*cell_properties['@dynamics:holding_current']

        h.dt = timeStep
        h.celsius = 34
        h.v_init = -65
        # h.cvode_active(0)
        h.tstop = timesimulation # ms
        h.run()

        time = np.array(recordings['time'])
        BBPTraces = [np.array(recordings['soma_'+str(i)+'(0.5)']) for i in range(len(protocols.keys()))]
        BBPTracesList = [list(array) for array in BBPTraces]
        
        return BBPTraces

    def runnetpyne(cellNumber,cfg):

        print(' --- Running NetPyNE implementation (pooled) --- ')

        cfg = CompareSimulators.configPaths(cfg)
        cell_properties = CompareSimulators.getCellProperties(cellNumber)

        
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
        
        allpops = [str(cell_properties.name)]
        # allpops = ['L1_1','L1_2','L1_3','L1_4']

        recordPops=[]
        for protocol_ind in protocols.keys(): recordPops.append(str(cell_properties.name)+'_'+protocol_ind)
        cfg.recordCells = recordPops
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
        cfg.saveDataInclude = ['simConfig', 'netParams', 'simData']# , 'simConfig', 'netParams'
        cfg.backupCfgFile = None 		##  
        cfg.gatherOnlySimData = False	##  
        cfg.saveCellSecs = False			##  
        cfg.saveCellConns = False		##  

        #------------------------------------------------------------------------------
        # Analysis and plotting 
        #------------------------------------------------------------------------------
        # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
        #------------------------------------------------------------------------------
        # Cells
        #------------------------------------------------------------------------------
        cfg.cellmod =  {str(cell_properties.name): 'HH_full'} # (?) is it necessary?

        #------------------------------------------------------------------------------
        # Current inputs 
        #------------------------------------------------------------------------------
        cfg.addIClamp = 1

        # --- Clamping to burst mode (~ -84 mV)
        cfg.IClamp_00 = {'pop': str(cell_properties.name)+'_'+str(0),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['0']['holding']/100)*cell_properties['@dynamics:holding_current']}
        # --- Clamping to tonic mode (~ -64 mV)
        cfg.IClamp_01 = {'pop': str(cell_properties.name)+'_'+str(1),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['1']['holding']/100)*cell_properties['@dynamics:holding_current']}
        cfg.IClamp_02 = {'pop': str(cell_properties.name)+'_'+str(2),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['1']['holding']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp_03 = {'pop': str(cell_properties.name)+'_'+str(3),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['1']['holding']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp_04 = {'pop': str(cell_properties.name)+'_'+str(4),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['1']['holding']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp_05 = {'pop': str(cell_properties.name)+'_'+str(5),   'sec': 'soma_0', 'loc': 0.5, 'start': delayhold, 'dur': durationhold, 'amp': (protocols['1']['holding']/100)*cell_properties['@dynamics:threshold_current']}
        
        # --- Stimulation
        cfg.IClamp0  = {'pop': str(cell_properties.name)+'_'+str(0),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['0']['threshold']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp1  = {'pop': str(cell_properties.name)+'_'+str(1),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['1']['threshold']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp2  = {'pop': str(cell_properties.name)+'_'+str(2),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['2']['threshold']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp3  = {'pop': str(cell_properties.name)+'_'+str(3),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['3']['threshold']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp4  = {'pop': str(cell_properties.name)+'_'+str(4),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['4']['threshold']/100)*cell_properties['@dynamics:threshold_current']}
        cfg.IClamp5  = {'pop': str(cell_properties.name)+'_'+str(5),  'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': (protocols['5']['threshold']/100)*cell_properties['@dynamics:threshold_current']}

        cfg.connType='' # to avoid running default pops on cfg/netParams file
        cfg.select_thal_gids=[] # to avoid running default pops on cfg/netParams file
        cfg.verbose=True

        # netParams = specs.NetParams()   # object of class NetParams to store the network parameters
        # from netParams import netParams

        #------------------------------------------------------------------------------
        # Cell parameters
        #------------------------------------------------------------------------------
        from cell_netParams import netParams

        gid = cell_properties.name
        MorphoName = cell_properties['morphology'] + '.swc'
        hocName = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc'
        cellName = hocName + '_' + str(cell_properties.name)

        cellRule = netParams.importCellParams(
            label           = cell_properties['mtype']+'__'+str(cell_properties.name), 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':str(cell_properties.name)},
            # conds           = {'cellType':str(thal_gid)},
            cellArgs        = [cell_properties.name,cfg.NetPyNE_exportedCells,cell_properties['morphology']+'.swc'], 
            importSynMechs  = True, 
            somaAtOrigin    = False, 
            cellInstance    = False,
            )
        # netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')

        # cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
        #     conds={'cellType': cellName, 'cellModel': 'HH_full'},
        #     fileName='validation_cellwrapper.py',
        #     cellName='loadCell',
        #     cellInstance = True,
        #     cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
        # # netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')
        #------------------------------------------------------------------------------
        # Population parameters
        #------------------------------------------------------------------------------
        netParams.popParams={}
        for protocol_ind in protocols.keys():
            netParams.popParams[str(cell_properties.name)+'_'+protocol_ind] = {'cellType': str(cell_properties.name), 'numCells': 1} 
        
        # netParams.popParams[str(cell_properties.name)] = {'cellType': str(cell_properties.name), 'numCells': 1} 
        
        #------------------------------------------------------------------------------
        # Current inputs (IClamp)
        #------------------------------------------------------------------------------
        if cfg.addIClamp:
            for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
                params = getattr(cfg, key, None)
                [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

                #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

                # add stim source
                netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}

                # connect stim source to target
                netParams.stimTargetParams[key+'_'+pop] =  {
                    'source': key, 
                    'conds': {'pop': pop},
                    'sec': sec, 
                    'loc': loc}

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
        for protocol_ind in protocols.keys():
        # for c in range(1):
        # for c in range(4):
            netpyneTraces.append(np.array(sim.simData['V_soma_0']['cell_'+protocol_ind]))
            netpyneTracesList.append(list(sim.simData['V_soma_0']['cell_'+protocol_ind]))        
    
        return netpyneTraces

    def plotOverlayedTraces(BBPTraces,netpyneTraces):
        
        timeRange = [0, timesimulation]
        cell_properties = CompareSimulators.getCellProperties(cellNumber)
        MEName = cell_properties['mtype'] + '_' + cell_properties['etype']
        
        for ind, protocol_ind in enumerate(protocols.keys()):
            figSize = (20,20)
            fig = plt.figure(figsize=figSize)  # Open a new figure
            plt.rcParams.update({'font.size': 40})
            t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
            # select_inds=[0,1,2,3,8,9,10,15]
            # if ind not in select_inds:continue
            netpyneTrace = netpyneTraces[ind]
            BBPTrace = BBPTraces[ind]
            # plt.subplot(len(ampPercent_threshold_current), 1, ind+1)
            plt.ylabel('V (mV)')
            plt.plot([t[0],t[-1]],[-84,-84],linewidth=4, color='b',linestyle=':')
            plt.plot([t[0],t[-1]],[-64,-64],linewidth=4, color='r',linestyle=':')
            plt.plot(t[:len(netpyneTrace)], netpyneTrace,   linewidth=10, color='k',        label='I = %.1f, NetPyNE' % (protocols[protocol_ind]['threshold']))
            plt.plot(t[:len(BBPTrace)],     BBPTrace,       linewidth=3, color='fuchsia',  label='I = %.1f, NEURON'  % (protocols[protocol_ind]['threshold']))  # linestyle=':'
            plt.xlabel('Time (ms)')
            # plt.xlim(0, timesimulation)
            plt.xlim(skipTime, timesimulation)
            plt.xticks([500,750,1000,1250,1500,1750])
            plt.ylim(-110, 40)
            plt.title('Percent of threshold current: '+ str(protocols[protocol_ind]['threshold'])+ ' + holding current: '+str(protocols[protocol_ind]['holding']))
            plt.grid(True)
            plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
            plt.ion()
            plt.tight_layout()

            # save_fig_path='./validation/nrn-7-8-2_'+str(cell_properties.name)+'_'+MEName+'_NetPyNEvsNEURON__thresh_'+str(protocols[protocol_ind]['threshold'])+'__holding_'+str(protocols[protocol_ind]['holding'])+'proposal_'+str(ind)+'.png'
            save_fig_path='../stats/validation/simulator_validation/traces_figures/nrn-7-8-2_'+str(cell_properties.name)+'_'+MEName+'_NetPyNEvsNEURON_'+str(ind)+'__thresh_'+str(protocols[protocol_ind]['threshold'])+'__holding_'+str(protocols[protocol_ind]['holding'])+'.png'
            plt.savefig(save_fig_path)
            print ('Figure Saved in ', save_fig_path)

    def calculatePearsonCorrelation(BBPTraces,netpyneTraces,time_range=[500,1500]):
        cell_properties = CompareSimulators.getCellProperties(cellNumber)
        MEName = cell_properties['mtype'] + '_' + cell_properties['etype']

        time_indexes = [int(t/timeStep) for t in time_range]

        import numpy as np
        correlation_coefficient_dict={}
        for ind, protocol_ind in enumerate(protocols.keys()):
            netpyneTrace = netpyneTraces[ind]
            BBPTrace = BBPTraces[ind]
            # --- Making sure both traces have the same number of samples
            if len(netpyneTrace)!=len(BBPTrace):
                netpyneTrace = netpyneTrace[0:min([len(netpyneTrace),len(BBPTrace)])]
                BBPTrace     = BBPTrace[    0:min([len(netpyneTrace),len(BBPTrace)])]
            # --- Windowing traces
            # print('traces length - before')
            # print(len(BBPTrace),len(netpyneTrace))
            netpyneTrace    = netpyneTrace[ time_indexes[0]:time_indexes[1]]
            BBPTrace        = BBPTrace[     time_indexes[0]:time_indexes[1]]
            # print('traces length - after')
            # print(len(BBPTrace),len(netpyneTrace))
            
            # Compute the Pearson correlation coefficient
            correlation_coefficient = np.corrcoef(netpyneTrace, BBPTrace)[0, 1]
            correlation_coefficient_dict.update({ind:correlation_coefficient})

            print("Pearson correlation coefficient - Trace # :", str(ind),' - ', correlation_coefficient)
        
        # Writing to json
        json_object = json.dumps(correlation_coefficient_dict, indent=4)
        save_json_path = '../stats/validation/simulator_validation/traces_correlation/nrn-7-8-2_'+str(cell_properties.name)+'_'+MEName+'_NetPyNEvsNEURON_pearson_correlation.json'
        print('saving traces correlation value to json in\n',save_json_path)
        with open(save_json_path, "w") as outfile: outfile.write(json_object)

        return correlation_coefficient_dict

########################################################################################################################################################################


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Loading a list of gids")
        
        # cellNumbers=[36636, 40408, 39690, 35628, 40957, 41864, 37916, 41105, 38048, 37184, 34832] # VPL_TC
        cellNumbers=[  
                          35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847, 
                        #   33798, 34368, 36219, 39232, 34389, 34242, 38102, 35703, 38487, 41067, 37463, 38468, 36711, 34932, 38346, 34503, 36248, 41454, 36721, 33741, 
                        #   40602, 34274, 41534, 33640, 36881, 34859, 36169, 38276, 37409, 34707, 38440, 41237, 38052, 36302, 33602, 41247, 38036, 39429, 38474, 35824, 
                        #   38651, 37968, 40213, 42177, 40168, 40215, 41723, 36655, 38134, 41695, 42422, 42460, 36521, 38775, 35220, 35162, 34349, 36440, 35739, 34954, 
                        #   37256, 41168, 39751, 38748, 33967, 35343, 40876, 39755, 36185, 41399, 39299, 38971, 37093, 37917, 37599, 34471, 39745, 39477, 42073, 36043, 
                        #   41388, 38169, 34773, 34401, 41379, 37475, 38090, 40659, 37782, 38709, 42405, 41353, 41307, 40641, 37685, 39390, 39239, 35684, 34363, 37548, 
                        #   34976, 35398, 34977, 34209, 37751, 39276, 38218, 41138, 37435, 37966, 42345, 35864, 34506, 40105, 38470, 34418, 37141, 39362, 33676, 36674, 
                        #   36748, 36059, 35158, 40735, 35483, 42198, 34433, 41390, 39229, 40044, 37740, 40122, 36364, 35113, 38793, 40560, 36857, 37553, 41271, 39981, 
                        #   41439, 38171, 39183, 41890, 37925, 37824, 38002, 35649, 41579, 38806, 37520, 40430, 33822, 39202, 37863, 41253, 33571, 35332, 35748, 39340, 
                        #   33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806
                          ] # VPL_TC

        # cellNumbers=[29876, 29504, 31571, 29227, 31493, 31775, 29925, 29704, 32718, 30882, 30203] # Rt_RC
        # cellNumbers=[29876, 29504, 31571, 29227, 31493, 31775, 29925, 29704, 32718, 30882, 30203, 33166, 32301, 29200, 30285, 29500, 30581, 32723, 30643, 29341, 30782, 28753, 33139, 32929, 30827, 28663, 31312, 32315, 32536, 29923, 31316, 29016, 30522, 32613, 28789, 29627, 31361, 29526, 31724, 30625, 33332, 33043, 29545, 32377, 28962, 29628, 29197, 32239, 29367, 28829, 32051, 31055, 31896, 30505, 29659, 32686, 29222, 31655, 29422, 31013, 33094, 30286, 33270, 29171, 30185, 29143, 33013, 31034, 29955, 33134, 29008, 32305, 30295, 30023, 33012, 33022, 29472, 29617, 28647, 32632, 29020, 30149, 29673, 32855, 31476, 33280, 31324, 29375, 30308, 32462, 32348, 33020, 29785, 32680, 30785, 30915, 29113, 30897, 33133, 30848] # Rt_RC
        print(cellNumbers)
    
    elif int(sys.argv[1])>=0 and int(sys.argv[1])==1:
        print('Processing TC cells')
        cellNumbers=[
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

    elif int(sys.argv[1])>=0 and int(sys.argv[1])==2:
        print('Processing RT cells')
        cellNumbers=[
                        30550, 30541, 29582, 29174, 29908, 33109, 31521, 32579, 32893, 32954, 32696, 32933, 33334, 31927, 30299, 29934, 30694, 31191, 31989, 32369, # sim 11
                        # 30242, 30823, 29379, 31241, 31793, 31492, 32974, 30653, 29993, 30022, 29770, 32501, 29195, 29892, 30730, 30655, 32740, 32640, 28671, 28831, # sim 12
                        # 28660, 29828, 31704, 28988, 29183, 29690, 31254, 30838, 31637, 30922, 30182, 33200, 28663, 31412, 31625, 31778, 29791, 31120, 30543, 29184, # sim 13
                        # 28612, 30652, 32453, 32047, 29522, 32049, 29342, 31907, 30072, 32729, 29735, 32221, 30986, 33224, 31309, 30551, 31296, 29803, 29007, 30947, # sim 14
                        # 28805, 30849, 33463, 29657, 30946, 32631, 31840, 30892, 31646, 31738, 31315, 29086, 29040, 28852, 29608, 30025, 31528, 32662, 32781, 31170, # sim 15
                        # 32479, 33190, 31420, 28785, 30084, 31972, 30225, 30872, 30506, 32036, 33089, 33362, 32299, 32620, 29371, 32292, 32978, 32313, 32267, 30174, # sim 16
                        # 33014, 30007, 31239, 28733, 32470, 31044, 28694, 29087, 29476, 29687, 30990, 29126, 31800, 28834, 31881, 28925, 30252, 29621, 29094, 29304, # sim 17
                        # 31400, 29526, 31674, 32147, 31113, 29861, 32413, 29052, 30152, 29731, 29205, 31864, 31393, 33031, 30772, 28731, 30090, 33325, 30891, 29863, # sim 18
                        # 30403, 31638, 32406, 33043, 30905, 32926, 30014, 30813, 30854, 29679, 29049, 31751, 31816, 29689, 32540, 29846, 28833, 32411, 32730, 29805, # sim 19
                        # 32846, 29328, 30216, 32641, 29663, 30936, 32371, 29722, 31923, 30609, 32591, 30670, 31012, 31181, 33204, 31924, 32040, 28873, 33230, 31602, # sim 20
                        ]
        
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=42510:
        print ("Comparing BBP and Netpyne Traces of:")
        cellNumber = int(sys.argv[1])

        cellNumbers=[cellNumber]

    ############################## SIM parameters ##############################
    
    # --- Run config
    skipTime        = 500
    timeStep        = 0.05
    recordStep      = timeStep
    durationstim    = 800
    delaystim       = 50  + skipTime
    afterStim       = 500
    timesimulation  = delaystim + durationstim + afterStim
    
    # skipTime        = 2000
    # timeStep        = 0.05
    # recordStep      = timeStep
    # durationstim    = 1350.0
    # delaystim       = 200.0  + skipTime
    # timesimulation  = 2250.0 + skipTime

    durationhold = timesimulation
    delayhold = 0

    # --- Amplitude of the currents to clamp the cell
    protocols={
                '0':{'holding': 100, 'threshold': 250},     # depolarizing from hyperpolarization
                '1':{'holding': 0,   'threshold': 250},     # depolarizing from RMP
                '2':{'holding': 0,   'threshold':-200},     # strong   hyperpolarization
                '3':{'holding': 0,   'threshold':-50},      # weak     hyperpolarization
                '4':{'holding': 0,   'threshold': 200},     # strong   depolarization
                '5':{'holding': 0,   'threshold': 300},     # stronger depolarization
                }
    tonic_hold_percent = 100

    correlation_coeffs={}
    for cellNumber in cellNumbers:

        # --- Moving Compare Traces outside of function to re-run plots iteractively
        # BBPTraces       = CompareSimulators.runneuron(cellNumber)
        BBPTraces       = CompareSimulators.runneuron2(cellNumber)

        from cell_cfg import cfg
        netpyneTraces   = CompareSimulators.runnetpyne(cellNumber,cfg)

        CompareSimulators.plotOverlayedTraces(BBPTraces,netpyneTraces)

        correlation_coeffs.update({cellNumber:CompareSimulators.calculatePearsonCorrelation(BBPTraces,netpyneTraces,time_range=[delaystim,durationstim+afterStim])})
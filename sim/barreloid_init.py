"""
barreloid_init.py

Starting script to run NetPyNE-based thalamus model for thesis.

Usage:
    python barreloid_init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 8 nrniv -python -mpi barreloid_init.py

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""


# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
from    matplotlib  import  pyplot  as plt
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

import matplotlib; matplotlib.use('agg')  # to avoid graphics error in servers
from netpyne import sim
import json
import sys
import os
import numpy as np
import pickle

############################################################################################################
# --- Running simulation
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='barreloid_cfg.py', netParamsDefault='barreloid_netParams.py')

createOnly          = False
updateNetwork       = False

updateGapConns      = False # overwrites createOnly=True    | rerun updateNetwork separately afterwards (new conns are not added directly in the network)  |   cfg.conn_TRNe_TRNe_connList must be set to FALSE in the cfg file


addNoise            = True; addHoldingCurrent = True; addThresholdCurrent = True; 
addConnVariability  = True; 
addWeightNorm       = False # (default: False) wasn't woking properly and made the model worse/would need retuning

# addConnVariability  = False; print('\n\n\n\n\n\nWARNING: addConnVariability set to false to speed up network conn debugging\n\n\n\n\n\n')

plotNetworkShape    = False

plotConnData        = False
quantifyGapConns    = False
pruneGapJunctions   = False
evaluateModel       = False
updateCurrents      = False
plotFromCreate      = False
runPostAnalysis     = False

updateRiCurrent     = False # (default: False) Runs a simulation to measure the resting membrane potential in each cell with the network connected, but not stimulation
addRiCurrent        = False  # overwrites createOnly=True     |   

addSEClamp          = False
addSECalibratedIClamp = False

# --- Flags to check if the gap conn variables are properly set
if updateGapConns and cfg.conn_TRNe_TRNe_connList:
    print('Cannot updated gap conns if the connection template is TRUE')
    sys.exit()
if not updateGapConns and not cfg.conn_TRNe_TRNe_connList:
    print('Gap conns are disabled')
    # sys.exit()

############################################################################################################
# --- Custom package imports
from PlotFromSimCreate      import PlotFromSimCreate    as pfsc
from ResampleParameters     import ResampleParameters
from CalculateProbability   import CalculateProbability as CalProb
from ProcessESyns           import SynapticPruning
from ProcessESyns           import QuantifyESyns
from CurrentStim            import CurrentStim          as CS
from EvaluateModel          import Evaluate             as Eval
from Network3D              import Network3D
from PlotHeatmap            import MembranePotentialHeatmap
from Build_Net              import BuildNetwork
from WeightNormalization    import WeightNormalization

############################################################################################################
# --- Creating network using NetPyNE high-level specifications
sim.initialize(
    simConfig = cfg,     
    netParams = netParams)                  # create network object and set cfg and net params
sim.net.createPops()                        # instantiate network populations
sim.net.createCells()                       # instantiate network cells based on defined populations
sim.net.connectCells()                      # create connections between cells based on params
sim.net.addStims()                          # add network stimulation
# print(' -------- gathering cells')
# allCells = sim._gatherCells()
# print(' -------- gathering cell tags')
# allCellTags = sim._gatherAllCellTags()
# sys.exit()

############################################################################################################
# --- Calculate the RMP for each cell when there is no spiking 
if updateRiCurrent: 
    print('\n\n\n ---- Running updateRiCurrent to measure cell rmp for individualized current clamping --- \n\n\n')
    createOnly=False; addNoise=False; addHoldingCurrent=True; addThresholdCurrent=False; addRiCurrent=False; addConnVariability=True
    # sim.cfg.duration = 1501
    sim.cfg.analysis.plotTraces['include']=['all']
    try:    del sim.cfg.recordTraces['I_soma'] # to save space
    except: print('-> No currents being recorded')

# --- Callibrate the currents using SEClamp 
if addSEClamp: 
    print('\n\n\n ---- Running SEClamp to measure cell i when returning to rmp with no stimulation --- \n\n\n')
    createOnly=False; addNoise=False; addHoldingCurrent=True; addThresholdCurrent=False; addRiCurrent=False; addConnVariability=True
    # sim.cfg.duration = 1501
    sim.cfg.analysis.plotTraces['include']=['all']
    sim = CS.addRmpVoltageToSEClamp(sim)
    sys.exit()

if addSECalibratedIClamp:
    createOnly=False; addNoise=False; addHoldingCurrent=True; addThresholdCurrent=False; addRiCurrent=False; addConnVariability=True
    print('\n\n\n ---- Adding callibrated currents to cells --- \n\n\n')
    # sim.cfg.duration = 1501
    sim.cfg.analysis.plotTraces['include']=['all']
    sim = CS.addSECalibratedIClamp(sim)
    # sys.exit()

############################################################################################################
# --- Update gap connectivity at the NEURON level
if updateGapConns: 
    print('Updating Gap junction connectivity')
    createOnly=True
    Network3D.runGapDetection(sim,              pops=['TRN__pop','TRN_ring__pop'],
                              max_distance=5,   max_appositions_per_section=6,
                              max_gaps_perCell=sim.cfg.max_gaps_perCell,
                            #   max_gaps_perCell=15,
                              skipSoma=True,
                              seed=sim.cfg.base_random_seed,      
                              plotFig=True,
                              filename='../conn/gap_connecivity/gap_connecivity_template_testingClass'
                              )

############################################################################################################
# --- Modifying network at the NEURON level
if quantifyGapConns: QuantifyESyns.quantifyGap(sim)

if pruneGapJunctions:
    CalProb.verifyConns(sim,1)
    # sim =   SynapticPruning.pruneGaps(sim, prune_pops=['TRN__pop','TRN_ring__pop'], prune_vals=[5,5], verbose=False)
    sim =   SynapticPruning.pruneGaps(sim, prune_pops=['TRN__pop','TRN_ring__pop'], prune_vals=[20,20], verbose=True)
    CalProb.verifyConns(sim,1)

if addNoise            and sim.cfg.addNoiseIClamp:      sim, vecs_dict = CS.addNoiseIClamp(sim)
if addHoldingCurrent   and sim.cfg.addHoldingCurrent:   sim = CS.addHoldingCurrent(sim)
# to test
if addThresholdCurrent and sim.cfg.addThresholdCurrent: 
    # sim = CS.addThresholdCurrent(sim)
    sim, vecs_dict1 = CS.addThresholdCurrent(sim)
if addRiCurrent        and sim.cfg.addRiCurrent:        sim = CS.addRiCurrent(sim)
if addConnVariability:
    if sim.nhosts>1 and updateNetwork==True:
        print('Network cannot be updated in parallel mode')
        sys.exit()
    print('Updating network')
    sim = ResampleParameters.processSimObj(sim,updateNetwork=updateNetwork,
                                        #    weightSorting=True,
                                           weightSorting=False,
                                           weightSorting_ExcludeConns=['syn|TRN|TRN|inh','syn|MLe|VPM|exc'],
                                        #    weightSorting=False,
                                           plotFig=False)
    # sim = ResampleParameters.processSimObj(sim,updateNetwork=updateNetwork,plotFig=False)
    
    if updateRiCurrent or addSEClamp:
        print('Removing impact of MLe imputs')
        for cell_ind,cell in enumerate(sim.net.cells):
            for conn_ind in range(len(sim.net.cells[cell_ind].conns)):
                if 'conn|MLe' in sim.net.cells[cell_ind].conns[conn_ind]['label']: 
                    sim.net.cells[cell_ind].conns[conn_ind]['weight']=0
                    sim.net.cells[cell_ind].conns[conn_ind]['hObj'].syn().gmax=0
                    # sim.net.cells[cell_ind].conns[conn_ind]['hObj'].weight=0

if addWeightNorm: 
    sim = WeightNormalization.normalizeWeights(sim)
    print('Debug <WeightNormalization.normalizeWeights> - check weight assignment to [hObj].weight and weight pointer _ref_weight[0]')
    print('Run some tests to see if code is properly applied')
    print('Result: Looks like it wasnt working properly')

# to test
if sim.cfg.addModifyConnWeight_vecPlay: 
    print('Adding vector play on conn weight')
    sim, vecs_dict2 = CS.addConnWeightVectorPlay_g(sim,
                                                  connLabels        = list(sim.cfg.modifyConnWeight_vecPlay_dict.keys()),
                                                  connWeight_vals   = list(sim.cfg.modifyConnWeight_vecPlay_dict.values()))

############################################################################################################
# --- Inspecting network properties
#  -  Temporary fix to allow plotting some figures without running the simulation
if plotConnData:                            pfsc.plotNetPyNEConnMatrix(sim)
#  -  Evaluate model properties, plots pathway ratios and prints Convergence/Divergence stats
if evaluateModel:
    Eval.printStats(sim,plotFigs=True,printExtendedStats=True)
    # Eval.verifyConns(sim, cellStats=False, popStats=True, plotFigs=True)
    # Eval.verifyPopConns(sim,'MLe','VPM') # --- Evaluate average connections between two isolated pops
    # --- Evaluate the average number of sources from MLe into VPM
    count_MLe_souce_cells=[]
    for cell_ind,cell in enumerate(sim.net.cells):
        if 'VPM' in sim.net.cells[cell_ind].tags['pop']:
            MLe_souce_cells=[]
            for conn_ind, conn in enumerate(sim.net.cells[cell_ind].conns):
                if sim.net.cells[cell_ind].conns[conn_ind]['synMech']=='syn|MLe|VPM|exc':
                    MLe_souce_cells.append(sim.net.cells[cell_ind].conns[conn_ind]['preGid'])
            MLe_souce_cells_set=list(set(MLe_souce_cells))
            count_MLe_souce_cells.append(len(MLe_souce_cells_set))
    mean_conn_num=np.mean(count_MLe_souce_cells)
    med_conn_num=np.median(count_MLe_souce_cells)
    print('Mean MLe->VPM cell sources:   mean = ', mean_conn_num, '\t| median = ', med_conn_num)
    no_driver_cell_percentage=100*sum([1 for i in count_MLe_souce_cells if i==0])/len(count_MLe_souce_cells)
    print('Cells with no driver input = ',no_driver_cell_percentage, ' %')
    from collections import Counter
    print(Counter(count_MLe_souce_cells))

    from CalculateConnectivityStatistics import CalculateConnectivityStatistics as CalcConnStats
    CalcConnStats.runConvergenceStatistics(sim)

    # --- Counting the total number of connections
    count_conns=0
    count_gaps_=0 # must be halved, because conns are bidirectional (netpyne creates an instance of the gap conn in the presynaptic and another in the postsynaptic cell)
    for cell_ind, cell in enumerate(sim.net.cells):
        for conn_ind, conn in enumerate(sim.net.cells[cell_ind].conns):
            if conn['synMech']=='esyn': count_gaps_+=1
            else:                       count_conns+=1
    count_gaps=count_gaps_/2
    print('Total number of connections: chem = ', count_conns,' | elec = ', count_gaps)


############################################################################################################
# --- Running the network
if createOnly: print('\n\n\n ---- Running CREATE ONLY mode to inspect for bugs during network creation --- \n\n\n')
else:

    gatherOnly=False
    if gatherOnly:
        sim.gatherDataFromFiles(saveMerged=True)
        sim.net.allCells = [cell.__getstate__() for cell in sim.net.cells]
        sim.net.allPops = {label: pop.__getstate__() for label, pop in sim.net.pops.items()}
        sim.analysis.plotData()
    else:

        distributedSaving=sim.cfg.distributedSaving
        if distributedSaving:
            sim.setupRecording()                          # setup variables to record for each cell (spikes, V traces, etc)
            print(' --- Running sim')
            sim.runSim()                                  # run parallel Neuron simulation  
            print(' --- Gathering sim data - distributed')
            # distributed saving (to avoid errors with large output data)
            sim.saveDataInNodes()
            sim.gatherDataFromFiles(saveMerged=True)
            sim.analysis.plotData()
        else:
            sim.setupRecording()                          # setup variables to record for each cell (spikes, V traces, etc)
            print(' --- Running sim')
            sim.runSim()                                  # run parallel Neuron simulation  
            print(' --- Gathering sim data')
            sim.gatherData()                              # gather spiking data and cell info from each node
            print(' --- Analyzing sim data')
            sim.analyze()


    if updateRiCurrent: 
        save_currents = '../stims/RiCurrent/rmp_values.json'
        figs, fig_data_dict = sim.analysis.plotTraces(include=sim.cfg.updateRiCurrentPops,
                                timeRange=[500,1500],
                                saveData=save_currents,
                                )
        tracesData = fig_data_dict['tracesData']
        store_mean_rmp={}
        for rec_ind in range(len(tracesData)):
            for trace in tracesData[rec_ind].keys():
                if '_V_soma' in trace:
                    cell_gid_str = trace.split('_V_soma')[0].split('cell_')[1]
                    store_mean_rmp.update({cell_gid_str:np.mean(tracesData[rec_ind][trace])})
        # save_mean_rmp = '../stims/RiCurrent/mean_rmp_values.json'
        try:
            with open(sim.cfg.storeRiCurrents, 'w') as file:
                json.dump(store_mean_rmp, file)
        except:
            save_mean_rmp =  '../stims/RiCurrent/mean_rmp_values.json'
            with open(save_mean_rmp, 'w') as file:
                json.dump(store_mean_rmp, file)
    
    if addSEClamp: 
        
        save_currents = '../stims/RiCurrent/seclamp_values.json'
        figs, fig_data_dict = sim.analysis.plotTraces(include=sim.cfg.addSEClampPops,
                                timeRange=[0,sim.cfg.duration],
                                saveData=save_currents,
                                )
        tracesData = fig_data_dict['tracesData']
        store_v={}
        store_i  ={}
        for rec_ind in range(len(tracesData)):
            for trace in tracesData[rec_ind].keys():
                if '_V_soma' in trace:
                    cell_gid_str = trace.split('_V_soma')[0].split('cell_')[1]
                    store_v.update({cell_gid_str:list(tracesData[rec_ind][trace])})
                elif '_I_soma' in trace:
                    cell_gid_str = trace.split('_I_soma')[0].split('cell_')[1]
                    store_i.update({cell_gid_str:list(tracesData[rec_ind][trace])})
        
        t_vector = list(tracesData[0]['t'])
        # t_vector = [i*sim.cfg.recordStep for i in range(len(store_i['0']))]

        with open(sim.cfg.storeSEClamp_v, 'w') as file:
            json.dump(store_v, file)
        
        with open(sim.cfg.storeSEClamp_i, 'w') as file:
            json.dump(store_i, file)

        pop_gids = {pop:sim.net.allPops[pop]['cellGids'] for pop in sim.net.allPops.keys()} # this structure works with MPI (sim.net.pops dont)
        
        store_current_val={}
        store_current_val_peak={}
        for cell_gid in store_i.keys():    
            for pop_ in pop_gids.keys(): 
                if int(cell_gid) in pop_gids[pop_]: pop = pop_
            
            store_current_val.update({str(cell_gid):{'ind':[],'t':[],'i':[],'neg_peak_i':0}})
            
            t1 = sim.cfg.addSEClampParameters[pop]['dur1']
            t2 = sim.cfg.addSEClampParameters[pop]['dur2']+t1
            t3 = sim.cfg.addSEClampParameters[pop]['dur3']+t2

            for i_index, i_val in enumerate(store_i[cell_gid]):
                t = sim.cfg.recordStep * i_index
                if (t<t3) or t>t3+5.0: continue # only picks a small interval of 5 ms after the stimulation ends
                store_current_val[str(cell_gid)]['ind'].append( i_index)
                store_current_val[str(cell_gid)]['t'].append(t_vector[i_index])
                store_current_val[str(cell_gid)]['i'].append(i_val)
            store_current_val[str(cell_gid)]['neg_peak_i']=min(store_current_val[str(cell_gid)]['i'])

            store_current_val_peak.update({str(cell_gid):min(store_current_val[str(cell_gid)]['i'])})

        with open('../stims/RiCurrent/_seclamp_i_vals.json', 'w') as file:
            json.dump(store_current_val, file)
        with open('../stims/RiCurrent/_seclamp_i_peak.json', 'w') as file:
            json.dump(store_current_val_peak, file)
       
        # pop_gids = {pop:sim.net.pops[pop].cellGids for pop in sim.net.pops}
        # pop_gids = {'VPM__pop':         [0,346],
        #             'TRN__pop':         [346,456],
        #             'TRN_ring__pop':    [456,1336]}
        
        # keep debugging here - save the current spike during SEClamp
        slice_currents_dict={}
        for cell_gid in store_i.keys():    
            for pop_ in pop_gids.keys(): 
                if int(cell_gid) in pop_gids[pop_]: pop = pop_
            
            # if      int(cell_gid)<346:                                          pop='VPM__pop'
            # elif    (int(cell_gid)>=346)     and (int(cell_gid)<346+110):       pop='TRN__pop'
            # elif    (int(cell_gid)>=346+110) and (int(cell_gid)<346+110+880):   pop='TRN_ring__pop'
            # else:continue

            # pop = sim.net.cells[int(cell_gid)].tags['pop']
            slice_currents_dict.update({str(cell_gid):{'t1':[],'i1':[],'t2':[],'i2':[]}})
            t1 = sim.cfg.addSEClampParameters[pop]['dur1']
            t2 = sim.cfg.addSEClampParameters[pop]['dur2']+t1

            for i_index, i_val in enumerate(store_i[cell_gid]):
                t = sim.cfg.recordStep * i_index
                if   (t>=t1-1.0) and (t<t1+4.0):
                    slice_currents_dict[str(cell_gid)]['t1'].append(t)
                    slice_currents_dict[str(cell_gid)]['i1'].append(i_val)
                elif (t>=t2-1.0) and (t<t2+4.0):
                    slice_currents_dict[str(cell_gid)]['t2'].append(t)
                    slice_currents_dict[str(cell_gid)]['i2'].append(i_val)
        with open('../stims/RiCurrent/seclamp_i_sliced.json', 'w') as file:
            json.dump(slice_currents_dict, file)

        # load Json
        with open('../stims/RiCurrent/seclamp_i_sliced.json', 'r') as file:
            slice_currents_dict=json.load(file)            
    
    ########################################################################################################################
    # --- Load saved network template - to make plotting compatible in parallel mode (MPI)
    try:    
        print('Loading updated network template')
        filename_json = sim.cfg.dump_cell_properties
        with open(filename_json, 'r') as file: allCells = json.load(file)
    except: 
        print('Failed to load updated network template - Loading default template')
        filename_json = '../network_template/barreloid_network_cell_properties.json'
        with open(filename_json, 'r') as file: allCells = json.load(file)
    # allCells = ResampleParameters._loadCellData(sim)
    ########################################################################################################################
    
    allpops = ['VPM__pop','TRN__pop','TRN_ring__pop']
    try:
        record_pops = [(pop,list(np.arange(0,netParams.popParams[pop]['numCells']))) for pop in allpops]
    except:
        record_pops = [(pop,list(np.arange(0,40))) for pop in allpops]

    mean_voltage_dict={}
    for pop_ind, pop in enumerate(allpops):
        print('\n\n',pop)
        # sim.analysis.plotTraces(
        figs, traces_dict = sim.analysis.plotTraces(
                                include=[pop], 
                                # include=[record_pops[pop_ind]], 
                                # timeRange=[490,550], 
                                overlay=True, oneFigPer='trace', 
                                # ylim=[-110,50], 
                                axis=True, 
                                figSize=(70, 15), 
                                # figSize=(40, 15), 
                                # figSize=(60, 18), 
                                fontSize=15, 
                                # saveFig=True,
                                # saveFig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_traces__'+pop+'.png',
                                )
        
        tracesData = traces_dict['tracesData']
        # store_v={}
        store_v=[]
        store_voltages={}
        for rec_ind in range(len(tracesData)):
            for trace in tracesData[rec_ind].keys():
                if '_V_soma' in trace:
                    cell_gid_str = trace.split('_V_soma')[0].split('cell_')[1]
                    # store_v.update({cell_gid_str:list(tracesData[rec_ind][trace])})
                    store_v.append(list(tracesData[rec_ind][trace]))
                    store_voltages.update({cell_gid_str:list(tracesData[rec_ind][trace])})

        t_vector = list(tracesData[0]['t'])
        mean_v = np.mean(store_v, axis=0)
        t_vector_=[t_vector[i] for i in range(len(mean_v))]
        plt.figure(figsize=(70, 15))
        plt.rcParams.update({'font.size': 30})
        for trace in store_v: plt.plot(t_vector_,trace,'gray',alpha=0.2)
        plt.plot(t_vector_,mean_v,'r')
        plt.ylim([-110,50])
        # plt.xlim([min(t_vector_),max(t_vector_)])
        try:    plt.xlim([1500,          max(t_vector_)])
        except: plt.xlim([min(t_vector_),max(t_vector_)])
        ax=plt.gca()
        ax.spines[['right', 'top']].set_visible(False)
        # plt.plot(mean_v,'k')
        # plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_mean_traces__'+pop+'.png')
        plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_mean_traces__'+pop+'.png', dpi=300)
        
        mean_voltage_dict.update({pop:{'mean_v':list(mean_v), 't':list(t_vector_)}})

        plotColoredTraces=False
        if plotColoredTraces:
            try:
                from matplotlib import cm
                cm_subsection = np.linspace(0, 1, len(store_v)) 
                colors = [ cm.jet(x) for x in cm_subsection ]
                plt.figure(figsize=(70, 15))
                for trace_ind,trace in enumerate(store_v): plt.plot(t_vector_,trace,colors[trace_ind],alpha=0.2)
                plt.ylim([-110,50])
                plt.xlim([min(t_vector_),max(t_vector_)])
                # plt.plot(mean_v,'k')
                plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_colored_traces__'+pop+'.png')
            except: print('Plotting ordered colored lines failed')

        print(pop, '\tlast val = ', mean_v[-1])
        print(pop, '\tmean val = ', np.mean(mean_v[-1000:-1]))

        plot_v = store_v
        indexing = '_byGID'
        plotHeatMap=False
        # plotHeatMap=True
        if plotHeatMap:
            if pop!='TRN_ring__pop':
                try:
                    heatmap = MembranePotentialHeatmap(plot_v)
                    # You can specify the colormap and vmin/vmax values as needed
                    heatmap.plot_heatmap(cmap='CMRmap', vmin=-90, vmax=-30, savefig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_heatmap'+indexing+'__'+pop+'.png', figsize=(45, 20))
                except:
                    print('Plotting Heatmap failed')
                    pass

    with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_meanVoltage__data.json', 'w') as fp:
        json.dump(mean_voltage_dict, fp)

    if ('deflectionAngle|None' in sim.cfg.deflection_dataset_path) or ('freq_10_Hz__downTimes__None__' in sim.cfg.deflection_dataset_path):
        deflection_times = [
                            [ 1500+sim.cfg.delayBS,  sim.cfg.duration],
                        ]
    else:
        # 2024_08_07 - changed all deflection times to 250 ms (like original paper) to allow statistical analysis
        deflection_times = [
                            [ 2000+sim.cfg.delayBS,  2250+sim.cfg.delayBS],
                            [ 3000+sim.cfg.delayBS,  3250+sim.cfg.delayBS],
                            [ 4000+sim.cfg.delayBS,  4250+sim.cfg.delayBS],
                            [ 5000+sim.cfg.delayBS,  5250+sim.cfg.delayBS],
                            [ 6000+sim.cfg.delayBS,  6250+sim.cfg.delayBS],
                            [ 7000+sim.cfg.delayBS,  7250+sim.cfg.delayBS],
                            [ 8000+sim.cfg.delayBS,  8250+sim.cfg.delayBS],
                            [ 9000+sim.cfg.delayBS,  9250+sim.cfg.delayBS],
                        ]
        if cfg.extend_spike_times:
            import math
            skip_time = 2000
            deflection_times.extend([([def_on-skip_time+((add_loop+1)*10000),def_off-skip_time+((add_loop+1)*10000)]) for add_loop in range(math.ceil(sim.cfg.duration/10000)-1) for [def_on,def_off] in deflection_times if def_off-skip_time+((add_loop+1)*10000)<sim.cfg.duration])

    popColors_dict={}
    for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
                'MLe__pop', 
                'L6A_activated__pop']:
        if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
        elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
        elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
        elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
        elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
        elif 'MLe' in           pop:    popColors_dict.update({pop:'olive'})
        else:                           popColors_dict.update({pop:'gray'})
    sim.analysis.plotSpikeHist(
                                include=['allCells', 'eachPop'], 
                                timeRange=[1500,cfg.duration],
                                # timeRange=[0,cfg.duration],
                                figSize=[70,15],
                                # binSize=1,      # bin size from Minnery paper PSTH (PrV and VPM)
                                binSize=3, 
                                # binSize=5, 
                                graphType='line', 
                                measure='rate', 
                                norm=False, 
                                smooth=None, 
                                filtFreq=None, 
                                filtOrder=3, 
                                axis=True, 
                                popColors=popColors_dict,
                                saveFig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeRate.png',
    )
    sim.analysis.plotSpikeHist(
                                include=['allCells', 'eachPop'], 
                                timeRange=[1500,cfg.duration],
                                # timeRange=[0,cfg.duration],
                                figSize=[70,15],
                                # binSize=1,      # bin size from Minnery paper PSTH (PrV and VPM)
                                binSize=3, 
                                # binSize=5, 
                                graphType='line', 
                                measure='', 
                                norm=False, 
                                smooth=None, 
                                filtFreq=None, 
                                filtOrder=3, 
                                axis=True, 
                                popColors=popColors_dict,
                                saveFig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist.png',
    )

    save_spikes_dict={}
    if ('deflectionAngle|None' in sim.cfg.deflection_dataset_path) or ('freq_10_Hz__downTimes__None__' in sim.cfg.deflection_dataset_path):
        print('Skipping spike histogram for no deflection group')

        store_cell_tags={}
        for cell_gid in allCells.keys():
            store_cell_tags.update({str(cell_gid):{'pop':allCells[cell_gid]['tags']['pop'],
                                                   'x':  allCells[cell_gid]['tags']['x'],
                                                   'y':  allCells[cell_gid]['tags']['y'],
                                                   'z':  allCells[cell_gid]['tags']['z'],}})

        save_spikes_dict.update({'cell_tags':store_cell_tags})

        pop_list=[
                    # 'L6A_activated__pop'
                    # 'TRN_ring__pop', 
                    'TRN__pop', 'VPM__pop', 'MLe__pop', ]
        for ipop,pop in enumerate(pop_list):
            fName=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__pop__'+pop
            fDesc='spike_data'
            fType='json'
            tRange=[0,sim.cfg.duration]
            binSz = 1.0

            histData = sim.analysis.prepareSpikeHist(   
                                                                sim=sim,
                                                                include=[pop],
                                                                timeRange=tRange, 
                                                                popRates=True,
                                                                saveData=False,
                                                                fileName=fName,
                                                                fileDesc=fDesc,
                                                                fileType=fType,
                                                                fileDir=None,
                                                                binSize=binSz,
                                                                )
            try:
                spkTimes = histData['spkTimes']
                spkGids  = histData['spkGids']
                spkInds  = histData['spkInds']
            except: # in case there are no spikes in the population
                spkTimes = []
                spkGids  = []
                spkInds  = []

            save_spikes_dict.update({pop:{  'spkTimes':list(spkTimes), 
                                            'spkGids': list(spkGids), 
                                            'spkInds': list(spkInds),}})
        # Spike times (w/ cell_tags, spkTimes, spkGids, spkInds)
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__spikes.json', 'w') as fp:
            json.dump(save_spikes_dict, fp)

    else:
        
        store_cell_tags={}
        for cell_gid in allCells.keys():
            store_cell_tags.update({str(cell_gid):{'pop':allCells[cell_gid]['tags']['pop'],
                                                   'x':  allCells[cell_gid]['tags']['x'],
                                                   'y':  allCells[cell_gid]['tags']['y'],
                                                   'z':  allCells[cell_gid]['tags']['z'],}})

        save_spikes_dict.update({'cell_tags':store_cell_tags})
        pop_list=[
                    # 'L6A_activated__pop'
                    # 'TRN_ring__pop', 
                    'TRN__pop', 'VPM__pop', 'MLe__pop', ]
        for ipop,pop in enumerate(pop_list):
            fName=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__pop__'+pop
            fDesc='spike_data'
            fType='json'
            tRange=[0,sim.cfg.duration]
            binSz = 1.0

            histData = sim.analysis.prepareSpikeHist(   
                                                                sim=sim,
                                                                include=[pop],
                                                                timeRange=tRange, 
                                                                popRates=True,
                                                                saveData=False,
                                                                fileName=fName,
                                                                fileDesc=fDesc,
                                                                fileType=fType,
                                                                fileDir=None,
                                                                binSize=binSz,
                                                                )
            try:
                spkTimes = histData['spkTimes']
                spkGids  = histData['spkGids']
                spkInds  = histData['spkInds']
            except: # in case there are no spikes in the population
                spkTimes = []
                spkGids  = []
                spkInds  = []

            save_spikes_dict.update({pop:{  'spkTimes':list(spkTimes), 
                                            'spkGids': list(spkGids), 
                                            'spkInds': list(spkInds),}})
        # Spike times (w/ cell_tags, spkTimes, spkGids, spkInds)
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__spikes.json', 'w') as fp:
            json.dump(save_spikes_dict, fp)
        
        # --- Dictionary to save deflection data
        save_deflection_histogram={}
        # --- Spike histogram for each single deflection event
        for deflection_index,deflection_time in enumerate(deflection_times):
            
            # --- Dictionary to save deflection data
            save_deflection_histogram.update({str(deflection_index):{}})
            
            if deflection_time[1]>sim.cfg.duration:continue # --- doesnt attemp to plot if the deflection time is outside of the sim duration
            sim.analysis.plotSpikeHist(
                                        include=['allCells', 'eachPop'], 
                                        timeRange=[deflection_time[0]-100,deflection_time[1]+200],
                                        figSize=[20,20],
                                        # binSize=1,      # bin size from Minnery paper PSTH (PrV and VPM)
                                        binSize=3, 
                                        # binSize=5, 
                                        graphType='line', 
                                        measure='', 
                                        norm=False, 
                                        smooth=None, 
                                        filtFreq=None, 
                                        filtOrder=3, 
                                        axis=True, 
                                        popColors=popColors_dict,
                                        saveFig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__deflection_'+str(deflection_index)+'.png',
            )

            plt.figure(figsize=(20,10))
            figName=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__deflection_'+str(deflection_index)
            pop_list=[
                        # 'L6A_activated__pop'
                        # 'TRN_ring__pop', 
                        'TRN__pop', 'VPM__pop', 'MLe__pop', ]
            for ipop,pop in enumerate(pop_list):
                fName=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__deflection_'+str(deflection_index)+'__pop__'+pop
                fDesc='spike_data'
                fType='json'
                # tRange=[deflection_time[0]-99,deflection_time[1]+200]
                # binSz=3
                tRange=[deflection_time[0]-100,deflection_time[1]+200]
                # binSz=2.5
                # binSz=sim.cfg.recordStep
                binSz=1 # ms
                histData = sim.analysis.prepareSpikeHist(   
                                                            sim=sim,
                                                            include=[pop],
                                                            timeRange=tRange, 
                                                            popRates=True,
                                                            saveData=False,
                                                            fileName=fName,
                                                            fileDesc=fDesc,
                                                            fileType=fType,
                                                            fileDir=None,
                                                            binSize=binSz,
                                                            )

                try:
                    spkTimes = histData['spkTimes']
                    spkGids  = histData['spkGids']
                    spkInds  = histData['spkInds']
                except: # in case there are no spikes in the population
                    spkTimes = []
                    spkGids  = []
                    spkInds  = []

                histoData = np.histogram(spkTimes, bins=np.arange(tRange[0], tRange[1], binSz))

                plt.subplot(len(pop_list),1,ipop+1)

                plt.rcParams.update({'font.size': 30})
                hist_data_all = plt.hist(spkTimes,bins=histoData[1],color='k',histtype='step')
                deflection_bins_baseline = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[0]-100) and (timestamp< deflection_time[0])]
                deflection_bins_ON       = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[0])     and (timestamp< deflection_time[1])]
                deflection_bins_OFF      = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[1])     and (timestamp<=deflection_time[1]+100)]
                # deflection_bins_baseline = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[0]-100) and (timestamp<=deflection_time[0])]
                # deflection_bins_ON       = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[0])     and (timestamp<=deflection_time[1])]
                # deflection_bins_OFF      = [timestamp for timestamp in histoData[1] if (timestamp>=deflection_time[1])     and (timestamp<=deflection_time[1]+100)]
                hist_data_baseline       = plt.hist(spkTimes,bins=deflection_bins_baseline,  color='lightgrey')
                hist_data_ON             = plt.hist(spkTimes,bins=deflection_bins_ON,        color='k')
                hist_data_OFF            = plt.hist(spkTimes,bins=deflection_bins_OFF,       color='grey')
                # plt.hist(spkTimes,bins=[deflection_time[0],deflection_time[0]+binSz],color='r')
                plt.ylim([0,150])
                plt.yticks([0,50,100,150])
                xticks=[deflection_time[0]-100,deflection_time[0],deflection_time[1],deflection_time[1]+100,deflection_time[1]+200]
                xticks_labels=[xtick-deflection_time[0] for xtick in xticks]
                # xticks_labels=[xtick-min(xticks) for xtick in xticks]
                plt.xticks(xticks,labels=xticks_labels)
                ax=plt.gca()
                ax.spines[['right', 'top']].set_visible(False)
                print('Deflection Histogram plot worked - ' + str(deflection_index))

                # Parsing histogram data
                hist_data_all_bins      = list(hist_data_all[1])
                hist_data_all_vals      = list(hist_data_all[0])
                
                hist_data_baseline_bins = list(hist_data_baseline[1])
                hist_data_baseline_vals = list(hist_data_baseline[0])
                
                hist_data_ON_bins       = list(hist_data_ON[1])
                hist_data_ON_vals       = list(hist_data_ON[0])

                hist_data_OFF_bins      = list(hist_data_OFF[1])
                hist_data_OFF_vals      = list(hist_data_OFF[0])

                deflection_data_dict={  
                                        'all':      { 'bins':hist_data_all_bins,      'data':hist_data_all_vals      },
                                        'baseline': { 'bins':hist_data_baseline_bins, 'data':hist_data_baseline_vals },
                                        'ON':       { 'bins':hist_data_ON_bins,       'data':hist_data_ON_vals       },
                                        'OFF':      { 'bins':hist_data_OFF_bins,      'data':hist_data_OFF_vals      }}

                save_deflection_histogram[str(deflection_index)].update({pop:deflection_data_dict})

            plt.savefig(figName+'_hist.png',dpi=500)

        # Histogram data from spike times
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist__data.json', 'w') as fp:
            json.dump(save_deflection_histogram, fp)

        # --- Analysis of deflection events
        from AnalyzeDeflection import AnalyzeDeflection

        mean_voltage_dict_ = {pop:mean_voltage_dict[pop]['mean_v'] for pop in mean_voltage_dict.keys()}
        mean_voltage_dict_.update({'t':mean_voltage_dict['VPM__pop']['t']})
        try:
            AnalyzeDeflection.analyzeEvents(
                                                # deflection_histogram=deflection_data_dict, 
                                                deflection_histogram=save_deflection_histogram, 
                                                recordStep=sim.cfg.recordStep, 
                                                # dt=sim.cfg.dt, 
                                                grouped_mean_voltages=mean_voltage_dict_,
                                                savefig_name=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist_mean.png',
                                                savefig_dpi=500,
                                                )
        except:    
            AnalyzeDeflection.analyzeEvents(
                                                # deflection_histogram=deflection_data_dict, 
                                                deflection_histogram=save_deflection_histogram, 
                                                recordStep=sim.cfg.recordStep, 
                                                # dt=sim.cfg.dt, 
                                                savefig_name=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikeHist_mean.png',
                                                savefig_dpi=500,
                                                )

    # --- Plotting angluar tuning
    import json
    from Build_Net import BuildNetwork
    import math
    from ProcessAngles import AngleProcessor
    angle_processor = AngleProcessor()

    # deflection_indexes=[[int(t1/0.1),int(t2/0.1)] for [t1,t2] in deflection_times]
    deflection_indexes=[[int(t1/sim.cfg.recordStep),int(t2/sim.cfg.recordStep)] for [t1,t2] in deflection_times]
    #
    # sliced rasters
    select_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop', 'L6A_activated__pop']
    spkTimes_   = []
    spkInds_    = []
    spkGids_    = []
    spkTimes_byDeflection = {}
    spkInds_byDeflection  = {}
    spkGids_byDeflection  = {}
    deflection_flag=''

    for deflection_time_ind,deflection_time in enumerate(deflection_times):
        print(deflection_time_ind,'\t',deflection_time)
        # --- Plotting raster tuning
        rasterData = sim.analysis.prepareRaster(include=select_pops, timeRange=deflection_time)
        try:
            spkTimes_.extend(rasterData['spkTimes'])
            spkInds_.extend(rasterData['spkInds'])
            spkGids_.extend(rasterData['spkGids'])

            spkTimes_byDeflection.update({str(deflection_time_ind):rasterData['spkTimes']})
            spkInds_byDeflection.update( {str(deflection_time_ind):rasterData['spkInds']})
            spkGids_byDeflection.update( {str(deflection_time_ind):rasterData['spkGids']})

        except:continue

    spkTimes_byDeflection.update({'all':spkTimes_})
    spkInds_byDeflection.update( {'all':spkInds_})
    spkGids_byDeflection.update( {'all':spkGids_})

    # Saving _all angular tuning spikes
    try:
        angular_tuning_dict = {'spkTimes':spkTimes_byDeflection, 'spkInds':spkInds_byDeflection, 'spkGids':spkGids_byDeflection}
        import json
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3_spikes_all.json', 'w') as fp: json.dump(angular_tuning_dict, fp)
    except:
        print('Failed to save angular_tuning_dict')
        pass

    # --- Plot angular tuning
    for deflection_flag in spkTimes_byDeflection.keys():
        
        # creates a dictionary with spike times for each cell
        spkGids_int = [int(spkGid_) for spkGid_ in spkGids_byDeflection[deflection_flag]]
        from collections import Counter
        count_spks = Counter(spkGids_int)

        # dict with all gids for the selected pops and their spike counts
        store_count_spks={}
        store_spikes_grouped={}
        count_spikes_grouped={}
        for cell_ind in allCells.keys():
            cell_gid = allCells[str(cell_ind)]['gid']
            cell_pop = allCells[str(cell_ind)]['tags']['pop']
            if cell_pop not in select_pops:continue
            if cell_pop not in store_count_spks.keys():     store_count_spks.update({cell_pop:{}})
            if int(cell_gid) in count_spks.keys():          store_count_spks[cell_pop].update({int(cell_gid):count_spks[int(cell_gid)]})
            else:                                           store_count_spks[cell_pop].update({int(cell_gid):0})
            #
            if cell_pop not in store_spikes_grouped.keys(): 
                angles_dict={angle:[] for angle in angle_processor.reference_angles_}
                store_spikes_grouped.update({cell_pop:angles_dict}) 
                angles_dict_counter={angle:0 for angle in angle_processor.reference_angles_}
                count_spikes_grouped.update({cell_pop:angles_dict_counter}) 
            #
            if ('MLe' in cell_pop) or ('VPM' in cell_pop) or ('TRN' in cell_pop):
                # diam,height         = BuildNetwork.getNetworkSize(sim.cfg.center_point)
                diam,height         = BuildNetwork.getNetworkSize(500)
                # --- debug boundaries --- seem inverted
                cell_upper_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_lower_boundary = height[cell_pop.split('__pop')[0]][0]
                # # --- debug boundaries --- seem inverted
                # cell_upper_boundary = height[cell_pop.split('__pop')[0]][0]
                # cell_lower_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_position   = allCells[str(cell_ind)]['tags']['y']
            else: # needs to be tested, because L6 is deactivated for now
                cell_upper_boundary=360
                cell_lower_boundary=0
                cell_position_x   = allCells[str(cell_ind)]['tags']['x']
                cell_position_z   = allCells[str(cell_ind)]['tags']['z']
                cell_position=(np.remainder((((np.arctan2(cell_position_x-500,cell_position_z-500))*(180/np.pi))+360),360)/360)
            cell_position_relative = (cell_position-cell_lower_boundary)/(cell_upper_boundary-cell_lower_boundary)
            matching_angle = angle_processor.get_closest_angle(cell_position_relative)

            count_spikes_grouped[cell_pop][matching_angle]+=store_count_spks[cell_pop][cell_gid]

        # --- calculates the normalized values
        normalized_spiking_grouped  = {cell_pop:{} for cell_pop in count_spikes_grouped.keys()}
        pop_spiking                 = {cell_pop:0  for cell_pop in count_spikes_grouped.keys()}
        for cell_pop in count_spikes_grouped.keys():
            if cell_pop in select_pops:
                # Counts the number of spikes for a given pop
                pop_spiking[cell_pop]=sum(count_spikes_grouped[cell_pop].values())
                for matching_angle in count_spikes_grouped[cell_pop].keys():
                    print(matching_angle)
                    if max(count_spikes_grouped[cell_pop].values())<=0:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:0})
                    else:normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped[cell_pop][matching_angle]/max(count_spikes_grouped[cell_pop].values())})
        
        popColors_dict={}
        for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
                    'MLe__pop', 
                    'L6A_activated__pop']:
            if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
            elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
            elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
            elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
            elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
            elif 'MLe' in           pop:    popColors_dict.update({pop:'k'})
            else:                           popColors_dict.update({pop:'gray'})
        
        print('Plotting normalized and non-normalized angular tuning')

        try:
            # --- Try plotting individual angular tuning figures with visuals that match PlotAngularTuning.py            
            plt.figure(figsize=(20,50))
            plt.rcParams.update({'font.size': 80})
            plt.subplot(2,1,1,projection='polar')
            for cell_pop in normalized_spiking_grouped.keys():
                polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            
            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine

            plt.subplot(2,1,2,projection='polar')
            for cell_pop in count_spikes_grouped.keys():
                polar_list_values=list(count_spikes_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(count_spikes_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))
            plt.suptitle('Angular tuning plot - '+deflection_flag)


            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning2__'+deflection_flag+'.png')
        except:
            # --- Keep debugging angular tuning figure
            plt.figure(figsize=(20,40))
            plt.rcParams.update({'font.size': 40})
            plt.subplot(2,1,1,projection='polar')
            for cell_pop in normalized_spiking_grouped.keys():
                polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.subplot(2,1,2,projection='polar')
            for cell_pop in count_spikes_grouped.keys():
                polar_list_values=list(count_spikes_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(count_spikes_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))

            plt.suptitle('Angular tuning plot - '+deflection_flag)
            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning2__'+deflection_flag+'.png')

    ################################################################################################################################################
    # Store data from angular tuning plot of 20 ms of activity
    select_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop', 'L6A_activated__pop']
    spkTimes_20ms   = []
    spkInds_20ms    = []
    spkGids_20ms    = []
    spkTimes_byDeflection_20ms = {}
    spkInds_byDeflection_20ms  = {}
    spkGids_byDeflection_20ms  = {}

    count_spikes_barPlot={}

    for deflection_time_ind,deflection_time in enumerate(deflection_times):
        print(deflection_time_ind,'\t',deflection_time,'\t 20 ms analysis')
        # --- Slicing a 20 ms raster after ON deflection for new angular tuning plot
        rasterData2 = sim.analysis.prepareRaster(include=select_pops, timeRange=[deflection_time[0],deflection_time[0]+20])
        try:
            spkTimes_20ms.extend(rasterData2['spkTimes'])
            spkInds_20ms.extend(rasterData2['spkInds'])
            spkGids_20ms.extend(rasterData2['spkGids'])

            spkTimes_byDeflection_20ms.update({str(deflection_time_ind):rasterData2['spkTimes']})
            spkInds_byDeflection_20ms.update( {str(deflection_time_ind):rasterData2['spkInds']})
            spkGids_byDeflection_20ms.update( {str(deflection_time_ind):rasterData2['spkGids']})

        except:continue

    spkTimes_byDeflection_20ms.update({'all':spkTimes_20ms})
    spkInds_byDeflection_20ms.update( {'all':spkInds_20ms})
    spkGids_byDeflection_20ms.update( {'all':spkGids_20ms})

    # Saving 20 ms spikes
    try:
        angular_tuning_20ms_dict = {'spkTimes':spkTimes_byDeflection_20ms, 'spkInds':spkInds_byDeflection_20ms, 'spkGids':spkGids_byDeflection_20ms}
        import json
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3_spikes_ON.json', 'w') as fp: json.dump(angular_tuning_20ms_dict, fp)
    except:pass

    # --- Plot angular tuning
    for deflection_flag in spkTimes_byDeflection_20ms.keys():
        
        # creates a dictionary with spike times for each cell
        spkGids_int = [int(spkGid_) for spkGid_ in spkGids_byDeflection_20ms[deflection_flag]]
        from collections import Counter
        count_spks = Counter(spkGids_int)

        # dict with all gids for the selected pops and their spike counts
        store_count_spks_ON={}
        store_spikes_grouped_ON={}
        count_spikes_grouped_ON={}
        for cell_ind in allCells.keys():
            cell_gid = allCells[str(cell_ind)]['gid']
            cell_pop = allCells[str(cell_ind)]['tags']['pop']
            if cell_pop not in select_pops:continue
            if cell_pop not in store_count_spks_ON.keys():      store_count_spks_ON.update({cell_pop:{}})
            if int(cell_gid) in count_spks.keys():              store_count_spks_ON[cell_pop].update({int(cell_gid):count_spks[int(cell_gid)]})
            else:                                               store_count_spks_ON[cell_pop].update({int(cell_gid):0})
            #
            if cell_pop not in store_spikes_grouped_ON.keys(): 
                angles_dict={angle:[] for angle in angle_processor.reference_angles_}
                store_spikes_grouped_ON.update({cell_pop:angles_dict}) 
                angles_dict_counter={angle:0 for angle in angle_processor.reference_angles_}
                count_spikes_grouped_ON.update({cell_pop:angles_dict_counter}) 
            #
            if ('MLe' in cell_pop) or ('VPM' in cell_pop) or ('TRN' in cell_pop):
                # diam,height         = BuildNetwork.getNetworkSize(sim.cfg.center_point)
                diam,height         = BuildNetwork.getNetworkSize(500)
                # --- debug boundaries --- seem inverted
                cell_upper_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_lower_boundary = height[cell_pop.split('__pop')[0]][0]
                # # --- debug boundaries --- seem inverted
                # cell_upper_boundary = height[cell_pop.split('__pop')[0]][0]
                # cell_lower_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_position   = allCells[str(cell_ind)]['tags']['y']
            else: # needs to be tested, because L6 is deactivated for now
                cell_upper_boundary=360
                cell_lower_boundary=0
                cell_position_x   = allCells[str(cell_ind)]['tags']['x']
                cell_position_z   = allCells[str(cell_ind)]['tags']['z']
                cell_position=(np.remainder((((np.arctan2(cell_position_x-500,cell_position_z-500))*(180/np.pi))+360),360)/360)
            cell_position_relative = (cell_position-cell_lower_boundary)/(cell_upper_boundary-cell_lower_boundary)
            matching_angle = angle_processor.get_closest_angle(cell_position_relative)

            count_spikes_grouped_ON[cell_pop][matching_angle]+=store_count_spks_ON[cell_pop][cell_gid]

        # --- calculates the normalized values
        normalized_spiking_grouped  = {cell_pop:{} for cell_pop in count_spikes_grouped_ON.keys()}
        pop_spiking                 = {cell_pop:0  for cell_pop in count_spikes_grouped_ON.keys()}
        for cell_pop in count_spikes_grouped_ON.keys():
            if cell_pop in select_pops:
                # Counts the number of spikes for a given pop
                pop_spiking[cell_pop]=sum(count_spikes_grouped_ON[cell_pop].values())
                for matching_angle in count_spikes_grouped_ON[cell_pop].keys():
                    print(matching_angle)
                    if max(count_spikes_grouped_ON[cell_pop].values())<=0:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:0})
                    else:normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped_ON[cell_pop][matching_angle]/max(count_spikes_grouped_ON[cell_pop].values())})
        
        popColors_dict={}
        for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
                    'MLe__pop', 
                    'L6A_activated__pop']:
            if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
            elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
            elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
            elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
            elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
            elif 'MLe' in           pop:    popColors_dict.update({pop:'k'})
            else:                           popColors_dict.update({pop:'gray'})
        
        print('Plotting normalized and non-normalized angular tuning')
        try:
            # --- Try plotting individual angular tuning figures with visuals that match PlotAngularTuning.py            
            plt.figure(figsize=(20,50))
            plt.rcParams.update({'font.size': 80})
            plt.subplot(2,1,1,projection='polar')
            for cell_pop in normalized_spiking_grouped.keys():
                polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            
            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine

            plt.subplot(2,1,2,projection='polar')
            for cell_pop in count_spikes_grouped_ON.keys():
                polar_list_values=list(count_spikes_grouped_ON[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(count_spikes_grouped_ON[cell_pop].keys())
                polar_list_keys.extend([360])
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))
            plt.suptitle('Angular tuning plot - '+deflection_flag)
            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3__'+deflection_flag+'_0_ON.png')
        except:
            # --- Keep debugging angular tuning figure
            plt.figure(figsize=(20,40))
            plt.rcParams.update({'font.size': 40})
            plt.subplot(2,1,1,projection='polar')
            for cell_pop in normalized_spiking_grouped.keys():
                polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
                polar_list_keys.extend([360])
                # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.subplot(2,1,2,projection='polar')
            for cell_pop in count_spikes_grouped_ON.keys():
                polar_list_values=list(count_spikes_grouped_ON[cell_pop].values())
                polar_list_values.extend([polar_list_values[0]])
                polar_list_keys=list(count_spikes_grouped_ON[cell_pop].keys())
                polar_list_keys.extend([360])
                try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
                except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
            plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))
            plt.suptitle('Angular tuning plot - '+deflection_flag)
            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3__'+deflection_flag+'_0_ON.png')

        if deflection_flag == 'all': count_spikes_barPlot.update({'ON':count_spikes_grouped_ON})

    ################################################################################################################################################
    # Store data from angular tuning plot of 20 ms of activity
    select_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop', 'L6A_activated__pop']
    spkTimes_20ms_OFF   = []
    spkInds_20ms_OFF    = []
    spkGids_20ms_OFF    = []
    spkTimes_byDeflection_20ms_OFF = {}
    spkInds_byDeflection_20ms_OFF  = {}
    spkGids_byDeflection_20ms_OFF  = {}

    for deflection_time_ind,deflection_time in enumerate(deflection_times):
        print(deflection_time_ind,'\t',deflection_time,'\t 20 ms analysis')
        # --- Slicing a 20 ms raster after ON deflection for new angular tuning plot
        rasterData3 = sim.analysis.prepareRaster(include=select_pops, timeRange=[deflection_time[1],deflection_time[1]+20])
        try:
            spkTimes_20ms_OFF.extend(rasterData3['spkTimes'])
            spkInds_20ms_OFF.extend(rasterData3['spkInds'])
            spkGids_20ms_OFF.extend(rasterData3['spkGids'])

            spkTimes_byDeflection_20ms_OFF.update({str(deflection_time_ind):rasterData3['spkTimes']})
            spkInds_byDeflection_20ms_OFF.update( {str(deflection_time_ind):rasterData3['spkInds']})
            spkGids_byDeflection_20ms_OFF.update( {str(deflection_time_ind):rasterData3['spkGids']})

        except:continue

    spkTimes_byDeflection_20ms_OFF.update({'all':spkTimes_20ms_OFF})
    spkInds_byDeflection_20ms_OFF.update( {'all':spkInds_20ms_OFF})
    spkGids_byDeflection_20ms_OFF.update( {'all':spkGids_20ms_OFF})

    # Saving 20 ms spikes
    try:
        angular_tuning_20ms_OFF_dict = {'spkTimes':spkTimes_byDeflection_20ms_OFF, 'spkInds':spkInds_byDeflection_20ms_OFF, 'spkGids':spkGids_byDeflection_20ms_OFF}
        import json
        with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3_spikes_OFF.json', 'w') as fp: json.dump(angular_tuning_20ms_OFF_dict, fp)
    except:pass

    # --- Plot angular tuning
    for deflection_flag in spkTimes_byDeflection_20ms_OFF.keys():
        
        # creates a dictionary with spike times for each cell
        spkGids_int = [int(spkGid_) for spkGid_ in spkGids_byDeflection_20ms_OFF[deflection_flag]]
        from collections import Counter
        count_spks = Counter(spkGids_int)

        # dict with all gids for the selected pops and their spike counts
        store_count_spks_OFF={}
        store_spikes_grouped_OFF={}
        count_spikes_grouped_OFF={}
        for cell_ind in allCells.keys():
            cell_gid = allCells[str(cell_ind)]['gid']
            cell_pop = allCells[str(cell_ind)]['tags']['pop']
            if cell_pop not in select_pops:continue
            if cell_pop not in store_count_spks_OFF.keys():     store_count_spks_OFF.update({cell_pop:{}})
            if int(cell_gid) in count_spks.keys():              store_count_spks_OFF[cell_pop].update({int(cell_gid):count_spks[int(cell_gid)]})
            else:                                               store_count_spks_OFF[cell_pop].update({int(cell_gid):0})
            #
            if cell_pop not in store_spikes_grouped_OFF.keys(): 
                angles_dict={angle:[] for angle in angle_processor.reference_angles_}
                store_spikes_grouped_OFF.update({cell_pop:angles_dict}) 
                angles_dict_counter={angle:0 for angle in angle_processor.reference_angles_}
                count_spikes_grouped_OFF.update({cell_pop:angles_dict_counter}) 
            #
            if ('MLe' in cell_pop) or ('VPM' in cell_pop) or ('TRN' in cell_pop):
                # diam,height         = BuildNetwork.getNetworkSize(sim.cfg.center_point)
                diam,height         = BuildNetwork.getNetworkSize(500)
                # --- debug boundaries --- seem inverted
                cell_upper_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_lower_boundary = height[cell_pop.split('__pop')[0]][0]
                # # --- debug boundaries --- seem inverted
                # cell_upper_boundary = height[cell_pop.split('__pop')[0]][0]
                # cell_lower_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_position   = allCells[str(cell_ind)]['tags']['y']
            else: # needs to be tested, because L6 is deactivated for now
                cell_upper_boundary=360
                cell_lower_boundary=0
                cell_position_x   = allCells[str(cell_ind)]['tags']['x']
                cell_position_z   = allCells[str(cell_ind)]['tags']['z']
                cell_position=(np.remainder((((np.arctan2(cell_position_x-500,cell_position_z-500))*(180/np.pi))+360),360)/360)
            cell_position_relative = (cell_position-cell_lower_boundary)/(cell_upper_boundary-cell_lower_boundary)
            matching_angle = angle_processor.get_closest_angle(cell_position_relative)

            count_spikes_grouped_OFF[cell_pop][matching_angle]+=store_count_spks_OFF[cell_pop][cell_gid]

        # --- calculates the normalized values
        normalized_spiking_grouped  = {cell_pop:{} for cell_pop in count_spikes_grouped_OFF.keys()}
        pop_spiking                 = {cell_pop:0  for cell_pop in count_spikes_grouped_OFF.keys()}
        for cell_pop in count_spikes_grouped_OFF.keys():
            if cell_pop in select_pops:
                # Counts the number of spikes for a given pop
                pop_spiking[cell_pop]=sum(count_spikes_grouped_OFF[cell_pop].values())
                for matching_angle in count_spikes_grouped_OFF[cell_pop].keys():
                    print(matching_angle)
                    if max(count_spikes_grouped_OFF[cell_pop].values())<=0:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:0})
                    else:normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped_OFF[cell_pop][matching_angle]/max(count_spikes_grouped_OFF[cell_pop].values())})
        
        popColors_dict={}
        for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
                    'MLe__pop', 
                    'L6A_activated__pop']:
            if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
            elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
            elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
            elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
            elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
            elif 'MLe' in           pop:    popColors_dict.update({pop:'k'})
            else:                           popColors_dict.update({pop:'gray'})
        
        print('Plotting normalized and non-normalized angular tuning')
        # --- Keep debugging angular tuning figure
        plt.figure(figsize=(20,40))
        plt.rcParams.update({'font.size': 40})
        plt.subplot(2,1,1,projection='polar')
        for cell_pop in normalized_spiking_grouped.keys():
            polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
            polar_list_values.extend([polar_list_values[0]])
            polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
            polar_list_keys.extend([360])
            # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
        plt.subplot(2,1,2,projection='polar')
        for cell_pop in count_spikes_grouped_OFF.keys():
            polar_list_values=list(count_spikes_grouped_OFF[cell_pop].values())
            polar_list_values.extend([polar_list_values[0]])
            polar_list_keys=list(count_spikes_grouped_OFF[cell_pop].keys())
            polar_list_keys.extend([360])
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
        plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))

        plt.suptitle('Angular tuning plot - '+deflection_flag)
        plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3__'+deflection_flag+'_1_OFF.png')

        if deflection_flag == 'all': count_spikes_barPlot.update({'OFF':count_spikes_grouped_OFF})

    ################################################################################################################################################
    try:
        for deflection_event_ind, deflection_event in enumerate(count_spikes_barPlot.keys()):
            print('Plotting firing in each pop - ', str(deflection_event))
            plt.figure(figsize=(20,20))
            plt.rcParams.update({'font.size': 40})
            for pop_ind,pop in enumerate(count_spikes_barPlot[deflection_event].keys()):
                plt.bar(pop_ind,sum(count_spikes_barPlot[deflection_event][pop].values()))
            plt.legend(list(count_spikes_barPlot[deflection_event].keys()))
            
            MLe_VPM_ratio=(sum(count_spikes_barPlot[deflection_event]['MLe__pop'].values())/sum(count_spikes_barPlot[deflection_event]['VPM__pop'].values()))
            plt.suptitle('Mean spike count - '+str(MLe_VPM_ratio))
            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_angularTuning3__'+deflection_event+'_3_bar.png')
    except:
        pass
    
    ################################################################################################################################################
    # plotLFP=True
    if sim.cfg.plotLFP:
        print('Preparing LFP timeseries manually to inspect bugs')
        LFP_timeRange=[1500,sim.cfg.duration]
        #             ]
        LFP_freqs=[ [1,     400],
                    [1,     50],
                    ]
        LFP_dict={}
        for LFP_ind,LFP_freq in enumerate(LFP_freqs):
            LFPData = sim.analysis.prepareLFP(
                    sim=sim,
                    timeRange=LFP_timeRange,
                    electrodes=['avg', 'all'],
                    # pop=pop,
                    LFPData=None,
                    normSignal=False,
                    filtFreq=LFP_freq,
                    # filtFreq=100,
                    # filtFreq=False,
                    filtOrder=3,
                    detrend=False,
                )
            print('LFPData:')
            print(LFPData)
            try:
                list_data = [list(data_) for data_ in LFPData['electrodes']['data']]
                LFP_dict.update({str(LFP_ind):{'data':list_data,'t':list(LFPData['t']),'names':list(LFPData['electrodes']['names'])}})
            except:
                print('Failed to calculate LFP data for frequencies ', LFP_freq)
                pass
        try:
            import json
            with open(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_LFP_data.json', 'w') as fp:
                json.dump(LFP_dict, fp)
        except:pass

        try:
            LFP_recordings=LFP_dict['0']['names']
            print('LFP_recordings worked')
        except:
            LFP_recordings=['avg', '0', '1']
            print('LFP_recordings failed')
        
        for LFP_recording_ind,LFP_recording in enumerate(LFP_recordings):
            plt.figure(figsize=(70,40))
            plt.rcParams.update({'font.size': 30})
            for LFP_freq_ind,LFP_freq_key in enumerate(LFP_dict.keys()):
                LFP_t     = list(LFP_dict[LFP_freq_key]['t'])           
                # LFP_trace = list(LFP_dict[LFP_freq_ind]['electrodes']['data'][LFP_recording_ind])    
                LFP_trace = list(LFP_dict[LFP_freq_key]['data'][LFP_recording_ind])    
                plt.subplot(len(list(LFP_dict.keys())),1,LFP_freq_ind+1)
                plt.plot(LFP_t,LFP_trace,'k',linewidth=5)
                if   LFP_recordings[LFP_recording_ind]=='0':
                    plt.ylim(-0.1,0.025)
                elif LFP_recordings[LFP_recording_ind]=='1':
                    plt.ylim(-0.025,0.1)
                else:    
                    plt.ylim(-0.125,0.125)
                try:    
                    plt.xlim(1500,sim.cfg.duration)
                    # plt.xlim(2000,sim.cfg.duration)
                except: pass
                ax=plt.gca()
                ax.spines[['right', 'top']].set_visible(False)

            plt.suptitle(LFP_recordings[LFP_recording_ind])
            plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_LFP_timeSeries_'+LFP_recordings[LFP_recording_ind]+'.png')

    if sim.cfg.plotLFP:
        print('Plotting Spectrogram')
        for pop in ['VPM__pop', 'TRN__pop']:
            print('Plotting Spectrogram - ', pop)
            
            pop_electrode = {'VPM__pop':[1], 'TRN__pop':[0]}
            sim.plotting.plotLFPSpectrogram(electrodes=pop_electrode[pop])
    
    # --- Analysis of ionic currents
    if sim.cfg.recordCurrents:
        from AnalyzeCurrents import AnalyzeCurrents
        timeRanges=[[10500,11500]]
        buffer_time=250 # ms
        
        deflection_times_withBeforeAfterTimes=[[def_on-buffer_time,def_off+buffer_time] for [def_on,def_off] in deflection_times]
        if ('deflectionAngle|None' not in sim.cfg.deflection_dataset_path): timeRanges.extend(deflection_times_withBeforeAfterTimes)

        for timeRange in timeRanges:
                
            print('\t---> Running current analysis on all currents')
            AnalyzeCurrents.runAnalysis(sim, 
                                        current_flag        = 'i__', 
                                        select_pops         = ['VPM__pop','TRN__pop'],
                                        select_currents     = None,
                                        time_range          = timeRange,
                                        singleTraces        = False, 
                                        allTraces           = False, 
                                        separateFigs        = True, 
                                        figSize             = (70,15), 
                                        figAlpha            = 0.5, 
                                        savefig_name        = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_currAnalysis', 
                                        # savefig_name        = 'pc_cell_curr', 
                                        savefig_format      = '.png')
            print('\t---> Running current analysis on selected currents - no HH')
            AnalyzeCurrents.runAnalysis(sim, 
                                        current_flag        = 'i__', 
                                        select_pops         = ['VPM__pop','TRN__pop'],
                                        select_currents = [ 'i__soma_0__SK_E2__ik', 
                                                            # 'i__soma_0__TC_HH__ina', 
                                                            # 'i__soma_0__TC_HH__ik', 
                                                            'i__soma_0__TC_Nap_Et2__ina', 
                                                            'i__soma_0__TC_iA__ik', 
                                                            'i__soma_0__TC_iL__ica', 
                                                            'i__soma_0__TC_iT_Des98__ica', 
                                                            'i__soma_0__TC_ih_Bud97__ih'],
                                        time_range          = timeRange,
                                        singleTraces        = False, 
                                        allTraces           = False, 
                                        separateFigs        = True, 
                                        figSize             = (70,15), 
                                        figAlpha            = 0.5, 
                                        savefig_name        = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_currAnalysis_noHH', 
                                        # savefig_name        = 'pc_cell_curr', 
                                        savefig_format      = '.png')
    '''
    All currents:
        'i__soma_0__SK_E2__ik'
        'i__soma_0__TC_HH__ina'
        'i__soma_0__TC_HH__ik'
        'i__soma_0__TC_Nap_Et2__ina'
        'i__soma_0__TC_iA__ik'
        'i__soma_0__TC_iL__ica'
        'i__soma_0__TC_iT_Des98__ica'
        'i__soma_0__TC_ih_Bud97__ih'
    '''
    if sim.cfg.plotLFP:
        try:
            # Saving LFP data

            # Full LFP data 
            try:
                LFP_data = sim.simData['LFP']
                cfg      = sim.cfg
            except:
                LFP_data = sim['simData']['LFP']
                cfg      = sim['simConfig']

            # Dictionary to store LFP data by electrode
            LFP_data_dict={'electrodes':{},'cfg':cfg}
            for i in range(len(LFP_data[0])):
                LFP_data_dict['electrodes'].update({i:[]})

            # Parsing dictionary data points into each electrode
            for j in range(len(LFP_data)):
                for electrode in LFP_data_dict['electrodes'].keys():
                    LFP_data_dict['electrodes'][electrode].append(LFP_data[j][electrode])

            # Adding a time vector for plotting
            LFP_data_dict.update({'t':[k*cfg['recordStep'] for k in range(len(LFP_data))]})

            save_lfp =  sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_lfp_values.json'
            with open(save_lfp, 'w') as file:
                json.dump(LFP_data_dict, file)

        except:pass

        print('--- Processing Rate PSD')
        try:
            data_ratePSD_dict={}
            for pop in ['MLe__pop','VPM__pop', 'TRN__pop']:
                print('Plotting custom Spike Rate PSD')
                # Function returns=> fig, {'allSignal': allSignal, 'allFreqs': allFreqs}
                fig, data_ratePSD_ = sim.analysis.plotRatePSD(   include         = [pop], 
                                            timeRange       = [1500,sim.cfg.duration], 
                                            binSize         = 1, 
                                            minFreq         = 1, 
                                            maxFreq         = 100, 
                                            transformMethod = 'morlet', 
                                            stepFreq        = 1, 
                                            NFFT            = 256, 
                                            noverlap        = 128, 
                                            smooth          = 0, 
                                            norm            = False, 
                                            overlay         = True, 
                                            # popColors       = popColors_dict[pop], 
                                            yLogScale       = True, 
                                            # ylim            = [0,1], 
                                            figSize         = (20, 20), 
                                            fontSize        = 25, 
                                            lineWidth       = 8, 
                                            saveData        = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikePSDData_'+pop+'_2.json', 
                                            saveFig         = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spikePSDPlot_'+pop+'_2.png', 
                                            showFig         = False,
                                            )
                data_ratePSD_dict.update({pop:{'signal':data_ratePSD_['allSignal'][0].tolist(),'freqs':data_ratePSD_['allFreqs'][0].tolist(),}})
                # data_ratePSD_dict.update({pop:data_ratePSD_})
            save_ratePSD =  sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_ratePSD_values.json'
            with open(save_ratePSD, 'w') as file:
                json.dump(data_ratePSD_dict, file)
        except: pass

        print('--- Processing LFP PSD')
        try:
            pop_electrode = {'VPM__pop':1, 'TRN__pop':0}
            lfp_pops = ['VPM__pop', 'TRN__pop']
            data_lfpPSD_dict={}
            data_lfpPSD_dict_windows={}
            data_spectrogram_dict_windows={}
            for ipop,pop in enumerate(lfp_pops):
                sim.plotting.plotTimeSeriesPSD.plotLFPPSD(  
                                                            PSDData=None, 
                                                            axis=None, timeRange=[1500,sim.cfg.duration], electrodes=[pop_electrode[pop]], 
                                                            pop=None, separation=1.0, roundOffset=True, 
                                                            NFFT=256, noverlap=128, nperseg=256, minFreq=1, maxFreq=100, stepFreq=1, smooth=0, 
                                                            logy=False, normSignal=False, normPSD=False, filtFreq=False, filtOrder=3, detrend=False, 
                                                            transformMethod='morlet', orderInverse=True, legend=True, colorList=None, returnPlotter=False,
                                                            saveFig = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_tracesPSDPlot_'+pop+'.png',
                                                            )

                NFFT        = 256
                noverlap    = 128
                nperseg     = 256
                minFreq     = 1
                maxFreq     = 100
                stepFreq    = 1
                smooth      = 0 
                PSDData_ = sim.analysis.preparePSD(
                    PSDData=None,
                    sim=sim,
                    timeRange=[1500,sim.cfg.duration],
                    electrodes=[pop_electrode[pop]],
                    pop=pop,
                    NFFT=NFFT,
                    noverlap=noverlap,
                    nperseg=nperseg,
                    minFreq=minFreq,
                    maxFreq=maxFreq,
                    stepFreq=stepFreq,
                    smooth=smooth,
                    logy=False,
                    normSignal=False,
                    normPSD=False,
                    filtFreq=False,
                    filtOrder=3,
                    detrend=False,
                    transformMethod='morlet',
                )
                data_lfpPSD_dict.update({pop:{'signal':PSDData_['psdSignal'][0].tolist(),'freqs':PSDData_['psdFreqs'][0].tolist(),'electrode':PSDData_['psdNames'][0]}})

                # processing windowed values to calculate PSD over time
                data_lfpPSD_dict_windows.update({pop:{}})
                LFP_PSD_window_size=250 # ms
                window_step = 50 # ms
                LFP_PSD_windows = [[time_val,time_val+LFP_PSD_window_size] for time_val in np.arange(1500,sim.cfg.duration-LFP_PSD_window_size+1,window_step)]
                for [time_val_0,time_val_1] in LFP_PSD_windows:
                    PSDData_window = sim.analysis.preparePSD(
                        PSDData=None,
                        sim=sim,
                        timeRange=[time_val_0,time_val_1],
                        electrodes=[pop_electrode[pop]],
                        pop=pop,
                        NFFT=NFFT,
                        noverlap=noverlap,
                        nperseg=nperseg,
                        minFreq=minFreq,
                        maxFreq=maxFreq,
                        stepFreq=stepFreq,
                        smooth=smooth,
                        logy=False,
                        normSignal=False,
                        normPSD=False,
                        filtFreq=False,
                        filtOrder=3,
                        detrend=False,
                        transformMethod='morlet',
                    )
                    data_lfpPSD_dict_windows[pop].update({str(int(time_val_0))+'__'+str(int(time_val_1)):{'signal':PSDData_window['psdSignal'][0].tolist(),'freqs':PSDData_window['psdFreqs'][0].tolist(),'electrode':PSDData_window['psdNames'][0]}})
                
                try:
                    if sim.cfg.record_windowed_spectrogram:
                        # --- up/down states together
                        data_spectrogram_dict_windows.update({pop:{}})
                        skip_up_time    = 0
                        skip_down_time  = 200
                        ct_times=[      [ 2500, 3000, 'u'],  [ 4500, 5000, 'u'],  [ 6500, 7000, 'u'],  [ 8500, 9000, 'u'],  [10500,11000, 'u'],  [12500,13000, 'u'],  [14500,15000, 'u'],
                                        [ 3000, 4500, 'd'],  [ 5000, 6500, 'd'],  [ 7000, 8500, 'd'],  [ 9000,10500, 'd'],  [11000,12500, 'd'],  [13000,14500, 'd'],  [15000,16000, 'd'],
                                        [ 1500, 16000,'a']
                                    ]
                        ct_times_=[]
                        for [state_start,state_end,state] in ct_times:
                            if   state=='u':    ct_times_.append([state_start+skip_up_time,   state_end,                state])
                            elif state=='d':    ct_times_.append([state_start+skip_down_time, state_end-skip_down_time, state])
                            else:               ct_times_.append([state_start+skip_up_time,   state_end,                state])
                        ct_times=ct_times_

                        window_number_States = len(ct_times)
                        for window_ind, [time_val_0,time_val_1,state] in enumerate(ct_times):
                            print('Processing window # ', window_ind, ' of ', window_number_States, ' for pop ', pop)
                            spectrogramData_window = sim.analysis.prepareSpectrogram(
                                sim=sim,
                                timeRange=[time_val_0,time_val_1],
                                electrodes=[pop_electrode[pop]],
                                pop=pop,
                                NFFT=NFFT,
                                noverlap=noverlap,
                                nperseg=nperseg,
                                minFreq=minFreq,
                                maxFreq=maxFreq,
                                stepFreq=stepFreq,
                                smooth=smooth,
                                logy=False,
                                normSignal=False,
                                normPSD=False,
                                normSpec=False,
                                filtFreq=False,
                                filtOrder=3,
                                detrend=False,
                                transformMethod='morlet',
                            )
                            data_spectrogram_dict_windows[pop].update({state+'__'+str(int(time_val_0))+'__'+str(int(time_val_1)):spectrogramData_window})
                    else:
                        print('Skipping saving windowed spectrogram dataset')
                except:
                    print('Failed to save windowed spectrogram')
                    pass

            save_lfpPSD =  sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_lfpPSD_values.json'
            with open(save_lfpPSD, 'w') as file:
                json.dump(data_lfpPSD_dict, file)

            # save windowed values to calculate PSD over time
            save_lfpPSD_windows =  sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_lfpPSD_values_windows.json'
            with open(save_lfpPSD_windows, 'w') as file:
                json.dump(data_lfpPSD_dict_windows, file)
            
            print('Saving windowed spectrogram data')
            # --- saving to .pkl because numpy arrays are not serializable and hard to work with, and too large
            save_spectrogram_windows_pkl =  sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_spectrogram_values_windows.pkl'
            # Open the file in binary write mode ('wb')
            with open(save_spectrogram_windows_pkl, 'wb') as file:
                # Use pickle.dump() to serialize and save the array
                pickle.dump(data_spectrogram_dict_windows, file)

        except:pass
    
    try:
        print('--- Processing Rate PSD')
        print(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_rasterCustom.png')
        sim.analysis.plotRaster(
            figSize         = [45, 25],
            include         = ['MLe__pop','VPM__pop', 'TRN__pop', 'CTvirtual_uniform__pop'],
            timeRange       = [1500,sim.cfg.duration],
            orderInverse    = True,
            orderBy         = 'y',
            spikeHist       = True,
            popColors       = popColors_dict,
            saveFig         = sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_rasterCustom_0.png',
            dpi             = 300,
            # tight_layout     = True,
            # tightLayout     = True,
        )

        import matplotlib.pyplot as plt
        fig, ax = plt.gcf(), plt.gca()
        fig.set_tight_layout(False)  # Ensure NetPyNE doesn't override this
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)  # Adjust margins manually
        plt.savefig(sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_rasterCustom_1.png')
        # plt.show()
    except:
        print('Failed to plot custom raster')
        pass

# --- Plotting 3D network shape
if plotNetworkShape:
    for cell_ind, cell in enumerate(sim.net.cells):
        try:
            if '_repr_mimebundle_' in sim.net.cells[cell_ind].secs.keys():
                del sim.net.cells[cell_ind].secs['_repr_mimebundle_']
        except:pass
    # Using the same default cmap as sim.analysis.plotShape
    cmap = plt.cm.hsv
    # cmap = plt.cm.viridis
    # Storing Y positions of each cell
    pop_yRanges = {pop:sim.net.pops[pop].tags['yRange'] for pop in sim.net.pops.keys()}

    pop_gids = {pop:sim.net.pops[pop].cellGids for pop in sim.net.pops.keys()}

    # density = 0.2
    density = 1.0
    sampled_pop_gids = {}
    for pop in pop_gids.keys():
        sampled_gids=[pop_gids[pop][gid_ind] for gid_ind in range(0,len(pop_gids[pop]),int(1/density))]
        sampled_pop_gids.update({pop:sampled_gids})

    # skips Vecstim pops
    morph_pops = [pop for pop in sim.net.pops.keys() if 'cellModel' not in sim.net.pops[pop].tags.keys()]
    
    includePre =morph_pops
    # includePost=morph_pops
    includePost=[(pop,sampled_pop_gids[pop]) for pop in morph_pops]
    includePost_gids = []
    for pop in sampled_pop_gids.keys():
        if pop in morph_pops: includePost_gids.extend(sampled_pop_gids[pop])
    
    cellsPreGids = [c.gid for c in sim.getCellsList(includePre)]  if includePre  else []
    cellsPost = sim.getCellsList(includePost_gids)

    secs=None
    includeAxon=True
    if not secs:
        secs = [s['hObj'] for cellPost in cellsPost for s in list(cellPost.secs.values())]
    if not includeAxon:
        secs = [sec for sec in secs if 'axon' not in sec.hname()]

    alpha=0.8
    # alpha_grey=0.2
    # grey_cells=5
    
    alpha_grey=0.1

    # grey_cells=15 # 1 out of every <number of grey_cells> cells will be colored
    # grey_gids = [gid for igid,gid in enumerate(includePost_gids) if (igid%grey_cells)!=0]
    
    grey_cells_VPM = 15 # 1 out of every <number of grey_cells> cells will be colored
    grey_cells_TRN = 10 # 1 out of every <number of grey_cells> cells will be colored

    grey_gids=[]
    for igid,gid in enumerate(includePost_gids):
        if gid in pop_gids['VPM__pop']:
            if (igid%grey_cells_VPM)!=0: grey_gids.append(gid)
        elif (gid in pop_gids['TRN__pop']) or (gid in pop_gids['TRN_ring__pop']):
            if (igid%grey_cells_TRN)!=0: grey_gids.append(gid)    

    cvals=[]
    color_strategy='relative_position'
    # color_strategy='center_distance'
    # color_strategy='random'
    # relative_normY = []

    if color_strategy=='random': 
        from neuron import h
        random_color = h.Random(100)
        cmap = plt.cm.turbo
    
    black_soma=False

    r_vals={}
    for sec in secs:
        gid = sec.cell().gid
        cell_pop=sim.net.cells[gid].tags['pop']
        print(cell_pop,gid,sec.name(),sec.nseg)

        if (color_strategy=='random') and (gid not in r_vals.keys()): r_vals.update({gid:random_color.uniform(0,1)})

        for seg in range(sec.nseg):
            cmap.set_bad(color='black')
            if sec.cell().gid in grey_gids:
                if ('soma' in sec.name()) and black_soma:
                    cvals.append((0,0,0,alpha_grey))
                else:
                    # cvals.append((0.5,0.5,0.5,alpha_grey))                                                # grey
                    cvals.append((0.9607843137254902, 0.9607843137254902, 0.9607843137254902,alpha_grey))   # whitesmoke (a bit dimmer than grey, to use with black background)
            else:
                if ('soma' in sec.name()) and black_soma:
                    cvals.append((0,0,0,alpha))
                else:
                    # relative_normY.append((sim.net.cells[gid].tags['y']-pop_yRanges[cell_pop][0])/(pop_yRanges[cell_pop][1]-pop_yRanges[cell_pop][0]))
                    if color_strategy == 'center_distance':
                        # --- Distance from pop center
                        pop_center  = (pop_yRanges[cell_pop][1]+pop_yRanges[cell_pop][0])/2
                        max_dist    = (pop_yRanges[cell_pop][1]-pop_yRanges[cell_pop][0])/2
                        relative_center_distance = abs(sim.net.cells[gid].tags['y']-pop_center)/max_dist
                        cval = 1-relative_center_distance
                    elif color_strategy == 'relative_position':
                        # --- Relative position in the pop
                        cval=cmap((sim.net.cells[gid].tags['y']-pop_yRanges[cell_pop][0])/(pop_yRanges[cell_pop][1]-pop_yRanges[cell_pop][0]))
                    elif color_strategy == 'random':
                        cval=cmap(r_vals[gid])
                    else:
                        # --- Relative position in the pop
                        cval=cmap((sim.net.cells[gid].tags['y']-pop_yRanges[cell_pop][0])/(pop_yRanges[cell_pop][1]-pop_yRanges[cell_pop][0]))
                    
                    # --- Reducing alpha
                    try:    cval_=(cval[0],cval[1],cval[2],alpha)
                    except: cval_=cval

                    cvals.append(cval_)
    
    # Plotting network shape for morphological cells
    sim.analysis.plotShape(
                            # includePre=includePre,includePost=includePost,
                            includePre=includePre,includePost=includePost_gids,
                            # includePre=['all'],includePost=['all'],
                            # includePre=['all'],includePost=morph_pops,
                            # includePre=morph_pops,includePost=['all'],
                            cvals=cvals,
                            # cvals=relative_normY,
                            clim=[0,1],
                            saveFig=sim.cfg.saveFolder+'/'+sim.cfg.simLabel+'_plotShape_'+color_strategy+'__cmap_'+cmap.name+'__every_'+str(int(1/density))+'_k4.png',
                            figSize=[60,60],
                            dpi=1000,
                            bkgColor='k',
                            axis=None,
                            # aspect='off',
                        )
    
    # Simplified network plot
    include_pops = {'MLe__pop':                 {'color':'w','contour':'k'},
                    'VPM__pop':                 {'color':'g','contour':'none'},
                    'TRN__pop':                 {'color':'b','contour':'none'},
                    'TRN_ring__pop':            {'color':'b','contour':'none'},
                    'CTvirtual_uniform__pop':   {'color':'w','contour':'r'},
                    'CTvirtual_activated__pop': {'color':'w','contour':'r'},
                    }
    pop_gids = {pop:sim.net.pops[pop].cellGids for pop in sim.net.pops.keys()}
    
    scramble_pops=[
                    'MLe__pop',
                    'CTvirtual_uniform__pop'
                    ]
    scramble_xz={}
    for pop in scramble_pops:
        if pop in include_pops.keys():

            max_x = max([sim.net.cells[gid].tags['x'] for gid in pop_gids[pop]])
            min_x = min([sim.net.cells[gid].tags['x'] for gid in pop_gids[pop]])
            max_z = max([sim.net.cells[gid].tags['z'] for gid in pop_gids[pop]])
            min_z = min([sim.net.cells[gid].tags['z'] for gid in pop_gids[pop]])

            size_x = max_x-min_x
            size_z = max_z-min_z
            from neuron import h
            random_x_vals_generator = h.Random(10004)
            random_x_vals = [random_x_vals_generator.uniform(min_x, max_x) for i in range(len(pop_gids[pop]))]
            random_z_vals_generator = h.Random(10005)
            random_z_vals = [random_z_vals_generator.uniform(min_z, max_z) for i in range(len(pop_gids[pop]))]

            scramble_xz.update({pop:{'x':random_x_vals,'z':random_z_vals}})

    scatter_points=[];contour_color=[];alpha_values=[]
    for pop in pop_gids.keys():
        if pop not in include_pops.keys():continue
        for gid_ind,gid in enumerate(pop_gids[pop]):
            if pop in scramble_xz.keys():
                scatter_points.append((scramble_xz[pop]['x'][gid_ind],sim.net.cells[gid].tags['y'],scramble_xz[pop]['z'][gid_ind],include_pops[pop]['color']))
            else: 
                scatter_points.append((sim.net.cells[gid].tags['x'],sim.net.cells[gid].tags['y'],sim.net.cells[gid].tags['z'],include_pops[pop]['color']))
            contour_color.append(include_pops[pop]['contour'])
            if pop=='TRN_ring__pop':    alpha_value=0.5
            else:                       alpha_value=1.0
            alpha_values.append(alpha_value)

    fig_size=(50,50)
    dpi = 400
    ax = plt.figure(figsize=fig_size).add_subplot(projection='3d')
    ys=[]
    for point_ind, (x,y,z,c) in enumerate(scatter_points):
        ax.scatter(x, z, y, c=c, alpha=alpha_values[point_ind], edgecolor=contour_color[point_ind]) #here I need to invert
        ys.append(y)
    # --- Invisible y point to fix aspect ratio of plot
    size_y    = (max(ys)-min(ys))
    x1_point = 500+(size_y/2)
    x2_point = 500-(size_y/2)
    z1_point = x1_point
    z2_point = x2_point
    ax.scatter(x1_point, z1_point, max(ys), c='k', alpha=0) #here I need to invert
    ax.scatter(x2_point, z2_point, min(ys), c='k', alpha=0) #here I need to invert

    # --- Adding scalebar
    ax.plot([400,400],[400,400],[5800,5600],'-',color='k',linewidth=5) # y-axis line
    ax.plot([400,600],[400,400],[5800,5800],'-',color='k',linewidth=5) # x-axis line
    ax.plot([400,400],[200,400],[5800,5800],'-',color='k',linewidth=5) # z-axis line

    # # make the panes transparent
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.set_axis_off()
    ax.invert_zaxis()   # here how you can invert --- obs: this is the PYPLOT axis reference, where Z=height

    # plt.savefig('test_2dnet__dpi_'+str(dpi)+'.png',dpi=dpi)
    # plt.savefig('debugCT_test_2dnet__dpi_'+str(dpi)+'.png',dpi=dpi)
    plt.savefig('../paper_figs/network_cell_positions'+str(dpi)+'.png',dpi=dpi)



# Code to plot the 3D distribution of the model network (and optinally plot the connections)
if plotFromCreate:
    # pfsc.plot2dNet(sim,select_sources=['VPM'],select_targets=['L6A'],plotEvery=1,singleConn=False)

    IN_pops=[pop.split('__pop')[0] for pop in netParams.popParams.keys() if 'L6IN' in pop]
    CC_pops=[pop.split('__pop')[0] for pop in netParams.popParams.keys() if 'L6CC' in pop]
    CT_pops=[pop.split('__pop')[0] for pop in netParams.popParams.keys() if 'L6A' in pop]
    pathways=[
                [CT_pops,CT_pops],
                [CT_pops,CC_pops],
                [CT_pops,IN_pops],
                
                [CC_pops,CT_pops],
                [CC_pops,CC_pops],
                [CC_pops,IN_pops],
                
                [IN_pops,CT_pops],
                [IN_pops,CC_pops],
                [IN_pops,IN_pops],
                ]
    for prePops,postPops in pathways:
        # for proj in ['3d','2d_top','2d']:
        # for proj in ['2d_top','2d']:
        # for proj in ['2d_top']:
        for proj in ['2d']:
            pfsc.plot2dNet(sim,select_sources=prePops,select_targets=postPops,plotEvery=1,singleConn=False,show_pops=prePops+postPops,fig_projection=proj)

    pop_pairs=[
                # ['VPM','L6A_activated'],
                # ['L6A_activated','VPM'],
                # ['VPM','TRN'],
                ['TRN','VPM'],
                # ['VPM','TRN_ring'],
                # ['TRN','TRN_ring'],
                # ['L6A_activated','VPM'],
                # ['L6A_sparse','VPM'],
                # ['L6IN_MC_cAC','L6A_activated'],
                # ['L6IN_MC_cAC','L6A_activated'],
                # ['L6A_activated','TRN'],
                ]
    pop_pairs=[
                ['MLe','VPM'],
                ['VPM','TRN'],
                ['TRN','VPM'],
                ['TRN','TRN'],
                ['VPM','L6A_activated'],
                ['L6A_activated','VPM'],
                ['L6A_activated','TRN'],
                ]
    show_pops=['L6A_activated','L6A_silent','L6A_sparse','L6A_suppressed','VPM','TRN','MLe']
    proj = '3d_special'
    angle_ranges=[[0,0.125],[0.125,0.25],[0.25,0.375],[0.375,0.5],[0.5,0.625],[0.625,0.75],[0.75,0.875],[0.875,1.0]]

    for [prePop,postPop] in pop_pairs:
        pfsc.plot2dNet(sim,select_sources=[prePop],select_targets=[postPop],plotEvery=1,singleConn=False,show_pops=show_pops,fig_projection=proj,plotConns=True, cells_alpha=0.25,fig_size=(50,50), dpi = 200)
        for angle_range in angle_ranges: 
            pfsc.plot2dNet(sim,select_sources=[prePop],select_targets=[postPop],plotEvery=1,singleConn=False,show_pops=show_pops,fig_projection=proj,plotConns=True, cells_alpha=0.25,fig_size=(50,50), dpi = 200, show_range=angle_range)

    proj = '3d_special'
    pfsc.plot2dNet(sim,select_sources=None,select_targets=None,plotEvery=600,singleConn=False,show_pops=['VPM', 'TRN', 'TRN_ring', 'CTvirtual_uniform', 'MLe',],fig_projection=proj,plotConns=False, cells_alpha=1.0,fig_size=(50,50), dpi = 200)
    # pfsc.plot2dNet(sim,select_sources=None,select_targets=None,plotEvery=600,singleConn=False,show_pops=['VPM', 'TRN', 'TRN_ring', 'CTvirtual_sparse', 'CTvirtual_suppressed', 'CTvirtual_activated', 'MLe',],fig_projection=proj,plotConns=False, cells_alpha=1.0,fig_size=(50,50), dpi = 200)

    for [prePop,postPop] in pop_pairs:
        for proj in ['3d','2d_top','2d','3d_special']:
        # for proj in ['3d','2d_top','2d']:
            pfsc.plot2dNet(sim,select_sources=None,select_targets=None,plotEvery=700,singleConn=False,show_pops='all',fig_projection=proj,plotConns=False, cells_alpha=1.0)
            pfsc.plot2dNet(sim,select_sources=None,select_targets=None,plotEvery=600,singleConn=False,show_pops=['L6A_activated','L6A_silent','L6A_sparse','L6A_suppressed','VPM','TRN','MLe'],fig_projection=proj,plotConns=False, cells_alpha=1.0)
            
            pfsc.plot2dNet(sim,select_sources=['VPM'],select_targets=['L6A_activated'],plotEvery=1,singleConn=False,show_pops='all',fig_projection=proj,plotConns=True, cells_alpha=1.0)
            pfsc.plot2dNet(sim,select_sources=['VPM'],select_targets=['TRN'],plotEvery=1,singleConn=False,show_pops='all',fig_projection=proj,plotConns=True, cells_alpha=0.25,fig_size=(15,15))
            
            pfsc.plot2dNet(sim,select_sources=[prePop],select_targets=[postPop],plotEvery=1,singleConn=False,show_pops=[prePop,postPop],fig_projection=proj)
            # angle_ranges=[[0,0.125],[0.125,0.25],[0.25,0.375],[0.375,0.5],[0.5,0.625],[0.625,0.75],[0.75,0.875],[0.875,1.0]]
            # for angle_range in angle_ranges: pfsc.plot2dNet(sim,select_sources=[prePop],select_targets=[postPop],plotEvery=1,singleConn=False,show_pops=[prePop,postPop],fig_projection=proj,show_range=angle_range)


############################################################################################################
'''
netParams.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

import NetPyNE_BBP
import numpy as np
from netpyne import specs
import pickle, json
import sys

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cell_cfg import cfg


#------------------------------------------------------------------------------
# --- OTHER IMPORTS
#------------------------------------------------------------------------------
import pandas as pd
import os

#------------------------------------------------------------------------------
# --- VERSION 
#------------------------------------------------------------------------------
netParams.version = 'validate_thalamus_v00'

#------------------------------------------------------------------------------
#
# --- NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# --- General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightNetStims = 1.0 #0.5  # scale conn weight factor for NetStims
# netParams.defaultThreshold = 0      # (threshold for BBP simulation_sonata.json)
netParams.defaultThreshold = -25.0  # spike threshold, 10 mV is NetCon default, lower it for all cells

### reevaluate these values
netParams.defaultDelay = 2.0 # default conn delay (ms) # DEFAULT
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)

### maybe add the edge effect parameter to compensate for the error form distance dependant conn
# netParams.correctBorder = {‘threshold’: [200, 40, 200]}

netParams.defineCellShapes = True # JV 2021-02-23 - Added to fix the lack of the pt3d term in the cells, which make it unable to record i_membrane_

#------------------------------------------------------------------------------
# --- Load BBP circuit
#------------------------------------------------------------------------------
if cfg.convertCellMorphologies: NetPyNE_BBP.Conversion.convertCellMorphology(inputFolder_h5=cfg.morphologyFolder_h5,outputFolder_swc=cfg.NetPyNE_exportedCells,inputFolder_asc=cfg.morphologyFolder_asc)

# --- Load dictionary with thalamic circuit properties
# circuit_dict   = NetPyNE_BBP.LoadBBPCircuit.getDataFrames ( cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.mc_number)
circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new( cfg_file=cfg.sonataConfigFile)

thal_gids = cfg.select_thal_gids

loading_failed=[]
for thal_gid in thal_gids:
    cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
    if cfg.loadCellModel:
        cfg.convertCellModel=False
        NetPyNE_BBP.Prompt.headerMsg('Loading cell '+str(thal_gid))
        try:        netParams.loadCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells_validation+'/netpyne_'+str(thal_gid)+'.json')
        except:     loading_failed.append(thal_gid);print('Loading failed - Cell '+str(thal_gid)+' will be created in the next step')
if len(loading_failed)>0:thal_gids=loading_failed;cfg.convertCellModel=True


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Cell Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
store_cell_properties={}
NetPyNE_BBP.Prompt.headerMsg('Processing cells')
for thal_gid_ind,thal_gid in enumerate(thal_gids):
    print('\t>>\t',str(thal_gid),'\t|\t',str(len(thal_gids)-thal_gid_ind),' cells left')
    cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
    if cfg.convertCellModel:
        thal_gid_conds={}
        # thal_gid_conds = df_thalamus_neurons.loc[thal_gid].to_dict()      # storing cell properties in NetPyNE model
        thal_gid_conds.update({'cell_index':thal_gid})
        netParams.importCellParams(
            label           = cell_properties['mtype']+'__'+str(thal_gid), 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':str(thal_gid)},
            # conds         = {'cellType':str(thal_gid),'mtype':cell_properties['mtype']},
            # conds         = {'cellType':cell_properties['mtype']},
            cellArgs        = [thal_gid,cfg.NetPyNE_exportedCells,cell_properties['morphology']+'.swc'], 
            importSynMechs  = True, 
            somaAtOrigin    = False, 
            cellInstance    = False,
        )

        # --- Fixing up split soma section
        if 'soma_1' in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys(): 
            netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)] = NetPyNE_BBP.ConvertMorphology.mergeSomaSections(cell = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)])
            single_soma_sec=True
        else: single_soma_sec=False

        # --- Adding cell properties in the conds dictionary
        store_cell_properties.update({cell_properties['mtype']+'__'+str(thal_gid):cell_properties.to_dict()})
        # netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds'] = cell_properties.to_dict()

        # --- Fixing bug in sections that don't have a pt3d list (or have an empty dictionary instead)
        for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys():
            if type(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']) is not list:
                netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']=[]

        # --- Creating secLists based on the path distance to soma during network setup
        cell_dict = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]
        soma_pathDist = NetPyNE_BBP.Utils.pathDistance(cell_dict=cell_dict,root_sec='soma_0',ignoreAxon=True)
        basalDends    = NetPyNE_BBP.Utils.findBasalDends(cell_dict,root_sec='soma_0',ignoreAxon=True)
        secLists_dict = NetPyNE_BBP.Utils.secListFromPathDistance(soma_pathDist,basalDendrites=basalDends,repeats=40)
        for secList in secLists_dict.keys(): netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secLists'].update({secList:secLists_dict[secList]})

        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if cfg.saveCellModel: netParams.saveCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells_validation+'/netpyne_'+str(thal_gid)+'.json')
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Pop Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NetPyNE_BBP.Prompt.headerMsg('Processing pops and conns')

if cfg.connType == 'original':
    # --- Declaring variables that store information at the network level
    store_edges_dict={}
    # --- Defining variables to point to the noise files
    noise_filePaths=[('CorticoThalamic_projections',cfg.ct_virtual_noise),('MedialLemniscus_projections',cfg.ml_virtual_noise)]
    
    # --- Adding cell diversity rule
    for thal_gid in thal_gids:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        # --- Creates a separate pop for each cell
        # if not cfg.th_singleTcPop: netParams.popParams[cell_properties['mtype']+'|'+str(thal_gid)+'__pop']={'cellType':cell_properties['mtype'], 'numCells': 1}
        if cfg.th_singleTcPop:  netParams.cellParams[cell_properties['mtype']+'__pop'].update({'diversityFraction':1/len(thal_gids)})
        else:                   netParams.popParams[cell_properties['mtype']+'__'+str(thal_gid)+'__pop']={'cellType':str(thal_gid), 'numCells': 1}

    # # --- Dictionary to store unique references to the presynaptic cells
    # edge_sources_dict_temp={}
    # for pathway in cfg.th_select_pathways: edge_sources_dict_temp.update({pathway:[]})
    # store_source_gids={}
    store_source_gids_byCell={}
    # for pathway in cfg.th_select_pathways: store_source_gids.update({pathway:[]})

    NetPyNE_BBP.Prompt.headerMsg('Loading conns')
    for thal_gid_ind,thal_gid in enumerate(thal_gids):
        print('\n\t>>\t',str(thal_gid),'\t|\t',str(len(thal_gids)-thal_gid_ind),' cells left')

        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        
        # --- Loads the cell directly from the SWC file to recreate the indexing done using 'afferent_section_pos'
        morphologyFullPath = cfg.NetPyNE_exportedCells+'/'+cell_properties['morphology']+'.swc'
        full_ids_list = NetPyNE_BBP.Conversion.getAfferentSectionIds(morphologyFullPath,single_soma_sec=single_soma_sec)
        # continue
        
        if   cfg.select_microcircuit is None:       cell_conns_file = str(thal_gid)+'_edges_full.json'
        elif type(cfg.select_microcircuit) == int:  cell_conns_file = str(thal_gid)+'_edges_mc'+str(cfg.select_microcircuit)+'.json'
        else: sys.exit('Select a microcircuit or choose <None> ')
            
        edges_fileName  = cfg.NetPyNE_node_pathway_edges_validation+'/'+cell_conns_file
        
        # --- Overwrites cfg.th_updateConns if conn file does not exits
        if (cfg.th_updateConns) or (cell_conns_file not in os.listdir(cfg.NetPyNE_node_pathway_edges_validation)):
            
            if (cell_conns_file not in os.listdir(cfg.NetPyNE_node_pathway_edges_validation)): print('\t\t--\tFile does not exist... Creating conns for cell '+str(thal_gid))
            else: print('\t\t--\tUpdating conns for cell '+str(thal_gid))

            # --- Generates a Dictionary of Dataframes with the EDGES (or conns) information for each pathway targeting the selected NODE
            #     Read mode about edges: https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md
            node_pathway_edges = NetPyNE_BBP.LoadBBPCircuit.getNodeEdges(       cfg.sonataConfigFile,
                                                                                thal_gid,
                                                                                special_conditions={'inter_sources':cfg.th_inter_sources,'inter_targets':[],'select_microcircuit':cfg.select_microcircuit}
                                                                                )
            # --- Adds section information in NetPyNE format to the edge properties
            node_pathway_edges = NetPyNE_BBP.ConvertSynapses.getTargetSections( cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                                node_pathway_edges  = node_pathway_edges,
                                                                                full_ids_list       = full_ids_list
                                                                                )
            # --- Adds a remapped version of the 3d positions of each synapse, once the soma is set to (0,0,0) when a SWC morphology is loaded
            node_pathway_edges = NetPyNE_BBP.ConvertSynapses.convert3DLocation( cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                                cell_properties     = cell_properties,
                                                                                node_pathway_edges  = node_pathway_edges,
                                                                                )
            # --- Converting conns from DICT of DFs to DICT of DICTs and writing to JSON
            node_pathway_edges_dict = NetPyNE_BBP.LoadBBPCircuit.saveNodeEdges(node_pathway_edges,edges_fileName)
        else:
            print('\t\t--\tLoading conns/edges from stored dataset: ', edges_fileName)
            # --- Loads the connectivity information previously saved
            node_pathway_edges = NetPyNE_BBP.LoadBBPCircuit.loadNodeEdges(edges_fileName)

        if cfg.plotSynLocation:
            # --- Plots the location each synapse in the cell model
            NetPyNE_BBP.ConvertSynapses.plotSynapseLocation(    cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                cell_properties     = cell_properties,
                                                                node_pathway_edges  = node_pathway_edges,
                                                                full_ids_list       = full_ids_list,
                                                                extra_flag          = '',
                                                                plotSimplified      = False
                                                                )
        if cfg.plotCellShape:
            NetPyNE_BBP.ConvertSynapses.plotCellShape(          cell=netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                extra_flag='_cellShape'
                                                                )

        edges = node_pathway_edges

        # --- Selects which pathways should be included from the original model (e.g.: ct, ml, rt, in)
        if cfg.th_select_pathways is not None: edges = NetPyNE_BBP.LoadBBPCircuit.selectPathways(cfg.th_select_pathways,edges)

        # --- Switches the Synaptic parameters by the average values reported in the Table S2 in the original paper instead of the model stored values
        if cfg.th_useTableVals: edges = NetPyNE_BBP.ConvertSynapses.replaceSynParams(cell_properties,edges)

        # --- Dictionary store the presynaptic cells that project to the evaluated cell
        store_source_gids_byCell.update({thal_gid:{}})
        # edge_sources_dict_temp={}
        for pathway in cfg.th_select_pathways: 
            store_source_gids_byCell[thal_gid].update({pathway:[]})
            # store_source_gids.update({pathway:[]})

        # --- Updates the dictionary to store the presynaptic cells that project to the evaluated cell
        for edge_name in edges.keys():
            a,b,edge_source_pop,c = edge_name.split('__')
            edge_sources = list(set(edges[edge_name]['@source_node']));edge_sources.sort()
            store_source_gids_byCell[thal_gid][edge_source_pop]+=edge_sources
            store_source_gids_byCell[thal_gid][edge_source_pop].sort()
            # --- Stores a list of presynaptic gids that will be later filtered
            # store_source_gids[edge_source_pop]+=edge_sources

        # --- Saves a dictionary of noise input to each postsynaptic cell
        if cfg.saveIndividualNoise:
            for (key,noise_filePath) in noise_filePaths: 
                if key in store_source_gids_byCell[thal_gid]: 
                    cell_source_gids = list(set(store_source_gids_byCell[thal_gid][key])); cell_source_gids.sort()
                    NetPyNE_BBP.LoadSimResults.getVirtualInputSpikes(noise_filePath, filterGids=cell_source_gids, saveName= cfg.NetPyNE_input_noise_savePath+'/post_cell_inputs/'+'tgt|'+str(thal_gid)+'__'+'src|'+key, loadFromSourceFile=False)
        
        # --- Stores the edge datasets so that it can be re-iterated over without overwriting
        for edge_name in edges.keys(): store_edges_dict.update({edge_name:edges[edge_name]})

    ##########################################################################################################################################################
    # --- Combines the presynaptic cells into a single dictionary to create the Vectstim populations only with connected cells
    #  -  Dictionary for each post cell, with all pre cells grouped by pathway
    edge_sources_dict_={}
    for edge_source_pop in cfg.th_select_pathways: edge_sources_dict_.update({edge_source_pop:[]})
    for thal_gid in store_source_gids_byCell.keys():
        for edge_source_pop in store_source_gids_byCell[thal_gid].keys():
            source_gids = list(set(store_source_gids_byCell[thal_gid][edge_source_pop]))
            source_gids.sort()
            edge_sources_dict_[edge_source_pop]+=source_gids
    
    #  -  Dictionary with all pre cells grouped by pathway
    edge_sources_dict={}
    for edge_source_pop in edge_sources_dict_.keys():
        edge_sources_dict.update({edge_source_pop:list(set(edge_sources_dict_[edge_source_pop]))})
        edge_sources_dict[edge_source_pop].sort()

    ##########################################################################################################################################################
    NetPyNE_BBP.Prompt.headerMsg('Loading internal inputs')
    # --- Loading dictionary of spike times
    th_spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=cfg.th_spikes_file, cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.select_microcircuit, showFig=False) # --- intrathalamic cells (biophysiscal)
    # --- Loading dictionary of spike times for connected cells only -- all cells together, regardless of pre - post conns
    th_connected_cells_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromConnectedCells(th_spk_dict, thal_gids, edge_sources_dict, showFig=False) # --- intrathalamic cells (biophysiscal)

    NetPyNE_BBP.Prompt.headerMsg('Loading external inputs')
    # --- Loading dictionaries of input noise
    noise_dict={}
    for (key,noise_filePath) in noise_filePaths: 
        if key in edge_sources_dict: noise_dict.update({key:NetPyNE_BBP.LoadSimResults.getVirtualInputSpikes(noise_filePath, filterGids=edge_sources_dict[key], saveName= cfg.NetPyNE_input_noise_savePath+'/'+key, loadFromSourceFile=False)})

    ##########################################################################################################################################################
    NetPyNE_BBP.Prompt.headerMsg('Creating presynaptic pops - Vecstims - one cell per pop')
    # --- Creating VecStims pops --- one cell per pop (pops based on preCell bbp_gid)
    for edge_source_pop in edge_sources_dict.keys():
        source_nodes = edge_sources_dict[edge_source_pop]
        if ('CorticoThalamic' in edge_source_pop) or ('MedialLemniscus' in edge_source_pop): 
            spkids = list(noise_dict[edge_source_pop].keys())
            spkts  = list(noise_dict[edge_source_pop].values())
        elif 'thalamus_neurons' in edge_source_pop:
            edge_source_pop_=edge_source_pop.split('|')
            spkids = list(th_connected_cells_dict[edge_source_pop_[1]].keys())
            spkts  = list(th_connected_cells_dict[edge_source_pop_[1]].values())

        # --- Adds a 'skipTime' to shift value of spike times
        if cfg.skipTime is not None: 
            for spkid_ind, spkid in enumerate(spkids):spkts[spkid_ind] = [spkt+cfg.skipTime for spkt in spkts[spkid_ind]]

        # --- Creates NetPyNE Vecstim Pops --- one cell per pop
        for spkid_ind, spkid in enumerate(spkids):
            if int(spkid) in edge_sources_dict[edge_source_pop]:
                try:
                    if len(spkts[spkid_ind])<1: spkts[spkid_ind]=[cfg.duration+1000]
                except: 
                    if type(spkts[spkid_ind])==None: spkts[spkid_ind]=[cfg.duration+1000] # adds a dummy spike at t>simDuration to prevent code from breaking because there were no spikes at the cell

                # netParams.popParams.update({ edge_source_pop+'__'+str(spkid)+'__vecstim_pop':{  'cellModel':'VecStim',
                #                                                                                 'numCells':1,
                #                                                                                 'rate': 0.01,
                #                                                                                 # 'spkTimes': [0,1,2,3],
                #                                                                                 # 'spkTimes': spkts[spkid_ind],
                #                                                                                 'pulses':[]}})
                netParams.popParams.update({ edge_source_pop+'__'+str(spkid)+'__vecstim_pop':{  'cellModel': 'VecStim',
                                                                                                'cellsList':[
                                                                                                                {
                                                                                                                    'x':0,
                                                                                                                    'y':0,
                                                                                                                    'z':0,
                                                                                                                    'spkTimes':  spkts[spkid_ind],   # spike outside simDuration
                                                                                                                    }
                                                                                                                ],
                                                                                                # 'noise':     1,
                                                                                                # 'numCells':  1, 
                                                                                                # 'rate':  0,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                                                                # 'xRange':   [0,1], 
                                                                                                # 'yRange':   [0,1], 
                                                                                                # 'zRange':   [0,1],
                                                                                                }})
    
    ##########################################################################################################################################################

    # --- Selecting type of MOD file to be used
    if cfg.modType == 'Prob_original':  modAMPANMDA = 'ProbAMPANMDA_EMS_original';         modGABA     = 'ProbGABAAB_EMS_original'  # original MOD from BBP model
    else:                               modAMPANMDA = 'ProbAMPANMDA_EMS';                  modGABA     = 'ProbGABAAB_EMS'           # Original Thalamus implementation of the BBP mod files
    NetPyNE_BBP.Prompt.headerMsg('MOD templateAMPA: '+modAMPANMDA+'   GABA: '+modGABA)
    # print('\n\t>>\tMOD template\tAMPA: ', modAMPANMDA, '\tGABA: ', modGABA)
    
    ##########################################################################################################################################################
    # --- Processing edges to create conns
    NetPyNE_BBP.Prompt.headerMsg('Processing edges to create conns')
    edge_modif_dict = {}
    for edge_name in store_edges_dict.keys():
        # print('\t>>\tProcessing ', edge_name)
        edge_type,edge_mech,edge_source_pop,edge_target = edge_name.split('__')
        print('\n\t>>\tProcessing ', edge_target)
        print('\t\t--\t',edge_source_pop,'\t|\t',edge_type,'\t|\t',edge_mech)
        
        if (cfg.removeChemical)     and edge_mech == 'chemical':            continue
        if (cfg.removeElectrical)   and edge_mech == 'electrical_synapse':  continue

        selected_edge_properties_list=[ '@target_node','@source_node','conductance','sec','afferent_section_pos',
                                        'efferent_section_pos',
                                        'u_syn','depression_time','facilitation_time','decay_time','n_rrp_vesicles','delay',
                                        #'NMDA_ratio', # <--- updated later
                                       ]
        
        # - Obs: (synsPerConn == 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
        if       edge_mech       == 'electrical_synapse':           mod_file = 'Gap'
        else:
            if   edge_source_pop == 'CorticoThalamic_projections':  mod_file = modAMPANMDA
            elif edge_source_pop == 'MedialLemniscus_projections':  mod_file = modAMPANMDA
            elif edge_source_pop == 'thalamus_neurons|VPL_TC':      mod_file = modAMPANMDA
            else:                                                   mod_file = modGABA

            TableS1         = NetPyNE_BBP.StoreParameters.getTableS1()
            cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,int(edge_target))
            try:
                syn_values      = TableS1[edge_source_pop][cell_properties['mtype']]
            except:
                print('Skipping invalid pathway:\t',edge_source_pop,' -> ', cell_properties['mtype'])
                continue
            syn_values['paper_reference_values'].update({'n_rrp_vesicles':1}) # only present in the Prob MOD, not in the Det MOD

        # --- Creating a dictionary with information referenced by edge gids
        bbp_mechs_dict={}

        for edge_index in store_edges_dict[edge_name].index:
            edge_data = store_edges_dict[edge_name].loc[edge_index]
            bbp_mechs_dict.update({edge_data.name:{}})
            edge_properties = edge_data.index
            edge_vals       = edge_data.values
            for edge_property_ind,edge_property in enumerate(edge_properties):
                if edge_property in selected_edge_properties_list: bbp_mechs_dict[edge_data.name].update({edge_property:edge_vals[edge_property_ind]})
            # --- Updates the NMDA_ratio value based on data from Table S2 on the paper
            if edge_mech == 'chemical': bbp_mechs_dict[edge_data.name].update({'NMDA_ratio':syn_values['paper_reference_values']['NMDA_ratio']}) # --- set to None in GABA mechanisms

        # # --- Connection Overrides - alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model
        if cfg.th_connOverrides: edge_modif_dict.update({edge_name:NetPyNE_BBP.CreateNetPyNE.modifySyn(circuit_dict, edge_name)})
        # if cfg.th_connOverrides: syn_mod = NetPyNE_BBP.CreateNetPyNE.modifySyn(circuit_dict, edge_name)

        # --- Creating NetPyNE syn mechs
        edge_dicts={}
        for mech_id in bbp_mechs_dict.keys():
            edge_dict = NetPyNE_BBP.CreateNetPyNE.modTemplate(bbp_mechs_dict[mech_id],mod_file)
            # --- Connection Overrides - alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model
            if cfg.th_connOverrides: edge_dict.update(edge_modif_dict[edge_name]) # --- updates the edge dict with the modifications from BlueConfig file
            edge_dicts.update({mech_id:edge_dict})

        for mech_id in bbp_mechs_dict.keys():
            netParams.synMechParams.update({edge_source_pop+'__'+str(mech_id)+'__vecstim_mech': edge_dicts[mech_id]})
        
        # --- Creating NetPyNE conns
        for mech_id in bbp_mechs_dict.keys():
            # print('conn: \t',edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node']),' --> ',bbp_mechs_dict[mech_id]['@target_node'])
            

            '''
            IMPORTANT NOTE:

            THIS CODE WAS CHANGED BACK FROM USING THE CONDUCTANCE AS THE WEIGHT, TO HAVING WEIGHT == 1 AND CONDUCTANCE SET AT THE SYNAPTIC MECHANISM, IN THE NetPyNEBBP.modTemplate method

            '''



            netParams.connParams.update({edge_source_pop+'__'+str(mech_id)+'__vecstim_conn': {  'preConds':  {'pop': edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node'])+'__vecstim_pop'}, 
                                                                                                'postConds': {'pop': cell_properties['mtype']+'__'+str(bbp_mechs_dict[mech_id]['@target_node'])+'__pop'},
                                                                                                'probability':  1,
                                                                                                'weight':       1,
                                                                                                # 'weight':       bbp_mechs_dict[mech_id]['conductance'], # test: 2023_11_23
                                                                                                'synsPerConn':  1, # - Obs: (synsPerConn = 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
                                                                                                'synMech':      edge_source_pop+'__'+str(mech_id)+'__vecstim_mech',
                                                                                                'sec':          bbp_mechs_dict[mech_id]['sec'],
                                                                                                'loc':          bbp_mechs_dict[mech_id]['afferent_section_pos'],
                                                                                                }})
            if   edge_mech == 'chemical':                                   netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'delay':       bbp_mechs_dict[mech_id]['delay']})
            # elif edge_mech == 'electrical_synapse':                         netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'gapJunction': True})
            if   'efferent_section_pos' in bbp_mechs_dict[mech_id].keys():  netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'preLoc':      bbp_mechs_dict[mech_id]['efferent_section_pos']})
            
            # # --- test - add minis using <netstim_inhpoisson.mod>
            # netParams.connParams.update({edge_source_pop+'__'+str(mech_id)+'__minis_conn': {    'preConds':  {'pop': edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node'])+'__minis_pop'}, 
            #                                                                                     'postConds': {'pop': cell_properties['mtype']+'__'+str(bbp_mechs_dict[mech_id]['@target_node'])+'__pop'},
            #                                                                                     'probability':  1,
            #                                                                                     'weight':       bbp_mechs_dict[mech_id]['conductance'],
            #                                                                                     'synsPerConn':  1, # - Obs: (synsPerConn = 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
            #                                                                                     'synMech':      edge_source_pop+'__'+str(mech_id)+'__vecstim_mech',
            #                                                                                     # 'synMech':      'minis_mech',
            #                                                                                     # 'delay':        bbp_mechs_dict[mech_id]['delay'],
            #                                                                                     'sec':          bbp_mechs_dict[mech_id]['sec'],
            #                                                                                     'loc':          bbp_mechs_dict[mech_id]['afferent_section_pos'],
            #                                                                                     }})
            # if   edge_mech == 'chemical':                                   netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'delay':       bbp_mechs_dict[mech_id]['delay']})
            # elif edge_mech == 'electrical_synapse':                         netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'gapJunction': True})
            # if   'efferent_section_pos' in bbp_mechs_dict[mech_id].keys():  netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'preLoc':      bbp_mechs_dict[mech_id]['efferent_section_pos']})

    # --- Rescale conn weights
    for conn_name in netParams.connParams.keys(): netParams.connParams[conn_name]['weight']*=cfg.rescale_conn_weight[conn_name.split('__')[0]]

    # --- Rescale USE parameter (probability of synapse activation)
    if cfg.rescaleUSE is not None:
        for mech in netParams.synMechParams.keys():
            try:    netParams.synMechParams[mech]['Use']*=cfg.rescaleUSE
            except: continue

    # --- Identifies which thalamic populations are present in the instance of the model
    pops_dict={}
    for pop_mtype in list(set(circuit_dict['thalamus_neurons']['mtype'])):
        pops_dict.update({pop_mtype:[]})
        pops_dict[pop_mtype]=[pop_ for pop_ in netParams.popParams.keys() if (pop_mtype in pop_) and ('vecstim' not in pop_)]

    ######################################################################################################################################################################################################
    # --- Adds Ornstein Uhlenbeck Noise current
    if cfg.addNoiseIClamp:
        for thal_gid in cfg.select_thal_gids:
        
            print('             ----                    ADDING STIMS')

        # for pop in cfg.NoiseIClampParams.keys():
            cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
            netParams.stimSourceParams['NoiseIClamp_source__'+str(thal_gid)] = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 
                                                                        #    'amp': 0,
                                                                           'amp': store_cell_properties[cell_properties['mtype']+'__'+str(thal_gid)]['@dynamics:holding_current']+(cfg.scale_threshold_current*store_cell_properties[cell_properties['mtype']+'__'+str(thal_gid)]['@dynamics:threshold_current']),
                                                                        #    'amp': cfg.NoiseIClampParams[pop]['amp'],
                                                                           }
            netParams.stimTargetParams['NoiseIClamp_target__'+str(thal_gid)] = {'source': 'NoiseIClamp_source__'+str(thal_gid), 'sec':'soma_0', 'loc': 0.5, 
                                                                                'conds': {
                                                                                            'cellType':str(thal_gid),
                                                                                            # 'pop':cell_properties['mtype']+'__'+str(thal_gid)+'__pop',
                                                                                            }
                                                                                }


    # if cfg.addNoiseIClamp:
    #     for pop in cfg.NoiseIClampParams.keys():
    #         netParams.stimSourceParams['NoiseIClamp_source__'+pop] = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': cfg.NoiseIClampParams[pop]['amp']}
    #         netParams.stimTargetParams['NoiseIClamp_target__'+pop] = {'source': 'NoiseIClamp_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}

    # ######################################################################################################################################################################################################
    # # --- Adds a threshold current stim to thalamic populations
    # if cfg.use_BBP_thresholdCurrent:
    #     NetPyNE_BBP.Prompt.headerMsg('Adding Threshold Current stim (IClamp)')
    #     try: 
    #         threshold_current_amp = cell_properties['@dynamics:threshold_current']
    #         print('\t--- Adding a Threshold current: ',threshold_current_amp,' (nA)')
    #     except:
    #         print('\t--- Failed to get Threshold current')
    #         threshold_current_amp = 0

    #     for ind_th_pop,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
    #     # for ind_th_pop,th_pop in enumerate(['VPL_TC']):
    #         if len(pops_dict[th_pop])>0:
    #             netParams.stimSourceParams['threshold_IClamp_'+str(ind_th_pop)] = {'type': 'IClamp', 'del': 0, 'dur': cfg.duration, 'amp': threshold_current_amp}
    #             netParams.stimTargetParams['threshold_IClamp_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'threshold_IClamp_'+str(ind_th_pop), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}

    # ######################################################################################################################################################################################################
    # # --- Adds a holding current stim to thalamic populations
    # if cfg.use_BBP_holdingCurrent:
    #     NetPyNE_BBP.Prompt.headerMsg('Adding Holding Current stim (IClamp)')
    #     try: 
    #         holding_current_amp = cell_properties['@dynamics:holding_current']
    #         print('\t--- Adding a Holding current: ',holding_current_amp,' (nA)')
    #     except:
    #         print('\t--- Failed to get Holding current')
    #         holding_current_amp = 0

    #     for ind_th_pop,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
    #     # for ind_th_pop,th_pop in enumerate(['VPL_TC']):
    #         if len(pops_dict[th_pop])>0:
    #             netParams.stimSourceParams['holding_IClamp_'+str(ind_th_pop)] = {'type': 'IClamp', 'del': 0, 'dur': cfg.duration, 'amp': holding_current_amp}
    #             netParams.stimTargetParams['holding_IClamp_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'holding_IClamp_'+str(ind_th_pop), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}
    
    # ######################################################################################################################################################################################################
    # # --- Adds a current stim to thalamic populations
    # if cfg.add_current_stims:
    #     NetPyNE_BBP.Prompt.headerMsg('Adding Current stim (IClamp)')
    #     for ind_th_pop,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
    #     # for ind_th_pop,th_pop in enumerate(['VPL_TC']):
    #         if len(pops_dict[th_pop])>0:
    #             netParams.stimSourceParams['IClamp_'+str(ind_th_pop)] = {'type': 'IClamp', 'del': cfg.current_stim_start, 'dur': cfg.duration, 'amp': cfg.current_stim_amp}
    #             netParams.stimTargetParams['IClamp_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'IClamp_'+str(ind_th_pop), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}
    
    # # --- Adds a noise current stim to thalamic populations
    # if cfg.add_current_stims_noise:
        
    #     # netParams.stimSourceParams['IClamp2'] = {'type': 'IClamp', 'del': 0, 'dur': cfg.skipTime, 'amp': cfg.noise_string}
    #     netParams.stimSourceParams['IClamp2'] = {'type': 'IClamp', 'del': cfg.current_stim_start, 'dur': cfg.current_stim_duration, 'amp': 'normal(0,0.001)'}
        
    #     for cell in netParams.cellParams.keys():
    #         netParams.stimTargetParams['IClamp2__'+cell+'__'+sec] = {'source': 'IClamp2', 'sec':'soma_0', 'loc': 0.5, 'conds': {'cellType':cell.split('__')[1]}}
            
    # ######################################################################################################################################################################################################

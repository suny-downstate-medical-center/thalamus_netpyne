'''
barreloid_netParams.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

#------------------------------------------------------------------------------
# --- LIBRARY IMPORTS
#------------------------------------------------------------------------------
from matplotlib.pyplot import plot
import NetPyNE_BBP
import numpy as np
from netpyne import specs
import pickle, json
import sys
import math
from GenerateStimDataset import SampleData

import Build_Net as BN

import pandas as pd
import os

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from barreloid_cfg import cfg

#------------------------------------------------------------------------------
# --- VERSION 
#------------------------------------------------------------------------------
netParams.version = 'thalamus_v00'

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
netParams.defaultThreshold = -20.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 100000.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
netParams.shape='cylinder'
netParams.defineCellShapes = True # JV 2021-02-23 - Added to fix the lack of the pt3d term in the cells, which make it unable to record i_membrane_
netParams.store_cell_properties={}

#------------------------------------------------------------------------------
# --- Load BBP circuit
#------------------------------------------------------------------------------
if cfg.convertCellMorphologies: NetPyNE_BBP.Conversion.convertCellMorphology(inputFolder_h5=cfg.morphologyFolder_h5,outputFolder_swc=cfg.NetPyNE_exportedCells,inputFolder_asc=cfg.morphologyFolder_asc)

# --- Load dictionary with thalamic circuit properties
# circuit_dict = NetPyNE_BBP.LoadBBPCircuit.getDataFrames ( cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.mc_number)

if cfg.loadCircuitProperties:
    NetPyNE_BBP.Prompt.headerMsg('Loading converted circuit from pickle file in \t\t'+cfg.stored_circuit_path)
    # with open(cfg.stored_circuit_path, 'rb') as f: circuit_dict = pickle.load(f)
    with open(cfg.stored_circuit_path, 'rb') as f: circuit_dict = pd.read_pickle(f)
else:
    NetPyNE_BBP.Prompt.headerMsg('Loading circuit from original project')
    circuit_dict = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new( cfg_file=cfg.sonataConfigFile)
    if cfg.saveCircuitProperties: 
        with open(cfg.stored_circuit_path, 'wb') as f: pd.to_pickle(circuit_dict, f)
        # with open(cfg.stored_circuit_path, 'wb') as f: pickle.dump(circuit_dict, f)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Cell Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
# --- Load cell models from converted NetPyNE morpohologies
loading_failed=[]
if cfg.loadCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Loading stored cells from \t\t'+cfg.NetPyNE_JSON_cells)
    for thal_gid in list(set(cfg.select_thal_gids)): # to avoid repetition if the same cell template is used multiple times
    # for thal_gid in cfg.select_thal_gids:
        print('Loading cell # ', str(thal_gid))
        cfg.convertCellModel=False
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        
        # --- Store the cell properties in a dictionary in netParams
        netParams.store_cell_properties.update(NetPyNE_BBP.LoadBBPCircuit.cell_properties_to_dict(cell_properties))
        
        # NetPyNE_BBP.Prompt.headerMsg('Loading cell '+str(thal_gid))
        try:        netParams.loadCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        except:     loading_failed.append(thal_gid);print('Loading failed - Cell '+str(thal_gid)+' will be created in the next step')
        # --- Check if the converted model has the necessary secLists to create the connections, otherwise, converts the cell again
        secList_names = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secLists'].keys()
        if ('inputs__proximal' or 'inputs__intermediate' or 'inputs__distal') not in secList_names:loading_failed.append(thal_gid);print('SecList missing - Cell '+str(thal_gid)+' will be created in the next step')

        # if cfg.reset_cell_rotation:
        #     desired_angle_y_deg=90
        #     angle = np.deg2rad(desired_angle_y_deg)-cell_properties.orientation_y
        #     rotated_secs=NetPyNE_BBP.RotateMorphology.rotate(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],angle=angle)
        #     # rotated_secs=NetPyNE_BBP.RotateMorphology.rotate(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],cell_properties)
        #     for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys(): netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']=rotated_secs[sec]
        #     plotCellShape=True
        #     if plotCellShape: 
        #         try: NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],'_rotated_cellShape',cell_label=cell_properties['mtype']+'__'+str(thal_gid),savefolder='../validation/rotated_morphologies',axis_off=False)
        #         except:pass # in case it tries to save the fig in the HPC

if len(loading_failed)>0:cfg.select_thal_gids=loading_failed;cfg.convertCellModel=True

# --- Load cell models from morpohology templates
if cfg.convertCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Processing cells')
    for thal_gid_ind,thal_gid in enumerate(list(set(cfg.select_thal_gids))): # to avoid repetition if the same cell template is used multiple times
    # for thal_gid_ind,thal_gid in enumerate(cfg.select_thal_gids):
        print('\t>>\t',str(thal_gid),'\t|\t',str(len(cfg.select_thal_gids)-thal_gid_ind),' cells left')
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)

        # --- Store the cell properties in a dictionary in netParams
        netParams.store_cell_properties.update(NetPyNE_BBP.LoadBBPCircuit.cell_properties_to_dict(cell_properties))

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
            print('\t--- Check behavior of model after fixing extra soma_1 section issue')
            netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)] = NetPyNE_BBP.ConvertMorphology.mergeSomaSections(cell = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)])
            single_soma_sec=True
        else: single_soma_sec=False

        # # --- Adding cell properties in the conds dictionary
        # netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds'].update(cell_properties.to_dict())
        # # cell_properties_dict=cell_properties.to_dict()
        # # for key in cell_properties_dict:
        # #     if '@dynamics:' in key:
        # #         key_ = key.split('@dynamics:')[1]
        # #         netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds'].update({key_:cell_properties_dict[key]})
        # #     else:
        # #         netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds'].update({key:cell_properties_dict[key]})
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        # --- Fixing bug in sections that don't have a pt3d list (or have an empty dictionary instead)
        for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys():
            if type(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']) is not list:
                netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']=[]
        
        # --- Creating secLists based on the path distance to soma during network setup
        cell_dict = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]
        soma_pathDist = NetPyNE_BBP.Utils.pathDistance(cell_dict=cell_dict,root_sec='soma_0',ignoreAxon=True)
        # basalDends    = NetPyNE_BBP.Utils.findBasalDends(cell_dict,root_sec='soma_0',ignoreAxon=True)
        basalDends    = [] # testing removing basal dendrites from 'inputs__proximal' list because some were too far from the soma
        secLists_dict = NetPyNE_BBP.Utils.secListFromPathDistance(soma_pathDist,basalDendrites=basalDends,repeats=40)
        for secList in secLists_dict.keys(): netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secLists'].update({secList:secLists_dict[secList]})
        
        if cfg.reset_cell_rotation:
            plotCellShape=False
            
            if plotCellShape: NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],'_beforeRotation_cellShape',cell_label=cell_properties['mtype']+'__'+str(thal_gid),savefolder='../validation/rotated_morphologies')
            
            print('Reseting cell rotation for cell ', cell_properties['mtype']+'__'+str(thal_gid))
            if cell_properties['mtype']+'__'+str(thal_gid) in cfg.rotation_angle.keys():    rotation_angle = cfg.rotation_angle[cell_properties['mtype']+'__'+str(thal_gid)]
            else:                                                                           rotation_angle = 0

            rotated_secs=NetPyNE_BBP.RotateMorphology.rotate(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],angle=rotation_angle)
            # rotated_secs=NetPyNE_BBP.RotateMorphology.rotate(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],cell_properties)
            for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys(): netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']=rotated_secs[sec]
            if plotCellShape: NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],'_rotated_cellShape',cell_label=cell_properties['mtype']+'__'+str(thal_gid),savefolder='../validation/rotated_morphologies')


        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if cfg.saveCellModel: netParams.saveCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for gids_list in cfg.gids_lists:
    # --- Adding cell diversity rule
    for thal_gid in gids_list:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        # netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(gids_list)})
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(list(set(gids_list)))})
        # --- Changes cellType to the same value, so that the different morphologies are added to the same pop can be identified by connParams
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds']['cellType']=cell_properties['mtype']

plotCellShape=False
if plotCellShape:
    for thal_gid in list(set(cfg.select_thal_gids)): # to avoid repetition if the same cell template is used multiple times
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],'_cellShape',cell_label=cell_properties['mtype']+'__'+str(thal_gid),savefolder='../validation/barreloid_morphologies',axis_off=False)

plotPaperCellShape = False
if plotPaperCellShape:
    for thal_gid in list(set(cfg.select_thal_gids)): # to avoid repetition if the same cell template is used multiple times
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        if   cell_properties['mtype']=='VPL_TC':    morph_color = 'g'
        elif cell_properties['mtype']=='Rt_RC':     morph_color = 'b'
        else:                                       morph_color = 'k'
        NetPyNE_BBP.ConvertSynapses.plotCellShape(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],'_cellShape',cell_label=cell_properties['mtype']+'__'+str(thal_gid),savefolder='../paper_figs/paper_barreloid_morphologies',axis_off=True, morph_color=morph_color)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- VPM and TRN pops 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if cfg.includeTH:
    thal_pops=['VPM','TRN']
    netParams.popParams  = BN.BuildNetwork.getPopTemplate(pops=thal_pops,center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape=cfg.pop_shape)
    netParams.popParams['VPM__pop']['cellType']='VPL_TC'
    netParams.popParams['TRN__pop']['cellType']='Rt_RC'

    if cfg.align_thal_y:
        diam,height         = BN.BuildNetwork.getNetworkSize(cfg.center_point)
        for pop in thal_pops:
            print('Aligning '+pop+' cells in the y-axis')
            # pop='VPM'
            # --- Creates cells in a random distribution in the x-z planes, and in a linear distribution in the y-plane
            from neuron import h
            # Initialize NEURON random number generators
            if cfg.base_random_seed is not None: 
                rand1 = h.Random(cfg.base_random_seed+10)
                rand2 = h.Random(cfg.base_random_seed+20)
            else:
                rand1 = h.Random(10)
                rand2 = h.Random(20)

            # new logic - extends the list from [0 .. 1] to [0 .. 1*<number of cell templates>] (e.g. [0 .. 5]), and then takes the remainder of the relative position in respect to 1, to get a value between 0 and 1
            # e.g.: when there are 5 templates, cell template 2 will be between 0.2 and 0.4, but with the new logic they are between [1 and 2], and the remainder of 1.xxx is .xxx, resulting in an even distribution between 0 and 1
            #       1.8%1 = 0.8, which would be in the relative 0.8 position
            pop_name = pop+'__pop'
            pop_cellType = netParams.popParams[pop_name]['cellType']
            count_cell_templates = len([cell_template for cell_template in netParams.cellParams.keys() if pop_cellType in cell_template])
            pop_cellsList=[{
                            'x':rand1.uniform(diam[pop][0], diam[pop][1]),
                            'z':rand2.uniform(diam[pop][0], diam[pop][1]),
                            'y':(height[pop][0])+((height[pop][1]-height[pop][0])*((y_pos/netParams.popParams[pop+'__pop']['numCells'])%1)) # (Timofeeva et al., 2003),
                            } for y_pos in np.arange(0,netParams.popParams[pop+'__pop']['numCells']*count_cell_templates,count_cell_templates)
                            ]

            # pop_cellsList = [{
            #                 'x':rand1.uniform(diam[pop][0], diam[pop][1]),
            #                 'z':rand2.uniform(diam[pop][0], diam[pop][1]),
            #                 'y':(height[pop][0])+((height[pop][1]-height[pop][0])*(y_pos/netParams.popParams[pop+'__pop']['numCells'])) # (Timofeeva et al., 2003),
            #                 } for y_pos in np.arange(0,netParams.popParams[pop+'__pop']['numCells'],1)
            #                 ]
            netParams.popParams[pop+'__pop'].update({'cellsList':pop_cellsList})

    # --- Adjacent TRN ring of cells
    if cfg.addTRNadjacentpops is not None: netParams.popParams.update( BN.BuildNetwork.getL6PopTemplateCtCcIn( pops=cfg.TRN_adacent_pops, center_point=cfg.center_point,pop_flag='__pop',cell_flag='',diversity=True,volume_shape=cfg.pop_shape))
    
    # --- Prints information about the cell types in the thalamus
    if cfg.printCellProperties:
        store_etypes={'dAD_ltb':[],'dNAD_ltb':[],'cNAD_noscltb':[],'other':[]}
        for gid in cfg.select_thal_gids:
            cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,gid)
            try:    store_etypes[cell_properties['etype']].append(gid)
            except: store_etypes['other'].append(gid)
            print(gid, cell_properties['mtype'], cell_properties['etype'])
        [print(key, len(store_etypes[key])) for key in store_etypes.keys()]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- L6A cells and pop
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if cfg.includeL6:
    if cfg.addL6Adetailed:
        # --- Based on (Dash ... Crandall, 2022) and (Kwegyir-Afful, 2009)
        print('\n\t>>>\tAdding L6A CT pops - ', cfg.L6ACTpops)
        netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells,file_format='pkl',cell_pop='L6A'))
        
        # --- testing using a new method that creates L6A and CC/IN pops together
        # netParams.popParams.update( BN.BuildNetwork.getL6ASubtypesPopTemplate(pops=cfg.L6ACTpops,cell_pop='L6A',  center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape=cfg.pop_shape,group_vecstims=True))
        netParams.popParams.update( BN.BuildNetwork.getL6PopTemplateCtCcIn(   pops=cfg.L6ACTpops,                 center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape=cfg.pop_shape,group_vecstims=True))
        # netParams.popParams.update( BN.BuildNetwork.getL6PopTemplateCtCcIn(   pops=cfg.L6ACTpops,                 center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape=cfg.pop_shape,group_vecstims=True))
        
        secTarget = 'basal'     # secList
        count_cells=0
        for pop_name in netParams.popParams.keys():
            if 'L6A' in pop_name:
                for CT_pop in cfg.L6ACTpops: 
                    if CT_pop in pop_name:
                        try:
                            try:    count_cells+=netParams.popParams[pop_name]['numCells']
                            except: count_cells+=len(netParams.popParams[pop_name]['cellsList'])
                            print('\t>>>\t', pop_name, ' count_cells: ', count_cells)
                        except:
                            print('\t>>>\tFailed to count cells in ', pop_name, ' pop')
                            continue
        CT_cells= count_cells
        # CT_cells = sum([netParams.popParams[ct_pop+'__pop']['numCells'] for ct_pop in cfg.L6ACTpops])
        
        CC_cells = round((CT_cells/0.563)*0.437)
        IN_cells = round((CT_cells+CC_cells)*0.082737031498839)
        if cfg.addL6Asubpops: 
            print('\n\t>>>\tAdding L6A subpops - ', cfg.L6Asubpops)
            for subpop in cfg.L6Asubpops:
                # print('\n\t>>>     Adding L6A subpop - ', subpop)
                netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells,file_format='pkl',cell_pop=subpop))
                netParams.popParams.update( BN.BuildNetwork.getL6PopTemplateCtCcIn(     pops=[subpop],center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape=cfg.pop_shape,))
            # netParams.popParams.update( BN.BuildNetwork.getL6ASubtypesPopTemplate(pops=['L6IN'],cell_pop='L6IN',center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape='cube',group_vecstims=True))
    
    elif cfg.addCTvirtual:
        try:
            for virtual_pop in cfg.ct_spk_files_dict.keys():
                print('\n\t>>>\tAdding virtual CT population - ', virtual_pop)

                if '.pkl' in cfg.ct_spk_files_dict[virtual_pop]:
                    print('Loading CT spike times from Pickle')
                    ct_spikes_dict = SampleData.LoadPickle(file_path=cfg.ct_spk_files_dict[virtual_pop])
                    ct_spk_dict = ct_spikes_dict['spkts']
                else:
                    print('Loading CT spike times from JSON')
                    with open(cfg.ct_spk_files_dict[virtual_pop], "r") as json_file: ct_spk_dict_=json.load(json_file)
                    print('Converting JSON file spkts into list')
                    ct_spk_dict = {int(c_ind):ct_spk_dict_[c_ind] for c_ind in ct_spk_dict_.keys()} # converting JSON dictionary of lists into list of lists

                if cfg.replace_spk_file is not None:
                    import copy
                    ct_spk_dict_=copy.deepcopy(ct_spk_dict)

                    print('Loading REPLACE spike times from Pickle')
                    angular_tuned_spikes = SampleData.LoadPickle(file_path=cfg.replace_spk_file)
                    replace_spk_dict = angular_tuned_spikes['spkts']

                    ct_keys = list(ct_spk_dict.keys())
                    replace_values = list(replace_spk_dict.values())
                    
                    # Calculate step size for dispersing replace_spk_dict values
                    step = len(ct_spk_dict) // len(replace_spk_dict)
                    
                    # Generate indices to spread replace values evenly across ct_spk_dict keys
                    replace_indices = np.linspace(0, len(ct_keys) - 1, len(replace_values), dtype=int)
                    
                    # Create a new dictionary with replaced values
                    new_dict = {key: ct_spk_dict[key] for key in ct_keys}  # Start with a copy of ct_spk_dict
                    
                    for i, index in enumerate(replace_indices): 
                        # print('Replacing spikes from cell ', i, ' | ', index, ' | ', i*(818/180), ' | ', i*(818/180)-index)
                        new_dict[ct_keys[index]] = replace_values[i]
                    
                    del ct_spk_dict
                    ct_spk_dict = copy.deepcopy(new_dict)

                    plotReplacedRaster=False
                    if plotReplacedRaster:
                        import matplotlib.pyplot as plt
                        plt.figure(figsize=(40,15))
                        plt.subplot(2,1,1)
                        for cell_spk_ind in ct_spk_dict_.keys(): 
                            plt.scatter(ct_spk_dict_[cell_spk_ind],[int(cell_spk_ind) for ind in range(len(ct_spk_dict_[cell_spk_ind]))], c='k',s=3)
                        plt.subplot(2,1,2)
                        for cell_spk_ind in ct_spk_dict.keys(): 
                            plt.scatter(ct_spk_dict[cell_spk_ind],[int(cell_spk_ind) for ind in range(len(ct_spk_dict[cell_spk_ind]))], c='k',s=3)
                        plt.savefig('replaced_ct_spikes_2.png')


                # --- Adding code to extend spike times
                pops_dict = BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='CTvirtual_'+virtual_pop,center_point=500,pop_flag='__pop',seed=100000,align_cells=cfg.align_virtual_cells,spkts_dict=ct_spk_dict)
                netParams.popParams.update(BN.BuildNetwork.extendSpikeTimes(pops_dict=pops_dict,
                                                                            max_spkt=10000,
                                                                            sim_duration_ms=cfg.duration,
                                                                            remove_spkt=[10001, 100001],
                                                                            skip_time=2000,
                                                                            ))

                # netParams.popParams.update( BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='CTvirtual_'+virtual_pop,center_point=500,pop_flag='__pop',seed=100000,align_cells=cfg.align_virtual_cells,spkts_dict=ct_spk_dict))
                if cfg.delayCTvirtual is not None:
                    print('\t>\tAdding ', cfg.delayCTvirtual,' ms shift in spike times')
                    if (type(cfg.delayCTvirtual) is int) or (type(cfg.delayCTvirtual) is float): 
                        for i in range(len(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'])):
                            netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes'] = list(map(lambda x: x + cfg.delayCTvirtual, netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes']))
                    delayCTvirtual=cfg.delayCTvirtual
                else: delayCTvirtual=0

                if cfg.ct_downTimes is not None:
                    try:
                        print('\t>\tRemoving selected spikes')
                        for i in range(len(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'])):
                            # print(i, sum(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes']))
                            spike_times = [spkt for spkt in netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes'] if not any(start <= spkt <= end for start, end in cfg.ct_downTimes)]
                            netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes']=spike_times
                            # print(i, sum(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes']))
                    except: print('spike removing failed')


                pop_mean_firing_rate = np.mean([len(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'][i]['spkTimes'])/((cfg.duration-delayCTvirtual)/1000) for i in range(len(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList']))])
                print('\t>\tMean Firing rate - ', virtual_pop, ' - ', str(pop_mean_firing_rate), ' Hz')
        except:        
            print('\n\t>>>\tAdding virtual CT population')
            with open(cfg.ct_spk_file, "r") as json_file: ct_spk_dict=json.load(json_file)
            netParams.popParams.update( BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='CTvirtual',center_point=500,pop_flag='__pop',seed=100000,align_cells=cfg.align_virtual_cells,spkts_dict=ct_spk_dict))
        cfg.conn_VPM_L6A=False
        cfg.conn_L6A_VPM=False
    else:
        if cfg.simplifyL6A: 
            netParams.cellParams.update(BN.BuildNetwork.getCellTemplate(    template='izhi',pops=['L6A'],cell_flag='__cell'))
            netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell'))
            secTarget = 'soma_0'    # single section
        else:
            netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells,file_format='pkl'))
            # netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells))
            netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape='cube'))
            secTarget = 'basal'     # secList
            # print('\t>>\tWarning: Add code to remove AXON segment from L6A cells')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- MLe cells and pops
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if cfg.includeBS:
    print('\n\t>>>\tAdding virtual Brainstem population')
    if '.pkl' in cfg.deflection_dataset_path:
        spikes_dict = SampleData.LoadPickle(file_path=cfg.deflection_dataset_path)
    else:
        spikes_dict={'spkts':{}}
        with open(cfg.deflection_dataset_path, "r") as json_file: spikes_dict['spkts']=json.load(json_file)

    if cfg.extend_spike_times:
        netParams.popParams.update(BN.BuildNetwork.extendSpikeTimes(
                                                                    BN.BuildNetwork.getMLePopTemplate_VecStim_cellsList(pop='MLe',step=2,center_point=cfg.center_point,pop_flag='__pop',seed=cfg.base_random_seed,align_cells=cfg.align_virtual_cells,spkts_dict=spikes_dict['spkts'],plotFig=False),
                                                                    max_spkt=10000,
                                                                    sim_duration_ms=cfg.duration,
                                                                    remove_spkt=[10001, 100001],
                                                                    skip_time=2000,
                                                                    ))
    else:
        netParams.popParams.update( BN.BuildNetwork.getMLePopTemplate_VecStim_cellsList(pop='MLe',step=2,center_point=cfg.center_point,pop_flag='__pop',seed=cfg.base_random_seed,align_cells=cfg.align_virtual_cells,spkts_dict=spikes_dict['spkts'],plotFig=False))

    if cfg.delayBS is not None:
        print('\t>\tAdding ', cfg.delayBS,' ms shift in spike times')
        if (type(cfg.delayBS) is int) or (type(cfg.delayBS) is float): 
            for i in range(len(netParams.popParams['MLe__pop']['cellsList'])):
                netParams.popParams['MLe__pop']['cellsList'][i]['spkTimes'] = list(map(lambda x: x + cfg.delayBS, netParams.popParams['MLe__pop']['cellsList'][i]['spkTimes']))
        delayBS=cfg.delayBS
    else: delayBS=0
    # Calculates the mean firing frequency of the stimulation dataset (frequency for None: ~22Hz) - 25Hz for (Iavarone,2023)
    mean_freq=np.mean([(len(netParams.popParams['MLe__pop']['cellsList'][cell_number]['spkTimes'])/((cfg.duration-delayBS)/1000)) for cell_number in range(len(netParams.popParams['MLe__pop']['cellsList']))])
    print('\t>\tMean Firing rate - MLe__pop - ', str(mean_freq), ' Hz')

plot_MLe_CT_spikingHist=False
if plot_MLe_CT_spikingHist:
    from PlotSpikingHist import PlotSpikingHist
    PlotSpikingHist.plot_hist(pop_name='MLe__pop',              netParams=netParams, spk_window=None, pop_color='k')
    PlotSpikingHist.plot_hist(pop_name='CTvirtual_uniform__pop',netParams=netParams, spk_window=None, pop_color='r')


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- VPM_ring population to drive the TRN_ring with baseline inputs (same dataset as the <None> deflection model because the firing rate is the same)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if cfg.includeBaselineDriver:
    baseline_dict = SampleData.LoadPickle(file_path=cfg.baseline_driver_path)
    baseline_spikes = {key:baseline_dict['spkts'][key] for key in range(25)} # --- reducing the number of noise sources to speed up computation
    print('\n\t>>>\tAdding virtual Baseline Driver population')
    if cfg.extend_spike_times:
        netParams.popParams.update(BN.BuildNetwork.extendSpikeTimes(
                                                                    BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='BaselineDriver',center_point=500,pop_flag='__pop',seed=100000,align_cells=True,spkts_dict=baseline_spikes),
                                                                    max_spkt=10000,
                                                                    sim_duration_ms=cfg.duration,
                                                                    remove_spkt=[10001, 100001],
                                                                    skip_time=2000,
                                                                    ))
    else:
        netParams.popParams.update( BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='BaselineDriver',center_point=500,pop_flag='__pop',seed=100000,align_cells=True,spkts_dict=baseline_spikes))
    # netParams.popParams.update( BN.BuildNetwork.getVirtualCT_Vecstim_cellsList(pop='BaselineDriver',center_point=500,pop_flag='__pop',seed=100000,align_cells=True,spkts_dict=baseline_spikes))

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Synaptic Mechanisms
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Dummy exc and inh synMechs for testing
netParams.synMechParams = BN.BuildNetwork.getSynMechParams()

# --- Selecting type of MOD file to be used
if   cfg.modType == 'Det':              modAMPANMDA = 'DetAMPANMDA';                modGABA     = 'DetGABAAB'              # S1 Deterministic  implementation of the BBP mod files
elif cfg.modType == 'Prob_S1':          modAMPANMDA = 'ProbAMPANMDA_EMS_S1';        modGABA     = 'ProbGABAAB_EMS_S1'      # S1 Probabilistic  implementation of the BBP mod files
elif cfg.modType == 'Prob_original':    modAMPANMDA = 'ProbAMPANMDA_EMS_original';  modGABA     = 'ProbGABAAB_EMS_original' # original MOD from BBP model
else:                                   modAMPANMDA = 'ProbAMPANMDA_EMS';           modGABA     = 'ProbGABAAB_EMS'         # Original Thalamus implementation of the BBP mod files
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
        
        # Post-creating modifications - from simulation_sonata.json
        if ('TRN' in pre_pop): 
            if 'VPM' in post_pop:   
                print('\t\tAdding BBP runtime modifications: \tmodifying TRN->VPM:   ', 'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag)
                edge_dict.update({'e_GABAA': -94.0, 'e_GABAB': -97.0, 'tau_d_GABAB': 77})
            else:                   
                print('\t\tAdding BBP runtime modifications: \tmodifying TRN->other: ', 'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag)
                edge_dict.update({'e_GABAA': -82.0, 'e_GABAB': -97.0, 'tau_d_GABAB': 77})
        
        netParams.synMechParams.update({'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag:edge_dict})

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Model adjustments for in vitro condition
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Changes extracellular Ca2+ concentration for all sections in the biophysical cell models
if cfg.cao_secs is not None:
    print('\t>>\tChanging extracellular Ca2+ concentration to ', str(cfg.cao_secs))
    for biophys_cell in netParams.cellParams.keys():
        for sec in netParams.cellParams[biophys_cell]['secs'].keys():
            if 'ions' in netParams.cellParams[biophys_cell]['secs'][sec].keys():
                if 'ca' in netParams.cellParams[biophys_cell]['secs'][sec]['ions'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['ions']['ca']['o'] = cfg.cao_secs

# --- Rescale USE parameter (probability of synapse activation)
if cfg.rescaleUSE is not None:
    print('\t>>\tRescaling synaptic USE to ', str(cfg.rescaleUSE))
    for mech in netParams.synMechParams.keys():
        try:    netParams.synMechParams[mech]['Use']*=cfg.rescaleUSE
        except: continue

# --- Modify NMDA ratio L6A->VPM
if cfg.modifyNMDAratio_L6A_VPM is not None:
    print('\t>>\tModifying NMDA ratio of L6A->VPM to ', str(cfg.modifyNMDAratio_L6A_VPM))
    for mech in netParams.synMechParams.keys():
        if mech == 'syn|L6A|VPM|exc':
            try:    netParams.synMechParams[mech]['NMDA_ratio']=cfg.modifyNMDAratio_L6A_VPM
            except: continue

# --- Modify NMDA ratio
if cfg.modifyNMDAratio_L6A_TRN is not None:
    print('\t>>\tModifying NMDA ratio of L6A->TRN to ', str(cfg.modifyNMDAratio_L6A_TRN))
    for mech in netParams.synMechParams.keys():
        if mech == 'syn|L6A|TRN|exc':
            try:    netParams.synMechParams[mech]['NMDA_ratio']=cfg.modifyNMDAratio_L6A_TRN
            except: continue

####################################################################################################################################################################################
# --- Modify currents for all Thalamus
####################################################################################################################################################################################

# --- Rescaling conductances

print('\t>>\tRescaling SK_E2 gmax by x', str(cfg.rescale_SK_E2))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'SK_E2' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['SK_E2']['gSK_E2bar']*=cfg.rescale_SK_E2

print('\t>>\tRescaling pas gmax by x', str(cfg.rescale_pas))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['g']*=cfg.rescale_pas

print('\t>>\tRescaling iT pcabar by x', str(cfg.rescale_iT))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_iT_Des98' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): 
                # print('iT before: ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['pcabar']))
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['pcabar']*=cfg.rescale_iT
                # print('iT after : ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['pcabar']))

print('\t>>\tRescaling iH gh_max by x', str(cfg.rescale_ih))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_ih_Bud97' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): 
                # print('iH before: ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['gh_max']))
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['gh_max']*=cfg.rescale_ih
                # print('iH after : ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['gh_max']))

print('\t>>\tRescaling iA gk_max by x', str(cfg.rescale_iA))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_iA' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): 
                # print('iA before: ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iA']['gk_max']))
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iA']['gk_max']*=cfg.rescale_iA
                # print('iA after : ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iA']['gk_max']))

print('\t>>\tRescaling iNap gNap_Et2bar by x', str(cfg.rescale_iNap))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_Nap_Et2' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): 
                # print('iNap before: ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Nap_Et2']['gNap_Et2bar']))
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Nap_Et2']['gNap_Et2bar']*=cfg.rescale_iNap
                # print('iNap after : ', str(netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Nap_Et2']['gNap_Et2bar']))

####################    IT shift - Thalamus   ###############
print('\t>>\tAdding a Thalamic shift of', str(cfg.add_iT_shift__Thalamus), '(mV)')
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            # if sec=='soma_0': print('\t\t- before:\t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])
            if 'TC_iT_Des98' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift']+=cfg.add_iT_shift__Thalamus
            # if sec=='soma_0': print('\t\t- after: \t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])

####################################################################################################################################################################################
# --- Modify currents by Thalamic pop
####################################################################################################################################################################################

print('\t>>\tRescaling VPM SK_E2 gmax by', str(cfg.rescale_SK_E2__VPM))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'SK_E2' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['SK_E2']['gSK_E2bar']*=cfg.rescale_SK_E2__VPM

print('\t>>\tRescaling VPM pas gmax by', str(cfg.rescale_pas__VPM))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['g']*=cfg.rescale_pas__VPM

print('\t>>\tRescaling TRN SK_E2 gmax by', str(cfg.rescale_SK_E2__TRN))
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'SK_E2' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['SK_E2']['gSK_E2bar']*=cfg.rescale_SK_E2__TRN

print('\t>>\tRescaling TRN pas gmax by', str(cfg.rescale_pas__TRN))
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['g']*=cfg.rescale_pas__TRN

####################    IH gh_max   ####################
print('\t>>\tRescaling VPM ih gmax by', str(cfg.rescale_ih_gmax__VPM))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_ih_Bud97' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['gh_max']*=cfg.rescale_ih_gmax__VPM

print('\t>>\tRescaling TRN ih gmax by', str(cfg.rescale_ih_gmax__TRN))
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'TC_ih_Bud97' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['gh_max']*=cfg.rescale_ih_gmax__TRN

####################    IH shift   ####################
print('\t>>\tAdding a VPM iH shift of', str(cfg.add_ih_shift__VPM), '(mV)')
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            # if sec=='soma_0': print('\t\t- before:\t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['shift'])
            if 'TC_ih_Bud97' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['shift']+=cfg.add_ih_shift__VPM
            # if sec=='soma_0': print('\t\t- after: \t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_ih_Bud97']['shift'])

#############################################################
####################    IT shift - VPM   ####################
print('\t>>\tAdding a VPM iT shift of', str(cfg.add_iT_shift__VPM), '(mV)')
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            # if sec=='soma_0': print('\t\t- before:\t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])
            if 'TC_iT_Des98' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift']+=cfg.add_iT_shift__VPM
            # if sec=='soma_0': print('\t\t- after: \t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])

####################    IT shift - TRN   ####################
print('\t>>\tAdding a TRN iT shift of', str(cfg.add_iT_shift__TRN), '(mV)')
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            # if sec=='soma_0': print('\t\t- before:\t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])
            if 'TC_iT_Des98' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift']+=cfg.add_iT_shift__TRN
            # if sec=='soma_0': print('\t\t- after: \t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_iT_Des98']['shift'])

#############################################################
print('\t>>\tChanging VPM pas e to ', str(cfg.modify_pas_e__VPM))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['e']=cfg.modify_pas_e__VPM
print('\t>>\tChanging VPM pas g to ', str(cfg.modify_pas_g__VPM))
for cell_name in netParams.cellParams.keys():
    if ('VPL_TC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['g']=cfg.modify_pas_g__VPM

print('\t>>\tChanging TRN pas e to ', str(cfg.modify_pas_e__TRN))
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['e']=cfg.modify_pas_e__TRN
print('\t>>\tChanging TRN pas g to ', str(cfg.modify_pas_g__TRN))
for cell_name in netParams.cellParams.keys():
    if ('Rt_RC' in cell_name):
        for sec in netParams.cellParams[cell_name]['secs'].keys():
            if 'pas' in netParams.cellParams[cell_name]['secs'][sec]['mechs'].keys(): netParams.cellParams[cell_name]['secs'][sec]['mechs']['pas']['g']=cfg.modify_pas_g__TRN

if cfg.add_KLeak__VPM > 0:
    print('\t>>\tAdding a KLeak current into VPM neurons ', str(cfg.add_KLeak__VPM), ' * (1.0e-5)')
    # pas_g=3.4702549429081374e-05
    for cell_name in netParams.cellParams.keys():
        if ('VPL_TC' in cell_name):
            for sec in netParams.cellParams[cell_name]['secs'].keys():
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Kleak']={'g':(1.0e-5)*cfg.add_KLeak__VPM}

if cfg.add_KLeak__TRN > 0:
    print('\t>>\tAdding a KLeak current into TRN neurons ', str(cfg.add_KLeak__TRN), ' * (1.0e-5)')
    # pas_g=8.617446501142974e-05
    for cell_name in netParams.cellParams.keys():
        if ('Rt_RC' in cell_name):
            for sec in netParams.cellParams[cell_name]['secs'].keys():
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Kleak']={'g':(1.0e-5)*cfg.add_KLeak__TRN}

if cfg.add_KLeak > 0:
    print('\t>>\tAdding a KLeak current into VPM and TRN neurons ', str(cfg.add_KLeak), ' * (1.0e-5)')
    # pas_g=3.4702549429081374e-05
    for cell_name in netParams.cellParams.keys():
        if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name):
            for sec in netParams.cellParams[cell_name]['secs'].keys():
                # if sec=='soma_0': print('\t\t- before:\t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Kleak']['g'])
                netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Kleak']={'g':(1.0e-5)*cfg.add_KLeak}
                # if sec=='soma_0': print('\t\t- after: \t', cell_name, sec, netParams.cellParams[cell_name]['secs'][sec]['mechs']['TC_Kleak']['g'])

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Connectivity
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# syns per conn - obs: conn_data is organized as conn_data[post_pop]['chem'][pre_pop][synsPerConn/convergence][MEAN,STD]

# --- MLe to VPM connections

#       Obs: Projections are organized so that:
#               the top    of the MLe projects to the top    of VPM (absolute top = 0, pia surface)
#               the bottom of the MLe projects to the bottom of VPM
#            Based on the data from (Timofeeva et al., 2003)

if cfg.conn_MLe_VPM:

    # using laminar projection based on y-axis relative positon to sharpen MLe->VPM inputs
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.laminarProjection(pre_pop='MLe',post_pop='VPM',conn_prob=1,
                                                                    #  y_thresh=0.0065, # test 2024_08_26
                                                                     y_thresh=cfg.MLe_VPM__y_tresh, # test 2024_08_26
                                                                     center_point=500)

    # # --- New conn method: distanceBasedProbability_1D_exponential
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.MLe_VPM__pre_pop,
                                                post_pop        = cfg.MLe_VPM__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.MLe_VPM__syn_mech,
                                                syns_per_conn   = cfg.MLe_VPM__syns_per_conn,
                                                conn_type       = cfg.MLe_VPM__conn_type,
                                                weight          = cfg.MLe_VPM__weight,
                                                target_secs     = cfg.MLe_VPM__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_L6A_VPM:
    if cfg.L6A_VPM__pre_pop == 'all_L6A_subpops':   L6A_VPM__pre_pops=[key.split('__')[0] for key in netParams.popParams.keys() if 'L6A' in key]
    else:                                           L6A_VPM__pre_pops = [cfg.L6A_VPM__pre_pop]
    print('CONTINUE DEBUGGING HERE!!! L6A->TC CONNS')
    conv_dict_L6A_VPM={'sparse': 113.79559748427673, 'silent': 75.27358490566037}
    for L6A_VPM__pre_pop in L6A_VPM__pre_pops:
        if ('sparse' in L6A_VPM__pre_pop) or ('silent' in L6A_VPM__pre_pop):
            if      'sparse' in L6A_VPM__pre_pop:   conv_val = conv_dict_L6A_VPM['sparse']
            elif    'silent' in L6A_VPM__pre_pop:   conv_val = conv_dict_L6A_VPM['silent']
            else:                                   conv_val = 1
            conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = L6A_VPM__pre_pop,
                                                        post_pop        = cfg.L6A_VPM__post_pop,
                                                        conn_method     = 'convergence',
                                                        conn_rule       = conv_val,
                                                        syn_mech        = cfg.L6A_VPM__syn_mech,
                                                        syns_per_conn   = cfg.L6A_VPM__syns_per_conn,
                                                        conn_type       = cfg.L6A_VPM__conn_type,
                                                        weight          = cfg.L6A_VPM__weight,
                                                        target_secs     = cfg.L6A_VPM__target_secs)
            print('\t>>\t Convergence L6A->TC',L6A_VPM__pre_pop, conv_val) 
        else:
            print('Connectivity not tested - must be retuned')
            conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop   = 'L6A',     post_pop     = cfg.L6A_VPM__post_pop,
                                                                                                    pre_pop_axis = cfg.L6A_VPM__pre_pop_axis,post_pop_axis= cfg.L6A_VPM__post_pop_axis,
                                                                                                    center_point = cfg.center_point,
                                                                                                    conn_prob=1,    baseline_prob=0,
                                                                                                    k=-2,            shift=0.13,          scale=25,
                                                                                                    plot_conn=True
                                                                                                )

            conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = L6A_VPM__pre_pop,
                                                        post_pop        = cfg.L6A_VPM__post_pop,
                                                        conn_method     = conn_method_1D,
                                                        conn_rule       = conn_rule_1D,
                                                        syn_mech        = cfg.L6A_VPM__syn_mech,
                                                        syns_per_conn   = cfg.L6A_VPM__syns_per_conn,
                                                        conn_type       = cfg.L6A_VPM__conn_type,
                                                        weight          = cfg.L6A_VPM__weight,
                                                        target_secs     = cfg.L6A_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)

if cfg.conn_L6A_TRN:
    if cfg.L6A_TRN__pre_pop == 'all_L6A_subpops':   L6A_TRN__pre_pops=[key.split('__')[0] for key in netParams.popParams.keys() if 'L6A' in key]
    else:                                           L6A_TRN__pre_pops = [cfg.L6A_TRN__pre_pop]
    print('CONTINUE DEBUGGING HERE!!! L6A->TRN CONNS')
    conv_dict_L6A_TRN={'sparse': 53.2970297029703, 'silent': 35.23762376237624}
    for L6A_TRN__pre_pop in L6A_TRN__pre_pops:
        if ('sparse' or 'silent') in L6A_TRN__pre_pop:
            if      'sparse' in L6A_VPM__pre_pop:   conv_val = conv_dict_L6A_TRN['sparse']
            elif    'silent' in L6A_VPM__pre_pop:   conv_val = conv_dict_L6A_TRN['silent']
            else:                                   conv_val = 1
            conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = L6A_TRN__pre_pop,
                                                        post_pop        = cfg.L6A_TRN__post_pop,
                                                        conn_method     = 'convergence',
                                                        conn_rule       = conv_val,
                                                        syn_mech        = cfg.L6A_TRN__syn_mech,
                                                        syns_per_conn   = cfg.L6A_TRN__syns_per_conn,
                                                        conn_type       = cfg.L6A_TRN__conn_type,
                                                        weight          = cfg.L6A_TRN__weight,
                                                        target_secs     = cfg.L6A_TRN__target_secs)
        else:
            print('Connectivity not tested - must be retuned')
            conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop   = 'L6A',     post_pop     = cfg.L6A_TRN__post_pop,
                                                                                                    pre_pop_axis = cfg.L6A_TRN__pre_pop_axis,post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                                    center_point = cfg.center_point,
                                                                                                    conn_prob=1,    baseline_prob=0,
                                                                                                    k=-2,            shift=0.13,          scale=25,
                                                                                                    plot_conn=True
                                                                                                )

            conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = L6A_TRN__pre_pop,
                                                        post_pop        = cfg.L6A_TRN__post_pop,
                                                        conn_method     = conn_method_1D,
                                                        conn_rule       = conn_rule_1D,
                                                        syn_mech        = cfg.L6A_TRN__syn_mech,
                                                        syns_per_conn   = cfg.L6A_TRN__syns_per_conn,
                                                        conn_type       = cfg.L6A_TRN__conn_type,
                                                        weight          = cfg.L6A_TRN__weight,
                                                        target_secs     = cfg.L6A_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)

if cfg.conn_TRN_VPM:
    print('\t --- Connectivity schema = ', cfg.dump_cell_properties)
    if '__TRN_fb_uniform' in cfg.dump_cell_properties:
        conn_method_1D='convergence'    
        conn_rule_1D = 43 # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 35 # 35 * 5 = 175 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)
    # Topological feedback (off-center)
    elif ('__TRN_fb_topolog' in cfg.dump_cell_properties) or ('__TRN_fb_open' in cfg.dump_cell_properties):
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)
    
    # Topological feedback (on-center)
    elif '__TRN_fb_closed' in cfg.dump_cell_properties:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)

    # Topological feedback (on-center and off-center)
    elif '__TRN_fb_mixed' in cfg.dump_cell_properties:
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0,    baseline_prob=0.5,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=0.5,  baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        conn_dict_1D__ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D__rename = {k+'2':conn_dict_1D__[k] for k in conn_dict_1D__.keys()}
        netParams.connParams.update(conn_dict_1D__rename)
    
    elif '__TRN_fb_remixed' in cfg.dump_cell_properties:
        
        print('\n\n ==== Running REMIXED TRN FB rule ==== \n\n')

        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        # 2025_01_09
        conn_rule_1D = int(43*topology_ratio[1]) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = int(43*(2.5/5)) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 35 # 35 * 5 = 175 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0*topology_ratio[0],  baseline_prob=0, # 2025_01_09
                                                                                                # conn_prob=1.0*(2.5/5),  baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        conn_dict_1D__ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D__rename = {k+'2':conn_dict_1D__[k] for k in conn_dict_1D__.keys()}
        netParams.connParams.update(conn_dict_1D__rename)
    
    elif '__TRN_fb_openRemixed' in cfg.dump_cell_properties:
        
        print('\n\n ==== Running OPEN REMIXED TRN FB rule ==== \n\n')

        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        # 2025_01_09
        conn_rule_1D = int(43*topology_ratio[1]) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = int(43*(2.5/5)) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 35 # 35 * 5 = 175 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=1-topology_ratio[0], # complement, because the probability equation is subtracted from 1 afterwards
                                                                                                # conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D__ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D__rename = {k+'2':conn_dict_1D__[k] for k in conn_dict_1D__.keys()}
        netParams.connParams.update(conn_dict_1D__rename)
    
    else:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0.48,
                                                                                                k=-2,           shift=0.13,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D)

    print('\t --- Connectivity schema = ', cfg.dump_cell_properties)
    if '__TRN_fb_uniform' in cfg.dump_cell_properties:
        conn_method_1D='convergence'    
        conn_rule_1D = 43 # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 9 # 9 * 5 = 45 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

    # Topological feedback (off-center)
    elif '__TRN_fb_topolog' in cfg.dump_cell_properties:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

    # Topological feedback (on-center)
    elif '__TRN_fb_closed' in cfg.dump_cell_properties:    
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'TRN_ring',                  post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

    # Topological feedback (on-center and off-center)
    elif '__TRN_fb_mixed' in cfg.dump_cell_properties:

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'TRN_ring',                  post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0,    baseline_prob=0.5,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)


        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'TRN_ring',                  post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=0.5,  baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)

    elif '__TRN_fb_remixed' in cfg.dump_cell_properties:
        
        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        # 2025_01_09
        conn_rule_1D = int(43*topology_ratio[1]) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = int(43*(2.5/5)) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 9 # 9 * 5 = 45 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'TRN_ring',                  post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0*topology_ratio[0],  baseline_prob=0, # 2025_01_09
                                                                                                # conn_prob=1.0*(2.5/5),  baseline_prob=0,
                                                                                                k=-2,           shift=0.225,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=False,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)

    elif '__TRN_fb_openRemixed' in cfg.dump_cell_properties:
        
        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        # 2025_01_09
        conn_rule_1D = int(43*topology_ratio[1]) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = int(43*(2.5/5)) # 43 * 5 = 215 conns received per cell - total 215 (215|215 ~= 50|50 ratio)
        # conn_rule_1D = 9 # 9 * 5 = 45 conns received per cell - total 220 (175|45 ~= 80|20 ratio)
        # conn_rule_1D = 22 # 22 * 5 = 110 conns received per cell - total 220 (110|110 ~= 50|50 ratio)
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.TRN_VPM__pre_pop,        post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=1-topology_ratio[0], # complement, because the probability equation is subtracted from 1 afterwards
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)

    else:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'TRN_ring',                  post_pop     = cfg.TRN_VPM__post_pop,
                                                                                                pre_pop_axis = cfg.TRN_VPM__pre_pop_axis,   post_pop_axis= cfg.TRN_VPM__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0.48,
                                                                                                k=-2,           shift=0.13,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )    
        print('\t>>\tconnecting adjacent TRN sectors ')
        conn_dict_1D_ = BN.BuildNetwork.getConnDict( pre_pop        = 'TRN_ring',
                                                    post_pop        = cfg.TRN_VPM__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.TRN_VPM__syn_mech,
                                                    syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                    conn_type       = cfg.TRN_VPM__conn_type,
                                                    weight          = cfg.TRN_VPM__weight,
                                                    target_secs     = cfg.TRN_VPM__target_secs)
        netParams.connParams.update(conn_dict_1D_)

if cfg.conn_TRN_TRN:
    if cfg.addTRNadjacentpops:
        print('\t>>\tIncluding adjacent TRN<->TRN connections - Chemical')
        if cfg.addTRNadjacentpops == 'biophysical':
            print('\t> \t\tAll TRN -> All TRN - biophisical adjacent pops')
            pre_pops__TRN_conn      = [pop for pop in netParams.popParams.keys() if cfg.TRN_TRN__pre_pop  in pop]
            post_pops__TRN_conn     = [pop for pop in netParams.popParams.keys() if cfg.TRN_TRN__post_pop in pop]
            pre_pops_flag           = cfg.TRN_TRN__pre_pop+ '@biophysical'
            post_pops_flag          = cfg.TRN_TRN__post_pop+'@biophysical'
        else:
            print('\t> \t\tAll TRN -> main TRN - vecstim adjacent pops')
            pre_pops__TRN_conn      = [pop for pop in netParams.popParams.keys() if cfg.TRN_TRN__pre_pop  in pop]
            post_pops__TRN_conn     = ['TRN__pop']
            pre_pops_flag           = cfg.TRN_TRN__pre_pop+ '@vecstim'
            post_pops_flag          = cfg.TRN_TRN__post_pop+'@biophysical'
        
        # --- testing topographical TRN interconnectivity - with minimum distance
        conn_method_3DDist, conn_rule_3DDist = BN.BuildNetwork.distanceBasedProbability_3DDist_minDist( pre_pop      = cfg.TRN_TRN__pre_pop,     post_pop     = cfg.TRN_TRN__post_pop,
                                                                                                        conn_prob    = cfg.TRN_TRN__conn_params_dict['conn_prob'],   
                                                                                                        decay_factor = cfg.TRN_TRN__conn_params_dict['decay_factor'],
                                                                                                        baseline_prob= cfg.TRN_TRN__conn_params_dict['baseline_prob'], 
                                                                                                        min_dist = 25,
                                                                                                        plot_conn=True
                                                                                                        )

        conn_dict_1D = BN.BuildNetwork.getConnDictCellType(  
                                                            pre_cellType    = 'Rt_RC',
                                                            post_cellType   = 'Rt_RC',
                                                            # post_pops       = [pop for pop in netParams.popParams.keys() if cfg.TRN_TRN__post_pop in pop],
                                                            pre_pops_flag   = pre_pops_flag,
                                                            post_pops_flag  = post_pops_flag,
                                                            # post_pops_flag  = cfg.TRN_TRN__post_pop+'@allPops',
                                                            conn_method     = conn_method_3DDist,
                                                            conn_rule       = conn_rule_3DDist,
                                                            syn_mech        = cfg.TRN_TRN__syn_mech,
                                                            syns_per_conn   = cfg.TRN_TRN__syns_per_conn,
                                                            conn_type       = cfg.TRN_TRN__conn_type,
                                                            weight          = cfg.TRN_TRN__weight,
                                                            target_secs     = cfg.TRN_TRN__target_secs)

    else:
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_TRN__pre_pop,
                                                    post_pop        = cfg.TRN_TRN__post_pop,
                                                    conn_method     = cfg.TRN_TRN__conn_method,
                                                    conn_rule       = cfg.TRN_TRN__conn_rule,
                                                    syn_mech        = cfg.TRN_TRN__syn_mech,
                                                    syns_per_conn   = cfg.TRN_TRN__syns_per_conn,
                                                    conn_type       = cfg.TRN_TRN__conn_type,
                                                    weight          = cfg.TRN_TRN__weight,
                                                    target_secs     = cfg.TRN_TRN__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_TRNe_TRNe_connList:
    # Connectivity parameters
    if cfg.singleCellPops:
        gap_preSec='soma_0'
        gap_preLoc=0.5
        gap_sec='soma_0'
        gap_loc=0.5
        gap_connList_relative=[[0,1]]
        # gap_connList_relative=[[0,1],[1,0]]
    else:
        # Read the JSON file and load it into a Python list
        if cfg.max_gaps_perCell is not None:    gap_conn_template='../conn/gap_connecivity/gap_connecivity_template_testingClass_unidirectional_maxGaps.json' # netpyne handles the creation of bidirectional connectivity
        else:                                   gap_conn_template='../conn/gap_connecivity/gap_connecivity_template_testingClass_unidirectional.json' # netpyne handles the creation of bidirectional connectivity
        
        print('Loading Gap template: ', gap_conn_template)

        with open(gap_conn_template, "r") as json_file: final_list = json.load(json_file)
        # final_list=final_list[0:1000]        
        
        gap_connList=[];gap_preGid=[];gap_postGid=[];gap_preSec=[];gap_preLoc=[];gap_sec=[];gap_loc=[]
        for ((post_gid,post_sec,post_pt3d_ind,post_loc),(pre_gid,pre_sec,pre_pt3d_ind,pre_loc)) in final_list:
            gap_connList.append([pre_gid,post_gid])
            if cfg.simplify_gap_conns:
                gap_preSec='soma_0'
                gap_preLoc=0.5
                gap_sec='soma_0'
                gap_loc=0.5
            else:
                gap_preSec.append(pre_sec)
                gap_preLoc.append(pre_loc)
                gap_sec.append(post_sec)
                gap_loc.append(post_loc)
        # --- Creates a list of relative GIDs, since the connParams orders the cells from 0 to n, regardless of the cell GID.
        relative_reference=min(min(gap_connList))
        gap_connList_relative = [[int(pre_gid-relative_reference),int(post_gid-relative_reference)] for [pre_gid,post_gid] in gap_connList]

    if cfg.individual_gap_conns:
        for ind, [pre_gid_r,post_gid_r] in enumerate(gap_connList_relative):
            netParams.connParams['conn@'+str(ind)+'|TRN|TRN|gap'] = {
                                                            'preConds':     {'pop': ['TRN__pop','TRN_ring__pop',]},     # pre pops
                                                            'postConds':    {'pop': ['TRN__pop','TRN_ring__pop',]},     # post pops
                                                            'synMech':      'esyn',                                   # gap junction synaptic mechanism
                                                            # 'probability' :    0.6,                               # (pre,post) cell gids
                                                            'connList' :    [gap_connList_relative[ind]],               # (pre,post) cell gids
                                                            'preSec':       gap_preSec[ind],                            # (pre)  cell sec
                                                            # 'preLoc':       0.3444,                            # (pre)  cell loc
                                                            'preLoc':       gap_preLoc[ind],                            # (pre)  cell loc
                                                            'sec':          gap_sec[ind],                               # (post) cell sec
                                                            # 'loc':          0.1222,                               # (post) cell loc
                                                            'loc':          gap_loc[ind],                               # (post) cell loc
                                                            'weight': cfg.TRNe_TRNe__weight,                                              # weight of each connection
                                                            # 'weight': 1.0,                                              # weight of each connection
                                                            'synsPerConn': 1,                                           # single syn per conn
                                                            # 'delay': 2,               # delay
                                                            'delay': '0.001+dist_3D/propVelocity',               # delay
                                                            # 'delay': ['defaultDelay+dist_3D/propVelocity' for i in range(len(gap_connList))],               # delay
                                                            # 'delay': [2 for i in range(len(gap_connList))],               # delay
                                                            
                                                            # 'distributeSynsUniformly': False,
                                                            # 'connRandomSecFromList': False, 
                                                            
                                                            }
    else:
        netParams.connParams['conn|TRN|TRN|gap'] = {
                                                        'preConds':     {'pop': ['TRN__pop','TRN_ring__pop',]},     # pre pops
                                                        'postConds':    {'pop': ['TRN__pop','TRN_ring__pop',]},     # post pops
                                                        'synMech':      'esyn',                                   # gap junction synaptic mechanism
                                                        # 'probability' :    0.6,                               # (pre,post) cell gids
                                                        'connList' :    gap_connList_relative,                               # (pre,post) cell gids
                                                        'preSec':       gap_preSec,                                 # (pre)  cell sec
                                                        'preLoc':       gap_preLoc,                                 # (pre)  cell loc
                                                        'sec':          gap_sec,                                    # (post) cell sec
                                                        'loc':          gap_loc,                                    # (post) cell loc
                                                        'weight': cfg.TRNe_TRNe__weight,                                              # weight of each connection
                                                        # 'weight': 1.0,                                              # weight of each connection
                                                        'synsPerConn': 1,                                           # single syn per conn
                                                        # 'delay': 2,               # delay
                                                        'delay': '0.001+dist_3D/propVelocity',               # delay
                                                        # 'delay': ['defaultDelay+dist_3D/propVelocity' for i in range(len(gap_connList))],               # delay
                                                        # 'delay': [2 for i in range(len(gap_connList))],               # delay
                                                        }

if cfg.conn_VPM_TRN:
    print('\t --- Connectivity schema = ', cfg.dump_cell_properties)
    if '__VPM_ff_uniform' in cfg.dump_cell_properties:
        conn_method_1D='convergence'    
        conn_rule_1D = 83 # 83 * 2 = 166 conns received per cell
        # conn_rule_1D = 36 # 36 * 2 = 74 conns received per cell
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
    
    elif '__VPM_ff_closed' in cfg.dump_cell_properties:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop=cfg.VPM_TRN__post_pop,
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
        
    elif '__VPM_ff_mixed' in cfg.dump_cell_properties:

        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]
        # approximation
        if   topology_ratio[1] == 0.2:  baseline_prob_equivalent = 0.96
        elif topology_ratio[1] == 0.3:  baseline_prob_equivalent = 0.945
        elif topology_ratio[1] == 0.4:  baseline_prob_equivalent = 0.93
        else:                           baseline_prob_equivalent = 1-(0.2*(topology_ratio[1]))

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.VPM_TRN__pre_pop,        post_pop     = cfg.VPM_TRN__post_pop,
                                                                                                pre_pop_axis = cfg.VPM_TRN__pre_pop_axis,   post_pop_axis= cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0,    baseline_prob=baseline_prob_equivalent,
                                                                                                # conn_prob=1.0,    baseline_prob=0.75,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
        
        #####

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop=cfg.VPM_TRN__post_pop,
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=1*topology_ratio[0],    baseline_prob=0,
                                                                                                # conn_prob=1*(4/5),    baseline_prob=0,
                                                                                                # conn_prob=1*(3/5),    baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)
    elif '__VPM_ff_remixed' in cfg.dump_cell_properties:
        
        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        conn_rule_1D = int(83*topology_ratio[1]) # 83 * 2 = 166 conns received per cell
        # conn_rule_1D = int(83*(2/5)) # 83 * 2 = 166 conns received per cell
        # conn_rule_1D = 36 # 36 * 2 = 74 conns received per cell
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)

        #####

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop=cfg.VPM_TRN__post_pop,
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=1*topology_ratio[0],    baseline_prob=0,
                                                                                                # conn_prob=1*(4/5),    baseline_prob=0,
                                                                                                # conn_prob=1*(3/5),    baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)
        
    else:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop=cfg.VPM_TRN__post_pop,
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=1,    baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = cfg.VPM_TRN__post_pop,
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
    
    print('\t --- Connectivity schema = ', cfg.dump_cell_properties)
    if '__VPM_ff_uniform' in cfg.dump_cell_properties:
        conn_method_1D='convergence'    
        conn_rule_1D = 20 # 20 * 2 = 40 conns received per cell (25% of main TRN)
        ''' adding VPM->TRN_ring projection to test '''
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
    elif '__VPM_ff_closed' in cfg.dump_cell_properties:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop='TRN_ring',
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=0.25,  baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)
    elif '__VPM_ff_mixed' in cfg.dump_cell_properties:
        
        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]
        # approximation
        if   topology_ratio[1] == 0.2:  baseline_prob_equivalent_ring = 0.96
        elif topology_ratio[1] == 0.3:  baseline_prob_equivalent_ring = 0.945
        elif topology_ratio[1] == 0.4:  baseline_prob_equivalent_ring = 0.93
        else:                           baseline_prob_equivalent_ring = 1-(0.2*(topology_ratio[1]))

        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = cfg.VPM_TRN__pre_pop,        post_pop     = 'TRN_ring',
                                                                                                pre_pop_axis = cfg.VPM_TRN__pre_pop_axis,   post_pop_axis= cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point = cfg.center_point,
                                                                                                conn_prob=1.0,    baseline_prob=baseline_prob_equivalent_ring,
                                                                                                # conn_prob=1.0,    baseline_prob=0.93,
                                                                                                k=-2,           shift=0.37,          scale=25,
                                                                                                plot_conn=True,
                                                                                                inverseDecay=True,
                                                                                            )
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)

        ####
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop='TRN_ring',
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=0.25*topology_ratio[0],  baseline_prob=0,
                                                                                                # conn_prob=0.25*(3/5),  baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)
    elif '__VPM_ff_remixed' in cfg.dump_cell_properties:

        # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
        topology_ratio = [
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
            int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
        ]

        conn_method_1D='convergence'    
        conn_rule_1D = int(20*topology_ratio[1]) # 20 * 2 = 40 conns received per cell (25% of main TRN)
        # conn_rule_1D = int(20*(2/5)) # 20 * 2 = 40 conns received per cell (25% of main TRN)
        ''' adding VPM->TRN_ring projection to test '''
        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)

        ####
        
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop='TRN_ring',
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=0.25*topology_ratio[0],  baseline_prob=0,
                                                                                                # conn_prob=0.25*(3/5),  baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D___ = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        conn_dict_1D___rename = {k+'2':conn_dict_1D___[k] for k in conn_dict_1D___.keys()}
        netParams.connParams.update(conn_dict_1D___rename)
    else:
        conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop=cfg.VPM_TRN__pre_pop,           post_pop='TRN_ring',
                                                                                                pre_pop_axis=cfg.VPM_TRN__pre_pop_axis, post_pop_axis=cfg.VPM_TRN__post_pop_axis,
                                                                                                center_point=cfg.center_point,
                                                                                                conn_prob=0.25,  baseline_prob=0,
                                                                                                k=-2,            shift=0.13,          scale=25,
                                                                                                # k=cfg.VPM_TRN__conn_params_dict['k'],            # Add this line if k is in conn_params_dict
                                                                                                # shift=cfg.VPM_TRN__conn_params_dict['shift'],    # Add this line if shift is in conn_params_dict
                                                                                                # scale=cfg.VPM_TRN__conn_params_dict['scale'],    # Add this line if scale is in conn_params_dict
                                                                                                plot_conn=True
                                                                                            )

        conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                    post_pop        = 'TRN_ring',
                                                    conn_method     = conn_method_1D,
                                                    conn_rule       = conn_rule_1D,
                                                    syn_mech        = cfg.VPM_TRN__syn_mech,
                                                    syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                    conn_type       = cfg.VPM_TRN__conn_type,
                                                    weight          = cfg.VPM_TRN__weight,
                                                    target_secs     = cfg.VPM_TRN__target_secs)
        netParams.connParams.update(conn_dict_1D)


if cfg.conn_VPM_L6A:

    print('Replace by distanceBasedProbability_1D_sigmoid if detailed L6A gets added')
    # --- TEST: 2024_01_11 - converting projections so that post synaptic pop receives angular projections across the circunference of the network
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential( pre_pop      = cfg.VPM_L6A__pre_pop,     post_pop     = cfg.VPM_L6A__post_pop,
                                                                                            pre_pop_axis = cfg.VPM_L6A__pre_pop_axis,post_pop_axis= 'theta',
                                                                                            center_point = cfg.center_point,
                                                                                            conn_prob    = cfg.VPM_L6A__conn_params_dict['conn_prob'],   
                                                                                            decay_factor = cfg.VPM_L6A__conn_params_dict['decay_factor'],
                                                                                            baseline_prob= cfg.VPM_L6A__conn_params_dict['baseline_prob'],
                                                                                            plot_conn=True,
                                                                                            )

    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_L6A__pre_pop,
                                                post_pop        = cfg.VPM_L6A__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.VPM_L6A__syn_mech,
                                                syns_per_conn   = cfg.VPM_L6A__TC_syns_per_conn,
                                                conn_type       = cfg.VPM_L6A__conn_type,
                                                weight          = cfg.VPM_L6A__weight,
                                                target_secs     = secTarget)
    netParams.connParams.update(conn_dict_1D)
    
    # 133*666/
    # (318*9)
    # see (Meyer 2010)

if cfg.addCTvirtual:

    ct_cell_num = sum([len(netParams.popParams['CTvirtual_'+pop+'__pop']['cellsList']) for pop in cfg.ct_spk_files_dict.keys()])

    '''
    silent 235  sparse 356  suppressed 92   activated 135
    Assuming (silent = L6B) and (sparse = L6A and L6B) -> 92+135+(356/2)
    pop_ratios={'L6A':{ 'L6A_activated':    0.165,
                        'L6A_suppressed':   0.113,
                        'L6A_sparse':       0.435/2,
                        # 'L6A_sparse':       0.435,
                        # 'L6A_silent':       0.287,
                        },}

    New values
    silent      0  
    sparse      818*(0.435/2)/(0.165+0.113+(0.435/2))   ~= 359
    suppressed  818*(0.113)/(0.165+0.113+(0.435/2))     ~= 187
    activated   818*(0.165)/(0.165+0.113+(0.435/2))     ~= 272
    '''

    for virtual_pop in cfg.ct_spk_files_dict.keys():
        pop_cell_num            = len(netParams.popParams['CTvirtual_'+virtual_pop+'__pop']['cellsList'])
        
        VPM_estimated_syns_per_conn = 15    # (Robson, 1983) - The paper shows that a single CT axon has multiple branches, and that around 15% of those branches make contact in a single TC cell
        TRN_estimated_syns_per_conn = 5     # (Robson, 1983) - The paper shows that a single CT axon has multiple branches, and that around 15% of those branches make contact in a single TC cell
        
        # convergence full
        convergence_VPM         = cfg.conn_data['VPM']['chem']['L6A']['synsPerConn'][0]*cfg.conn_data['VPM']['chem']['L6A']['convergence'][0]/VPM_estimated_syns_per_conn
        convergence_TRN         = cfg.conn_data['TRN']['chem']['L6A']['synsPerConn'][0]*cfg.conn_data['TRN']['chem']['L6A']['convergence'][0]/TRN_estimated_syns_per_conn
        # # convergence reduced
        # convergence_VPM         = 278 # from (MLe convergence * (L6A/VPM / MLe/VPM)) = 58.77167630057804 * 4.728101509424807  = 277.8784514281892
        # convergence_TRN         = 134 # from (convergence_VPM / (L6A/VPM / L6A/TRN)) = 277.8784514281892 * 2.0754747354032723 = 133.8866943009048
        # # convergence_TRN         = 217 # from (MLe convergence * (TRN/VPM / MLe/VPM)) = 58.77167630057804 * 3.6931049237127365 = 217.0499671205159

        scaled_convergence_VPM  = int(convergence_VPM*(pop_cell_num/ct_cell_num))
        scaled_convergence_TRN  = int(convergence_TRN*(pop_cell_num/ct_cell_num))

        if '__L6A_fb_closed' in cfg.dump_cell_properties:
            print('\t --- CT Connectivity schema = closed')
            
            # --- Global parameters            
            #  -  shift controls where in the x-axis the sigmoid crosses prob=0.5 in the y-axis
            # sigmoid_shift = 0.13 # default for thalamic conns
            sigmoid_shift_VPM = 0.04255 # default for thalamic conns
            sigmoid_shift_TRN = 0.056 # default for thalamic conns
            
            # sigmoid_k = -2 # default for thalamic conns
            sigmoid_k = -5 # default for thalamic conns

            # ---------------------
            # L6A CT virtual -> VPM
            # ---------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_VPM, conn_rule_1D_L6A_VPM = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop      = cfg.L6A_VPM__post_pop,
                                                                                                                    pre_pop_axis = cfg.L6A_VPM__pre_pop_axis, post_pop_axis = cfg.L6A_VPM__post_pop_axis,
                                                                                                                    center_point = cfg.center_point,
                                                                                                                    conn_prob=1,    baseline_prob=0,
                                                                                                                    k=sigmoid_k,            shift=sigmoid_shift_VPM,          scale=25,
                                                                                                                    plot_conn=True
                                                                                                                )

            conn_dict_1D_L6A_VPM = BN.BuildNetwork.getConnDict( pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = cfg.L6A_VPM__post_pop,
                                                                conn_method     = conn_method_1D_L6A_VPM,
                                                                conn_rule       = conn_rule_1D_L6A_VPM,
                                                                syn_mech        = cfg.L6A_VPM__syn_mech,
                                                                syns_per_conn   = VPM_estimated_syns_per_conn,
                                                                conn_type       = cfg.L6A_VPM__conn_type,
                                                                weight          = cfg.L6A_VPM__weight,
                                                                target_secs     = cfg.L6A_VPM__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_VPM)

            # ---------------------
            # L6A CT virtual -> TRN
            # ---------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_TRN, conn_rule_1D_L6A_TRN = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop     = cfg.L6A_TRN__post_pop,
                                                                                                                    pre_pop_axis = cfg.L6A_TRN__pre_pop_axis, post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                                                    center_point = cfg.center_point,
                                                                                                                    conn_prob=1,    baseline_prob=0,
                                                                                                                    k=sigmoid_k,            shift=sigmoid_shift_TRN,          scale=25,
                                                                                                                    plot_conn=True
                                                                                                                )

            conn_dict_1D_L6A_TRN = BN.BuildNetwork.getConnDict( pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = cfg.L6A_TRN__post_pop,
                                                                conn_method     = conn_method_1D_L6A_TRN,
                                                                conn_rule       = conn_rule_1D_L6A_TRN,
                                                                syn_mech        = cfg.L6A_TRN__syn_mech,
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                conn_type       = cfg.L6A_TRN__conn_type,
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = cfg.L6A_TRN__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_TRN)

            # --------------------------
            # L6A CT virtual -> TRN_ring
            # --------------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_TRN_ring, conn_rule_1D_L6A_TRN_ring = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(   pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop     = 'TRN_ring',
                                                                                                                            pre_pop_axis = cfg.L6A_TRN__pre_pop_axis, post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                                                            center_point = cfg.center_point,
                                                                                                                            conn_prob=1,    baseline_prob=0,
                                                                                                                            k=sigmoid_k,            shift=sigmoid_shift_TRN,          scale=25,
                                                                                                                            plot_conn=True
                                                                                                                        )

            conn_dict_1D_L6A_TRN_ring = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                        post_pop        = 'TRN_ring',
                                                                        conn_method     = conn_method_1D_L6A_TRN_ring,
                                                                        conn_rule       = conn_rule_1D_L6A_TRN_ring,
                                                                        syn_mech        = cfg.L6A_TRN__syn_mech,
                                                                        syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                        conn_type       = cfg.L6A_TRN__conn_type,
                                                                        weight          = cfg.L6A_TRN__weight,
                                                                        target_secs     = cfg.L6A_TRN__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_TRN_ring)

        elif '__L6A_fb_remixed' in cfg.dump_cell_properties:
            print('\t --- CT Connectivity schema = closed + uniform')
            
            # --- Global parameters            
            #  -  shift controls where in the x-axis the sigmoid crosses prob=0.5 in the y-axis
            # sigmoid_shift = 0.13 # default for thalamic conns
            sigmoid_shift_VPM = 0.04255 # default for thalamic conns
            sigmoid_shift_TRN = 0.056 # default for thalamic conns
            
            # sigmoid_k = -2 # default for thalamic conns
            sigmoid_k = -5 # default for thalamic conns

            # --- Grabbing the values of cfg.topology_ratio directly from the flag in cfg.dump_cell_properties to make it batch compatible
            topology_ratio = [
                int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[0])/100,
                int(cfg.dump_cell_properties.split('_tr_')[1].split('_')[1])/100
            ]

            # ---------------------
            # L6A CT virtual -> VPM
            # ---------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_VPM, conn_rule_1D_L6A_VPM = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop      = cfg.L6A_VPM__post_pop,
                                                                                                                    pre_pop_axis = cfg.L6A_VPM__pre_pop_axis, post_pop_axis = cfg.L6A_VPM__post_pop_axis,
                                                                                                                    center_point = cfg.center_point,
                                                                                                                    conn_prob=1*topology_ratio[0],    baseline_prob=0,
                                                                                                                    # conn_prob=1*(4/5),    baseline_prob=0,
                                                                                                                    k=sigmoid_k,            shift=sigmoid_shift_VPM,          scale=25,
                                                                                                                    plot_conn=True
                                                                                                                )

            conn_dict_1D_L6A_VPM = BN.BuildNetwork.getConnDict( pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = cfg.L6A_VPM__post_pop,
                                                                conn_method     = conn_method_1D_L6A_VPM,
                                                                conn_rule       = conn_rule_1D_L6A_VPM,
                                                                syn_mech        = cfg.L6A_VPM__syn_mech,
                                                                syns_per_conn   = VPM_estimated_syns_per_conn,
                                                                conn_type       = cfg.L6A_VPM__conn_type,
                                                                weight          = cfg.L6A_VPM__weight,
                                                                target_secs     = cfg.L6A_VPM__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_VPM)

            # vCT->VPM
            conn_dict_vCT_VPM = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'VPM',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_VPM*topology_ratio[1],
                                                                # conn_rule       = scaled_convergence_VPM*(1/5),
                                                                syn_mech        = 'syn|L6A|VPM|exc',
                                                                syns_per_conn   = VPM_estimated_syns_per_conn,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_VPM__weight,
                                                                target_secs     = 'inputs__distal')
            
            conn_dict_vCT_VPM_rename = {k+'2':conn_dict_vCT_VPM[k] for k in conn_dict_vCT_VPM.keys()}
            netParams.connParams.update(conn_dict_vCT_VPM_rename)

            # ---------------------
            # L6A CT virtual -> TRN
            # ---------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_TRN, conn_rule_1D_L6A_TRN = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(     pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop     = cfg.L6A_TRN__post_pop,
                                                                                                                    pre_pop_axis = cfg.L6A_TRN__pre_pop_axis, post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                                                    center_point = cfg.center_point,
                                                                                                                    conn_prob=1*topology_ratio[0],    baseline_prob=0,
                                                                                                                    # conn_prob=1*(4/5),    baseline_prob=0,
                                                                                                                    k=sigmoid_k,            shift=sigmoid_shift_TRN,          scale=25,
                                                                                                                    plot_conn=True
                                                                                                                )

            conn_dict_1D_L6A_TRN = BN.BuildNetwork.getConnDict( pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = cfg.L6A_TRN__post_pop,
                                                                conn_method     = conn_method_1D_L6A_TRN,
                                                                conn_rule       = conn_rule_1D_L6A_TRN,
                                                                syn_mech        = cfg.L6A_TRN__syn_mech,
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                conn_type       = cfg.L6A_TRN__conn_type,
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = cfg.L6A_TRN__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_TRN)

            # vCT->TRN
            conn_dict_vCT_TRN = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN*topology_ratio[1],
                                                                # conn_rule       = scaled_convergence_TRN*(1/5),
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            conn_dict_vCT_TRN_rename = {k+'2':conn_dict_vCT_TRN[k] for k in conn_dict_vCT_TRN.keys()}
            netParams.connParams.update(conn_dict_vCT_TRN_rename)

            # --------------------------
            # L6A CT virtual -> TRN_ring
            # --------------------------
            print('Connectivity not tested - must be retuned')
            conn_method_1D_L6A_TRN_ring, conn_rule_1D_L6A_TRN_ring = BN.BuildNetwork.distanceBasedProbability_1D_sigmoid(   pre_pop      = 'CTvirtual_'+virtual_pop,  post_pop     = 'TRN_ring',
                                                                                                                            pre_pop_axis = cfg.L6A_TRN__pre_pop_axis, post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                                                            center_point = cfg.center_point,
                                                                                                                            conn_prob=1*topology_ratio[0],    baseline_prob=0,
                                                                                                                            # conn_prob=1*(4/5),    baseline_prob=0,
                                                                                                                            k=sigmoid_k,            shift=sigmoid_shift_TRN,          scale=25,
                                                                                                                            plot_conn=True
                                                                                                                        )

            conn_dict_1D_L6A_TRN_ring = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                        post_pop        = 'TRN_ring',
                                                                        conn_method     = conn_method_1D_L6A_TRN_ring,
                                                                        conn_rule       = conn_rule_1D_L6A_TRN_ring,
                                                                        syn_mech        = cfg.L6A_TRN__syn_mech,
                                                                        syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                        conn_type       = cfg.L6A_TRN__conn_type,
                                                                        weight          = cfg.L6A_TRN__weight,
                                                                        target_secs     = cfg.L6A_TRN__target_secs)
            netParams.connParams.update(conn_dict_1D_L6A_TRN_ring)

            # vCT->TRN_ring
            conn_dict_vCT_TRNring = BN.BuildNetwork.getConnDict(pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN_ring',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN*topology_ratio[1],
                                                                # conn_rule       = scaled_convergence_TRN*(1/5),
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            conn_dict_vCT_TRNring_rename = {k+'2':conn_dict_vCT_TRNring[k] for k in conn_dict_vCT_TRNring.keys()}
            netParams.connParams.update(conn_dict_vCT_TRNring_rename)

        elif '__L6A_fb_uniform' in cfg.dump_cell_properties:
            print('\t --- CT Connectivity schema = uniform (convergence)')
            # vCT->VPM
            conn_dict_vCT_VPM = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'VPM',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_VPM,
                                                                syn_mech        = 'syn|L6A|VPM|exc',
                                                                syns_per_conn   = VPM_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_VPM__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_VPM)

            # vCT->TRN
            conn_dict_vCT_TRN = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN,
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_TRN)

            # vCT->TRN_ring
            conn_dict_vCT_TRNring = BN.BuildNetwork.getConnDict(pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN_ring',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN,
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_TRNring)
    
        else:
            print('\t --- CT Connectivity schema = uniform (convergence)')
            # vCT->VPM
            conn_dict_vCT_VPM = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'VPM',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_VPM,
                                                                syn_mech        = 'syn|L6A|VPM|exc',
                                                                syns_per_conn   = VPM_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_VPM__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_VPM)

            # vCT->TRN
            conn_dict_vCT_TRN = BN.BuildNetwork.getConnDict(    pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN,
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_TRN)

            # vCT->TRN_ring
            conn_dict_vCT_TRNring = BN.BuildNetwork.getConnDict(pre_pop         = 'CTvirtual_'+virtual_pop,
                                                                post_pop        = 'TRN_ring',
                                                                conn_method     = 'convergence',
                                                                conn_rule       = scaled_convergence_TRN,
                                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                                syns_per_conn   = TRN_estimated_syns_per_conn,
                                                                # syns_per_conn   = 1,
                                                                conn_type       = 'chem',
                                                                weight          = cfg.L6A_TRN__weight,
                                                                target_secs     = 'inputs__distal')
            netParams.connParams.update(conn_dict_vCT_TRNring)

if cfg.includeBaselineDriver:
    conn_dict_vBS_TRN = BN.BuildNetwork.getConnDict(    pre_pop         = 'BaselineDriver',
                                                        post_pop        = 'TRN_ring',
                                                        conn_method     = 'convergence',
                                                        conn_rule       = 161, # from BBP VPM->TRN convergence 161.43545983636517
                                                        syn_mech        = 'syn|VPM|TRN|exc',
                                                        syns_per_conn   = 1,
                                                        conn_type       = 'chem',
                                                        weight          = cfg.BaselineDriver_weights_VPM_TRN_ring,
                                                        target_secs     = 'inputs__proximal')
    netParams.connParams.update(conn_dict_vBS_TRN)
    

if cfg.addL6Ainterconns:
    print('\n\t>>Adding L6A interconnections')
    print('\n\t\t>Loading L6A syn mechs')
    S1L6SynMechs_dict       = NetPyNE_BBP.StoreParameters.getS1L6SynMechs()
    netParams.synMechParams.update(S1L6SynMechs_dict)

    print('\t\t>Loading L6A conns')
    S1L6Connectivity_dict   = NetPyNE_BBP.StoreParameters.getS1L6Connectivity()
    netParams.connParams.update(S1L6Connectivity_dict)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Stimulation
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if cfg.addNoiseIClamp:
    for pop in cfg.NoiseIClampParams.keys():
        netParams.stimSourceParams['NoiseIClamp_source__'+pop] = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': cfg.NoiseIClampParams[pop]['amp']}
        netParams.stimTargetParams['NoiseIClamp_target__'+pop] = {'source': 'NoiseIClamp_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}

if cfg.addHoldingCurrent:
    for pop in cfg.addHoldingCurrentPops:
        netParams.stimSourceParams['HoldingCurrent_source__'+pop]   = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': 0}
        netParams.stimTargetParams['HoldingCurrent_target__'+pop]   = {'source': 'HoldingCurrent_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}
if cfg.addThresholdCurrent:
    for pop in cfg.addThresholdCurrentPops:
        netParams.stimSourceParams['ThresholdCurrent_source__'+pop] = {'type': 'IClamp', 'del': 0, 'dur': 1e9, 'amp': 0}
        netParams.stimTargetParams['ThresholdCurrent_target__'+pop] = {'source': 'ThresholdCurrent_source__'+pop, 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Network mods - Section for specific changes after network creation
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if cfg.scale_density is not None: 
    if cfg.scale_density!= 1.0: 
        NetPyNE_BBP.Prompt.headerMsg('WARNING: SCALING CELL DENSITY TO DEBUG NETWORK - '+str(cfg.scale_density))
        verbose=False
        for pop_name in netParams.popParams.keys():
            if 'MLe' in pop_name: continue
            if verbose: print('Pop name: ', pop_name)
            if 'cellsList' in netParams.popParams[pop_name].keys():             # 'cellsList' must come before 'numCells' because numCells is created internally when a cellsList is passed
                if verbose: print('scalling cellsList ',netParams.popParams[pop_name]['cellsList'])
                cellsList_=netParams.popParams[pop_name]['cellsList']
                netParams.popParams[pop_name]['cellsList'] = netParams.popParams[pop_name]['cellsList'][0:int(len(netParams.popParams[pop_name]['cellsList'])*cfg.scale_density)]
                if len(netParams.popParams[pop_name]['cellsList'])<1:netParams.popParams[pop_name]['cellsList'] = [cellsList_[0]]
                if verbose: print('         cellsList ',netParams.popParams[pop_name]['cellsList'])
            elif   'density'   in netParams.popParams[pop_name].keys(): 
                if verbose: print('scalling density ',netParams.popParams[pop_name]['density'])
                netParams.popParams[pop_name]['density']   = math.ceil(netParams.popParams[pop_name]['density']*cfg.scale_density)
                if verbose: print('         density ',netParams.popParams[pop_name]['density'])
            elif 'numCells'  in netParams.popParams[pop_name].keys(): 
                if verbose: print('scalling numCells ', netParams.popParams[pop_name]['numCells'])
                netParams.popParams[pop_name]['numCells']  = math.ceil(netParams.popParams[pop_name]['numCells']*cfg.scale_density)
                if verbose: print('         numCells ', netParams.popParams[pop_name]['numCells'])

if cfg.removeConns:
    if cfg.removeConns == True:
        print('\n\n\n ---- REMOVING ALL CONNS FOR DEBUGGING --- \n\n\n')
        netParams.connParams={}
    elif cfg.removeConns == 'chem':
        conns=list(netParams.connParams.keys())
        for conn in conns: 
            if 'chem' in conn: del netParams.connParams[conn]
    elif cfg.removeConns == 'elec':
        conns=list(netParams.connParams.keys())
        for conn in conns: 
            if 'elec' in conn: del netParams.connParams[conn]

if cfg.removeESyns:
    conn_names = list(netParams.connParams.keys())
    for conn_name in conn_names:
        if '|elec' in conn_name: 
            print('\n\t>>Removing Electrical conn: ', conn_name)
            del netParams.connParams[conn_name]

if cfg.singleCellPops:
    print('\n\n\n ---- SINGLE CELL POPS FOR DEBUGGING --- \n\n\n')
    if type(cfg.singleCellPops)==list:  pops=cfg.singleCellPops
    else:                               pops=netParams.popParams.keys()


    for pop in pops:
        if   ('VPM' in pop):                                        numCells=1
        elif ('TRN' in pop):                                        numCells=10
        elif ('L6A' in pop) or ('L6CC' in pop) or ('L6IN' in pop):  numCells=5
        else:                                                       numCells=1

        if 'cellsList' in netParams.popParams[pop].keys():
            if len(netParams.popParams[pop]['cellsList'])>=numCells: 
                netParams.popParams[pop]['cellsList']=netParams.popParams[pop]['cellsList'][0:numCells]
                netParams.popParams[pop]['numCells'] = numCells
                print('singleCellPops ', pop,' cellsList: ',netParams.popParams[pop]['cellsList'])
            else:
                print('nope: ', pop)
                print(pop,netParams.popParams[pop]['cellsList'])
                # netParams.popParams[pop]['cellsList']=[netParams.popParams[pop]['cellsList'][0]]
            # netParams.popParams[pop]['numCells']=numCells
        elif 'density' in netParams.popParams[pop].keys():
            del netParams.popParams[pop]['density']
            netParams.popParams[pop]['numCells']=numCells
        elif 'numCells' in netParams.popParams[pop].keys():
            netParams.popParams[pop]['numCells']=numCells
        else:
            netParams.popParams[pop]['numCells']=numCells
            continue

    plotCellShape=False
    if plotCellShape:
        for cell_name in netParams.cellParams.keys():
            NetPyNE_BBP.ConvertSynapses.plotCellShape(          cell=netParams.cellParams[cell_name],
                                                                extra_flag='_debuggingCellShape_'+cell_name.split('__cell')[0]
                                                                )

if      cfg.simulateThalonly:
    print('\n\n\n ---- Running THALAMUS ONLY mode to inspect tune Thalamic connectivity --- \n\n\n')
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name)}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('VPM'    in pop_name)  or ('TRN'   in pop_name) }
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if (('VPM' in conn_name) and ('TRN' in conn_name)) or ('TRN' in conn_name)}
elif      cfg.simulateThalMle:
    print('\n\n\n ---- Running THALAMUS+BRAINSTEM ONLY mode to inspect tune Thalamic connectivity --- \n\n\n')
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name)}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('VPM'    in pop_name)  or ('TRN'   in pop_name) or ('MLe'   in pop_name)}
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if (('VPM' in conn_name)  and ('TRN' in conn_name)) or ('TRN' in conn_name) or ('MLe' in conn_name)}
elif      cfg.simulateThalMle_vCT:
    print('\n\n\n ---- Running THALAMUS+BRAINSTEM ONLY mode to inspect tune Thalamic connectivity --- \n\n\n')
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name)}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('VPM'    in pop_name)  or ('TRN'   in pop_name) or ('MLe'   in pop_name) or ('CTvirtual'   in pop_name) or ('BaselineDriver'   in pop_name)}
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if (('VPM' in conn_name)  and ('TRN' in conn_name)) or ('TRN' in conn_name) or ('MLe' in conn_name) or ('CTvirtual' in conn_name) or ('BaselineDriver' in conn_name)}
elif    cfg.simulateL6only:
    print('\n\n\n ---- Running L6A ONLY mode to inspect tune L6 connectivity --- \n\n\n')
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if 'L6'   in cell_name}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if 'L6'   in pop_name }
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if '_L6_' in conn_name}
elif      cfg.simulateMleonly:
    print('\n\n\n ---- Running BRAINSTEM ONLY mode to inspect Brainstem population creation --- \n\n\n')
    netParams.cellParams    = {}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('MLe'   in pop_name)}
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if ('MLe' in conn_name)}
elif      cfg.simulateMleVPM:
    print('\n\n\n ---- Running VPM+BRAINSTEM ONLY mode to inspect driver connectivity --- \n\n\n')
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if ('VPL_TC' in cell_name)}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('VPM'    in pop_name) or ('MLe'   in pop_name)}
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if ('VPM' in conn_name) or ('MLe' in conn_name)}
elif cfg.simulateTRNgap:
    netParams.cellParams    = {cell_name: cell_params for cell_name, cell_params in netParams.cellParams.items() if ('VPL_TC' in cell_name) or ('Rt_RC' in cell_name) or ('SECell' in cell_name)}
    netParams.popParams     = {pop_name:  pop_params  for pop_name,  pop_params  in netParams.popParams.items()  if ('VPM'    in pop_name)  or ('TRN'   in pop_name) or ('MLe'   in pop_name) or ('SEPop' in pop_name)}
    netParams.connParams    = {conn_name: conn_params for conn_name, conn_params in netParams.connParams.items() if (('gap' in conn_name) and ('TRN' in conn_name)) or ('SEConn' in conn_name)}

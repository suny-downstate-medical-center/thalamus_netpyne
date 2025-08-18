'''
Class to load the parameters to build the network model

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

########################################################################################################################
import numpy as np
import json
import math
import os
import pandas as pd
import pickle
########################################################################################################################

class BuildNetwork():
    def getCellDensity():
        # --- Network dimensions
        cell_density_dict={
                            'L6A':              {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN':             {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'TRN':              {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'VPM':              {'mean':62924.10, 'sd':2886.00,     'e':1,                    'i':0},
                            'S1':               {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},

                            'L1':               {'mean':19154.90,  'sd':14132.90,   'e':0,                    'i':1},
                            'L23':              {'mean':70348.50,  'sd':2308.00,    'e':0.8405069049,         'i':0.1594930951},
                            'L4':               {'mean':110972.40, 'sd':7612.00,    'e':0.8950099304,         'i':0.1049900696},
                            'L5A':              {'mean':64202.70,  'sd':9478.90,    'e':0.7463595768,         'i':0.2536404232},
                            'L5B':              {'mean':64202.70,  'sd':9478.90,    'e':0.7463595768,         'i':0.2536404232},
                            'L6B':              {'mean':48133.20,  'sd':4068.10,    'e':0.8585924061,         'i':0.1414075939},
                            
                            'TRN_ring':         {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_1':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_2':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_3':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_4':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_5':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_6':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_7':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},
                            'TRN_8':            {'mean':56160.40, 'sd':2278.30,     'e':0,                    'i':1},

                            'L6A_activated':    {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6A_suppressed':   {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6A_sparse':       {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6A_silent':       {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},

                            'L6CC_TPC_L1_cAD':  {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6CC_UTPC_cAD':    {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6CC_BPC_cAD':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6CC_IPC_cAD':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            
                            'L6IN_LBC_bAC':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN_LBC_bNA':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN_LBC_cNA':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN_MC_bAC':      {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN_MC_bNA':      {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'L6IN_MC_cAC':      {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},

                            'CTvirtual':        {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            
                            'CTvirtual_activated':  {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'CTvirtual_suppressed': {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'CTvirtual_sparse':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            'CTvirtual_silent':     {'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            
                            'CTvirtual_uniform':{'mean':82257.00, 'sd':8995.60,     'e':0.9172629685011611,   'i':0.082737031498839},
                            
                            'BaselineDriver':   {'mean':62924.10, 'sd':2886.00,     'e':1,                    'i':0},

                            }
        return cell_density_dict
    
    def getNetworkSize(center_point,
                       pops_dict={
                            'VPM':              {'diam':100,   'y':[3300,4000],},
                            'TRN':              {'diam':100,   'y':[2550,2800],},
                            'MLe':              {'diam':55,    'y':[4500,5700],},
                            # 'TRN_ring':         {'diam':200,   'y':[2550,2800],},
                            'TRN_ring':         {'diam':140,   'y':[2550,2800],},
                           
                            'TRN_1':            {'diam':100,   'y':[2550,2800],},
                            'TRN_2':            {'diam':100,   'y':[2550,2800],},
                            'TRN_3':            {'diam':100,   'y':[2550,2800],},
                            'TRN_4':            {'diam':100,   'y':[2550,2800],},
                            'TRN_5':            {'diam':100,   'y':[2550,2800],},
                            'TRN_6':            {'diam':100,   'y':[2550,2800],},
                            'TRN_7':            {'diam':100,   'y':[2550,2800],},
                            'TRN_8':            {'diam':100,   'y':[2550,2800],},

                            # 'L6A':              {'diam':280,   'y':[ 890,1090],}, # --- old dimension - 2024-01-26
                            'L6A':              {'diam':350,   'y':[ 890,1090],},
                            'L6A_activated':    {'diam':350,   'y':[ 890,1090],}, # --- L6 TC neurons
                            'L6A_suppressed':   {'diam':350,   'y':[ 890,1090],}, # --- L6 TC neurons
                            'L6A_silent':       {'diam':350,   'y':[ 890,1090],}, # --- L6 TC neurons
                            'L6A_sparse':       {'diam':350,   'y':[ 890,1090],}, # --- L6 TC neurons
                            'L6IN':             {'diam':350,   'y':[ 890,1090],}, # --- old L6 IN neurons

                            'L6CC_TPC_L1_cAD':  {'diam':350,   'y':[ 890,1090],}, # --- L6 CC neurons
                            'L6CC_UTPC_cAD':    {'diam':350,   'y':[ 890,1090],}, # --- L6 CC neurons
                            'L6CC_BPC_cAD':     {'diam':350,   'y':[ 890,1090],}, # --- L6 CC neurons
                            'L6CC_IPC_cAD':     {'diam':350,   'y':[ 890,1090],}, # --- L6 CC neurons
                           
                            'L6IN_LBC_bAC':     {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            'L6IN_LBC_bNA':     {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            'L6IN_LBC_cNA':     {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            'L6IN_MC_bAC':      {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            'L6IN_MC_bNA':      {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            'L6IN_MC_cAC':      {'diam':350,   'y':[ 890,1090],}, # --- L6 IN neurons
                            
                            'CTvirtual':        {'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            
                            'CTvirtual_activated':  {'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            'CTvirtual_suppressed': {'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            'CTvirtual_sparse':     {'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            'CTvirtual_silent':     {'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            
                            'CTvirtual_uniform':{'diam':350,   'y':[ 890,1090],}, # --- L6 CT neurons
                            
                            'BaselineDriver':   {'diam':200,   'y':[2500,2550],}, # --- Virtual VPM neurons to drive the TRN_ring pop

                           }
                           ):
        
        diam={};height={}

        for pop in pops_dict.keys():diam.update({pop:[center_point-pops_dict[pop]['diam']/2,center_point+pops_dict[pop]['diam']/2]})
        for pop in pops_dict.keys():height.update({pop:pops_dict[pop]['y']})

        # diam.update({key:[center_point-pops_dict[key]['diam']/2,center_point+pops_dict[key]['diam']/2]} for key in pops_dict.keys())
        # height.update({key:pops_dict[key]['y']} for key in pops_dict.keys())

        return diam,height
    
    def getCellTemplate(template='izhi',pops=['VPM','TRN','MLe','L6A'],cell_flag='__cell'):
        if template == 'izhi':
            # --- Cell parameters
            izhi_template = {'secs': {}}
            izhi_template['secs']['soma_0'] = {'geom': {}, 'pointps': {}}                        # soma params dict
            izhi_template['secs']['soma_0']['geom'] = {'diam': 0.1, 'L': 0.1, 'cm': 31.831}      # soma geometry

            cells_dict={}
            for pop in pops: cells_dict.update({pop+cell_flag:izhi_template})
        
        return cells_dict
    
    def getL6ACellTemplate(cellsFolder, file_format='json',cell_pop='L6A',verbose=False):
        if verbose: print('\t-\tLoading L6A detailed cells')
        cells_dict={}
        list_files = os.listdir(cellsFolder)
        L6A_fileNames = [fileName for fileName in list_files if ('.json' in fileName) and (cell_pop in fileName)]
        # L6A_fileNames = [fileName for fileName in list_files if '.json' in fileName]
        for L6A_cell_ind,L6A_fileName in enumerate(L6A_fileNames):
            theta1  = str(L6A_cell_ind)
            theta2  = theta1.zfill(3)
            cellMe  = cell_pop+'@'+str(theta2)+'__cell'
            if file_format == 'json':
                with open(cellsFolder+cellMe+'_cellParams.json', 'r') as openfile: cell_dict_ = json.load(openfile)
            elif file_format == 'pkl':
                if verbose: print('\t-\tLoading cell ', cellMe, ' in pkl format')
                # with open(cellsFolder+cellMe+'_cellParams.pkl', 'rb') as openfile: cell_dict_ = pickle.load(openfile)
                with open(cellsFolder+cellMe+'_cellParams.pkl', 'rb') as openfile: cell_dict_ = pd.read_pickle(openfile)
            else:
                with open(cellsFolder+cellMe+'_cellParams.json', 'r') as openfile: cell_dict_ = json.load(openfile)
            cells_dict.update({cellMe:cell_dict_})
        return cells_dict

    def getPopTemplate(pops=['VPM','TRN','L6A'],center_point=500,pop_flag='__pop',cell_flag='__cell',diversity=False,volume_shape='cylinder',verbose=False):
        cell_density_dict   = BuildNetwork.getCellDensity()
        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict={}

        for pop in pops:
            if pop=='L6A':
                if verbose: print('\t-\tL6A density scaled to reflect the ratio of CT-projecting L6A cells (56.3%) - Crandall, 2017')
                scale_density = 0.563   # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
            else:scale_density = 1

            if   'TRN' in pop:  mech='i'
            else:               mech='e'

            # pops_dict.update({pop+pop_flag:{'cellType': pop+cell_flag, 
            #                                 'density':  cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density,
            #                                 'xRange':   diam[pop], 
            #                                 'yRange':   height[pop], 
            #                                 'zRange':   diam[pop]}})
            
            if diversity: 
                if verbose: print('\t-\tAdding Cell diversity')
                sizeY   = (height[pop][1] - height[pop][0])
                sizeX   = (diam[pop][1]   - diam[pop][0]  )
                sizeZ   = (diam[pop][1]   - diam[pop][0]  )
                if   volume_shape=='cylinder':  vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder
                elif volume_shape=='cube':      vol = sizeY * (sizeX)   * (sizeZ)                # cube
                else:                           vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder

                if verbose: print('\t-\tBuilding network shape as a ', volume_shape)

                pops_dict.update({pop+pop_flag:{'cellType': pop+cell_flag, 
                                                'numCells': round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*1e-9)*vol),
                                                'xRange':   diam[pop], 
                                                'yRange':   height[pop], 
                                                'zRange':   diam[pop],
                                                'diversity': True,
                                                }})
            else:
                pops_dict.update({pop+pop_flag:{'cellType': pop+cell_flag, 
                                                'density':  cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density,
                                                'xRange':   diam[pop], 
                                                'yRange':   height[pop], 
                                                'zRange':   diam[pop],
                                                }})

        return pops_dict

    def getL6ASubtypesPopTemplate(pops=['L6A_activated','L6A_suppressed','L6A_sparse','L6A_silent'],cell_pop='L6A',center_point=500,pop_flag='__pop',cell_flag='__cell',diversity=False,volume_shape='cylinder',group_vecstims=True,random_seed=100000,verbose=False):
        from neuron import h
        cell_density_dict   = BuildNetwork.getCellDensity()
        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict={}

        pop_ratios={'L6A_activated':    0.165,
                    'L6A_suppressed':   0.113,
                    'L6A_sparse':       0.435,
                    'L6A_silent':       0.287}
                    
        for pop_name in pops:
            if verbose: print('\t-\tL6A density scaled to reflect the ratio of CT-projecting L6A cells (56.3%) - Crandall, 2017')
            scale_density = 0.563 # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
            
            mech='e'
            # --- sparse vecstim pop
            if ('sparse' in pop_name):
                if group_vecstims:
                # --- grouped vecstim pops
                    pops_dict.update({pop_name+pop_flag:{   'cellModel': 'VecStim', 
                                                            'numCells':  round((cell_density_dict[cell_pop][mech]*cell_density_dict[cell_pop]['mean']*scale_density*1e-9)*vol*pop_ratios[pop_name]), 
                                                            'rate':  0.22,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                            # 'spkTimes':  [100000000],   # spike outside simDuration
                                                            # 'spkTimes':  cell_spks,   # spike outside simDuration
                                                            'noise':     1,
                                                            'xRange':   diam[cell_pop], 
                                                            'yRange':   height[cell_pop], 
                                                            'zRange':   diam[cell_pop],
                                            }})
                else:
                    # --- individual vecstim pops - to allow a different firing rate for each cell
                    numCells = round((cell_density_dict[cell_pop][mech]*cell_density_dict[cell_pop]['mean']*scale_density*1e-9)*vol*pop_ratios[pop_name])
                    r = h.Random(random_seed)
                    mean=0.22; var =0.55**2
                    for i in range(numCells):
                        rate_neuron=r.normal(mean, var)
                        if rate_neuron<0:rate_neuron=mean/10
                        pops_dict.update({pop_name+'@'+str(i)+pop_flag:{   'cellModel': 'VecStim', 
                                                                'numCells':  1, 
                                                                'rate':  rate_neuron,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                                # 'spkTimes':  [100000000],   # spike outside simDuration
                                                                # 'spkTimes':  cell_spks,   # spike outside simDuration
                                                                'noise':     1,
                                                                'xRange':   diam[cell_pop], 
                                                                'yRange':   height[cell_pop], 
                                                                'zRange':   diam[cell_pop],
                                                }})
            # --- silent vecstim pop
            elif ('silent' in pop_name):
                pops_dict.update({pop_name+pop_flag:{   'cellModel': 'VecStim', 
                                                        'numCells':  round((cell_density_dict[cell_pop][mech]*cell_density_dict[cell_pop]['mean']*scale_density*1e-9)*vol*pop_ratios[pop_name]), 
                                                        # 'rate':  0,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                        'spkTimes':  [100000000],   # spike outside simDuration
                                                        'noise':     1,
                                                        'xRange':   diam[cell_pop], 
                                                        'yRange':   height[cell_pop], 
                                                        'zRange':   diam[cell_pop],
                                        }})
            # ---  compartmental pops
            else:
                if diversity: 
                    if verbose: print('\t-\tAdding Cell diversity')
                    sizeY   = (height[cell_pop][1] - height[cell_pop][0])
                    sizeX   = (diam[cell_pop][1]   - diam[cell_pop][0]  )
                    sizeZ   = (diam[cell_pop][1]   - diam[cell_pop][0]  )
                    if   volume_shape=='cylinder':  vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder
                    elif volume_shape=='cube':      vol = sizeY * (sizeX) * (sizeZ)                # cube
                    else:                           vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder

                    if verbose: print('\t-\tBuilding network shape as a ', volume_shape)

                    pops_dict.update({pop_name+pop_flag:{'cellType': 'L6A__cell', 
                                                    'numCells': round((cell_density_dict[cell_pop][mech]*cell_density_dict[cell_pop]['mean']*scale_density*1e-9)*vol*pop_ratios[pop_name]),
                                                    'xRange':   diam[cell_pop], 
                                                    'yRange':   height[cell_pop], 
                                                    'zRange':   diam[cell_pop],
                                                    'diversity': True,
                                                    }})
                else:
                    pops_dict.update({cell_pop+pop_flag:{'cellType': cell_pop+cell_flag, 
                                                    'density':  cell_density_dict[cell_pop][mech]*cell_density_dict[cell_pop]['mean']*scale_density*pop_ratios[cell_pop],
                                                    'xRange':   diam[cell_pop], 
                                                    'yRange':   height[cell_pop], 
                                                    'zRange':   diam[cell_pop],
                                                    }})
        return pops_dict

    def getL6PopTemplateCtCcIn(pops=['L6A_activated','L6A_suppressed','L6A_sparse','L6A_silent'],center_point=500,pop_flag='__pop',cell_flag='__cell',diversity=False,volume_shape='cylinder',set_cell_positions=True,base_seed=100000,verbose=False,vecstim_pops=[],group_vecstims=True):
        from CreateCellPositions import CreateCellPositions as CCP
        from neuron import h     

        shift_seeds         = 10000
        cell_density_dict   = BuildNetwork.getCellDensity()
        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict={}
        # --------------------------------------------------------------------------------------------------------------
        # --- Ratios of each cell type in the CT (L6A), CC (L6CC) and IN (L6IN) pops of L6A
        '''
        cell_count = {  'L6IN':{'L6IN_LBC_bAC':     24,'L6IN_LBC_bNA':     13,'L6IN_LBC_cNA':     12,'L6IN_MC_bAC':      19,'L6IN_MC_bNA':      5,'L6IN_MC_cAC':      11,}, 'L6CC':{'L6CC_TPC_L1_cAD':  66, 'L6CC_UTPC_cAD':    70, 'L6CC_BPC_cAD':     129, 'L6CC_IPC_cAD':     141,},'L6A':{'L6A_activated':    0.165,'L6A_suppressed':   0.113,'L6A_sparse':       0.435,'L6A_silent':       0.287,},}
        pop_ratios = {key1: {key2: value / sum(inner_dict.values()) for key2, value in inner_dict.items()} for key1, inner_dict in cell_count.items()}
        # cell_count = {'L6IN':{'L6IN_LBC_bAC': 18,'L6IN_LBC_bNA': 9,'L6IN_LBC_cNA': 9,'L6IN_MC_bAC': 9,'L6IN_MC_bNA': 6,'L6IN_MC_cAC': 20,}, 'L6CC':{'L6CC_TPC_L1_cAD':150, 'L6CC_UTPC_cAD': 150, 'L6CC_BPC_cAD': 250, 'L6CC_IPC_cAD': 300,}}
        # pop_ratios = {key1: {key2: value / sum(inner_dict.values()) for key2, value in inner_dict.items()} for key1, inner_dict in cell_count.items()}
        '''
        
        pop_ratios={   
            'L6IN': {  
                'L6IN_LBC_bAC':     0.2857142857142857,
                'L6IN_LBC_bNA':     0.15476190476190477,
                'L6IN_LBC_cNA':     0.14285714285714285,
                'L6IN_MC_bAC':      0.2261904761904762,
                'L6IN_MC_bNA':      0.05952380952380952,
                'L6IN_MC_cAC':      0.13095238095238096,
                },
            'L6CC': {
                'L6CC_TPC_L1_cAD':  0.1625615763546798,
                'L6CC_UTPC_cAD':    0.1724137931034483,
                'L6CC_BPC_cAD':     0.31773399014778325,
                'L6CC_IPC_cAD':     0.3472906403940887,
                },
            'L6A':{
                'L6A_activated':    0.165,
                'L6A_suppressed':   0.113,
                'L6A_sparse':       0.435,
                'L6A_silent':       0.287,
                },
            'TRN':{
                'TRN_ring':         1,
                },
            }
        # --------------------------------------------------------------------------------------------------------------
        # --- Seeds for random number generation using NEURON methods
        pop_seeds={};seed_=1
        for k1 in pop_ratios.keys():
            pop_seeds.update({k1:{}})
            for k2 in pop_ratios[k1].keys():
                pop_seeds[k1].update({k2:base_seed+shift_seeds+seed_})
                seed_+=1
        # --------------------------------------------------------------------------------------------------------------
        # --- Regions and Relative Density
        # Obs: region is defined as the radius value (not diameter)
        cell_dist_dict={
            'L6IN': {
                'L6IN_LBC_bAC':     {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                'L6IN_LBC_bNA':     {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                'L6IN_LBC_cNA':     {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                'L6IN_MC_bAC':      {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                'L6IN_MC_bNA':      {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                'L6IN_MC_cAC':      {'regions':[[0.0, 175.0]], 'relative_density':[1.0],},
                },
            'L6CC': {
                'L6CC_TPC_L1_cAD':  {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.01, 0.025, 0.065, 0.3, 0.6],},
                'L6CC_UTPC_cAD':    {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.01, 0.025, 0.065, 0.3, 0.6],},
                'L6CC_BPC_cAD':     {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.01, 0.025, 0.065, 0.3, 0.6],},
                'L6CC_IPC_cAD':     {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.01, 0.025, 0.065, 0.3, 0.6],},
                },
            'L6A':{
                'L6A_activated':    {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.10, 0.15, 0.25, 0.25, 0.25],},
                'L6A_suppressed':   {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.10, 0.15, 0.25, 0.25, 0.25],},
                'L6A_sparse':       {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.10, 0.15, 0.25, 0.25, 0.25],},
                'L6A_silent':       {'regions':[[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]], 'relative_density':[0.10, 0.15, 0.25, 0.25, 0.25],},
                },
            'TRN': {
                'TRN_ring':         {'regions':[[0.0, (diam['TRN'][1]-diam['TRN'][0])/2],[(diam['TRN'][1]-diam['TRN'][0])/2, (diam['TRN_ring'][1]-diam['TRN_ring'][0])/2]], 'relative_density':[0.0, 1.0],},
                # 'TRN_ring':         {'regions':[[0.0, 50.0],[50.0, 100.0]], 'relative_density':[0.0, 1.0],},
                # 'TRN_ring':         {'regions':[[0.0, 50.0],[50.0, 75.0]], 'relative_density':[0.0, 1.0],},
                # 'TRN_ring':         {'regions':[[0.0, 50.0],[50.0, 70.0]], 'relative_density':[0.0, 1.0],},
                }
            }
        # --------------------------------------------------------------------------------------------------------------
        for pop in pops:
            if 'L6A' in pop:
                if verbose: print('\t-\tL6A density scaled to reflect the ratio of CT-projecting L6A cells (56.3%) - Crandall, 2017')
                scale_density   = 0.563   # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
                pop_key         = 'L6A'
                mech='e'
            elif 'L6CC' in pop: 
                if verbose: print('\t-\tL6A density scaled to reflect the ratio of CC-projecting L6A cells (43.7%) - Crandall, 2017')
                scale_density   = 0.437   # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
                pop_key         = 'L6CC'
                mech='e'
            elif 'L6IN' in pop: 
                scale_density   = 1       # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
                pop_key         = 'L6IN'
                mech='i'
            elif 'TRN' in pop: 
                
                # calculates the volume of the inner and outer rings
                inner_volume = math.pi*(cell_dist_dict['TRN']['TRN_ring']['regions'][0][1])**2
                outer_volume = math.pi*(cell_dist_dict['TRN']['TRN_ring']['regions'][1][1])**2

                # rescales the density term, because the class that calculates the ring population only displaces the inner ring cells to the outside, instead of removing them
                density_rescaling = (outer_volume-inner_volume)/outer_volume

                scale_density   = 1 * density_rescaling       # Ring of TRN cells around the main TRN projection
                pop_key         = 'TRN'
                mech='i'
            else:               scale_density = 1

            # --- Scale proportion of cells in pop
            try: 
                scale_proportion = pop_ratios[pop_key][pop]
                if verbose: print('\t-\tScaling ', pop, ' proportion of cells - ', scale_proportion)
            except:
                if verbose: print('\t-\tFailed to find '+pop_key+' pop: '+pop)
                scale_proportion = 1

            if   'L6A' in pop:  cellType_name = 'L6A'    
            elif 'VPM' in pop:  cellType_name = 'VPL_TC'
            elif 'TRN' in pop:  cellType_name = 'Rt_RC'
            else:               cellType_name = pop
            # print('\t-\t',pop, cellType_name)

            if pop in vecstim_pops:
                # --- sparse vecstim pop
                if ('sparse' in pop):
                    if group_vecstims:
                        print('\t-\t\tCreating grouped sparse population: ', pop)
                        # --- grouped vecstim pops
                        numCells = round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*scale_proportion*1e-9)*vol)
                        pops_dict.update({pop+pop_flag:{   'cellModel': 'VecStim', 
                                                                'numCells':  numCells, 
                                                                'rate':  0.22,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                                # 'spkTimes':  [100000000],   # spike outside simDuration
                                                                # 'spkTimes':  cell_spks,   # spike outside simDuration
                                                                'noise':     1,
                                                                'xRange':   diam[pop], 
                                                                'yRange':   height[pop], 
                                                                'zRange':   diam[pop],
                                                }})
                    else:
                        print('\t-\t\tCreating individual sparse population: ', pop)
                        # --- individual vecstim pops - to allow a different firing rate for each cell
                        numCells = round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*scale_proportion*1e-9)*vol)
                        # r = h.Random(random_seed)
                        r = h.Random(pop_seeds[pop_key][pop])
                        mean=0.22; var =0.55**2
                        for i in range(numCells):
                            rate_neuron=r.normal(mean, var)
                            if rate_neuron<0:rate_neuron=mean/10
                            pops_dict.update({pop+'@'+str(i)+pop_flag:{   'cellModel': 'VecStim', 
                                                                    'numCells':  1, 
                                                                    'rate':  rate_neuron,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                                    # 'spkTimes':  [100000000],   # spike outside simDuration
                                                                    # 'spkTimes':  cell_spks,   # spike outside simDuration
                                                                    'noise':     1,
                                                                    'xRange':   diam[pop], 
                                                                    'yRange':   height[pop], 
                                                                    'zRange':   diam[pop],
                                                    }})
                # --- silent vecstim pop
                elif ('silent' in pop):
                    print('\t-\t\tCreating silent population: ', pop)
                    numCells = round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*scale_proportion*1e-9)*vol)
                    pops_dict.update({pop+pop_flag:{   'cellModel': 'VecStim', 
                                                            'numCells':  numCells, 
                                                            # 'rate':  0,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                            'spkTimes':  [100000000],   # spike outside simDuration
                                                            'noise':     1,
                                                            'xRange':   diam[pop], 
                                                            'yRange':   height[pop], 
                                                            'zRange':   diam[pop],
                                            }})
            # ---  compartmental pops
            else:
                print('\t-\t\tCreating compartmental population: ', pop)
                if diversity: 
                    if verbose: print('\t-\tAdding Cell diversity')
                    sizeY   = (height[pop][1] - height[pop][0])
                    sizeX   = (diam[pop][1]   - diam[pop][0]  )
                    sizeZ   = (diam[pop][1]   - diam[pop][0]  )
                    if   volume_shape=='cylinder':  vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder
                    elif volume_shape=='cube':      vol = sizeY * (sizeX)   * (sizeZ)                # cube
                    else:                           vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder

                    if verbose: print('\t-\tBuilding network shape as a ', volume_shape)
                    numCells = round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*scale_proportion*1e-9)*vol)
                    # pops_dict.update({pop+pop_flag:{'cellType': pop+cell_flag, 
                    pops_dict.update({pop+pop_flag:{'cellType': cellType_name+cell_flag, 
                                                    'numCells': numCells,
                                                    'xRange':   diam[pop], 
                                                    'yRange':   height[pop], 
                                                    'zRange':   diam[pop],
                                                    'diversity': True,
                                                    }})
                    if set_cell_positions:
                        cell_positions_3D = CCP.ring_population(num_cells=numCells, regions=cell_dist_dict[pop_key][pop]['regions'], relative_density=cell_dist_dict[pop_key][pop]['relative_density'], center_point=center_point, y_interval=height[pop], seed=pop_seeds[pop_key][pop], 
                                                                save_fig=None,
                                                                # save_fig='figs_cell_positions/cell_positions__'+pop+'_numCells'+str(numCells)+'_seed'+str(pop_seeds[pop_key][pop])+'.png',
                                                                )
                        # Convert to list of dictionaries
                        cellsList = [{'x': x, 'y': y, 'z': z} for x, y, z in cell_positions_3D]
                        pops_dict[pop+pop_flag].update({'cellsList':cellsList})
                        # print(pop)
                        # print('numCells: ', numCells)
                        # print('cellsList: ', len(cellsList),len(cell_positions_3D), '\n\n\n')
                        # del pops_dict[pop+pop_flag]['numCells']
                else:
                    # pops_dict.update({pop+pop_flag:{'cellType': pop+cell_flag, 
                    pops_dict.update({pop+pop_flag:{'cellType': cellType_name+cell_flag, 
                                                    'density':  cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*scale_proportion,
                                                    'xRange':   diam[pop], 
                                                    'yRange':   height[pop], 
                                                    'zRange':   diam[pop],
                                                    }})

        return pops_dict

    def getTRNadjacentPopTemplate(pops=['TRN'],pop_type='biophysical',center_point=500,pop_flag='__pop',cell_flag='__cell',volume_shape='cylinder',numPops=8,diversity=True,verbose=False):
        from neuron import h
        cell_density_dict   = BuildNetwork.getCellDensity()
        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict={}

        for pop_name in pops:
            pop='TRN'
            mech='i'
            scale_density=1
            
            if verbose: print('\t-\tAdding Cell diversity')
            sizeY   = (height[pop][1] - height[pop][0])
            sizeX   = (diam[pop][1]   - diam[pop][0]  )
            sizeZ   = (diam[pop][1]   - diam[pop][0]  )
            if   volume_shape=='cylinder':  vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder
            elif volume_shape=='cube':      vol = sizeY * (sizeX) * (sizeZ)                # cube
            else:                           vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder

            if verbose: print('\t-\tBuilding network shape as a ', volume_shape)
            # --- grouped vecstim pops
            # print(pop_name)
            # reference: (x increases from l to r)(z increases moving 'into the screen')
            #   [6 7 8]
            #   [4 0 5]     (0 is the center TRN pop, aligned with the thalamic barreloid)
            #   [1 2 3]
            # positons  =[  1,  2,  3,  4,  5,  6,  7,  8]
            shift_x     =[ -1,  0,  1, -1,  1, -1,  0,  1]
            shift_z     =[ -1, -1, -1,  0,  0,  1,  1,  1]
            
            for ind in range(numPops):
                ind_=ind+1 # so that the center pop is the 0th

                if pop_type == 'biophysical':
                    if diversity:
                        if verbose: print('\t-\tAdding biophisical TRN pop # ', str(ind))
                        pops_dict.update({pop_name+'_'+str(ind_)+pop_flag:{ 'cellType': pop+cell_flag, 
                                                                            'numCells': round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*1e-9)*vol),
                                                                            'xRange':   [d+(shift_x[ind]*sizeX) for d in diam[pop]], 
                                                                            'yRange':   height[pop], 
                                                                            'zRange':   [d+(shift_z[ind]*sizeZ) for d in diam[pop]], 
                                                                            'diversity': True,
                                                                            }})
                    else:
                        pops_dict.update({pop_name+'_'+str(ind_)+pop_flag:{ 'cellType': pop+cell_flag, 
                                                                            'density':  cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density,
                                                                            'xRange':   [d+(shift_x[ind]*sizeX) for d in diam[pop]], 
                                                                            'yRange':   height[pop], 
                                                                            'zRange':   [d+(shift_z[ind]*sizeZ) for d in diam[pop]], 
                                                                            }})
                else:
                    pops_dict.update({pop_name+'_'+str(ind_)+pop_flag:{ 'cellModel': 'VecStim', 
                                                                        'numCells':  round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*1e-9)*vol), 
                                                                        'rate':  0.22,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
                                                                        # 'spkTimes':  cell_spks,   # spike outside simDuration
                                                                        'noise':     1,
                                                                        'xRange':   [d+(shift_x[ind]*sizeX) for d in diam[pop]], 
                                                                        'yRange':   height[pop], 
                                                                        'zRange':   [d+(shift_z[ind]*sizeZ) for d in diam[pop]], 
                                                                        # 'xRange':   diam[pop], 
                                                                        # 'yRange':   height[pop], 
                                                                        # 'zRange':   diam[pop],
                                            }})
        return pops_dict

    # def getTRNadjacentPopTemplate(pops=['TRN'],center_point=500,pop_flag='__pop',volume_shape='cylinder',numPops=8):
    #     from neuron import h
    #     cell_density_dict   = BuildNetwork.getCellDensity()
    #     diam,height         = BuildNetwork.getNetworkSize(center_point)

    #     pops_dict={}

    #     for pop_name in pops:
    #         pop='TRN'
    #         mech='i'
    #         scale_density=1
            
    #         print('\t>>>\tAdding Cell diversity')
    #         sizeY   = (height[pop][1] - height[pop][0])
    #         sizeX   = (diam[pop][1]   - diam[pop][0]  )
    #         sizeZ   = (diam[pop][1]   - diam[pop][0]  )
    #         if   volume_shape=='cylinder':  vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder
    #         elif volume_shape=='cube':      vol = sizeY * (sizeX) * (sizeZ)                # cube
    #         else:                           vol = sizeY * (sizeX/2) * (sizeZ/2) * math.pi  # cylinder

    #         print('\t>>>\tBuilding network shape as a ', volume_shape)
    #         # --- grouped vecstim pops
    #         print(pop_name)
    #         # reference: (x increases from l to r)(z increases moving 'into the screen')
    #         #   [6 7 8]
    #         #   [4 0 5]     (0 is the center TRN pop, aligned with the thalamic barreloid)
    #         #   [1 2 3]
    #         # positons  =[  1,  2,  3,  4,  5,  6,  7,  8]
    #         shift_x     =[ -1,  0,  1, -1,  1, -1,  0,  1]
    #         shift_z     =[ -1, -1, -1,  0,  0,  1,  1,  1]
            
    #         for ind in range(numPops):
    #             ind_=ind+1 # so that the center pop is the 0th
    #             pops_dict.update({pop_name+'_'+str(ind_)+pop_flag:{   'cellModel': 'VecStim', 
    #                                                     'numCells':  round((cell_density_dict[pop][mech]*cell_density_dict[pop]['mean']*scale_density*1e-9)*vol), 
    #                                                     'rate':  0.22,   # (Dash ... Crandall, 2022) (0.22 +- 0.55 spikes/s)
    #                                                     # 'spkTimes':  cell_spks,   # spike outside simDuration
    #                                                     'noise':     1,
    #                                                     'xRange':   [d+(shift_x[ind]*sizeX) for d in diam[pop]], 
    #                                                     'yRange':   height[pop], 
    #                                                     'zRange':   [d+(shift_z[ind]*sizeZ) for d in diam[pop]], 
    #                                                     # 'xRange':   diam[pop], 
    #                                                     # 'yRange':   height[pop], 
    #                                                     # 'zRange':   diam[pop],
    #                                     }})
    #     return pops_dict

    def getMLePopTemplate(step=2,center_point=500,align_cells=True):
        # Default: 
        #   step=2              # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (160~200 cells per barrelette)
        #                       # (200+160)/2 = 180 fibers --> (360 degrees / 180 fibers) = 2 degrees per fiber
        
        #   Barrelettes size    # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (diam ~55um / length ~1.2mm)
        #   height = 1.2 mm
        #   diam   = 55 um

        diam,height         = BuildNetwork.getNetworkSize(center_point)
        
        pops_dict   = {}
        thetas      = 360
        theta_range = np.arange(0,thetas,step)
        
        if align_cells: x_positions = np.arange(diam['MLe'][0],diam['MLe'][1],(diam['MLe'][1]-diam['MLe'][0])/int(thetas/step))

        for ind,theta in enumerate(theta_range):
            theta1=str(theta)
            theta2=theta1.zfill(3)

            if align_cells: x_range = [x_positions[ind],x_positions[ind]]
            else:           x_range = diam['MLe']
            if align_cells: z_range = [x_positions[ind],x_positions[ind]]
            else:           z_range = diam['MLe']

            y_position = (height['MLe'][0])+((height['MLe'][1]-height['MLe'][0])*(theta/thetas))
            y_range = [y_position,y_position]

            pops_dict['MLe@'+theta2+'__pop']={  'cellType': 'MLe__cell', 'numCells': 1, 'xRange': x_range, 'yRange': y_range, 'zRange': z_range}
            # pops_dict['MLe@'+theta2+'__pop']={  'cellType': 'MLe__cell', 'numCells': 1, 'xRange': x_range, 'yRange': y_range, 'zRange': [center_point,center_point]}

        return pops_dict
    
    def getMLePopTemplate_VecStim(step=2,center_point=500,align_cells=True,spkts_dict=None):
        # Default: 
        #   step=2              # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (160~200 cells per barrelette)
        #                       # (200+160)/2 = 180 fibers --> (360 degrees / 180 fibers) = 2 degrees per fiber
        
        #   Barrelettes size    # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (diam ~55um / length ~1.2mm)
        #   height = 1.2 mm
        #   diam   = 55 um

        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict   = {}
        thetas      = 360
        theta_range = np.arange(0,thetas,step)

        if align_cells: x_positions = np.arange(diam['MLe'][0],diam['MLe'][1],(diam['MLe'][1]-diam['MLe'][0])/int(thetas/step))

        for ind,theta in enumerate(theta_range):
            theta1=str(theta)
            theta2=theta1.zfill(3)

            if align_cells: x_range = [x_positions[ind],x_positions[ind]]
            else:           x_range = diam['MLe']
            if align_cells: z_range = [x_positions[ind],x_positions[ind]]
            else:           z_range = diam['MLe']

            cell_spks = []  # creating empty list to avoid error
            if spkts_dict is not None:
                try: 
                    try:    cell_spks = spkts_dict[int(theta/step)]         # --- from pkl format
                    except: cell_spks = spkts_dict[str(int(theta/step))]    # --- from json format
                    # print('position: ', theta,'\t cell: ',int(theta/step),'\t spk num: ',len(cell_spks))
                except:
                    print('\t-\tspike assigning failed')
                    cell_spks = [100000]
            
            cell_spks+= [100001] # adding a spike outside of simulation range to prevent empty 'spkTimes' arrays

            # y_position = (height['MLe'][0])+((height['MLe'][1]-height['MLe'][0])*((thetas-theta)/thetas)) # testing reversal of positions so that (bottom: 0 degree -> top: 360 degree)
            y_position = (height['MLe'][0])+((height['MLe'][1]-height['MLe'][0])*(theta/thetas)) # (Timofeeva et al., 2003)
            y_range = [y_position,y_position]

            pops_dict['MLe@'+theta2+'__pop']={  'cellModel': 'VecStim', 
                                                'numCells':  1, 
                                                'spkTimes':  cell_spks,   # spike outside simDuration
                                                'noise':     0,
                                                'xRange': x_range, 'yRange': y_range, 'zRange': z_range,
                                                # 'xRange': x_range, 'yRange': y_range, 'zRange': [center_point,center_point],
                                                }
        
        return pops_dict
    
    def getVirtualCT_Vecstim_cellsList(pop='CTvirtual',center_point=500,pop_flag='__pop',seed=100000,align_cells=True,spkts_dict=None):
        
        diam,height         = BuildNetwork.getNetworkSize(center_point)
        
        pops_dict={}
        numCells=len(spkts_dict.keys())
        # step=1
        step=360/numCells           # 2024_11_05 - reformulating the assignment of positions, so that cells are evenly distributed between [0-360] degrees, instead of [1-numCells], resulting in more than one lap around the unit circle for polar coordinates, in order to use position based connectivity.
        theta_range = np.arange(0,numCells,step)

        theta_array_ = np.arange(0,360,step)
        theta_array  = [theta_array_[v] for v in range(len(spkts_dict.keys()))]
        

        if align_cells: 
            # --- Creates cells in a ring-like distribution in the x-z planes, and in a linear distribution in the y-plane
            cellsList = [{
                            'x':((diam[pop][1]-diam[pop][0])/2)*math.cos(np.deg2rad(theta))+((diam[pop][1]+diam[pop][0])/2),
                            'z':((diam[pop][1]-diam[pop][0])/2)*math.sin(np.deg2rad(theta))+((diam[pop][1]+diam[pop][0])/2),
                            'y':(height[pop][0])+((height[pop][1]-height[pop][0])*(theta/numCells))
                            } 
                            # for theta in np.arange(0,numCells,step)
                            # for theta in np.arange(0,360,step)          # 2024_11_05 - reformulating the assignment of positions, so that cells are evenly distributed between [0-360] degrees, instead of [1-numCells], resulting in more than one lap around the unit circle for polar coordinates, in order to use position based connectivity.
                            for theta in theta_array                      # 2024_11_05 - fixing rounding error - reformulating the assignment of positions, so that cells are evenly distributed between [0-360] degrees, instead of [1-numCells], resulting in more than one lap around the unit circle for polar coordinates, in order to use position based connectivity.
                            ]
        else:
            # --- Creates cells in a random distribution in the x-z planes, and in a linear distribution in the y-plane
            from neuron import h
            # Initialize NEURON random number generators
            if seed is not None: 
                rand1 = h.Random(seed+1)
                rand2 = h.Random(seed+2)

            cellsList = [{
                            'x':rand1.uniform(diam[pop][0], diam[pop][1]),
                            'z':rand2.uniform(diam[pop][0], diam[pop][1]),
                            'y':(height[pop][0])+((height[pop][1]-height[pop][0])*(theta/numCells))
                            } 
                            # for theta in np.arange(0,numCells,step)
                            # for theta in np.arange(0,360,step)          # 2024_11_05 - reformulating the assignment of positions, so that cells are evenly distributed between [0-360] degrees, instead of [1-numCells], resulting in more than one lap around the unit circle for polar coordinates, in order to use position based connectivity.
                            for theta in theta_array                      # 2024_11_05 - fixing rounding error - reformulating the assignment of positions, so that cells are evenly distributed between [0-360] degrees, instead of [1-numCells], resulting in more than one lap around the unit circle for polar coordinates, in order to use position based connectivity.
                            ]

        # --- Adds spike times for each cell in the population
        for ind,ct_gid in enumerate(spkts_dict.keys()):
            cell_spks = []  # creating empty list to avoid error
            if spkts_dict is not None:
                try: 
                    cell_spks = spkts_dict[ct_gid]
                    # print('\t-\tWORKED:\tposition: ', str(ct_gid),'\t cell: ',str(ct_gid),'\t spk num: ',str(len(cell_spks)))
                    # try:    cell_spks = spkts_dict[int(ct_gid)]         # --- from pkl format
                    # except: cell_spks = spkts_dict[str(int(ct_gid))]    # --- from json format
                except:
                    print('\t-\tWARNING: CTvirtual spike assigning failed')
                    cell_spks = [100000]
            cell_spks+= [100001] # adding a spike outside of simulation range to prevent empty 'spkTimes' arrays
            cellsList[ind].update({'spkTimes':cell_spks})

        pops_dict[pop+pop_flag] = { 'cellModel':    'VecStim', 
                                    'noise':        0,
                                    'xRange':       diam[pop], 
                                    'yRange':       height[pop], 
                                    'zRange':       diam[pop],
                                }
        pops_dict[pop+pop_flag].update({'cellsList':cellsList})

        # print('spkts_dict keys: ', spkts_dict.keys())
        # print('len cells list = '+ str(len(cellsList)))

        return pops_dict

    def getMLePopTemplate_VecStim_cellsList(pop='MLe',step=2,center_point=500,pop_flag='__pop',seed=100000,align_cells=True,spkts_dict=None,plotFig=False):
        # Default: 
        #   step=2              # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (160~200 cells per barrelette)
        #                       # (200+160)/2 = 180 fibers --> (360 degrees / 180 fibers) = 2 degrees per fiber
        
        #   Barrelettes size    # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (diam ~55um / length ~1.2mm)
        #   height = 1.2 mm
        #   diam   = 55 um

        diam,height         = BuildNetwork.getNetworkSize(center_point)

        pops_dict   = {}
        thetas      = 360
        theta_range = np.arange(0,thetas,step)

        if align_cells: 
            # --- Creates cells in a ring-like distribution in the x-z planes, and in a linear distribution in the y-plane
            cellsList = [{
                            'x':((diam[pop][1]-diam[pop][0])/2)*math.cos(np.deg2rad(theta))+((diam[pop][1]+diam[pop][0])/2),
                            'z':((diam[pop][1]-diam[pop][0])/2)*math.sin(np.deg2rad(theta))+((diam[pop][1]+diam[pop][0])/2),
                            'y':(height[pop][0])+((height[pop][1]-height[pop][0])*(theta/thetas)) # (Timofeeva et al., 2003),
                            } for theta in np.arange(0,thetas,step)
                            ]
        else:
            # --- Creates cells in a random distribution in the x-z planes, and in a linear distribution in the y-plane
            from neuron import h
            # Initialize NEURON random number generators
            if seed is not None: 
                rand1 = h.Random(seed+1)
                rand2 = h.Random(seed+2)

            cellsList = [{
                            'x':rand1.uniform(diam[pop][0], diam[pop][1]),
                            'z':rand2.uniform(diam[pop][0], diam[pop][1]),
                            'y':(height[pop][0])+((height[pop][1]-height[pop][0])*(theta/thetas)) # (Timofeeva et al., 2003),
                            } for theta in np.arange(0,thetas,step)
                            ]
        
        # --- Adds spike times for each cell in the population
        for ind,theta in enumerate(theta_range):
            cell_spks = []  # creating empty list to avoid error
            if spkts_dict is not None:
                try: 
                    try:    cell_spks = spkts_dict[int(theta/step)]         # --- from pkl format
                    except: cell_spks = spkts_dict[str(int(theta/step))]    # --- from json format
                    # print('position: ', theta,'\t cell: ',int(theta/step),'\t spk num: ',len(cell_spks))
                except:
                    print('\t-\tWARNING: MLe spike assigning failed')
                    cell_spks = [100000]
            cell_spks+= [100001] # adding a spike outside of simulation range to prevent empty 'spkTimes' arrays
            cellsList[ind].update({'spkTimes':cell_spks})
        
        pops_dict[pop+pop_flag] = { 'cellModel':    'VecStim', 
                                    'noise':        0,
                                    'xRange':       diam[pop], 
                                    'yRange':       height[pop], 
                                    'zRange':       diam[pop],
                                }
        pops_dict[pop+pop_flag].update({'cellsList':cellsList})
        

        if plotFig:
            import matplotlib.pyplot as plt
            # Extract x, y, z values
            y_values = [item['x'] for item in cellsList]
            x_values = [item['y'] for item in cellsList]
            z_values = [item['z'] for item in cellsList]

            # Create figure and 3D axis
            fig = plt.figure(8,8)
            ax = fig.add_subplot(111, projection='3d')
            # Plot points
            ax.scatter(x_values, y_values, z_values, c='b', marker='o')
            # Set labels and title
            ax.set_xlabel('Y')
            ax.set_ylabel('X')
            ax.set_zlabel('')
            ax.set_zticks([])
            ax.set_title('Scatter plot of coordinates')
            ax.view_init(elev=90, azim=0)
            # Show plot
            saveFigName = 'test_plot_mle_coords'
            if align_cells: saveFigName+='_aligned'
            plt.savefig(saveFigName+'.png',dpi=200)

        return pops_dict
    
    # Updating method to skip silent period and start spiking at t=10s instead of waiting until t=12s
    def extendSpikeTimes(pops_dict,max_spkt=10000,sim_duration_ms=40000,remove_spkt=[10001,100001],skip_time=2000):
        import copy
        if sim_duration_ms<max_spkt:
            print('\t\t - Skipping adding spikes: (sim_duration<max_spkt) | ',sim_duration_ms,' < ',max_spkt)
            pass
        else:
            for pop_name in pops_dict.keys():
                print('\t\t - Extending spike times for population - ', pop_name)
                if 'cellsList' in pops_dict[pop_name].keys():
                    for i in range(len(pops_dict[pop_name]['cellsList'])):
                        try:    
                            spk_times = pops_dict[pop_name]['cellsList'][str(i)]['spkTimes']
                            # print(str(i)+' str')
                        except: 
                            # print(pop_name, i, pops_dict[pop_name]['cellsList'][i].keys())
                            spk_times = pops_dict[pop_name]['cellsList'][i]['spkTimes']
                            # print(str(i)+' int')
                        # spk_times = pops_dict[pop_name]['cellsList'][i]['spkTimes']
                        spk_times_new = [spkt for spkt in spk_times if spkt<max_spkt]
                        # spk_times_maxTime = [spkt for spkt in spk_times if spkt<max_spkt]
                        if skip_time is not None:   spk_loops=int(sim_duration_ms/(max_spkt-skip_time))
                        else:                       spk_loops=int(sim_duration_ms/max_spkt)
                        if spk_loops<1: pass
                        else:
                            for spk_loop in range(spk_loops):                                               # skips 1 lap to account for the 0th loop, and adds +1 to the 0th index
                                spk_times_maxTime = [spkt for spkt in spk_times if spkt<max_spkt]           # slices the spike times vector at max_spkt
                                
                                # adds a shift in the appended spikes
                                if skip_time is not None: 
                                    spk_times_maxTime_  = [spkt_-skip_time for spkt_ in spk_times_maxTime]
                                    spk_times_maxTime__ = [spkt_ for spkt_ in spk_times_maxTime_ if spkt_>0]

                                
                                # spks_=[(max_spkt*(spk_loop+1))+spk_val for spk_val in spk_times_maxTime]    # appends a copy of the vector spike times * a loop index
                                spks_=[(max_spkt*(spk_loop+1))+spk_val for spk_val in spk_times_maxTime__]    # appends a copy of the vector spike times * a loop index
                                
                                spks__ = [spk_ for spk_ in spks_ if spk_<sim_duration_ms]
                                spk_times_new.extend(spks__)
                                # spk_times_new.extend(spks_)
                            spk_times_new.append(sim_duration_ms+1) # adds a single spike at (t = sim_duration_ms + 1)
                            # print('appending spike at time: ',sim_duration_ms+1 )
                        pops_dict[pop_name]['cellsList'][i]['spkTimes']=spk_times_new
        return pops_dict
    # def extendSpikeTimes(pops_dict,max_spkt=10000,sim_duration_ms=40000,remove_spkt=[10001,100001]):
    #     import copy
    #     if sim_duration_ms<max_spkt:
    #         print('\t\t - Skipping adding spikes: (sim_duration<max_spkt) | ',sim_duration_ms,' < ',max_spkt)
    #         pass
    #     else:
    #         for pop_name in pops_dict.keys():
    #             print('\t\t - Extending spike times for population - ', pop_name)
    #             if 'cellsList' in pops_dict[pop_name].keys():
    #                 for i in range(len(pops_dict[pop_name]['cellsList'])):
    #                     try:    
    #                         spk_times = pops_dict[pop_name]['cellsList'][str(i)]['spkTimes']
    #                         # print(str(i)+' str')
    #                     except: 
    #                         # print(pop_name, i, pops_dict[pop_name]['cellsList'][i].keys())
    #                         spk_times = pops_dict[pop_name]['cellsList'][i]['spkTimes']
    #                         # print(str(i)+' int')
    #                     # spk_times = pops_dict[pop_name]['cellsList'][i]['spkTimes']
    #                     spk_times_new = [spkt for spkt in spk_times if spkt<max_spkt]
    #                     # spk_times_maxTime = [spkt for spkt in spk_times if spkt<max_spkt]
    #                     spk_loops=int(sim_duration_ms/max_spkt)
    #                     if spk_loops<1: pass
    #                     else:
    #                         for spk_loop in range(spk_loops):                                               # skips 1 lap to account for the 0th loop, and adds +1 to the 0th index
    #                             spk_times_maxTime = [spkt for spkt in spk_times if spkt<max_spkt]           # slices the spike times vector at max_spkt
    #                             spks_=[(max_spkt*(spk_loop+1))+spk_val for spk_val in spk_times_maxTime]    # appends a copy of the vector spike times * a loop index
    #                             spks__ = [spk_ for spk_ in spks_ if spk_<sim_duration_ms]
    #                             spk_times_new.extend(spks__)
    #                             # spk_times_new.extend(spks_)
    #                         spk_times_new.append(sim_duration_ms+1) # adds a single spike at (t = sim_duration_ms + 1)
    #                         # print('appending spike at time: ',sim_duration_ms+1 )
    #                     pops_dict[pop_name]['cellsList'][i]['spkTimes']=spk_times_new
    #     return pops_dict
    ''' 
    test values:
    max_spkt=10000
    remove_spkt=[10001,100001]
    sim_duration_ms=40000
    spk_times=[82.0, 124.275, 201.97500000000002, 370.95000000000005, 382.15000000000003, 426.6, 496.52500000000003, 585.075, 682.6750000000001, 702.6500000000001, 733.1, 774.825, 831.7, 873.2, 897.0250000000001, 903.9250000000001, 997.4250000000001, 1056.55, 1247.15, 1297.6750000000002, 1405.8000000000002, 1462.65, 1657.775, 1663.7250000000001, 1686.6750000000002, 1723.4, 1736.1000000000001, 1781.15, 1801.9, 1837.8500000000001, 1879.15, 1969.875, 1985.1000000000001, 1988.4, 2000.9250000000034, 2004.2750000000156, 2008.3750000000305, 2023.375000000085, 2034.650000000126, 2045.3000000001648, 2051.3000000001866, 2058.525000000213, 2084.600000000308, 2093.40000000034, 2098.3250000003577, 2098.97500000036, 2104.950000000382, 2161.1750000000407, 2161.4250000000416, 2253.9750000000145, 2282.100000000117, 2284.925000000127, 2381.450000000478, 2409.725000000581, 2432.6750000006646, 2442.0750000006988, 2505.57500000093, 2505.60000000093, 2509.6750000009447, 2528.5500000010134, 2565.4000000011474, 2575.400000001184, 2592.0250000012443, 2599.7250000012723, 2609.275000001307, 2644.950000001437, 2654.2500000014707, 2679.650000001563, 2714.850000001691, 2721.475000001715, 2752.0500000018264, 2811.900000002044, 2821.1250000020777, 2851.92500000219, 2864.075000002234, 2898.6250000023597, 2909.5000000023992, 2969.5000000026175, 3008.4000000000306, 3011.325000000041, 3012.8000000000466, 3016.3750000000596, 3028.4500000001035, 3116.050000000422, 3119.900000000436, 3142.1750000005172, 3153.4750000000126, 3164.4000000000524, 3167.6750000000643, 3176.8750000000978, 3275.1000000000913, 3296.225000000168, 3316.8750000002433, 3356.400000000387, 3367.0250000004257, 3400.100000000546, 3403.90000000056, 3451.9500000007347, 3475.2250000008194, 3563.0000000011387, 3569.3250000011617, 3616.950000001335, 3631.4500000013877, 3633.6250000013956, 3642.8750000014293, 3666.450000001515, 3678.6500000015594, 3765.375000001875, 3833.900000002124, 3854.3750000021987, 3881.950000002299, 3984.225000002671, 4003.325000000012, 4039.4000000001433, 4041.500000000151, 4049.1000000001786, 4049.3250000001794, 4057.0000000002074, 4088.1000000003205, 4096.450000000351, 4099.00000000036, 4109.825000000399, 4132.575000000483, 4161.225000000586, 4176.225000000641, 4254.12499999994, 4345.6249999986085, 4384.224999999502, 4393.249999999371, 4417.524999999017, 4427.124999998878, 4428.9249999988515, 4449.274999998555, 4503.3499999977685, 4550.67499999708, 4562.9999999969, 4691.674999995028, 4910.324999991846, 5004.624999999933, 5009.59999999986, 5034.724999999495, 5072.899999998939, 5087.224999998731, 5091.049999998675, 5100.474999998538, 5106.649999998448, 5126.724999998156, 5142.24999999793, 5150.47499999781, 5166.7499999975735, 5190.499999997228, 5215.09999999978, 5234.499999999498, 5269.374999999718, 5377.774999999596, 5463.799999998344, 5465.724999998316, 5481.92499999808, 5555.0249999970165, 5580.024999996653, 5631.599999995902, 5671.449999995322, 5698.874999994923, 5708.72499999478, 5854.324999992661, 5898.699999992015, 5905.174999991921, 5925.349999991628, 5938.574999991435, 5997.424999990579, 6005.074999999926, 6031.549999999541, 6043.424999999368, 6073.449999998931, 6080.099999998834, 6154.074999997758, 6168.1249999975535, 6196.0249999971475, 6214.42499999979, 6251.474999999251, 6261.17499999911, 6266.299999999035, 6272.424999998946, 6356.674999999903, 6357.374999999893, 6360.199999999852, 6394.874999999347, 6396.124999999329, 6409.174999999139, 6423.249999998934, 6470.0999999997075, 6568.374999998277, 6597.199999997858, 6670.0749999967975, 6804.949999994835, 6808.674999994781, 6841.549999994302, 6879.0499999937565, 6929.174999993027, 7000.824999999988, 7001.774999999974, 7003.924999999943, 7033.2249999995165, 7106.799999998446, 7126.349999998161, 7164.9249999976, 7302.899999998503, 7317.1999999982945, 7360.399999999849, 7412.999999999083, 7423.374999998932, 7442.07499999866, 7449.549999998551, 7486.674999999466, 7493.8749999993615, 7502.874999999231, 7541.949999998662, 7544.374999998627, 7547.999999998574, 7573.874999998197, 7731.449999995904, 7734.299999995863, 7762.924999995446, 7778.574999995219, 7855.949999994093, 7881.074999993727, 7897.724999993485, 8009.199999991863, 8039.774999991418, 8084.674999990764, 8184.449999989312, 8201.999999989057, 8295.399999987698, 8326.049999987252, 8373.124999986567, 8451.19999998543, 8551.149999983976, 8564.299999983785, 8573.874999983645, 8581.449999983535, 8685.474999982021, 8686.949999982, 8775.299999980714, 8791.499999980479, 8898.124999978927, 8911.62499997873, 9127.024999975596, 9133.799999975497, 9218.09999997427, 9237.474999973989, 9284.7749999733, 9393.474999971719, 9421.249999971315, 9603.849999968657, 9680.324999967544, 9728.299999966846, 9759.449999966393, 9781.149999966077, 9809.274999965668, 9856.324999964983, 9937.024999963809, 9961.67499996345, 9989.12499996305, 10000.774999999989, 10027.424999999892, 10037.89999999974, 10058.924999999434, 10094.22499999892, 100001, 100001, 100001]
    '''

    # def getMLePopTemplate_VecStim(step=2):
    #     pops_dict   = {}
    #     thetas      = 360
    #     theta_range = np.arange(0,thetas,step)

    #     for ind,theta in enumerate(theta_range):
    #         theta1=str(theta)
    #         theta2=theta1.zfill(3)

    #         pops_dict['MLe@'+theta2+'__pop']={  'cellModel': 'VecStim', 
    #                                             'numCells':  1, 
    #                                             'spkTimes':  [100000],   # spike outside simDuration
    #                                             'noise':     0}
        
    #     return pops_dict

    def getSynMechParams():
        mechs_dict={}
        # --- Synaptic mechanism parameters
        mechs_dict['exc']   = {'mod': 'Exp2Syn', 'tau1': 0.8, 'tau2': 5.3, 'e': 0}   # NMDA synaptic mechanism
        mechs_dict['inh']   = {'mod': 'Exp2Syn', 'tau1': 0.6, 'tau2': 8.5, 'e': -75} # GABA synaptic mechanism
        mechs_dict['gap']   = {'mod': 'Gap', 'g': 0.2} # g(nS) - # Gap synaptic mechanism
        mechs_dict['gap_nr']= {'mod': 'Gap_NR', 'g': 0.2} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap_nr2']= {'mod': 'Gap_NR', 'g': 1} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap_nr_testConductance']= {'mod': 'Gap_NR', 'g': 0.05} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['ggap_mod']= {'mod': 'gGapPar', 'g': 0.05} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap_gpt']= {'mod': 'gapGPT', 'g': 0.05} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap_gpt2']= {'mod': 'gapGPT2', 'g': 0.05} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap_gpt3']= {'mod': 'gapGPT3', 'g': 0.05} # g(nS) - # Gap synaptic mechanism
        # mechs_dict['gap2']  = {'mod': 'HalfGap', 'g': 0.2} # g(nS) - # Gap synaptic mechanism from https://github.com/BlueBrain/neuron_reduce/blob/08bea3520c0f535cdba27ef0c3e4d8f970d08604/tests/TestsFiles/Test_4_LBC_amsalem/mod/halfgap.mod#L4 
        mechs_dict['esyn']  = { 'mod': 
                                        'ElectSyn','g': 0.2,
                                        # 'ElectSyn','g': 0.2*0.001,
                                        'pointerParams': {
                                                            'target_var':  'vpeer',
                                                            'source_var': 'v', # already there by default:
                                                            'bidirectional': True # already there by default:
                                                            }
                                    }
        return mechs_dict

    def connectPops():
        # Connection radius - (axonal footprints)
        #                PRE        POST
        conn_radius = { 'VPM':   {  'TRN': 103.57,              # (Lam 2011)
                                    'L6A': 100,     },          # needs data
                        'L6A':  {   'VPM': 155,                 # (Bourassa, 1994) - [... rods in the VPm nucleus. Rod  dimension was fairly constant:  width = 100 um: rostrocaudal  length  = 1.2 mm: thickness ~200 um.]
                                    'TRN': 100,     },          # (Hoogland, 1987) - 
                        'TRN':  {   'VPM': 64.33,               # (Lam 2007)
                                    'TRN': 264.63,  },          # (Lam 2006 - 350 x 200 um)
                        }
        # Connection probability - (adjustable to modify final values of convergence or number of conns)
        conn_prob   = { 'VPM':   {  'TRN': 0.15,                # needs data
                                    'L6A': 0.15,    },          # needs data
                        'L6A':  {   'VPM': 1,                   # needs data
                                    'TRN': 1,       },          # needs data
                        'TRN':  {   'VPM': 1 ,                # needs data
                                    'TRN': 0.036,   },          # needs data
                        }
        return conn_radius,conn_prob

    def cylinderProjection_expDecay(pre_pop,post_pop,center_point=500,y_spread=0.2,set_prob=None,set_radius=None):
        '''
        # probabilistic variables
            conn_prob                   # connection probability (given that the conditions are met)
            l                           # arbitrary value for capping the distance-based probability decay
            conn_radius                 # axonal footprint / radius projected onto the postsynaptic population
            y_spread                    # percentage of conn_radius projected into y-axis / defines the height of the cylinder
        # conditional variables:
            pre_y                       # position of the presynaptic cell
            height[pre_pop][0]          # min boundary on pre_pop
            height[post_pop][0]         # min boundary on post_pop
            height[pre_pop][1]          # max boundary on pre_pop
            height[post_pop][1]         # max boundary on post_pop

            - pre_y must fit within the normalized radius that is calculated using the min and max boundaries of the pre_pop and post_pop

        '''
        conn_method = 'probability'

        diam,height           = BuildNetwork.getNetworkSize(center_point)
        conn_radius,conn_prob = BuildNetwork.connectPops()

        if set_prob:    prob = set_prob
        else:           prob = conn_prob[pre_pop][post_pop]

        if set_radius:  radius = set_radius
        else:           radius = conn_radius[pre_pop][post_pop]

        if 'TRN' in pre_pop: l = radius*10000   # arbitrary value for capping the distance-based probability decay
        else:                l = radius         # arbitrary value for capping the distance-based probability decay
        
        y_thresh    = radius*y_spread           # y-axis radius – proportional to the x-z radius

        conn_rule   = '%f * exp(-dist_2D/%f)*(dist_2D<%f)*(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)'\
                %  (prob, l, radius,
                    height[pre_pop][0],height[post_pop][1],height[post_pop][0],
                    height[pre_pop][1],height[pre_pop][0] , height[post_pop][0],
                    y_thresh
                    )
        return (conn_method,conn_rule)
    
    def cylinderProjection_uniform(pre_pop,post_pop,center_point=500,y_spread=0.2,set_prob=None,set_radius=None):
        '''
        # probabilistic variables
            conn_prob                   # connection probability (given that the conditions are met)
            conn_radius                 # axonal footprint / radius projected onto the postsynaptic population
            y_spread                    # percentage of conn_radius projected into y-axis / defines the height of the cylinder
        # conditional variables:
            pre_y                       # position of the presynaptic cell
            height[pre_pop][0]          # min boundary on pre_pop
            height[post_pop][0]         # min boundary on post_pop
            height[pre_pop][1]          # max boundary on pre_pop
            height[post_pop][1]         # max boundary on post_pop

            - pre_y must fit within the normalized radius that is calculated using the min and max boundaries of the pre_pop and post_pop

        '''
        conn_method = 'probability'

        diam,height           = BuildNetwork.getNetworkSize(center_point)
        conn_radius,conn_prob = BuildNetwork.connectPops()

        if set_prob:    prob = set_prob
        else:           prob = conn_prob[pre_pop][post_pop]

        if set_radius:  radius = set_radius
        else:           radius = conn_radius[pre_pop][post_pop]
        
        y_thresh    = radius*y_spread           # y-axis radius – proportional to the x-z radius

        conn_rule   = '%f *(dist_2D<%f)*(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)'\
                %  (prob, radius,
                    height[pre_pop][0],height[post_pop][1],height[post_pop][0],
                    height[pre_pop][1],height[pre_pop][0] , height[post_pop][0],
                    y_thresh
                    )
        return (conn_method,conn_rule)

    def laminarProjection(pre_pop='MLe',post_pop='VPM',conn_prob=0.8,y_thresh=0.005,center_point=500,thetas=360,step=2,plot_conn=True):

        conn_rules={}
        diam,height           = BuildNetwork.getNetworkSize(center_point)

        # topographycal connectivity
        conn_method = 'probability'
        # y_thresh    = int((height[post_pop][1]-height[post_pop][0])/int(thetas/step))
        # l           = conn_radius[pre_conn][post_conn]      # arbitrary value for capping the distance-based probability decay
        # y_thresh    = conn_radius[pre_conn][post_conn]/5    # y-axis radius – proportional to the x-z radius
        
        conn_rule = '%f *   (\
                                (\
                                    abs(((pre_y-%f)/(%f-%f))-((post_y-%f)/(%f-%f)))\
                                )<%f\
                            )'\
                        % ( conn_prob,
                            height[pre_pop][0], height[pre_pop][1], height[pre_pop][0],
                            height[post_pop][0],height[post_pop][1],height[post_pop][0],
                            y_thresh)

        ############################################################################################################################################

        if plot_conn:
            import math
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20,20))
            print('\t-\tPlotting conn figure')
            
            # dist_list = [0.01*(d-100) for d in range(200)]
            dist_list_  = [(d*0.01)<y_thresh for d in range(100)]
            dist_list = [dist_.real for dist_ in dist_list_]
            dist_list_inds  = [(d*0.01) for d in range(100)]

            plt.subplot(2,1,1)
            exp_decay_list=[]
            for dist in dist_list: exp_decay_list.append(dist)
            plt.plot(dist_list_inds,exp_decay_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.ylabel('Conn probability')
            plt.ylim([0,1.1])
            
            plt.subplot(2,1,2)
            conn_prob_list=[]
            for dist in dist_list: conn_prob_list.append(conn_prob*dist)
            plt.plot(dist_list_inds,conn_prob_list, color='k', linestyle='-', linewidth=4)
            plt.grid()
            plt.xlabel('Absolute distance')
            plt.ylabel('Probability')
            # plt.xlim([-max_dist,max_dist])
            # plt.yticks([0,0.25,0.5,0.75,1.0])
            plt.savefig('../conn/conn_prob_figs/conn_prob_laminarProjection_'+pre_pop+'_to_'+post_pop+'__connProb_'+str(conn_prob)+'_y_thresh_'+str(y_thresh)+'.png')
        
        return conn_method, conn_rule

    def laminarProjection_VecStim(pre_pop='MLe',post_pop='VPM',pre_angle=45,
                                  conn_prob=1,   # test 2023_11_03
                                  y_thresh=0.004,  # test 2023_11_03
                                #   y_thresh=0.0045,  # test 2023_11_03
                                
                                #   conn_prob=0.65,   # test 2023_11_03
                                #   y_thresh=0.0075,  # test 2023_11_03

                                #   conn_prob=0.8,
                                #   y_thresh=0.005,
                                  center_point=500,thetas=360,step=2):

        conn_rules={}
        diam,height           = BuildNetwork.getNetworkSize(center_point)

        # topographycal connectivity
        conn_method = 'probability'
        # y_thresh    = int((height[post_pop][1]-height[post_pop][0])/int(thetas/step))
        # l           = conn_radius[pre_conn][post_conn]      # arbitrary value for capping the distance-based probability decay
        # y_thresh    = conn_radius[pre_conn][post_conn]/5    # y-axis radius – proportional to the x-z radius
        
        # pre_angle_relative = # relative position in the 0-360 degrees space -> pre_angle_relative = pre_angle/360
        # pre_angle_relative = 1-(pre_angle/360) # testing reversal of positions so that (bottom: 0 degree -> top: 360 degree)
        pre_angle_relative = pre_angle/360 # top: 0 degree -> bottom: 360 degree (Timofeeva et al., 2003)

        conn_rule = '%f *   (\
                                (\
                                    abs(%f-((post_y-%f)/(%f-%f)))\
                                )<%f\
                            )'\
                        % ( conn_prob,
                            pre_angle_relative,
                            height[post_pop][0],height[post_pop][1],height[post_pop][0],
                            y_thresh)
        
        return conn_method, conn_rule
    
    def getBoudaries(pre_pop='L6A',post_pop='VPM',pre_pop_axis='y',post_pop_axis='y',center_point=500):
        
        pre_pop_boundaries, post_pop_boundaries = None, None
        
        diam,height           = BuildNetwork.getNetworkSize(center_point)
        if pre_pop is not None:
            # --- Pre pop reference position
            if   (pre_pop_axis == 'x'):         pre_pop_boundaries  = ('pre_x',             diam[pre_pop][0],    diam[pre_pop][1],    diam[pre_pop][0])
            elif (pre_pop_axis == 'z'):         pre_pop_boundaries  = ('pre_z',             diam[pre_pop][0],    diam[pre_pop][1],    diam[pre_pop][0])
            elif (pre_pop_axis == 'y'):         pre_pop_boundaries  = ('pre_y',             height[pre_pop][0],  height[pre_pop][1],  height[pre_pop][0])
            elif (pre_pop_axis == 'theta'):     pre_pop_boundaries  = (['pre_x','pre_z'],   center_point,        diam[pre_pop][1],    center_point)
            else:                               pre_pop_boundaries  = ('pre_y',             height[pre_pop][0],  height[pre_pop][1],  height[pre_pop][0])
        if post_pop is not None:
            # --- Post pop reference position
            if   (post_pop_axis == 'x'):        post_pop_boundaries = ('post_x',            diam[post_pop][0],   diam[post_pop][1],   diam[post_pop][0])
            elif (post_pop_axis == 'z'):        post_pop_boundaries = ('post_z',            diam[post_pop][0],   diam[post_pop][1],   diam[post_pop][0])
            elif (post_pop_axis == 'y'):        post_pop_boundaries = ('post_y',            height[post_pop][0], height[post_pop][1], height[post_pop][0])
            elif (post_pop_axis == 'theta'):    post_pop_boundaries = (['post_x','post_z'], center_point,        diam[post_pop][1],   center_point)
            else:                               post_pop_boundaries = ('post_y',            height[post_pop][0], height[post_pop][1], height[post_pop][0])

        return pre_pop_boundaries,post_pop_boundaries
    
    def find_angle(point, center, unit='radians'):
        """
        Calculate the angle of a point in polar coordinates with respect to a center point.

        Parameters:
        - point:    Tuple (x, y) representing the coordinates of the point.
        - center:   Tuple (x, y) representing the coordinates of the center point.
        - unit:     String, either 'radians' or 'degrees' (default is 'radians').

        Returns:
        - The angle in radians or degrees.
        """
        x, y = point
        x_center, y_center = center
        # Calculate the angle in radians
        angle_rad = math.atan2(y - y_center, x - x_center)
        # Convert to degrees if required
        if unit == 'degrees':   return math.degrees(angle_rad)
        else:                   return angle_rad

    def distanceBasedProbability_1D_uniform(pre_pop='L6A',post_pop='VPM',pre_pop_axis='x',post_pop_axis='y',conn_prob=0.8,center_point=500):

        diam,height           = BuildNetwork.getNetworkSize(center_point)
        # conn_radius,conn_prob = BuildNetwork.connectPops()

        pre_pop_boundaries,post_pop_boundaries = BuildNetwork.getBoudaries(pre_pop=pre_pop,post_pop=post_pop,pre_pop_axis=pre_pop_axis,post_pop_axis=post_pop_axis,center_point=center_point)

        # 1D distance-based probability connectivity
        conn_method = 'probability'

        conn_rule = '%f *   (\
                                    1-abs(((pre_x-%f)/(%f-%f))-((post_y-%f)/(%f-%f)))\
                            )'\
                        % ( conn_prob,
                            diam[pre_pop][0], diam[pre_pop][1], diam[pre_pop][0],
                            height[post_pop][0],height[post_pop][1],height[post_pop][0])


        # conn_rule = '%f *   (1-abs(pre_xnorm-post_ynorm))'% (conn_prob)
        return conn_method, conn_rule

    def distanceBasedProbability_3DDist(    pre_pop='L6A',post_pop='VPM',conn_prob=0.8,baseline_prob=0.05,
                                            center_point=500,
                                            decay_factor=10, # changes the slope of the exponential decay - the larger the number, the narrower the peak and the sharper the decay
                                            plot_conn=False):
        
        # 3D distance-based probability connectivity
        conn_method = 'probability'
        conn_rule = '   %f* (\
                                %f+ exp(\
                                    -%f* dist_3D\
                                        )\
                                )\
                                '% (conn_prob,baseline_prob,
                                    decay_factor)

        ############################################################################################################################################

        if plot_conn:
            import math
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20,20))
            print('\t-\tPlotting conn figure')
            
            # dist_list = [0.01*(d-100) for d in range(200)]
            dist_list = [(d) for d in range(500)]

            plt.subplot(2,1,1)
            exp_decay_list=[]
            for dist in dist_list: exp_decay_list.append(baseline_prob+math.exp(-conn_prob*decay_factor*abs(dist)))
            plt.plot(dist_list,exp_decay_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.ylabel('Exponential factor decay')
            plt.ylim([0,1.1])
            
            plt.subplot(2,1,2)
            conn_prob_list=[]
            for dist in dist_list: conn_prob_list.append(baseline_prob*conn_prob+conn_prob*math.exp(-conn_prob*decay_factor*abs(dist)))
            
            plt.plot(dist_list,conn_prob_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.xlabel('Absolute distance')
            plt.ylabel('Probability')
            # plt.xlim([-max_dist,max_dist])
            # plt.yticks([0,0.25,0.5,0.75,1.0])
            plt.savefig('../conn/conn_prob_figs/conn_prob_decay_3DDist_'+pre_pop+'_to_'+post_pop+'__connProb_'+str(conn_prob)+'_decayFactor_'+str(decay_factor)+'.png')

        return conn_method, conn_rule

    def distanceBasedProbability_3DDist_minDist(    pre_pop='L6A',post_pop='VPM',conn_prob=0.8,baseline_prob=0.05,
                                                    center_point=500,
                                                    decay_factor=10,    # changes the slope of the exponential decay - the larger the number, the narrower the peak and the sharper the decay
                                                    min_dist=10,        # adds a minimum distance between the cell somas to enable the connection
                                                    plot_conn=False):
        
        # 3D distance-based probability connectivity
        conn_method = 'probability'
        conn_rule = '   %f* (\
                                %f+ exp(\
                                    -%f* dist_3D\
                                        )\
                                )\
                                *(dist_3D>%f)\
                                '% (conn_prob,baseline_prob,
                                    decay_factor,
                                    min_dist,)

        ############################################################################################################################################

        if plot_conn:
            import math
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20,20))
            print('\t-\tPlotting conn figure')
            
            # dist_list = [0.01*(d-100) for d in range(200)]
            dist_list = [(d) for d in range(500)]

            plt.subplot(2,1,1)
            exp_decay_list=[]
            for dist in dist_list: 
                if dist>min_dist:   exp_decay_list.append(baseline_prob+math.exp(-conn_prob*decay_factor*abs(dist)))
                else:               exp_decay_list.append(0)
            plt.plot(dist_list,exp_decay_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.ylabel('Exponential factor decay')
            plt.ylim([-0.1,1.1])
            
            plt.subplot(2,1,2)
            conn_prob_list=[]
            for dist in dist_list: conn_prob_list.append(baseline_prob*conn_prob+conn_prob*math.exp(-conn_prob*decay_factor*abs(dist)))
            
            plt.plot(dist_list,conn_prob_list, color='k', linestyle='-', linewidth=4)
            # plt.grid()
            plt.xlabel('Absolute distance')
            plt.ylabel('Probability')
            plt.ylim([-0.1,1.1])
            plt.xlim([-10,510])
            # plt.xlim([-max_dist,max_dist])
            # plt.yticks([0,0.25,0.5,0.75,1.0])
            plt.savefig('../conn/conn_prob_figs/conn_prob_decay_3DDist_minDist_'+pre_pop+'_to_'+post_pop+'__connProb_'+str(conn_prob)+'_decayFactor_'+str(decay_factor)+'.png')

        return conn_method, conn_rule


    def distanceBasedProbability_1D_exponential(pre_pop='L6A',post_pop='VPM',pre_pop_axis='x',post_pop_axis='y',conn_prob=0.8,baseline_prob=0.05,
                                                decay_factor=10, # changes the slope of the exponential decay - the larger the number, the narrower the peak and the sharper the decay
                                                # decay_factor=250,
                                                center_point=500,plot_conn=False,
                                                inverseDecay=False
                                                ):

        diam,height           = BuildNetwork.getNetworkSize(center_point)
        pre_pop_boundaries,post_pop_boundaries = BuildNetwork.getBoudaries(pre_pop=pre_pop,post_pop=post_pop,pre_pop_axis=pre_pop_axis,post_pop_axis=post_pop_axis,center_point=center_point)

        # 1D distance-based probability connectivity
        conn_method = 'probability'

        
        ############################################################################################################################################
        if ('theta' in pre_pop_axis) or ('theta' in post_pop_axis):

            # --- Usign Arctangent function (atan2)
            if ('theta' in pre_pop_axis):
                # find the angle, convert from rad to deg, normalizes [0,360] to [0,1], subtract from normalized position in the post pop
                conn_rule = '   %f* (\
                                    %f+ exp(\
                                        -%f* abs(\
                                                (remainder((((arctan2(%s - %f, %s - %f))*(180/pi))+360),360)/360)-((%s-%f)/(%f-%f))\
                                                )\
                                            )\
                                    )\
                                    '%( conn_prob,baseline_prob,
                                        decay_factor,
                                        pre_pop_boundaries[0][1], pre_pop_boundaries[1],pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                                        post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3]
                                        )
                # conn_rule = '   %f* (\
                #                     %f+ exp(\
                #                         -%f* abs(\
                #                                 (((arctan2(%s - %f, %s - %f))*(180/pi))/360)-((%s-%f)/(%f-%f))\
                #                                 )\
                #                             )\
                #                     )\
                #                     '%( conn_prob,baseline_prob,
                #                         decay_factor,
                #                         pre_pop_boundaries[0][1], pre_pop_boundaries[1],pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                #                         post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3]
                #                         )
            elif ('theta' in post_pop_axis):
                # find the angle, convert from rad to deg, normalizes [0,360] to [0,1], subtract from normalized position in the pre pop
                conn_rule = '   %f* (\
                                    %f+ exp(\
                                        -%f* abs(\
                                                ((%s-%f)/(%f-%f))-(remainder((((arctan2(%s - %f, %s - %f))*(180/pi))+360),360)/360)\
                                                )\
                                            )\
                                    )\
                                    '%( conn_prob,baseline_prob,
                                        decay_factor,
                                        pre_pop_boundaries[0],  pre_pop_boundaries[1],  pre_pop_boundaries[2],  pre_pop_boundaries[3],
                                        post_pop_boundaries[0][1], post_pop_boundaries[1],post_pop_boundaries[0][0], post_pop_boundaries[1],
                                        )
            
            # # --- Usign Arcsin function (asin)
            # if ('theta' in pre_pop_axis):
            #     # find the angle, convert from rad to deg, normalizes [0,360] to [0,1], subtract from normalized position in the post pop
            #     conn_rule = '   %f* (\
            #                         %f+ exp(\
            #                             -%f* abs(\
            #                                     (\
            #                                         ((asin((%s - %f)/sqrt((%s - %f)**2-(%s - %f)**2)))*(180/3.1415))/360)-((%s-%f)/(%f-%f))\
            #                                     )\
            #                                 )\
            #                         )\
            #                         '%( conn_prob,baseline_prob,
            #                             decay_factor,
            #                             pre_pop_boundaries[0][1], pre_pop_boundaries[1],
            #                             pre_pop_boundaries[0][0], pre_pop_boundaries[1],
            #                             pre_pop_boundaries[0][1], pre_pop_boundaries[1],
            #                             post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3]
            #                             )
            # elif ('theta' in post_pop_axis):
            #     # find the angle, convert from rad to deg, normalizes [0,360] to [0,1], subtract from normalized position in the pre pop
            #     conn_rule = '   %f* (\
            #                         %f+ exp(\
            #                             -%f* abs(\
            #                                     (\
            #                                         ((%s-%f)/(%f-%f))-((cos((%s - %f)/sqrt((%s - %f)**2-(%s - %f)**2)))*(180/3.1415))/360)\
            #                                     )\
            #                                 )\
            #                         )\
            #                         '%( conn_prob,baseline_prob,
            #                             decay_factor,
            #                             pre_pop_boundaries[0],  pre_pop_boundaries[1],  pre_pop_boundaries[2],  pre_pop_boundaries[3],
            #                             post_pop_boundaries[0][1], post_pop_boundaries[1],
            #                             post_pop_boundaries[0][0], post_pop_boundaries[1],
            #                             post_pop_boundaries[0][1], post_pop_boundaries[1],
            #                             )

        elif inverseDecay:
            conn_rule = '   1-(%f+\
                            (%f-%f)* (\
                                    exp(\
                                    -%f* abs(\
                                            ((%s-%f)/(%f-%f))-((%s-%f)/(%f-%f))\
                                            )\
                                        )\
                                )\
                            )\
                            '% (baseline_prob,
                                conn_prob,baseline_prob,
                                decay_factor,
                                pre_pop_boundaries[0],  pre_pop_boundaries[1],  pre_pop_boundaries[2],  pre_pop_boundaries[3],
                                post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3])
        else:
            conn_rule = '   %f+\
                            (%f-%f)* (\
                                    exp(\
                                    -%f* abs(\
                                            ((%s-%f)/(%f-%f))-((%s-%f)/(%f-%f))\
                                            )\
                                        )\
                                )\
                                '% (baseline_prob,
                                    conn_prob,baseline_prob,
                                    decay_factor,
                                    pre_pop_boundaries[0],  pre_pop_boundaries[1],  pre_pop_boundaries[2],  pre_pop_boundaries[3],
                                    post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3])
            # conn_rule = '   %f* (\
            #                     %f+ exp(\
            #                         -%f* abs(\
            #                                 ((%s-%f)/(%f-%f))-((%s-%f)/(%f-%f))\
            #                                 )\
            #                             )\
            #                     )\
            #                     '% (conn_prob,baseline_prob,
            #                         decay_factor,
            #                         pre_pop_boundaries[0],  pre_pop_boundaries[1],  pre_pop_boundaries[2],  pre_pop_boundaries[3],
            #                         post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3])

        ############################################################################################################################################

        if plot_conn:
            import math
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20,20))
            print('\t-\tPlotting conn figure')
            
            dist_list = [0.01*(d-100) for d in range(200)]

            plt.subplot(2,1,1)
            exp_decay_list=[]
            for dist in dist_list: exp_decay_list.append(baseline_prob+math.exp(-conn_prob*decay_factor*abs(dist)))
            plt.plot(dist_list,exp_decay_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.ylabel('Exponential factor decay')
            plt.ylim([0,1.1])
            
            plt.subplot(2,1,2)
            conn_prob_list=[]
            if inverseDecay:
                for dist in dist_list: conn_prob_list.append(1-(baseline_prob+((conn_prob-baseline_prob)*math.exp(-decay_factor*abs(dist)))))
            else: 
                for dist in dist_list: conn_prob_list.append(baseline_prob+((conn_prob-baseline_prob)*math.exp(-decay_factor*abs(dist))))
            # for dist in dist_list: conn_prob_list.append(baseline_prob*conn_prob+conn_prob*math.exp(-conn_prob*decay_factor*abs(dist)))
            
            plt.plot(dist_list,conn_prob_list,color='red',linestyle=':',linewidth=1)
            plt.grid()
            plt.xlabel('Normalized distance')
            plt.ylabel('Probability')
            # plt.xlim([-max_dist,max_dist])
            # plt.yticks([0,0.25,0.5,0.75,1.0])
            plt.savefig('../conn/conn_prob_figs/conn_prob_decay_'+pre_pop+'_to_'+post_pop+'__connProb_'+str(conn_prob)+'_decayFactor_'+str(decay_factor)+'.png')
        # import sys;sys.exit()
        return conn_method, conn_rule

    def distanceBasedProbability_1D_sigmoid(pre_pop='L6A', post_pop='VPM', pre_pop_axis='x', post_pop_axis='y',
                                            conn_prob=0.8, baseline_prob=0.05,
                                            k=1, shift=0.5, scale=10,  # Parameters for the reverse sigmoid function
                                            center_point=500, plot_conn=False,
                                            inverseDecay=False):

        print('Running sigmoid function')
        # Assume BuildNetwork functions are available
        diam, height = BuildNetwork.getNetworkSize(center_point)
        pre_pop_boundaries, post_pop_boundaries = BuildNetwork.getBoudaries(pre_pop=pre_pop, post_pop=post_pop,
                                                                            pre_pop_axis=pre_pop_axis, post_pop_axis=post_pop_axis,
                                                                            center_point=center_point)

        # 1D distance-based probability connectivity
        conn_method = 'probability'

        if ('theta' in pre_pop_axis) or ('theta' in post_pop_axis):

            # Using Arctangent function (atan2)
            if ('theta' in pre_pop_axis):

                
                # (2024_11_07) - Updating theta conn rule 
                # theta_pre

                pre_cell_angle_radians_str  = 'arctan2(%s - %f, %s - %f)'
                rad_to_deg_str              = '180/pi'

                # dist_pre        = np.remainder((np.arctan2()*180/np.pi)+360,360)/360
                postion_pre_norm_string    = 'remainder(('+pre_cell_angle_radians_str+'*'+rad_to_deg_str+')+360,360)/360'
                postion_post_norm_string   = '(%s - %f) / (%f - %f)'

                distance_abs_string = 'abs(('+postion_pre_norm_string+') - ('+postion_post_norm_string+'))'

                # Equation reference string:
                #   y = 1 / (1 + np.exp(-k * (abs(x) - shift) * scale))
                #   sigmoid_decay_list.append(baseline_prob + (conn_prob - baseline_prob) * y)
                
                # adding k, shift and scale strings
                sigmoid_template_string = '1 / (1+exp(-1 * (%f) * ('+distance_abs_string+' - %f) * %f))'

                # adding baseline_prob, conn_prob, baseline_prob
                sigmoid_prob_string     = '%f + (%f - %f) * ('+sigmoid_template_string+')'

                print('Using Theta in pre pop rule: ', sigmoid_prob_string)

                conn_rule = sigmoid_prob_string%(   baseline_prob, conn_prob, baseline_prob, k,
                                                    pre_pop_boundaries[0][1], pre_pop_boundaries[1],pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                                                    post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                                                    shift, scale
                                                    )


                # dist_pre_string = 'remainder(\
                #                                 (\
                #                                     (\
                #                                         (\
                #                                             arctan2(%s - %f, %s - %f)\
                #                                         )\
                #                                         *\
                #                                         (\
                #                                             180/pi\
                #                                         )\
                #                                     )+360\
                #                                 ),360\
                #                             )/360'

                # conn_rule = '   %f + ( \
                #                         (%f - %f) * ( \
                #                                         1 / (\
                #                                                 1 + exp(\
                #                                                             -%f * (\
                #                                                                         (\
                #                                                                             abs(\
                #                                                                                     (\
                #                                                                                         remainder((((arctan2(%s - %f, %s - %f))*(180/pi))+360),360)/360\
                #                                                                                     )\
                #                                                                                     - (\
                #                                                                                         (%s - %f) / (%f - %f)\
                #                                                                                     )\
                #                                                                                 )\
                #                                                                         - %f)\
                #                                                                     * %f)\
                #                                                         )\
                #                                             ) \
                #                                     ) \
                #                     ) \
                #                     '%( baseline_prob, conn_prob, baseline_prob, k,
                #                         pre_pop_boundaries[0][1], pre_pop_boundaries[1],pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                #                         post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                #                         shift, scale
                #                         )

                # # conn_rule = '   %f* (\
                # #                 %f + ( \
                # #                     1 / (1 + exp(-%f * ((abs((remainder((((arctan2(%s - %f, %s - %f))*(180/pi)) + 360), 360) / 360) - ((%s - %f) / (%f - %f)))) - %f) * %f))) \
                # #                 ) \
                # #             ' % (conn_prob, baseline_prob, k,
                # #                 pre_pop_boundaries[0][1], pre_pop_boundaries[1], pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                # #                 post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                # #                 shift, scale)
                # # # # --- Fixing expression with a missing parenthesis
                # # # conn_rule = '   %f* (\
                # # #                 %f + ( \
                # # #                     1 / (1 + exp(-%f * ((abs((remainder((((arctan2(%s - %f, %s - %f))*(180/pi)) + 360), 360) / 360) - ((%s - %f) / (%f - %f))) - %f) * %f))) \
                # # #                 ) \
                # # #             ' % (conn_prob, baseline_prob, k,
                # # #                 pre_pop_boundaries[0][1], pre_pop_boundaries[1], pre_pop_boundaries[0][0], pre_pop_boundaries[1],
                # # #                 post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                # # #                 shift, scale)
                

            elif ('theta' in post_pop_axis):
                
                # (2024_11_07) - Updating theta conn rule 
                # theta_post
                conn_rule = '   %f + ( \
                                        (%f - %f) * ( \
                                                        1 / (\
                                                                1 + exp(\
                                                                            -%f * (\
                                                                                        (\
                                                                                            abs(\
                                                                                                    (\
                                                                                                        (%s - %f) / (%f - %f)\
                                                                                                    )\
                                                                                                    - (\
                                                                                                        remainder((((arctan2(%s - %f, %s - %f))*(180/pi))+360),360)/360\
                                                                                                    )\
                                                                                                )\
                                                                                        - %f)\
                                                                                    * %f)\
                                                                        )\
                                                            ) \
                                                    ) \
                                    ) \
                                    '%( baseline_prob, conn_prob, baseline_prob, k,
                                        pre_pop_boundaries[0], pre_pop_boundaries[1], pre_pop_boundaries[2], pre_pop_boundaries[3],
                                        post_pop_boundaries[0][1], post_pop_boundaries[1], post_pop_boundaries[0][0], post_pop_boundaries[1],
                                        shift, scale
                                        )

                # conn_rule = '   %f* (\
                #                 %f + ( \
                #                     1 / (1 + exp(-%f * ((abs(((%s - %f) / (%f - %f)) - (remainder((((arctan2(%s - %f, %s - %f))*(180/pi)) + 360), 360) / 360))) - %f) * %f))) \
                #                 ) \
                #             ' % (conn_prob, baseline_prob, k,
                #                 pre_pop_boundaries[0], pre_pop_boundaries[1], pre_pop_boundaries[2], pre_pop_boundaries[3],
                #                 post_pop_boundaries[0][1], post_pop_boundaries[1], post_pop_boundaries[0][0], post_pop_boundaries[1],
                #                 shift, scale)

        elif inverseDecay:
            conn_rule = '   1 - ( \
                            %f + ( \
                                (%f - %f) * ( \
                                    1 / (1 + exp(-%f * ((abs(((%s - %f) / (%f - %f)) - ((%s - %f) / (%f - %f))) - %f) * %f))) \
                                ) \
                            ) \
                        ) \
                        ' % (baseline_prob, conn_prob, baseline_prob, k,
                            pre_pop_boundaries[0], pre_pop_boundaries[1], pre_pop_boundaries[2], pre_pop_boundaries[3],
                            post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                            shift, scale)
        else:
            conn_rule = '   %f + ( \
                            (%f - %f) * ( \
                                1 / (1 + exp(-%f * ((abs(((%s - %f) / (%f - %f)) - ((%s - %f) / (%f - %f))) - %f) * %f))) \
                            ) \
                        ) \
                        ' % (baseline_prob, conn_prob, baseline_prob, k,
                            pre_pop_boundaries[0], pre_pop_boundaries[1], pre_pop_boundaries[2], pre_pop_boundaries[3],
                            post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                            shift, scale)

        if plot_conn:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20, 20))
            print('\t-\tPlotting conn figure')

            dist_list = [0.01 * (d - 100) for d in range(200)]
            print(pre_pop_boundaries,post_pop_boundaries)
            plt.subplot(2, 1, 1)
            sigmoid_decay_list = []
            for dist in dist_list:
                # Calculate sigmoid decay
                x = dist
                y = 1 / (1 + np.exp(-k * (abs(x) - shift) * scale))
                sigmoid_decay_list.append(baseline_prob + (conn_prob - baseline_prob) * y)
            plt.plot(dist_list, sigmoid_decay_list, color='red', linestyle=':', linewidth=1)
            plt.grid()
            plt.ylabel('Sigmoid factor decay')
            plt.ylim([0, 1.1])

            plt.subplot(2, 1, 2)
            conn_prob_list = []
            if inverseDecay:
                for dist in dist_list:
                    x = dist
                    y = 1 / (1 + np.exp(-k * (abs(x) - shift) * scale))
                    conn_prob_list.append(1 - (baseline_prob + (conn_prob - baseline_prob) * y))
            else:
                for dist in dist_list:
                    x = dist
                    y = 1 / (1 + np.exp(-k * (abs(x) - shift) * scale))
                    conn_prob_list.append(baseline_prob + (conn_prob - baseline_prob) * y)

            plt.plot(dist_list, conn_prob_list, color='k', linestyle='-', linewidth=8)
            plt.ylim([-0.1, 1.1])
            plt.yticks([0, 0.25, 0.5, 0.75, 1])
            plt.xlim([-1.05, 1.05])
            plt.xticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
            # plt.grid()
            plt.xlabel('Normalized distance')
            plt.ylabel('Probability')
            plt.savefig('../conn/conn_prob_figs/conn_prob_sigmoid_'+pre_pop+'_to_'+post_pop+'__connProb_'+str(conn_prob)+'_k_'+str(k)+'_shift_'+str(shift)+'_scale_'+str(scale)+'.png')

        return conn_method, conn_rule




    def positionBasedProbability_1D_linear( pre_pop='VPM',post_pop='L6A',pre_pop_axis='y',post_pop_axis='y',conn_prob=0.8,
                                            coef=[-185.45820847,  228.16359477], # from (Meyer,2010) m and b terms from Y=(mX+b) x=depth/y=#_of_boutons
                                            center_point=500,
                                            ):

        diam,height           = BuildNetwork.getNetworkSize(center_point)
        pre_pop_boundaries,post_pop_boundaries = BuildNetwork.getBoudaries(pre_pop=pre_pop,post_pop=post_pop,pre_pop_axis=pre_pop_axis,post_pop_axis=post_pop_axis,center_point=center_point)

        # 1D distance-based probability connectivity
        conn_method = 'probability'

        # post_y is the 'x' in the y=mx+b equation and 'y' is the #_of_boutons
        
        conn_rule = ' %f *  (\
                                %f * (\
                                        (%s-%f)/(%f-%f)\
                                    )+%f\
                            )/%f \
                    '\
                    % ( conn_prob,
                        coef[0],
                        post_pop_boundaries[0], post_pop_boundaries[1], post_pop_boundaries[2], post_pop_boundaries[3],
                        coef[1],
                        coef[1] # divided by coef[1] (or 'b') to normalize boutons between 0-1 range, so that the only variable is the probability and the average number of conns*syns_per_conn, which is calculated from the output network
                        )
                            # *(post_y<910)\
                            # *(post_y>1070)\
        
        return conn_method, conn_rule

    def getConnDict(pre_pop,post_pop,conn_method,conn_rule,syn_mech,syns_per_conn,conn_type,weight=1,target_secs='soma_0'):

        # diam,height = BuildNetwork.getNetworkSize(center_point=500)
        # if pre_pop != post_pop: # avoid case where distance between pops = 0
        #     pops_distance = abs(((height[post_pop][0]+height[post_pop][1])/2) - ((height[pre_pop][0]+height[pre_pop][1])/2))
        #     # delay_string = 'defaultDelay+'+str(pops_distance)+'/propVelocity'
        #     # delay_string = 'defaultDelay+lognormal('+str(pops_distance)+'/propVelocity,defaultDelay)'
        #     # delay_string = 'defaultDelay+('+str(pops_distance)+'/propVelocity)+poisson(1)'
        #     delay_string = 'lognormal(defaultDelay,defaultDelay/4)+('+str(pops_distance)+'/propVelocity)'
        # else:
        #     # delay_string = 'defaultDelay+dist_3D/propVelocity'
        #     delay_string = 'lognormal(defaultDelay,defaultDelay/4)+(dist_3D/propVelocity)'
        # print('\t\t\t\t',pre_pop, ' -> ', post_pop, ': ',delay_string,'\n')
        
        # # print('\t\t\t\t',pre_pop, ' -> ', post_pop, ': ',2.0+pops_distance/500,'\n')
        # print('\t\t\t\t',pre_pop, ' -> ', post_pop, ': ',2.0+pops_distance/100000,'\n')

        
        # delay_string = 'lognormal(defaultDelay,defaultDelay/4)'
        if ('MLe' in pre_pop):  
            print('Adding fixed delay to MLe->VPM conns')
            delay_string = 'defaultDelay'
        else:                   delay_string = 'defaultDelay+lognormal(defaultDelay,defaultDelay/4)'


        conn_dict={'conn|'+pre_pop+'|'+post_pop+'|'+conn_type: {    'preConds':  {'pop': pre_pop+'__pop'}, 
                                                                    'postConds': {'pop': post_pop+'__pop'},
                                                                    'synMech': syn_mech,
                                                                    conn_method:  conn_rule,
                                                                    'weight': weight, 
                                                                    'delay': delay_string,
                                                                    # 'delay': 'defaultDelay+uniform(0,0.25)',
                                                                    # 'delay': 'defaultDelay+dist_3D/propVelocity',
                                                                    'synsPerConn': syns_per_conn,
                                                                    'sec': target_secs}}
        return conn_dict
    
    # --- Not in use - Gap junctions are being connected using the general getConnDict method
    def getConnDictGap(pre_pop,post_pop,conn_method,conn_rule,syn_mech,syns_per_conn,target_secs='soma_0'):
        conn_dict={'conn|'+pre_pop+'|'+post_pop+'|elec': {  'preConds':  {'pop': pre_pop+'__pop'}, 
                                                    'postConds': {'pop': post_pop+'__pop'},
                                                    'synMech': syn_mech,
                                                    conn_method:  conn_rule,
                                                    'weight': 1, 
                                                    'delay': 0.000001,
                                                    # 'delay': 'defaultDelay+dist_3D/propVelocity',
                                                    'synsPerConn': syns_per_conn,
                                                    'sec': target_secs}}
        return conn_dict

    def getConnDictPopList(pre_pops,post_pops,pre_pops_flag,post_pops_flag,conn_method,conn_rule,syn_mech,syns_per_conn,conn_type,weight=1,target_secs='soma_0'):
        
        # delay_string = 'lognormal(defaultDelay,defaultDelay/4)'
        delay_string = 'defaultDelay+lognormal(defaultDelay,defaultDelay/4)'

        conn_dict={'conn|'+pre_pops_flag+'|'+post_pops_flag+'|'+conn_type: {    'preConds':  {'pop': pre_pops}, 
                                                                                'postConds': {'pop': post_pops},
                                                                                'synMech': syn_mech,
                                                                                conn_method:  conn_rule,
                                                                                'weight': weight, 
                                                                                'delay': delay_string,
                                                                                # 'delay': 'defaultDelay+uniform(0,0.25)',
                                                                                # 'delay': 'defaultDelay+dist_3D/propVelocity',
                                                                                'synsPerConn': syns_per_conn,
                                                                                'sec': target_secs}}
        return conn_dict
    
    def getConnDictCellType(pre_cellType,post_cellType,pre_pops_flag,post_pops_flag,conn_method,conn_rule,syn_mech,syns_per_conn,conn_type,weight=1,target_secs='soma_0'):
        
        # delay_string = 'lognormal(defaultDelay,defaultDelay/4)+(dist_3D/propVelocity)'
        # delay_string = 'lognormal(defaultDelay,defaultDelay/4)'
        delay_string = 'defaultDelay+lognormal(defaultDelay,defaultDelay/4)'
        
        conn_dict={'conn|'+pre_pops_flag+'|'+post_pops_flag+'|'+conn_type: {    'preConds':  {'cellType': pre_cellType}, 
                                                                                'postConds': {'cellType': post_cellType},
                                                                                'synMech': syn_mech,
                                                                                conn_method:  conn_rule,
                                                                                'weight': weight, 
                                                                                'delay': delay_string,
                                                                                # 'delay': 'defaultDelay+uniform(0,0.25)',
                                                                                # 'delay': 'defaultDelay+dist_3D/propVelocity',
                                                                                'synsPerConn': syns_per_conn,
                                                                                'sec': target_secs}}
        return conn_dict
    


########################################################################################################################

class NetworkConversion():
    
    # --- Recursive function to replace the keys of a dictionary at any nested level 
    # (source: https://stackoverflow.com/questions/38491318/replace-keys-in-a-nested-dictionary)
    # old_dict:     dictionary to be modified
    # key_dict:     dictionary with <current keys> as <keys>, and <new keys> as <values> (e.g.: {old_key1:new_key1, old_key2:new_key2})
    def replace_keys(old_dict, key_dict):
        new_dict = { }
        for key in old_dict.keys():
            new_key = key_dict.get(key, key)
            if isinstance(old_dict[key], dict): new_dict[new_key] = NetworkConversion.replace_keys(old_dict[key], key_dict)
            else:                               new_dict[new_key] = old_dict[key]
        return new_dict

    # def convertConnProperties(filePath='Users/joao/Research/Models/BBP/thalamus_netpyne/conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'):
    #     with open(filePath, 'r') as openfile: conn_data_BBP = json.load(openfile)
    #     replaces={'VPL':'VPM','Rt_RC':'TRN','CorticoThalamic':'L6A','MedialLemniscus':'MLe'}
    #     new_conn_data_BBP={}
    #     for key1 in conn_data_BBP.keys():
    #         for replace in replaces.keys():
    #             if replace in key1:
    #                 key1_=replaces[replace]
    #                 new_conn_data_BBP.update({key1_:{}})
            
    #         for key2 in conn_data_BBP[key1].keys():
    #             if ('Rt_RC' in key1) and ('MedialLemniscus' in key2):
    #                 print('\t>>\tRemoving innexistent connection between ', key1, ' and ', key2)
    #                 continue
    #             for replace in replaces.keys():
    #                 if replace in key2:
    #                     # print(key2,replace)
    #                     new_conn_data_BBP[key1_].update({replaces[replace]:conn_data_BBP[key1][key2]})

    #     return new_conn_data_BBP

    def convertConnProperties(filePath=None):
        if filePath==None:
            print('\t-\tplease, provide the path to the BBP conn properties file')
            return
        
        with open(filePath, 'r') as openfile: conn_data_BBP = json.load(openfile)
        replaces={'VPL':'VPM','Rt_RC':'TRN','CorticoThalamic':'L6A','MedialLemniscus':'MLe'}
        new_conn_data_BBP={}
        for key1 in conn_data_BBP.keys():
            for replace in replaces.keys():
                if replace in key1: 
                    key1_=replaces[replace]
                    new_conn_data_BBP.update({key1_:{'chem':{},'elec':{}}})

            for key2 in conn_data_BBP[key1].keys():
                if 'electrical' in key2:    syn_type='elec'
                elif 'chemical' in key2:    syn_type='chem'
                else:                       continue
            
                if ('Rt_RC' in key1) and ('MedialLemniscus' in key2):
                    print('\t-\tRemoving innexistent connection between ', key1, ' and ', key2)
                    continue
                for replace in replaces.keys():
                    if replace in key2:
                        # print(key2,replace)
                        new_conn_data_BBP[key1_][syn_type].update({replaces[replace]:conn_data_BBP[key1][key2]})

        return new_conn_data_BBP




        # # center_point = 500
        # VPM_size = 85   # (Claus, 2018) - barrel width of 50-120 um
        # TRN_size = 85   # needs data - should probably be increased to represent realistic axonal footprints
        # # L6A_size = 200  # (Crandall, 2017) - Layer 6A infrabarrels only - measured using digitizer tool
        # L6A_size = 280  # (Welker and Woolsey, 1974)
        # S1_size  = 280  # (Welker and Woolsey, 1974)
        # PR5_size = 55   # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (diam ~55um / length ~1.2mm)

        # VPM_range_xz = [center_point-(VPM_size/2),center_point+(VPM_size/2)]
        # TRN_range_xz = [center_point-(TRN_size/2),center_point+(TRN_size/2)]
        # L6A_range_xz = [center_point-(L6A_size/2),center_point+(L6A_size/2)]
        # L6A_scale_density = 0.563 # (Crandall, 2017) (CT = 56.3% ± 0.9%, non-CT = 43.7% ± 0.9%; mean ± SEM; n = 8 barrel columns, 4 hemispheres, and 3 mice)
        # S1_range_xz  = [center_point-(S1_size/2),center_point+(S1_size/2)]
        # PR5_range_xz = [center_point-(PR5_size/2),center_point+(PR5_size/2)]

        # xz_dimensions= {'VPM': VPM_range_xz,
        #                 'TRN': TRN_range_xz,
        #                 'L6A': L6A_range_xz,
        #                 'S1':S1_range_xz,
        #                 'PR5': PR5_range_xz,
        #                 }
        # y_dimensions = {
        #                 'VPM':   [3300,4000],   # (Claus, 2018) - barrel height of 700 um
        #                 'TRN':   [2550,2800],   # needs data
        #                 'TRNs1': [2800,2900],   # (Li, 2020) Core/shell structure in TRN
        #                 'TRNs2': [2450,2550],   # (Li, 2020) Core/shell structure in TRN
        #                 # 'L6A':   [1000,1250],   # (Crandall, 2017) - Layer 6A infrabarrels only - measured using digitizer tool
        #                 'L6A':   [890,1090],    # (Crandall, 2017) - Layer 6A infrabarrels only - measured using digitizer tool
        #                 'PR5':   [4500,5700],   # (Timofeeva, 2003) - Measurements of Barrelettes shown in the discussion (diam ~55um / length ~1.2mm)

        #                 'L1':    [0,128],       # (Lefort, 2009) Mouse Layer dimensions
        #                 'L23':    [128,418],     # (Lefort, 2009) Mouse Layer dimensions
        #                 # 'L2':    [128,269],     # (Lefort, 2009) Mouse Layer dimensions
        #                 # 'L3':    [269,418],     # (Lefort, 2009) Mouse Layer dimensions
        #                 'L4':    [418,588],     # (Lefort, 2009) Mouse Layer dimensions
        #                 'L5A':   [588,708],     # (Lefort, 2009) Mouse Layer dimensions
        #                 'L5B':   [708,890],     # (Lefort, 2009) Mouse Layer dimensions
        #                 'L6B':    [1090,1154],    # (Lefort, 2009) Mouse Layer dimensions
        #                 }

from neuron import h
import matplotlib.pyplot as plt
import numpy as np
import NetPyNE_BBP
import json
import pickle
import copy
import math 

########################################################################################################################################################################################################

class ResampleParameters:

    def truncated_normal_neuron(mean, std_dev, num_samples, rand=None, seed=100000, verbose=False):
        # Create a NEURON Random object and set the seed for reproducibility
        if rand is None:
            if seed is not None:    rand = h.Random(seed)
            else:                   rand = h.Random(1000)

        # Create a Vector to store the samples
        samples = h.Vector()

        # Generate random samples from the truncated normal distribution using h.Random().normal
        while len(samples) < num_samples:
            value = rand.normal(mean, std_dev**2)  # Pass variance instead of standard deviation
            if (value > 0) and value>(mean-(3*std_dev)) and value<(mean+(3*std_dev)):
                samples.append(value)
        result = list(samples)

        # Calculate and print the mean and standard deviation of the sampled data
        if verbose:
            sampled_mean    = np.mean(result)
            sampled_std_dev = np.std(result)
            print('\n\t Number of samples',num_samples)
            print(f"\t\t Sampled Mean: {sampled_mean:.4f},   Reference Mean: {mean:.4f}")
            print(f"\t\t Sampled Std:  {sampled_std_dev:.4f},   Reference Std: {std_dev:.4f}")

        # Return the samples as a Python list
        return result

    def plot_truncated_normal(result, mean, std_dev, savefig=None):
        # Plotting the histogram of sampled values
        plt.figure(figsize=(10,10))
        plt.hist(result, bins=30, density=True, alpha=0.7, color='blue', label='Sampled Values')

        # Plot the original normal distribution for comparison
        x = np.linspace(max(0,(mean-(5*std_dev))) , (mean+(5*std_dev)), 1000)
        y = (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-(x - mean)**2 / (2 * std_dev**2))
        plt.plot(x, y, color='red', label='Original Distribution')
        plt.title('Truncated Normal Distribution')
        plt.xlabel('Value')
        plt.ylabel('Probability Density')
        plt.legend()
        if savefig is not None: plt.savefig(savefig)
    
    def plot_conn_weights(sim,sorted_conn_dist,allCells,result,synMech_name,var,mean,std_dev,thresh=None,plotEvery=1,verbose=False):
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        # --- Identifies pre and post pops
        prePop      = synMech_name.split('|')[1]
        postPop     = synMech_name.split('|')[2]
        
        pre_cell_gids=[]
        post_cell_gids=[]
        for result_index,[(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)] in enumerate(sorted_conn_dist):
            pre_cell_gids.append(pre_gid)
            post_cell_gids.append(cell_gid)
        pre_cell_gids=list(set(pre_cell_gids))
        post_cell_gids=list(set(post_cell_gids))

        # --- Scatter plot of pre and post cell positions
        pre_cell_xyz=[]
        for pre_cell_gid in pre_cell_gids:      pre_cell_xyz.append((allCells[str(pre_cell_gid)]['tags']['x'],allCells[str(pre_cell_gid)]['tags']['y'],allCells[str(pre_cell_gid)]['tags']['z']))
        post_cell_xyz=[]
        for post_cell_gid in post_cell_gids:    post_cell_xyz.append((allCells[str(post_cell_gid)]['tags']['x'],allCells[str(post_cell_gid)]['tags']['y'],allCells[str(post_cell_gid)]['tags']['z']))

        # --- Extracting first and last elements of each tuple
        pre_x_values  = [item[0] for item in pre_cell_xyz]
        pre_y_values  = [item[1] for item in pre_cell_xyz]
        post_x_values = [item[0] for item in post_cell_xyz]
        post_y_values = [item[1] for item in post_cell_xyz]

        # --- Creating a scatter plot of cell positions
        plt.figure(figsize=(10, 10))
        
        # --- Histogram of all resampled results
        plt.subplot(2,2,4)
        plt.hist(result,bins=100)
        plt.title('Hist: '+prePop+'-'+postPop)
        # Plot the original normal distribution for comparison
        x = np.linspace(max(0,(mean-(5*std_dev))) , (mean+(5*std_dev)), len(result))
        y = (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-(x - mean)**2 / (2 * std_dev**2))
        plt.plot(x, y, color='red', label='Original Distribution')
        plt.xlim((mean-(7*std_dev),mean+(7*std_dev)))


        plt.subplot(1,2,1)
        plt.scatter(pre_x_values, pre_y_values,   color='grey', alpha=0.1, s=1)
        plt.scatter(post_x_values, post_y_values, color='blue', alpha=0.1, s=1)
        plt.legend(['pre pop: '+prePop,'post pop: '+postPop])
        plt.title('Scatter Plot of Connection weights')
        plt.xlabel('First Element')
        plt.ylabel('Last Element')

        

        absolute_distances=[]
        for result_index,[(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)] in enumerate(sorted_conn_dist): absolute_distances.append(absolute_distance)
        max_absolute_distance = max(absolute_distances)
        min_absolute_distance = min(absolute_distances)
        
        max_linewidth = 0.1
        min_linewidth = 0.01

        max_alpha = 1.0
        # max_alpha = 1.0
        min_alpha = 0.05
        
        # norm = Normalize(vmin=min(absolute_distances), vmax=max(absolute_distances))
        # sm = ScalarMappable(cmap=cmap_name, norm=norm)

        plotted_result=[]
        
        cmap_name = 'jet'

        # norm = Normalize(vmin=min(result), vmax=max(result))
        # sm = ScalarMappable(cmap=cmap_name, norm=norm)

        cmap = plt.get_cmap(cmap_name)

        max_result = max(result)
        min_result = min(result)
        norm_result=[((r-min_result)/(max_result-min_result)) for r in result]

        for result_index,[(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)] in enumerate(sorted_conn_dist):
            if result_index%plotEvery!=0:continue
            conns = [conn for conn in allCells[str(cell_ind)]['conns']]
            
            # --- Only tries to modify the conn if the stored dataset matches the conn properties (same mechanism and preGid) - otherwise breaks
            target_conn = conns[conn_index]
            if (target_conn['synMech']==synMech_name) and (target_conn['preGid']==pre_gid):
                if verbose: print('plotting conn: ',target_conn)
            else:
                if verbose: print('plotting failed - invalid conn: ', target_conn)
                continue

            # --- The larger the value, the larger the line size 
            linewidth = 0.05
            # linewidth = ((absolute_distance-min_absolute_distance)/(max_absolute_distance-min_absolute_distance))*(max_linewidth-min_linewidth)
            # --- The larger the distance, the smaller the line alpha
            # alpha = (((absolute_distance-min_absolute_distance)/(max_absolute_distance-min_absolute_distance)))*(max_alpha-min_alpha)
            alpha=0.35
            # --- The larger the distance, the smaller the line alpha
            # norm_result = ((result[result_index]-min_result)/(max_result-min_result))
            norm_res = norm_result[result_index]
            
            if thresh is not None:
                if norm_res<thresh[0] or norm_res>thresh[1]:
                    # print('skipping plotting conn # ', result_index)
                    continue
                else:
                    # --- To plot sliced histogram later
                    plotted_result.append(result[result_index])

            # color = sm.to_rgba(norm_res)
            color = cmap(norm_res)

            # print(linewidth,alpha)
            x_coords = [allCells[str(pre_gid)]['tags']['x'],allCells[str(cell_gid)]['tags']['x']]
            y_coords = [allCells[str(pre_gid)]['tags']['y'],allCells[str(cell_gid)]['tags']['y']]
            plt.plot(x_coords,y_coords,linewidth=linewidth,alpha=alpha,color=color)

            # plt.plot(result[result_index])
        
        plt.xlim([100,900])
        plt.ylim([0,6000])
        plt.gca().invert_yaxis()

        if thresh is not None:
            # --- Histogram of all resampled results
            plt.subplot(2,2,2)
            plt.hist(plotted_result,bins=int(100*(thresh[1]-thresh[0])))
            plt.title('Thresh Hist: '+prePop+'-'+postPop)
            plt.xlim((mean-(7*std_dev),mean+(7*std_dev)))

        savefigName = 'figs_resampleParameters/validate_conn/validate_connectivity_'+var+'_'+prePop+'_'+postPop
        if thresh is not None:savefigName+='_thresh_'+str(thresh[0])+'_'+str(thresh[1])
        plt.savefig(savefigName+'.png',dpi=1000)

    # --- Convert dictionary to serializable (turns objects into strings)
    def to_dict(obj):
        return json.loads(json.dumps(obj, default=lambda o: o.__dict__))

    # --- Method to remove 'hObj' objects from dictionary and make it serializable across different nodes
    def remove_key_recursive(dictionary, key_to_remove):
        if not isinstance(dictionary, dict): return

        for key, value in list(dictionary.items()):
            if key == key_to_remove:            del dictionary[key]
            elif isinstance(value, dict):       ResampleParameters.remove_key_recursive(value, key_to_remove)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):  ResampleParameters.remove_key_recursive(item, key_to_remove)

    # ------------------------------------------------------------------------------
    # Gather cells from different nodes - based on <_gatherAllCellTags> from netpyne/sim/gather.py
    # ------------------------------------------------------------------------------
    def _gatherParallelCells(sim):

        print('Running gather for sim rank ', sim.rank, ' in ', sim.nhosts, ' hosts')

        # data = [{cell.gid: ResampleParameters.remove_key_recursive(sim.net.cells[ind], 'hObj') for ind, cell in enumerate(sim.net.cells)}] * sim.nhosts  # send cells data to other nodes
        sim_copy = copy.deepcopy(sim)
        data = [{cell.gid:{ 
                                'gid':      sim_copy.gid, 
                                'tags':     sim_copy.tags, 
                                'conns':    ResampleParameters.remove_key_recursive(sim_copy.net.cells[cell.gid].conns,'hObj'), 
                                'secs':     ResampleParameters.remove_key_recursive(sim_copy.net.cells[cell.gid].secs,'hObj')
                                } 
                                for ind, cell in enumerate(sim_copy.net.cells)}] * sim.nhosts  # send cells data to other nodes

                # ResampleParameters.remove_key_recursive(sim.net.cells[ind], 'hObj')

        # data_conns={}
        # for cell in sim.net.cells:
        #     if cell.gid not in data_conns.keys(): data_conns.update({cell.gid:{}})
            
        #     # # --- Creates a version of the conns that excludes the hObj, because those can't be synchronized and break the code
        #     # conns_list=[]
        #     # for conn_ind,conn in enumerate(cell.conns):
        #     #     conns_list.append({key:cell.conns[conn_ind][key] for key in cell.conns[conn_ind].keys() if 'hObj' not in key})

        #     # --- Creates a version of the secs that excludes the hObj, because those can't be synchronized and break the code
        #     sec_list=[]
        #     for sec_ind,sec in enumerate(cell.secs):
        #         sec_list.append({key:cell.secs[sec][key] for key in cell.secs[sec].keys() if 'hObj' not in key})

        #     # conns_list = ResampleParameters.to_dict(cell.conns)

        #     data_conns[cell.gid]={'gid': cell.gid, 'tags': cell.tags, 'conns': conns_list}
        # data = [data_conns] * sim.nhosts  # send cells data to other nodes


        # data = [{cell.gid: ResampleParameters.to_dict(cell) for cell in sim.net.cells}] * sim.nhosts  # send cells data to other nodes

        # data = [{cell.gid: ResampleParameters.to_dict(cell) for cell in sim.net.cells}] * sim.nhosts  # send cells data to other nodes
        
        
        # for cell in sim.net.cells:
        #     serializable_data = json.dumps(cell, default=lambda o: o.__dict__)
        #     serializable_conns = json.dumps(cell.conns, default=lambda o: o.__dict__)

        # data_conns={}
        # for cell in sim.net.cells:
        #     if cell.gid not in data_conns.keys(): data_conns.update({cell.gid:{}})
        #     conns_list=[]
        #     # --- Creates a version of the conns that excludes the hObj, because those can't be synchronized and break the code
        #     for conn_ind,conn in enumerate(cell.conns):
        #         conns_list.append({key:cell.conns[conn_ind][key] for key in cell.conns[conn_ind].keys() if 'hObj' not in key})
        #     data_conns[cell.gid]=conns_list

        # data = [data_conns] * sim.nhosts  # send cells data to other nodes

        gather = sim.pc.py_alltoall(data)  # collect cells data from other nodes (required to generate connections)
        sim.pc.barrier()
        allCells = {}
        for dataNode in gather:
            allCells.update(dataNode)

        # clean to avoid mem leaks
        for node in gather:
            if node:
                node.clear()
                del node
        for item in data:
            if item:
                item.clear()
                del item

        return allCells

    def _gatherAllCellConns(sim):
    
        data_conns={}
        for cell in sim.net.cells:
            if cell.gid not in data_conns.keys(): data_conns.update({cell.gid:{}})
            conns_list=[]
            # --- Creates a version of the conns that excludes the hObj, because those can't be synchronized and break the code
            for conn_ind,conn in enumerate(cell.conns):
                conns_list.append({key:cell.conns[conn_ind][key] for key in cell.conns[conn_ind].keys() if ('hObj' not in key) or ('hObj' not in cell.conns[conn_ind][key].keys())})
            data_conns[cell.gid]=conns_list

        data = [data_conns] * sim.nhosts  # send cells data to other nodes
        
        gather = sim.pc.py_alltoall(data)  # collect cells data from other nodes (required to generate connections)
        sim.pc.barrier()
        allCellConns = {}
        for dataNode in gather:
            allCellConns.update(dataNode)

        # clean to avoid mem leaks
        for node in gather:
            if node:
                node.clear()
                del node
        for item in data:
            if item:
                item.clear()
                del item

        return allCellConns

    def _gatherAllCellProperties2(sim): # worked - but very slow because it gathers all information from the network
        
        cells = [c.__getstate__() for c in sim.net.cells]

        data = [{cell.gid: cell for cell in cells}] * sim.nhosts  # send cells data to other nodes
        gather = sim.pc.py_alltoall(data)  # collect cells data from other nodes (required to generate connections)
        sim.pc.barrier()
        allCellTags = {}
        for dataNode in gather:
            allCellTags.update(dataNode)

        # clean to avoid mem leaks
        for node in gather:
            if node:
                node.clear()
                del node
        for item in data:
            if item:
                item.clear()
                del item

        return allCellTags

    def divide_chunks(l, n): 
        
        # looping till length l 
        for i in range(0, len(l), n):  
            yield l[i:i + n] 

    

    def chunk_into_n(lst, n):
        size = math.ceil(len(lst) / n)
        return list(
            map(lambda x: lst[x * size:x * size + size],
            list(range(n)))
        )

    def _saveCellData(sim,num_conn_files=False):
        
        if sim.nhosts > 1: 
            print('Code cannot be executed in parallel mode\nRun network creation step first to enable synapse resampling')
            return
        try: 
            cfg_dict = sim.cfg.__dict__
            filename = cfg_dict['dump_cell_properties']
            filename_json = filename
            filename_secs = cfg_dict['dump_cell_properties'].split('.json')[0]+'__secs.pkl'
            filename_conns = cfg_dict['dump_cell_properties'].split('.json')[0]+'__conns.pkl'
            if sim.rank==0: print('Saving network template in cfg path: ', filename)
        except:
            try:
                filename = '../network_template/barreloid_network_cell_properties.json'
                filename_json = '../network_template/barreloid_network_cell_properties.json'
                filename_secs = '../network_template/barreloid_network_cell_properties__secs.pkl'
                filename_conns = '../network_template/barreloid_network_cell_properties__conns.pkl'
                if sim.rank==0: print('Saving network template in default path: ', filename)
            except:
                filename = 'barreloid_network_cell_properties.json'
                filename_json = 'barreloid_network_cell_properties.json'
                filename_secs = 'barreloid_network_cell_properties__secs.pkl'
                filename_conns = 'barreloid_network_cell_properties__conns.pkl'
                if sim.rank==0: print('Saving network template in current folder: ', filename)

        cell_data={}
        for cell in sim.net.cells:
            try:    conns_list = [{'preGid': conn['preGid'],'synMech': conn['synMech'],'sec': conn['sec']} for conn_ind,conn in enumerate(cell.conns)]
            except: conns_list = []

            if cell.gid not in cell_data.keys(): cell_data.update({str(cell.gid):{  'gid':     cell.gid,
                                                                                    'tags':    cell.tags,
                                                                                    'conns':   conns_list,
                                                                                    }})
        # with open(filename, 'wb') as file:
        #     pickle.dump(cell_data, file)
        with open(filename_json, 'w') as file:
            json.dump(cell_data, file)
        if sim.rank==0: print(f"Cell properties dumped into {filename} successfully.")

        return cell_data

    def _loadCellData(sim,fileFormat='json'):

        try: 
            cfg_dict = sim.cfg.__dict__
            filename = cfg_dict['dump_cell_properties']
            filename_json = filename
            # filename_conns = cfg_dict['dump_cell_properties'].split('.json')[0]+'__conns.pkl'
            if sim.rank==0: print('Loading network template in cfg path: ', filename)
        except:
            try:
                filename = '../network_template/barreloid_network_cell_properties.json'
                filename_json = '../network_template/barreloid_network_cell_properties.json'
                # filename_conns = '../network_template/barreloid_network_cell_properties__conns.pkl'
                if sim.rank==0: print('Loading network template in default path: ', filename)
            except:
                filename = 'barreloid_network_cell_properties.json'
                filename_json = 'barreloid_network_cell_properties.json'
                # filename_conns = 'barreloid_network_cell_properties__conns.pkl'
                if sim.rank==0: print('Loading network template in current folder: ', filename)

        if (fileFormat == 'json') or (fileFormat == '.json'):
            with open(filename_json, 'r') as file:
                cell_data = json.load(file)
        elif (fileFormat == 'pkl') or (fileFormat == '.pkl'):
            with open(filename, 'rb') as file:
                cell_data = pickle.load(file)
        else:
            try:
                with open(filename_json, 'r') as file:
                    cell_data = json.load(file)
            except:
                if sim.rank==0: print('Loading data failed')
                cell_data={}
        
        # with open(filename_conns, 'rb') as file:
        #     cell_data_conns = pickle.load(file)
        
        # # --- Add conns back to the returned dataset 
        # for gid in cell_data.keys(): cell_data[gid].update({'conns':cell_data_conns[gid]['conns']})

        print('\t --- Loaded network template file: ', filename_json)

        if sim.rank==0: print('Cell properties loaded successfully.')

        return cell_data
 
    def _gatherAllCellProperties(sim): # worked
    
        cell_data={}
        for cell in sim.net.cells:
            secs = cell.secs
            # --- compartment cells vs point cells
            try:    secs_dict  = {sec:{'geom':cell.secs[sec]['geom'],'topol':cell.secs[sec]['topol']}    for sec in secs.keys()}
            except: secs_dict  = {}
            try:    conns_list = [{'preGid': conn['preGid'],'synMech': conn['synMech'],'sec': conn['sec']} for conn_ind,conn in enumerate(cell.conns)]
            except: conns_list = []

            if cell.gid not in cell_data.keys(): cell_data.update({cell.gid:{'gid':     cell.gid,
                                                                             'tags':    cell.tags,
                                                                             'secs':    secs_dict,
                                                                             'conns':   conns_list,
                                                                             }})

        data = [cell_data] * sim.nhosts  # send cells data to other nodes
        
        gather = sim.pc.py_alltoall(data)  # collect cells data from other nodes (required to generate connections)
        sim.pc.barrier()
        allCellConns = {}
        for dataNode in gather:
            allCellConns.update(dataNode)

        # clean to avoid mem leaks
        for node in gather:
            if node:
                node.clear()
                del node
        for item in data:
            if item:
                item.clear()
                del item

        return allCellConns

    def processSimObj(sim,updateNetwork=True,weightSorting=True, weightSorting_ExcludeConns=None, verbose=False,plotFig=True):
        if sim.rank==0: print('\t>> Resampling Conn parameters to add mean+-std variability from literature ', sim.rank)

        if updateNetwork:
            if sim.nhosts==1:
                print('Updating network template')
                allCells = ResampleParameters._saveCellData(sim)
            else:
                if sim.rank==0: print('Must be run in a single core, but number of cores is ', sim.nhosts)
                return
        else:
            try: allCells = ResampleParameters._loadCellData(sim)
            except:
                if sim.rank==0: print('Loading cell properties failed\nAttempting to create the dictionary of cell properties')
                if sim.nhosts==1:
                    allCells = ResampleParameters._saveCellData(sim)
                else:
                    allCells = ResampleParameters._gatherAllCellProperties(sim)

        # --- Converting cfg to dict to process variables
        cfg_dict = sim.cfg.__dict__ 

        # --- Stores the GID of the cells that are actually present in the given node - in case of parallel simulation
        cells_in_node = [sim.net.cells[cell_i].gid for cell_i,cell in enumerate(sim.net.cells)]
        # --- Matches the GID:LID of the cells that are actually present in the given node - in case of parallel simulation
        cell_lids = {sim.net.cells[cell_i].gid:cell_i for cell_i,cell in enumerate(sim.net.cells)}

        # if sim.rank == 1:
        #     print(allCellTags.keys())

        TableS1_barreloidThalamus   = NetPyNE_BBP.StoreParameters.getTableS1_barreloidThalamus()
        resample_vars               = ['conductance','u_syn','depression_time','facilitation_time']
        count_synMechs              = {}
        seed                        = cfg_dict['base_random_seed']
        rand                        = h.Random(seed)

        store_conn_gids     = {}
        store_conn_distance = {}

        store_conn_dist = {} # testing a simpler implementation of data storage
        
        # --- Loops through all cells
        store_path_distance={}
        if sim.rank==0: print('Calculating path Distance for secs in cells')
        for cell_ind in range(len(allCells)):
            cell_gid = allCells[str(cell_ind)]['gid']

            if 'cellModel' in allCells[str(cell_ind)]['tags'].keys():  continue # skips Vecstims
            else:
                cell_dict = sim.net.params['cellParams'][allCells[str(cell_ind)]['tags']['label'][0]]   # --- finds the cell morphology based on the cell template, to avoid having to store the morphology of each cell
                store_path_distance.update({cell_gid:NetPyNE_BBP.Utils.pathDistance(cell_dict)})        # --- calculates the path distance for each cell/section

        # --- Loops through all cells
        if sim.rank==0: print('Processing connectivity')
        for cell_ind in range(len(allCells)):
            cell_gid = allCells[str(cell_ind)]['gid'] # --- IMPORTANT: <CELL_GID> IS THE GLOBAL GID FOR THE SIMULATION, AND NOT THE GID AT EACH RANK

            count_synMechs.update({cell_gid:{}})        # --- count the conns (each synapse) and adds an individual index to each one
            store_conn_gids.update({cell_gid:{}})
            # --- Loops through all conns
            for conn_ind in range(len(allCells[str(cell_ind)]['conns'])):
                # --- Skips changing NetStims for now
                preGid = allCells[str(cell_ind)]['conns'][conn_ind]['preGid']
                if preGid == 'NetStim': continue

                synMech_name = allCells[str(cell_ind)]['conns'][conn_ind]['synMech']
                
                # --- Condition to enforce that gap junctions are not modified
                if '|gap' in synMech_name: continue

                if weightSorting and (weightSorting_ExcludeConns is not None):
                    if synMech_name in weightSorting_ExcludeConns: 
                        # print('Excluding conn with synMech: ', cell_ind, ' | ', conn_ind, ' | ', allCells[str(cell_ind)]['conns'][conn_ind]['synMech'])
                        continue # --- skips the conn if the synMech is in the list of excluded conns
                
                # --- Modify only the chemical synaptic mechanisms from Thalamus<->Cortex 
                if 'syn|' in synMech_name:
                    # print('Modifying gid/conn ', allCells[str(cell_ind)]['gid'], '\t',  allCells[str(cell_ind)]['conns'][conn_ind], '\t mech: ', synMech_name)
                    # --- Stores the index of the connection the sim obj for referencing (so that it skips the conns from other mechs)
                    if synMech_name not in count_synMechs[cell_gid].keys(): count_synMechs[cell_gid].update({synMech_name:[conn_ind]})
                    else:                                                   count_synMechs[cell_gid][synMech_name].append(conn_ind)

                    if synMech_name not in store_conn_gids[cell_gid].keys(): store_conn_gids[cell_gid].update({synMech_name:[preGid]})
                    else:                                                    store_conn_gids[cell_gid][synMech_name].append(preGid)
                    

            # ---  Dictionary to store the distance between the cell somas for each synaptic contact
            store_conn_distance.update({cell_gid:{}})
            # --- Measure distance between cell somas to rank conn weights
            # --- Loops through the stored mechs
            for synMech_name in store_conn_gids[cell_gid].keys():

                if synMech_name not in store_conn_dist.keys():store_conn_dist.update({synMech_name:[]})
                
                store_conn_distance[cell_gid].update({synMech_name:{}})

                # --- Identifies pre and post pops
                prePop      = synMech_name.split('|')[1]
                postPop     = synMech_name.split('|')[2]
                try:
                    pre_pop_axis  = cfg_dict[prePop+'_'+postPop+'__pre_pop_axis']
                    post_pop_axis = cfg_dict[prePop+'_'+postPop+'__post_pop_axis']
                except:
                    pre_pop_axis  = 'y'
                    post_pop_axis = 'y'
                
                if verbose: print('CELL INDEX: ', cell_ind, ' - CELL GID ACROSS CORES: ', cell_gid,'\t',prePop,'-',pre_pop_axis,'\t',postPop,'-',post_pop_axis)               

                import math
                from Build_Net import BuildNetwork
                post_gid=cell_gid
                # post_gid=cell_ind

                # print('all conns: ', allCells[str(cell_ind)]['conns'])
                # print(sim.net.cells,'\n\n',store_conn_gids[cell_ind],'\n\n\n\n')
                for gid_index, pre_gid in enumerate(store_conn_gids[cell_gid][synMech_name]):

                    # print(pre_gid)
                    
                    conn_index = count_synMechs[cell_gid][synMech_name][gid_index]

                    if pre_gid not in store_conn_distance[cell_gid][synMech_name].keys():   store_conn_distance[cell_gid][synMech_name].update({pre_gid:[]})

                    diam,height           = BuildNetwork.getNetworkSize(center_point=cfg_dict['center_point'])
                    pre_pop_boundaries,post_pop_boundaries = BuildNetwork.getBoudaries(pre_pop=prePop,post_pop=postPop,pre_pop_axis=pre_pop_axis,post_pop_axis=post_pop_axis,center_point=cfg_dict['center_point'])

                    if 'theta' in pre_pop_axis:     pre_cell_norm_dist  = ((np.remainder((((np.arctan2(allCells[str(pre_gid)]['tags']['x']- pre_pop_boundaries[1], allCells[str(pre_gid)]['tags']['z']- pre_pop_boundaries[1]))*(180/np.pi))+360),360))/360)
                    else:                           pre_cell_norm_dist  = ((allCells[str(pre_gid)]['tags'][pre_pop_axis]-pre_pop_boundaries[1])/(pre_pop_boundaries[2]-pre_pop_boundaries[3]))

                    if 'theta' in post_pop_axis:    post_cell_norm_dist = ((np.remainder((((np.arctan2(allCells[str(post_gid)]['tags']['x']- post_pop_boundaries[1], allCells[str(post_gid)]['tags']['z']- post_pop_boundaries[1]))*(180/np.pi))+360),360))/360)
                    else:                           post_cell_norm_dist = ((allCells[str(post_gid)]['tags'][post_pop_axis]-post_pop_boundaries[1])/(post_pop_boundaries[2]-post_pop_boundaries[3]))

                    # if 'theta' in pre_pop_axis:     pre_cell_norm_dist  = ((np.remainder((((np.arctan2(sim.net.cells[pre_gid].tags['x']- pre_pop_boundaries[1], sim.net.cells[pre_gid].tags['z']- pre_pop_boundaries[1]))*(180/np.pi))+360),360))/360)
                    # else:                           pre_cell_norm_dist  = ((sim.net.cells[pre_gid].tags[pre_pop_axis]-pre_pop_boundaries[1])/(pre_pop_boundaries[2]-pre_pop_boundaries[3]))

                    # if 'theta' in post_pop_axis:    post_cell_norm_dist = ((np.remainder((((np.arctan2(sim.net.cells[post_gid].tags['x']- post_pop_boundaries[1], sim.net.cells[post_gid].tags['z']- post_pop_boundaries[1]))*(180/np.pi))+360),360))/360)
                    # else:                           post_cell_norm_dist = ((sim.net.cells[post_gid].tags[post_pop_axis]-post_pop_boundaries[1])/(post_pop_boundaries[2]-post_pop_boundaries[3]))

                    absolute_distance = abs(pre_cell_norm_dist-post_cell_norm_dist)

                    # print('conn_ind: ', conn_ind, ' | conn_index: ', conn_index, ' | len conns: ', len(allCells[str(cell_ind)]['conns']))
                    # print(allCells.keys())
                    conn = allCells[str(cell_ind)]['conns'][conn_index]
                    # conn = allCells[str(cell_ind)]['conns'][conn_ind][conn_index]
                    # print('passed condition: ', sim.rank)
                    post_sec_dist = store_path_distance[cell_gid][conn['sec']]
                    # print('conn # ',conn_index,': ',conn, ' - post sec dist: ', post_sec_dist)

                    if verbose: print('\n',synMech_name, conn_index, post_gid, '\t',prePop, pre_gid, pre_pop_axis, ': ', pre_cell_norm_dist, ' | ', postPop, post_gid, post_pop_axis, ': ', post_cell_norm_dist,   '\t Absolute distance: ' ,absolute_distance, ' - post sec dist: ', post_sec_dist)

                    post_sec_name = conn['sec']
                    store_conn_distance[cell_gid][synMech_name][pre_gid].append((absolute_distance,post_sec_dist,post_sec_name))
                    
                    store_conn_dist[synMech_name].append([(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)])

                    # pre_cell_dict = {'secs':sim.net.cells[cell_ind].secs}
                    # cell_dict['secs']

                    # NetPyNE_BBP.Utils.pathDistance({'secs':sim.net.cells[cell_ind].secs})

        if weightSorting:
            if sim.rank==0: print('Sorting stored conns based on cell alignment')
            sorted_store_conn_dist={}
            for synMech_name in store_conn_dist.keys():
                sorted_data = sorted(store_conn_dist[synMech_name], key=lambda x: (x[1][0:2])) # sorts based on (1st -> absolute_distance / 2nd -> post_sec_dist)
                sorted_store_conn_dist.update({synMech_name:sorted_data})
        else:
            sorted_store_conn_dist = copy.deepcopy(store_conn_dist)
        
        # --- Loops through the stored mechs
        if sim.rank==0: 
            print('Modifying syn Mechs')
            if weightSorting:
                if weightSorting_ExcludeConns is not None: print('\tExcept ', str(weightSorting_ExcludeConns))
        # print(sorted_store_conn_dist.keys())
        for synMech_name in sorted_store_conn_dist.keys():
            # --- Identifies the number of repeats of the same mech for unified resampling
            num_samples = len(sorted_store_conn_dist[synMech_name])
            # print('sim rank: ',sim.rank,'\tnum_samples: ', num_samples)
            # --- Identifies pre and post pops
            prePop      = synMech_name.split('|')[1]
            postPop     = synMech_name.split('|')[2]
            # --- Loops through the variables to be resampled
            for var in resample_vars:
                mean_value      = TableS1_barreloidThalamus[prePop][postPop]['paper_reference_values'][var]
                std_dev_value   = TableS1_barreloidThalamus[prePop][postPop]['paper_reference_values'][var+'_std']
                # print('Mean/std values - ',prePop,'->',postPop,': ', mean_value,'\t',std_dev_value)

                # --- Resampled values
                result = ResampleParameters.truncated_normal_neuron(mean_value, std_dev_value, num_samples, rand, seed=seed, verbose=False)


                ''' 
                Removing conductance sorting, since it introduces a bug, making the network biased in the y-axis, with cells receiving decreasing weights based on GID/LID
                '''

                # if var == 'conductance': 
                #     print('Sorting conductance values!!!!!!')
                #     result.sort(reverse = True) # --- sorts only the conductance, and allows the other variables to have a random distriution of values



                # print(sim.rank, ' | ', var, '\n result value: ',result[0:50])
                for result_index,[(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)] in enumerate(sorted_store_conn_dist[synMech_name]):
                    if verbose: print([(cell_gid,pre_gid,conn_index,cell_ind),(absolute_distance,post_sec_dist,post_sec_name)])

                    if cell_gid not in cells_in_node:continue

                    # conns = [conn for conn in sim.net.cells[cell_ind].conns if (conn['synMech']==synMech_name) and (conn['preGid']==pre_gid)]
                    conns = [conn for conn in allCells[str(cell_ind)]['conns']]
                    
                    # --- Only tries to modify the conn if the stored dataset matches the conn properties (same mechanism and preGid) - otherwise breaks
                    target_conn = conns[conn_index]
                    if (target_conn['synMech']==synMech_name) and (target_conn['preGid']==pre_gid):
                        if verbose: print(target_conn)
                    else:
                        if verbose: print('invalid conn: ', target_conn)
                        continue

                    cell_lid = cell_lids[cell_gid]

                    # print(cell_lid, '\t' ,var,'\t new value: ', result[result_index])

                    if   var == 'conductance':
                        # # print('before: ', sim.net.cells[cell_ind].conns[conn_index]['hObj'].syn().gmax, ' | distance: ', absolute_distance)
                        # print('conn index: ', conn_index, ' | Result index: ', result_index, ' | Result: ', result[result_index])
                        # print('Conn: ', sim.net.cells[cell_lid].conns[conn_index])
                        sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().gmax   = result[result_index]
                        # print('after : ', sim.net.cells[cell_ind].conns[conn_index]['hObj'].syn().gmax, ' | distance: ', absolute_distance)
                        # if 'MLe' in sim.net.cells[pre_gid].tags['pop']:print('pre pop = MLe | conductance = ', sim.net.cells[cell_ind].conns[conn_index]['hObj'].syn().gmax)
                    elif var == 'u_syn':                
                        if 'rescaleUSE' in cfg_dict.keys(): 
                            sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use    = result[result_index]*cfg_dict['rescaleUSE']
                            # print('rescaling USE after resampling from source Table: ', result[result_index], '\t',cfg_dict['rescaleUSE'], '\t',result[result_index]*cfg_dict['rescaleUSE'])
                        else:                               sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use    = result[result_index]

                        if 'addUSE_MLe' in cfg_dict.keys(): 
                            if 'syn|MLe' in sim.net.cells[cell_lid].conns[conn_index]['synMech']:
                                # print('adding fixed value to MLe syn mech USE: ', sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use)
                                sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use += cfg_dict['addUSE_MLe']
                                if sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use > 1: sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use=1
                                if sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use < 0: sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use=0
                                # print('added  fixed value to MLe syn mech USE: ', sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Use, '\n\n')

                    elif var == 'depression_time':      
                        if 'syn|MLe' in sim.net.cells[cell_lid].conns[conn_index]['synMech']:
                            if 'rescaleDep_MLe' in cfg_dict.keys(): sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Dep    = result[result_index] * cfg_dict['rescaleDep_MLe']
                            else:                                   sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Dep    = result[result_index]
                        else:                                       sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Dep    = result[result_index]
                    
                    elif var == 'facilitation_time':    
                        if 'syn|MLe' in sim.net.cells[cell_lid].conns[conn_index]['synMech']:
                            if 'rescaleFac_MLe' in cfg_dict.keys(): sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Fac    = result[result_index] * cfg_dict['rescaleFac_MLe']
                            else:                                   sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Fac    = result[result_index]
                        else:                                       sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().Fac    = result[result_index]
                    
                    else:                               print('invalid variable --> ', var)
                    
                    if verbose:
                        if var == 'conductance': print(cell_lid, ' - pre pop = ',prePop,'\t post pop = ',postPop,' | conductance = ', sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().gmax, '\t| absolute distance: ', absolute_distance)
                    
                    # if var == 'conductance': print(cell_lid, ' - pre pop = ',prePop,'\t post pop = ',postPop,' | conductance = ', sim.net.cells[cell_lid].conns[conn_index]['hObj'].syn().gmax, '\t| absolute distance: ', absolute_distance)

                    if verbose:
                        try:    print(conns[conn_index])
                        except: print(conn_index,' | ',conns)

                if plotFig:
                    if var == 'conductance':
                        for thresh in [[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4],[0.4,0.5],[0.5,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]:
                            # plot_conn_weights(sim,sorted_store_conn_dist[synMech_name],result,synMech_name,var,mean=mean_value,std_dev=std_dev_value,thresh=thresh,plotEvery=1)
                            ResampleParameters.plot_conn_weights(sim,sorted_store_conn_dist[synMech_name],allCells,result,synMech_name,var,mean=mean_value,std_dev=std_dev_value,thresh=thresh,plotEvery=1)

        return sim

########################################################################################################################################################################################################
if __name__ == '__main__':
    print('\t>> Running code to test resampling of synaptic parameters - see output in folder /figs_resampleParameters')
    TableS1_barreloidThalamus   = NetPyNE_BBP.StoreParameters.getTableS1_barreloidThalamus()
    resample_vars               = ['conductance','u_syn','depression_time','facilitation_time']
    for source in TableS1_barreloidThalamus.keys():
        # print(source)
        for target in TableS1_barreloidThalamus[source].keys():
            # print('\t', target)
            for var in resample_vars:
                mean_value      = TableS1_barreloidThalamus[source][target]['paper_reference_values'][var]
                std_dev_value   = TableS1_barreloidThalamus[source][target]['paper_reference_values'][var+'_std']
                # print('\t\t', var,': ',mean_value, '\t std: ', std_dev_value)

                num_samples_value = 100000  # Increase the number of samples for better visualization
                seed=100000
                result = ResampleParameters.truncated_normal_neuron(mean_value, std_dev_value, num_samples_value, rand=h.Random(seed), verbose=False)

                ResampleParameters.plot_truncated_normal(result, mean_value, std_dev_value, savefig='figs_resampleParameters/validate_resampling_'+source+'_'+target+'_'+var+'.png')




'''

# cell_data={}
# for cell in sim.net.cells:
#     if cell.gid not in cell_data.keys(): cell_data.update({cell.gid:{'gid':     cell.gid,
#                                                                      'tags':    cell.tags,
#                                                                      'secs':    sim.copyRemoveItemObj(cell.secs, keystart='hObj', exclude_list=['hebbwt']),
#                                                                      'conns':   sim.copyRemoveItemObj(cell.conns, keystart='hObj', exclude_list=['hebbwt'])}})
    


    # # --- Creates a version of the conns that excludes the hObj, because those can't be synchronized and break the code
    # secs_list=[]
    # for sec_ind,sec in enumerate(cell.secs): 
    #     secs_dict = sim.copyRemoveItemObj(cell.secs, keystart='hObj', exclude_list=['hebbwt'])  # , newval=None)  # replace h objects with None so can be pickled
    # cell_data[cell.gid]['secs']=secs_dict


    #     # for key in cell.secs[sec].keys():
    #     #     if ('hObj' not in key):
    #     #         if type(cell.secs[sec][key])==dict:
    #     #             if ('hObj' in cell.secs[sec][key].keys()):continue

    #     #         secs_list.append({key:cell.secs[sec][key]})
        
    #     # secs_list.append({key:cell.secs[sec][key] for key in cell.secs[sec].keys() if ('hObj' not in key) or ('hObj' not in cell.secs[sec][key].keys())})
        
    # # cell_data[cell.gid]['secs']=secs_list
    
    # # --- Creates a version of the conns that excludes the hObj, because those can't be synchronized and break the code
    # conns_list=[]
    # for conn_ind,conn in enumerate(cell.conns):
        
    #     for key in cell.conns[conn_ind].keys():
    #         if ('hObj' not in key): 
    #             conns_list.append({key:cell.conns[conn_ind][key]})
        
    #     # conns_list.append({key:cell.conns[conn_ind][key] for key in cell.conns[conn_ind].keys() if ('hObj' not in key) or ('hObj' not in cell.conns[conn_ind][key].keys())})

    # cell_data[cell.gid]['conns']=conns_list
'''
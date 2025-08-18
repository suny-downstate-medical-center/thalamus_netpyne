from matplotlib import pyplot as plt
import numpy as np
import json

class AnalyzeCurrents():

    def loadSimObj(filename):
        with open(filename, 'r') as fileObj: simFile=fileObj.read()
        sim = json.loads(simFile)
        return sim

    def runAnalysis(sim, 
                    current_flag='i__',
                    select_pops=None,
                    # select_currents=['i__soma_0__kap__ik', 'i__soma_0__kdr__ik'],
                    select_currents=None,
                    time_range=None,
                    singleTraces=False,
                    allTraces=False,
                    separateFigs=True,
                    figSize=(70,15),
                    figAlpha=0.5,
                    savefig_name = 'current_analysis',
                    savefig_format = '.png',
                    ):
        # Parse sim data from Obj or dict format
        try:    
            print('Loading data from OBJ')
            simData_all = sim.simData
            cfg         = sim.cfg
            net         = sim.net
            pops_list   = list(sim.net.pops.keys())
            pops_gids   = {pop:sim.net.pops[pop].cellGids for pop in pops_list}
            print('Loading data from OBJ - end')
        except: 
            print('Loading data from DICT')
            simData_all = sim['simData']
            try:    cfg = sim['simConfig']
            except: cfg = sim['cfg']
            net         = sim['net']
            pops_list   = list(sim['net']['pops'].keys())
            pops_gids   = {pop:sim['net']['pops'][pop]['cellGids'] for pop in pops_list}
            print('Loading data from DICT - end')

        # time_vec=simData['t']
        time_vec_=list(simData_all['t'])

        # time_step = cfg['dt']
        time_step = cfg['recordStep']
        
        try:    skipTime = cfg['skipTime']
        except: skipTime = 0
        
        # --- Checks if there are any currents that match the current flag
        any_currents=False
        for recordedTrace in simData_all.keys():
            if recordedTrace.startswith(current_flag):
                any_currents=True
                break
        if not any_currents: 
            print('No currents to plot')
            return
        
        # --- Lists the populations in <sim> that will be considered for plotting, as specified in select_pops
        if select_pops is not None: plot_pops = [pop for pop in pops_list if pop in select_pops]
        else:                       plot_pops = pops_list
        plot_pops_dict = {pop:pops_gids[pop] for pop in plot_pops}
        
        # --- Lists the recorded currents that will be considered for plotting, as specified in select_currents
        if select_currents is not None: plot_currents = [curr_name for curr_name in cfg['recordTraces'] if (current_flag in curr_name) and (curr_name in select_currents)]
        else:                           plot_currents = [curr_name for curr_name in cfg['recordTraces'] if (current_flag in curr_name)]

        # --- Filtering simData to include only current traces, and the ones that should be plotted
        simData_ = {trace:simData_all[trace] for trace in simData_all.keys() if (trace in plot_currents)}

        if time_range is not None:
            # try:
            time_range_indexes=[int(t_range/time_step) for t_range in time_range]
            simData={}
            for trace in simData_.keys():
                if trace not in simData.keys():simData.update({trace:{}})
                for cell in simData_[trace].keys():
                    simData[trace].update({cell:list(simData_[trace][cell])[time_range_indexes[0]:time_range_indexes[1]]})

            time_vec = time_vec_[time_range_indexes[0]:time_range_indexes[1]]
            # savefig_format='_timeRange_'+str(time_range[0])+'__'+str(time_range[1])+savefig_format
            savefig_name=savefig_name+'_timeRange_'+str(time_range[0])+'__'+str(time_range[1])

            # except: 
            #     print('Time range selection failed')
            #     simData  = simData_
            #     time_vec = time_vec_
        else: 
            simData  = simData_
            time_vec = time_vec_

        store_mean_currents_byPop_dict={} # dictionary to store mean and std values for the currents over time, across all cells in each pop

        for pop in plot_pops:
            pop_gids = plot_pops_dict[pop]

            current_colors_dict = AnalyzeCurrents.get_current_colors_dict(figAlpha = figAlpha)
            current_handles=[]
            for current_color in current_colors_dict.keys(): current_handles.append(plt.Line2D([], [], color=current_colors_dict[current_color], marker="o", linewidth=0, alpha=figAlpha, markeredgecolor='k'))
            
            # --- Plots each trace separately - using the mapped colors
            if singleTraces:    AnalyzeCurrents.plotSingleTraces(simData, time_vec, current_colors_dict, current_flag, pop_gids, pop, savefig_name, savefig_format)
            # --- Plots all traces together - using the mapped colors
            if allTraces:       AnalyzeCurrents.plotAllTraces(simData, time_vec, current_colors_dict, current_handles, current_flag, pop_gids, pop, savefig_name, savefig_format)
            # --- Generates a list of traces and label/color info for plotting
            store_currents, store_labels, store_colors  = AnalyzeCurrents.storePlotInfo(simData, current_colors_dict, current_flag, pop_gids)
            store_currents, store_colors                = AnalyzeCurrents.sortCurrents(store_currents, store_labels, store_colors)

            store_currents_np = np.array(store_currents)
            zeros_trace = [0 for i in range(len(time_vec))]

            if separateFigs: plt.figure(figsize=figSize)
            else:
                fig = plt.subplots(figsize=figSize)  # Open a new figure
                plt.subplot(2, 1, 1)
                
            for ind,current in enumerate(store_currents_np): plt.fill_between(time_vec,zeros_trace,current,color=store_colors[ind],alpha=figAlpha)
            plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
            plt.title('Fill plot - current traces')
            if time_range is not None:  plt.xlim(time_range)
            else:                       plt.xlim([0,cfg['duration']])
            if separateFigs: plt.savefig(savefig_name+'_0_current_traces__'+pop+savefig_format)
            else: plt.xticks([])

            # --- Fill plot of currents ratio at each timepoint - version 2
            ######
            # - Set of unique labels to group currents
            store_labels_set        = list(set(store_labels))
            # - Dictionary to store currents grouped by label
            store_currents_dict     = {current_label_:[] for current_label_ in store_labels_set}        
            # - Store currents grouped by label
            for current_label_ind, current_label in enumerate(store_labels): store_currents_dict[current_label].append(store_currents_np[current_label_ind])
            # - Calculate mean/std values for each current grouped by label
            store_mean_currents_list        = [np.mean(store_currents_dict[current_label],0) for current_label in store_currents_dict.keys()]
            store_std_currents_list         = [np.std( store_currents_dict[current_label],0) for current_label in store_currents_dict.keys()]
            # - Absolute value for each current (inward and outward currents are accounted the same in the current ratio plot)
            abs_store_mean_currents_list_np = np.abs(store_mean_currents_list)
            # - Ratio of currents - Normalizing current values in the [0-1] range
            ratio_mean_currents = abs_store_mean_currents_list_np / abs_store_mean_currents_list_np.sum(axis=0)
            
            # - Colors of each fill area
            fill_colors=[]
            for current_label_ in store_labels_set:
                try:    fill_colors.append(current_colors_dict[current_label_]) 
                except: fill_colors.append('cyan')
            
            # - Dictionary to store values for histogram plotting
            store_mean_currents_byPop_dict.update({pop:{'labels':store_labels_set,'mean':store_mean_currents_list,'std':store_std_currents_list,'colors':fill_colors}})
            ######
            
            # --- Creates an array of currents that adds the values of the previous plotted current to it, to make the fill_between plot
            store_ratio_mean_currents=[]
            for i, sublist in enumerate(ratio_mean_currents): store_ratio_mean_currents.append(sublist+np.sum(ratio_mean_currents[0:i], axis=0))
            if separateFigs: plt.figure(figsize=figSize)
            else:            plt.subplot(2, 1, 2)
            for ind,trace in enumerate(store_ratio_mean_currents):
                if ind ==0: plt.fill_between(time_vec,zeros_trace,                     store_ratio_mean_currents[ind],color=fill_colors[ind],alpha=figAlpha)
                else:       plt.fill_between(time_vec,store_ratio_mean_currents[ind-1],store_ratio_mean_currents[ind],color=fill_colors[ind],alpha=figAlpha)
            if separateFigs: plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
            plt.title('Fill plot - mean currents ratio')
            if time_range is not None:  plt.xlim(time_range)
            else:                       plt.xlim([0,cfg['duration']])
            plt.ylim([0,1])
            if separateFigs:    plt.savefig(savefig_name+'_1_instantaneous_current_ratio__'+pop+savefig_format)
            else:               plt.savefig(savefig_name+'_0_currents__'+pop+savefig_format)

            # --- Instantaneous current ratio plot
            plt.figure(figsize=figSize)
            for ind,trace in enumerate(ratio_mean_currents):
                if ind ==0: plt.fill_between(time_vec,zeros_trace,               ratio_mean_currents[ind],color=fill_colors[ind],alpha=0.1)
                else:       plt.fill_between(time_vec,ratio_mean_currents[ind-1],ratio_mean_currents[ind],color=fill_colors[ind],alpha=0.1)
                plt.plot(time_vec,ratio_mean_currents[ind],color=fill_colors[ind])
            plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
            plt.title('Fill plot - instantaneous current ratio')
            if time_range is not None:  plt.xlim(time_range)
            else:                       plt.xlim([0,cfg['duration']])
            plt.ylim([0,1])
            plt.savefig(savefig_name+'_2_instantaneous_current_normalized__'+pop+savefig_format)

        # --- Currents histogram
        plt.figure(figsize=figSize)
        plt.suptitle('Average of Chunks with Standard Deviation')
        npops=len(store_mean_currents_byPop_dict.keys())
        for pop_ind, pop in enumerate(store_mean_currents_byPop_dict.keys()):
            plt.subplot(npops,1,pop_ind+1)
            plt.xlabel('Chunks')
            plt.ylabel('Average Value')
            plt.title(pop)
            binSize=3 # ms
            binNumSamples=int(binSize/time_step)

            if time_range is not None:  time0_shift = time_range[0]
            else:                       time0_shift = 0

            # --- Calculates the average 
            for label_ind,label in enumerate(store_mean_currents_byPop_dict[pop]['labels']):
                mean_averages, mean_chunks = AnalyzeCurrents.chunked_averages(store_mean_currents_byPop_dict[pop]['mean'][label_ind], binNumSamples)
                std_averages,  std_chunks  = AnalyzeCurrents.chunked_averages(store_mean_currents_byPop_dict[pop]['std'][label_ind],  binNumSamples)

                mean_chunks_timebins    = [[(t0*time_step)+time0_shift,(t1*time_step)+time0_shift] for [t0,t1] in mean_chunks]
                mean_chunks_timebins_t0 = [ (t0*time_step)+time0_shift for [t0,t1] in mean_chunks]
                mean_chunks_timebins_t1 = [ (t1*time_step)+time0_shift for [t0,t1] in mean_chunks]
                mean_chunks_timebins_midpoint = [(t1+t0)/2 for [t0,t1] in mean_chunks_timebins]

                plt.step(       mean_chunks_timebins_midpoint,mean_averages,color=store_mean_currents_byPop_dict[pop]['colors'][label_ind],linewidth=3)
                plt.errorbar(   mean_chunks_timebins_t0, mean_averages, yerr=std_averages, 
                                fmt='', linestyle='', 
                                capsize=3,elinewidth=2,
                                color = store_mean_currents_byPop_dict[pop]['colors'][label_ind],
                                )
            if time_range is not None:  plt.xlim(time_range)
            else:                       plt.xlim([0,cfg['duration']])
        plt.savefig(savefig_name+'_3_hist'+savefig_format)

        # --- Currents histogram - Normalized
        plt.figure(figsize=figSize)
        plt.title('Normalized Average of Chunks with Standard Deviation')
        npops=len(store_mean_currents_byPop_dict.keys())
        for pop_ind, pop in enumerate(store_mean_currents_byPop_dict.keys()):
            plt.subplot(npops,1,pop_ind+1)
            plt.xlabel('Chunks')
            plt.ylabel('Average Value')
            plt.title(pop)
            
            binSize=3 # ms
            binNumSamples=int(binSize/time_step)
            
            if time_range is not None:  time0_shift = time_range[0]
            else:                       time0_shift = 0

            # --- Calculates the average 
            for label_ind,label in enumerate(store_mean_currents_byPop_dict[pop]['labels']):
                mean_averages, mean_chunks = AnalyzeCurrents.chunked_averages(store_mean_currents_byPop_dict[pop]['mean'][label_ind], binNumSamples)
                std_averages,  std_chunks  = AnalyzeCurrents.chunked_averages(store_mean_currents_byPop_dict[pop]['std'][label_ind],  binNumSamples)

                mean_chunks_timebins    = [[(t0*time_step)+time0_shift,(t1*time_step)+time0_shift] for [t0,t1] in mean_chunks]
                mean_chunks_timebins_t0 = [ (t0*time_step)+time0_shift for [t0,t1] in mean_chunks]
                mean_chunks_timebins_t1 = [ (t1*time_step)+time0_shift for [t0,t1] in mean_chunks]
                mean_chunks_timebins_midpoint = [(t1+t0)/2 for [t0,t1] in mean_chunks_timebins]

                if len(mean_averages)==0:continue # if no mean values

                if abs(min(mean_averages))>abs(max(mean_averages)):
                    mean_averages_norm_ = list((mean_averages-np.min(mean_averages))/(np.max(mean_averages)-np.min(mean_averages)))
                    mean_averages_norm = [val-1 for val in mean_averages_norm_]
                else:
                    mean_averages_norm = list((mean_averages-np.min(mean_averages))/(np.max(mean_averages)-np.min(mean_averages)))
                    
                # std_averages_norm = [abs((mean_averages_norm[std_ind]/mean_averages[std_ind])*std_averages[std_ind]) for std_ind in range(len(std_averages))]

                plt.step(       mean_chunks_timebins_midpoint,mean_averages_norm,color=store_mean_currents_byPop_dict[pop]['colors'][label_ind],linewidth=3)
                # plt.errorbar(   mean_chunks_timebins_t0, mean_averages_norm, yerr=std_averages_norm, 
                #                 fmt='', linestyle='', 
                #                 capsize=3,elinewidth=2,
                #                 color = store_mean_currents_byPop_dict[pop]['colors'][label_ind],
                #                 )
            if time_range is not None:  plt.xlim(time_range)
            else:                       plt.xlim([0,cfg['duration']])
        plt.savefig(savefig_name+'_4_histNorm'+savefig_format)

    ###############################################################################################################################################

    def get_current_colors_dict(figAlpha=0.5):
        current_colors_dict = { 
                                'SK_E2__ik'         : 'magenta',
                                'TC_HH__ina'        : 'brown',
                                'TC_HH__ik'         : 'teal',
                                'TC_Nap_Et2__ina'   : 'gold',
                                'TC_iA__ik'         : 'green',
                                'TC_iL__ica'        : 'r',
                                'TC_iT_Des98__ica'  : 'darkblue',
                                'TC_ih_Bud97__ih'   : 'orange',
                                
                                # test currents using example of netpyne/M1 pt cell 
                                'nax__ina'          : 'red',
                                'kap__ik'           : 'blue',
                                'kdr__ik'           : 'green',

                                'other'             : 'cyan'
                            }
        current_handles=[]
        for current_color in current_colors_dict.keys(): current_handles.append(plt.Line2D([], [], color=current_colors_dict[current_color], marker="o", linewidth=0, alpha=figAlpha, markeredgecolor='k'))
        return current_colors_dict

    def plotSingleTraces(simData, time_vec, current_colors_dict, current_flag, pop_gids, pop, savefig_name, savefig_format):
        print('Plotting single traces')
        defaultColor='cyan'

        figSize = (9,9)
        for recordedTrace in simData.keys():
            if recordedTrace.startswith(current_flag):
                for cell_name in simData[recordedTrace].keys():
                    gid = int(cell_name.split('cell_')[1])
                    # --- skips cells that are not in the population being plotted
                    if gid not in pop_gids:continue
                    c = defaultColor
                    for current_color in current_colors_dict.keys():
                        if current_color in recordedTrace: c=current_colors_dict[current_color]
                    plt.figure(figsize=figSize)
                    plt.plot(time_vec,simData[recordedTrace][cell_name],c=c)
                    # plt.plot(time_vec,simData[recordedTrace]['cell_'+str(gid)],c=c)
                    plt.title(pop+'__'+cell_name+'__'+recordedTrace)
                    plt.savefig(savefig_name+'_01_SingleTraces__'+pop+'__'+str(gid)+savefig_format)

    def plotAllTraces(simData, time_vec, current_colors_dict, current_handles, current_flag, pop_gids, pop, savefig_name, savefig_format):
        print('Plotting all traces')
        defaultColor='cyan'

        figSize = (9,9)
        fig2 = plt.subplots(figsize=figSize)  # Open a new figure
        # plt.suptitle('Currents of cell '+ str(sim['simConfig']['select_thal_gids'][gid]))
        num_traces = len([trace  for trace in simData.keys() if trace.startswith(current_flag)])
        trace_index_counter=0
        for recordedTrace_ind, recordedTrace in enumerate(simData.keys()):
            if recordedTrace.startswith(current_flag):
                trace_index_counter+=1
                plt.subplot(num_traces, 1, trace_index_counter)
                plt.title(recordedTrace)
                for cell_name in simData[recordedTrace].keys():
                    gid = int(cell_name.split('cell_')[1])
                    # --- skips cells that are not in the population being plotted
                    if gid not in pop_gids:continue
                    c = defaultColor
                    for current_color in current_colors_dict.keys():
                        if current_color in recordedTrace: c=current_colors_dict[current_color]
                    plt.plot(time_vec,list(simData[recordedTrace][cell_name]),c=c)
        plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
        plt.savefig(savefig_name+'_02_AllTraces__'+pop+savefig_format)

    def storePlotInfo(simData, current_colors_dict, current_flag, pop_gids):        
        defaultLabel = 'other'
        defaultColor = 'cyan'
        
        store_currents=[];store_labels=[];store_colors=[]
        for recordedTrace in simData.keys():
            if recordedTrace.startswith(current_flag):
                for cell_name in simData[recordedTrace].keys():
                    gid = int(cell_name.split('cell_')[1])
                    # --- skips cells that are not in the population being plotted
                    if gid not in pop_gids:continue
                    l = defaultLabel; c = defaultColor
                    for current_color in current_colors_dict.keys():
                        if current_color in recordedTrace: 
                            l = current_color
                            c = current_colors_dict[current_color]
                    store_currents.append(simData[recordedTrace][cell_name])
                    store_labels.append(l)
                    store_colors.append(c)
        return store_currents, store_labels, store_colors
    
    def sortCurrents(store_currents, store_labels, store_colors):
        sum_abs_current=[]
        for current_i,current in enumerate(store_currents):
            current_np  = np.array(current)
            abs_current = np.abs(current_np)
            sum_abs_current.append((sum(abs_current),current_i))
        sum_abs_current.sort(reverse=True)
        sorted_store_currents=[];sorted_store_labels=[];sorted_store_colors=[]
        for (sum_abs_curr,sorting_index) in sum_abs_current:
            sorted_store_currents.append(store_currents[sorting_index])
            sorted_store_labels.append(store_labels[sorting_index])
            sorted_store_colors.append(store_colors[sorting_index])

        del store_currents, store_labels, store_colors
            
        store_currents  = sorted_store_currents
        store_labels    = sorted_store_labels
        store_colors    = sorted_store_colors
        return store_currents, store_colors
    
    # --- Function to calculate the average value of a list with a defined chunck size
    def chunked_averages(input_list, chunk_size):
        num_chunks = len(input_list) // chunk_size
        averages = []
        indices = []

        for i in range(num_chunks):
            start = i * chunk_size
            end = start + chunk_size
            
            chunk_average = sum(input_list[start:end]) / chunk_size
            averages.append(chunk_average)
            indices.append([start, end])
        
        return averages, indices
    
    # def plot_chunked_averages(averages, std_devs, indices, savefig_name, savefig_format):
    #     chunk_labels = [f"{start}-{end-1}" for start, end in indices]
        
    #     plt.figure(figsize=(10, 6))
    #     plt.bar(chunk_labels, averages, yerr=std_devs, capsize=5, color='skyblue', ecolor='black')
    #     plt.xlabel('Chunks')
    #     plt.ylabel('Average Value')
    #     plt.title('Average of Chunks with Standard Deviation')
    #     plt.xticks(rotation=45)
    #     plt.savefig(savefig_name+'_AllTraces__'+pop+savefig_format)
    #     # plt.show()

    #########################################################################################################################

if __name__ == "__main__":

    '''
    Usage:
        -  to run AnalyzeCurrents.runAnalysis during runtime, add the code in the init.py script as follows:

    # --- Example of NetPyNE defaults
    from netpyne import sim
    cfg, netParams = sim.loadFromIndexFile('index.npjson')
    sim.createSimulateAnalyze(netParams, cfg)

    # --- Add the code as follows, passing the <sim> object/dictionary into the <AnalyzeCurrents.runAnalysis> method and configuring the other flags as desired
    from AnalyzeCurrents import AnalyzeCurrents
    AnalyzeCurrents.runAnalysis(sim, current_flag='i__', singleTraces=False, allTraces=False, separateFigs=True, figSize=(70,15), figAlpha=0.5)
    plt.show()

    example dataset path:
    /Users/joao/Desktop/bWgt_9227m_net_16sec/bWgt_9227m_net_16sec_currAn_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_data.json

    '''
    import sys
    if len(sys.argv)>1:
        try: 
            if sys.argv[1]=='-i':   arg_ind=2
            else:                   arg_ind=1
            filenames=[sys.argv[arg_ind]]
        except:
            try:
                # --- Loading a simulation output to plot the current analysis
                dataFolder='../data/'
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk()
                root.withdraw()
                filenames = filedialog.askopenfilenames(initialdir=dataFolder,typevariable='.json',defaultextension='.json')
                root.update()
            except:
                sys.exit()

    for filename in filenames:
        sim = AnalyzeCurrents.loadSimObj(filename)        
        AnalyzeCurrents.runAnalysis(sim, 
                                    current_flag='i__',
                                    select_pops=['VPM__pop','TRN__pop'],
                                    # select_currents=['i__soma_0__kap__ik', 'i__soma_0__kdr__ik'],
                                    select_currents=None,
                                    time_range=None,
                                    singleTraces=False,
                                    allTraces=False,
                                    separateFigs=True,
                                    figSize=(70,15),
                                    figAlpha=0.5,
                                    savefig_name = 'debug_current_analysis',
                                    savefig_format = '.png',
                                    )

    
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import json

class AnalyzeDeflection():

    def analyzeEvents(  
                        deflection_histogram, 
                        recordStep=0.025,
                        before_time_ms=100,
                        grouped_mean_voltages=None,
                        deflection_window = [-100, 450],
                        savefig_name='spikeHist_mean.png',
                        savefig_dpi=500,
                        minnery_units=True,
                        ):

        deflection_times = list(set([deflection_histogram[deflection_index]['VPM__pop']['ON']['bins'][0] for deflection_index in deflection_histogram.keys() if len(deflection_histogram[deflection_index].keys())>0]))
        deflection_times.sort()

        pop_colors={
            'MLe__pop':                 'k',
            'VPM__pop':                 'g',
            'TRN__pop':                 'b',
            'CTvirtual_uniform__pop':   'r',}

        normalize_spiking_activity_perCell=False
        num_cells={
            'MLe__pop':                 180,
            'VPM__pop':                 346,
            'TRN__pop':                 216,
            'CTvirtual_uniform__pop':   818,}

        # --- Store values used to plot the deflection histogram and calculcate the histogram difference
        store_deflection_mean={}

        # Create a bar plot
        plt.figure(figsize=(7, 10))
        npops = len(deflection_histogram['0'].keys())
        for pop_ind, pop in enumerate(deflection_histogram['0'].keys()):
            # def calculateDeflectionStatistics(deflection_histogram):
            # --- Analysis of deflection events
            #  -  PSTH values
            all_deflections_list        = [deflection_histogram[deflection_index][pop]['all']['data'] for deflection_index in deflection_histogram.keys() if pop in deflection_histogram[deflection_index].keys()]
            
            print('length of all deflections: '+ str(len(all_deflections_list)))

            #  -  all mean values
            all_deflections_list_mean      = list(np.mean(all_deflections_list, 0))
            all_deflections_list_std       = list(np.std( all_deflections_list, 0))

            bar_samples_mean = all_deflections_list_mean
            bar_samples_std  = all_deflections_list_std

            if pop in deflection_histogram["0"].keys(): time_vec_all_  = deflection_histogram["0"][pop]['all']['bins']

            if len(time_vec_all_)>len(all_deflections_list_mean):   del time_vec_all_[-1]
            # else:                                                   del time_vec_all_[0]

            # skip_samples=int(before_time_ms/recordStep)
            # index_zero = skip_samples+1
            # t_zero=time_vec_all_[index_zero]

            time_vec_all        = [t-min(time_vec_all_)-before_time_ms for t in time_vec_all_]
            time_vec_offset     = (time_vec_all[1]-time_vec_all[0])/2
            time_vec_all_offset = [t-time_vec_offset for t in time_vec_all]

            if minnery_units:
                binSize = abs(time_vec_all[1]-time_vec_all[0])

                print('Bin Size: ', binSize)
                bar_samples_mean = [(b*0.1)/binSize for b in bar_samples_mean]
                bar_samples_std  = [(b*0.1)/binSize for b in bar_samples_std]


            if normalize_spiking_activity_perCell:
                num_cells
                bar_samples_mean = [b/num_cells[pop] for b in bar_samples_mean]
                bar_samples_std  = [b/num_cells[pop] for b in bar_samples_std]


            ax1 = plt.subplot(npops+2, 1, pop_ind+1)
            # Set title
            # ax1.set_title(pop)

            # ax1.errorbar(   time_vec_all_offset, bar_samples_mean, yerr=bar_samples_std, 
            #                 fmt='', linestyle='', 
            #                 capsize=3,elinewidth=1,
            #                 color = 'k', alpha=0.5
            #                 )

            # ax1.fill_between(time_vec_all_offset, 
            #                  np.array(bar_samples_mean)-np.array(bar_samples_std), 
            #                  np.array(bar_samples_mean)+np.array(bar_samples_std),
            #                  )

            ax1.bar(        time_vec_all_offset, np.array(bar_samples_mean)+np.array(bar_samples_std), 
                            width=binSize,
                            color = pop_colors[pop], alpha=0.25
                            # color = 'k', alpha=0.25
                            )
            ax1.bar(        time_vec_all_offset, np.array(bar_samples_mean)-np.array(bar_samples_std), 
                            width=binSize,
                            color = 'w', alpha=1
                            )

            ax1.step(       time_vec_all,        bar_samples_mean, color=pop_colors[pop],linewidth=1)
            # ax1.step(       time_vec_all,        bar_samples_mean, color='k',linewidth=1)


            # ax1.legend(loc='upper left')
            
            # print(grouped_mean_voltages)

            if grouped_mean_voltages is not None:
                print('Plotting traces for ', pop)
                # print(mean_deflection_voltages.keys())
                if pop in grouped_mean_voltages.keys():
                    print('Plotting traces for ', pop)
                    ax2 = ax1.twinx()
                    # Overlay the time series plot on ax2
                    mean_deflection_voltages = AnalyzeDeflection.sliceMeanVoltages(grouped_mean_voltages, deflection_times=deflection_times, deflection_window = deflection_window, recordStep=recordStep)
                    ax2.plot(mean_deflection_voltages['t'], mean_deflection_voltages[pop]['mean'], 'k', alpha = 0.25)
                    # ax2.plot(mean_deflection_voltages['t'], mean_deflection_voltages[pop]['mean'], 'k', alpha = 0.25, label='Time Series')
                    ax2.set_ylabel('Mean Voltages')
                    ax2.legend(loc='upper right')
                    ax2.set_ylim([-60,-40])
                else: print('Population ', pop, ' not in mean_voltages')

            # Add labels and legends
            ax1.set_ylabel('0.1 * spikes/ms')
            if minnery_units:   
                # ax1.set_ylim([0,6.0])
                if max(bar_samples_mean)>2.5:   
                    ax1.set_ylim([0,6.0])
                    ax1.set_yticks([0,2.5,5])
                else:                           
                    ax1.set_ylim([0,3.0])
                    ax1.set_yticks([0,2.5])

            else:               ax1.set_ylim([0,110])

            # Formating plot frame and ticks
            # ax1.set_yticks([0,2.5,5])
            ax1.set_xticks([0,250])
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)

            store_deflection_mean.update({pop:{'t':time_vec_all,'vals': bar_samples_mean}})

        # plt.savefig(savefig_name,dpi=savefig_dpi)

        # Subtracting VPM from MLe histograms
        max_val_VPM=0
        max_ind_VPM=0
        for val_ind, val in enumerate(store_deflection_mean['VPM__pop']['vals']):
            if (store_deflection_mean['VPM__pop']['t'][val_ind]>-25) and (store_deflection_mean['VPM__pop']['t'][val_ind]<25):
                if store_deflection_mean['VPM__pop']['vals'][val_ind]>max_val_VPM:
                    max_val_VPM=val
                    max_ind_VPM=val_ind
                    print(max_val_VPM, max_ind_VPM)

        max_val_MLe=0
        max_ind_MLe=0
        for val_ind, val in enumerate(store_deflection_mean['MLe__pop']['vals']):
            if (store_deflection_mean['MLe__pop']['t'][val_ind]>-25) and (store_deflection_mean['MLe__pop']['t'][val_ind]<25):
                if store_deflection_mean['MLe__pop']['vals'][val_ind]>max_val_MLe:
                    max_val_MLe=val
                    max_ind_MLe=val_ind
                    print(max_val_MLe, max_ind_MLe)
        
        skip_indexes = max_ind_VPM-max_ind_MLe

        MLe_vals = store_deflection_mean['MLe__pop']['vals'][0:-skip_indexes]
        VPM_vals = store_deflection_mean['VPM__pop']['vals'][skip_indexes::]
        t_vals   = store_deflection_mean['MLe__pop']['t'][0:-skip_indexes]

        hist_difference = [VPM_vals[ind]-MLe_vals[ind] for ind in range(len(MLe_vals))]

        ax1 = plt.subplot(npops+2, 1, npops+1)
        for t_val_ind, t_val in enumerate(t_vals):
            if hist_difference[t_val_ind]>0:    bar_color='g'
            else:                               bar_color='k'
            ax1.bar(t_vals[t_val_ind], hist_difference[t_val_ind], width=2.5, color=bar_color, linewidth=2)
        # ax1.bar(t_vals, hist_difference, width=2.5, color='k', linewidth=2)

        # Formating plot frame and ticks
        ax1.set_yticks([-2.5,0,2.5])
        ax1.set_xticks([0,250])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        ax1.set_xlabel('Time (ms)')
        ax1.set_ylabel('0.1 * spikes/ms')

        ax1.set_ylim(-4.0,4.0)


        
        # Making a plot of TRN - VPM
        TRN_vals = store_deflection_mean['TRN__pop']['vals'][skip_indexes::]
        hist_difference2 = [TRN_vals[ind]-VPM_vals[ind] for ind in range(len(VPM_vals))]

        ax1 = plt.subplot(npops+2, 1, npops+2)
        for t_val_ind, t_val in enumerate(t_vals):
            if hist_difference2[t_val_ind]>0:   bar_color='b'
            else:                               bar_color='g'
            ax1.bar(t_vals[t_val_ind], hist_difference2[t_val_ind], width=2.5, color=bar_color, linewidth=2)
        # ax1.bar(t_vals, hist_difference, width=2.5, color='k', linewidth=2)

        # Formating plot frame and ticks
        ax1.set_yticks([-2.5,0,2.5])
        ax1.set_xticks([0,250])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        ax1.set_xlabel('Time (ms)')
        ax1.set_ylabel('0.1 * spikes/ms')

        ax1.set_ylim(-4.0,4.0)

        plt.savefig(savefig_name,dpi=savefig_dpi)
        
    ##  Normalized version of the figure


    
    
    
    ################################################################################################################################################################
            

    def poolDeflectionDicts(filesPath):
        import os
        files = [filesPath+file for file in os.listdir(filesPath) if 'spikeHist__data.json' in file]

        deflection_histogram={}
        deflection_index_counter=0
        for filePath in files:
            with open(filePath, 'r') as file:
                deflection_histogram_=json.load(file) 

            for ind_,i in enumerate(deflection_histogram_.keys()):
                if len(deflection_histogram_[i].keys())>0:
                    deflection_histogram.update({str(deflection_index_counter+ind_):deflection_histogram_[i]})
                else:print('Empty deflection dict on deflection # ', str(i))

            deflection_index_counter+=len(deflection_histogram_.keys())
        
        print('Loaded spike histogram data')
        return deflection_histogram

    def poolMeanVoltages_fromFiles(filesPath):
        import os
        files = [filesPath+file for file in os.listdir(filesPath) if 'meanVoltage__data.json' in file]

        group_mean_voltage_dicts={}
        update_t_vec=True

        for filePath_ind,filePath in enumerate(files):
            print('Processing file: ', filePath)
            with open(filePath, 'r') as file:
                mean_voltage_dict=json.load(file)
            
            if update_t_vec:
                mean_voltage_dict_keys = list(mean_voltage_dict.keys())
                t_vec = mean_voltage_dict[mean_voltage_dict_keys[0]]['t']
                update_t_vec=False
            
            # combines all membrane voltages into a list of lists, to calculate the mean of means, and compress the dataset
            for pop in mean_voltage_dict.keys():
                if pop not in group_mean_voltage_dicts.keys(): group_mean_voltage_dicts.update({pop:[]})
                
                group_mean_voltage_dicts[pop].append(mean_voltage_dict[pop]['mean_v'])
        
        # dictionary of lists with the mean value for each population (still needs to be sliced for each deflection event)
        grouped_mean_voltages = {pop:list(np.mean(group_mean_voltage_dicts[pop],0)) for pop in group_mean_voltage_dicts.keys()}
        grouped_mean_voltages.update({'t':t_vec})

        print('Loaded average membrane voltage data')
        return grouped_mean_voltages
    
    def sliceMeanVoltages(grouped_mean_voltages, deflection_times=[], deflection_window = [-100,350], recordStep=0.1):

        deflection_intervals = [[deflection_time+deflection_window[0],deflection_time+deflection_window[1]] for deflection_time in deflection_times]

        deflection_indexes = [[int(deflection_t0/recordStep), int(deflection_t1/recordStep)] for [deflection_t0, deflection_t1] in deflection_intervals]

        deflection_voltages={}
        for pop in grouped_mean_voltages.keys():
            if '__pop' in pop:
                if pop not in deflection_voltages.keys():deflection_voltages.update({pop:[]})
                for [t0,t1] in deflection_indexes:
                    deflection_voltages[pop].append(grouped_mean_voltages[pop][t0:t1])

        mean_deflection_voltages={pop:{'mean':np.mean(deflection_voltages[pop],0),'std':np.std(deflection_voltages[pop],0)} for pop in deflection_voltages.keys()}

        mean_deflection_voltages.update({'t':np.arange(deflection_window[0],deflection_window[1],recordStep)})

        return mean_deflection_voltages

        # # Gather indexes and values within the specified intervals
        # indexes={}
        # in_range_values = {}
        # for deflection_intervals_index, [start, end] in enumerate(deflection_intervals):
        #     mask = (grouped_mean_voltages['t'] >= start) & (grouped_mean_voltages['t'] <= end)
        #     # indexes.update(         {str(deflection_intervals_index):np.where(mask)[0]})
        #     indexes.update(         {str(deflection_intervals_index):[np.where(mask)[0][0],np.where(mask)[0][-1]]})
        #     in_range_values.update( {str(deflection_intervals_index):grouped_mean_voltages['t'][mask]})


        # for t in grouped_mean_voltages['t']:
        #     if 

if __name__ == "__main__":
    
    '''
    Example of command for pooled deflection analysis
        # ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9227n_net_16sec/'
        # ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9229a_net_16sec/'
        # ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9229d_net_16sec/'
        # ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9229f/'
        # ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9229g/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9230d/spike_hist_data/0/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9230d/spike_hist_data/1/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9230d/spike_hist_data/2/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9230d/spike_hist_data/3/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9235u_shift_iT_fixedDelay_1dot5ms/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9235v_fixedDelay_1dot5ms/spike_hist_data/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240r/spike_hist/uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240r/spike_hist/open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240r/spike_hist/closed/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240w/spike_hist/uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240w/spike_hist/open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/bWgt_9240w/spike_hist/closed/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241a/spike_hist/uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241a/spike_hist/open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241a/spike_hist/closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241a/spike_hist/mixed/'
        ipython -i AnalyzeDeflection.py '/Users/joao/Desktop/SfN_2024_data/bWgt_9241a/spike_hist/uniform/bWgt_9241a_connSchemas_testMixedFeedback_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_0_0_0_0_0_0_0_0_0_0_0_spikeHist__data.json'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241e/spike_hist/uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241e/spike_hist/open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241e/spike_hist/closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Desktop/SfN_2024_data/bWgt_9241e/spike_hist/mixed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/1_closed_open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/2_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/3_closed_mixed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/4_mixed_open/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/5_mixed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9247/bWgt_9247d_newFF/spike_hist/6_mixed_mixed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/1_uniform_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/2_uniform_remixed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/3_closed_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/4_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/5_closed_remixed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/6_remixed_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/7_remixed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9248/bWgt_9248d_newFF/spike_hist/8_remixed_remixed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9255/bWgt_9255d_newFF_newCtFb/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9255/bWgt_9255d_newFF_newCtFb/spike_hist/4_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9255/bWgt_9255d_newFF_newCtFb/spike_hist/8_remixed_remixed/'

        # bWgt_9256
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replaceNone/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replaceNone/spike_hist/4_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replaceNone/spike_hist/8_remixed_remixed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace180/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace180/spike_hist/4_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace180/spike_hist/8_remixed_remixed/'

        # bWgt_9256
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace360/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace360/spike_hist/4_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9256/bWgt_9256d_newFF_newCtFb_remixed_replace360/spike_hist/8_remixed_remixed/'

        # --- previous good sim
        # bWgt_9257
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/1_closed_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/2_remixed_remixed_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/3_remixed_remixed_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/4_remixed_remixed_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/5_remixed_remixed_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9257/bWgt_9257d_ct360/spike_hist/_remixed_remixed_90_10/'
    
        # bWgt_9266
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/0_uniform_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/1_remixed_remixed_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/2_remixed_remixed_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/3_remixed_remixed_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/4_remixed_remixed_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/5_remixed_remixed_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/6_remixed_remixed_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/7_remixed_remixed_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/8_remixed_remixed_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/9_remixed_remixed_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9266/bWgt_9266d_ct360_tuneMLe/spike_hist/10_closed_closed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d1_ct360_tuneMLe2/spike_hist/10_closed/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9270/bWgt_9270d2_ct360_tuneMLe2/spike_hist/10_closed/'
        
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9272/bWgt_9272d1_ct360_tuneMLe2/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9272/bWgt_9272d1_ct360_tuneMLe2/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9272/bWgt_9272d1_ct360_tuneMLe2/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9272/bWgt_9272d1_ct360_tuneMLe2/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9272/bWgt_9272d1_ct360_tuneMLe2/spike_hist/6_60_40/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_debug/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_debug3/spike_hist/'


        
        # With CT AT                
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2/spike_hist/6_60_40/'

        # No CT AT                
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat/spike_hist/6_60_40/'

        ##########################################################################################
        # no CT
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_noCTat_awake/spike_hist/10_closed/'

        # CT 10 / 25 / 50 / 100 %
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_CTat_10pct/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_CTat_25pct/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_CTat_50pct/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9273/bWgt_9273d1_ct360_tuneMLe2_CTat_100pct/spike_hist/4_40_60/'
        
        ##########################################################################################

        # bWgt_9275
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9275/bWgt_9275d1_MLe6_noCTat_awake/spike_hist/10_closed/'
        ##########################################################################################
        # bWgt_9278
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9278/bWgt_9278_GABAgExtra_015_RMP_4/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9278/bWgt_9278_GABAgExtra_015_RMP_2/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9278/bWgt_9278_GABAgExtra_015_RMP_0/spike_hist/'
        ##########################################################################################
        # bWgt_9279
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9279/bWgt_9279_RMP_2_MLe_2_0/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9279/bWgt_9279_RMP_2_MLe_4_0/spike_hist/'
        ##########################################################################################
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9280/bWgt_9280_GABAe_015_MLe_2_0/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9280/bWgt_9280_GABAe_015_MLe_3_0/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9280/bWgt_9280_GABAe_015_MLe_4_0/spike_hist/'
        ##########################################################################################
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9281/bWgt_9281a_GABAe_025_MLe_2_0/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9281/bWgt_9281b_GABAe_025_MLe_1_0/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9281/bWgt_9281c_GABAe_015_MLe_1_0/spike_hist/'
        ##########################################################################################
        # --- Closed-loop
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/10_closed/'
        
        # Plotting PSTH for a single deflection angle for the paper - 180 deg
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake/spike_hist/6_60_40_180deg/'

        # --- Open-loop
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/0_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/1_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/2_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/3_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/4_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/5_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/6_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/7_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/8_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/9_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/10_closed/'
        
        # Plotting PSTH for a single deflection angle for the paper - 180 deg
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9300/bWgt_9300e_awake_open/spike_hist/6_60_40_180deg/'
        ##########################################################################################

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9306/bWgt_9306e2_awake_closedUni/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9306/bWgt_9306e1_awake_closedUni/spike_hist/'
        ##########################################################################################
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/0_closed/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/1_closed_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/2_closed_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/3_closed_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/4_closed_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/5_closed_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/6_closed_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/7_closed_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/8_closed_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/9_closed_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/10_uniform/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/11_open_10_90/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/12_open_20_80/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/13_open_30_70/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/14_open_40_60/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/15_open_50_50/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/16_open_60_40/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/17_open_70_30/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/18_open_80_20/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/19_open_90_10/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9308/bWgt_9308e_awake/spike_hist/20_open/'
        ##########################################################################################

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_1_mixed5050_CTmodulation_0/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_1_mixed5050_CTmodulation_10/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_1_mixed5050_CTmodulation_25/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_1_mixed5050_CTmodulation_50/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_1_mixed5050_CTmodulation_100/sim_output/spike_hist/'

        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_2_closed_CTmodulation_0/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_2_closed_CTmodulation_10/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_2_closed_CTmodulation_25/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_2_closed_CTmodulation_50/sim_output/spike_hist/'
        ipython -i AnalyzeDeflection.py 'all_deflections' '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/c_0003/c_0003_300_awake_2_closed_CTmodulation_100/sim_output/spike_hist/'

        ##########################################################################################

    '''

    import sys
    from    matplotlib  import  pyplot  as plt
    print("Matplotlib backend (default): %s" %plt.get_backend())
    modules = []
    for module in sys.modules:
        if module.startswith('matplotlib'):
            modules.append(module)
    for module in modules:
        sys.modules.pop(module)
    import matplotlib
    matplotlib.use("MacOSX")
    from    matplotlib  import  pyplot  as plt
    print("Matplotlib backend (dynamic): %s" %plt.get_backend())

    import sys
    if len(sys.argv)>1:
        if 'all_deflections' in sys.argv: 
            print('Running pooled analysis')
            filesPath = sys.argv[-1]
            deflection_histogram  = AnalyzeDeflection.poolDeflectionDicts(filesPath=filesPath)

            try:    grouped_mean_voltages = AnalyzeDeflection.poolMeanVoltages_fromFiles(filesPath=filesPath)
            except:
                print('Failed to pool voltages')
                # try:    grouped_mean_voltages = {   
                #                                     # 'MLe__pop':    list(np.random.uniform(low=-75*1e-3, high=-57*1e-3, size=(int(16000/0.1),))),
                #                                     'VPM__pop':    list(np.random.uniform(low=-75*1e-3, high=-57*1e-3, size=(int(16000/0.1),))),
                #                                     'TRN__pop':    list(np.random.uniform(low=-75*1e-3, high=-57*1e-3, size=(int(16000/0.1),))),
                #                                     't':           np.arange(0,16000,0.1),
                #                                  }
                # except: 
                #     print('Failed to assign voltages')
                #     grouped_mean_voltages = None
                grouped_mean_voltages=None

        else: 
            if sys.argv[1]=='-i':   arg_ind=2
            else:                   arg_ind=1

            print(sys.argv)

            filePath=sys.argv[arg_ind]
            # try: 
            #     print(filePath)
            # except:
            #     try:

            #         from    matplotlib  import  pyplot  as plt
            #         import sys
            #         print("Matplotlib backend (default): %s" %plt.get_backend())
            #         modules = []
            #         for module in sys.modules:
            #             if module.startswith('matplotlib'):
            #                 modules.append(module)
            #         for module in modules:
            #             sys.modules.pop(module)
            #         import matplotlib
            #         matplotlib.use("MacOSX")
            #         from    matplotlib  import  pyplot  as plt
            #         print("Matplotlib backend (dynamic): %s" %plt.get_backend())

            #         filePath='/Users/joao/Desktop/bWgt_9227i_net_16sec_currAn_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_0_0_0_0_0_0_0_0_0_0_0_0_0_spikeHist__data.json'

            #     except:
            #         sys.exit()
            with open(filePath, 'r') as file:
                deflection_histogram=json.load(file)   
    else:sys.exit()

    print('Analyzing '+str(len(deflection_histogram.keys()))+' deflection events')
    recordStep=0.1
    try:
        try:
            figName = '__'.join(filesPath.split('/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/')[1].split('/'))
            AnalyzeDeflection.analyzeEvents(
                                                deflection_histogram,
                                                # dt=0.1, 
                                                recordStep=recordStep, 
                                                # dt=sim.cfg.recordStep, 
                                                # dt=sim.cfg.dt, 
                                                grouped_mean_voltages = grouped_mean_voltages,
                                                savefig_name=figName+'.png',
                                                savefig_dpi=500,
                )
        except:
            figName = '__'.join(filesPath.split('/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/')[1].split('/'))
            AnalyzeDeflection.analyzeEvents(
                                                deflection_histogram,
                                                # dt=0.1, 
                                                recordStep=recordStep, 
                                                # dt=sim.cfg.recordStep, 
                                                # dt=sim.cfg.dt, 
                                                grouped_mean_voltages = grouped_mean_voltages,
                                                savefig_name='debug_spikeHist_mean2.png',
                                                savefig_dpi=500,
                )
    except:
        AnalyzeDeflection.analyzeEvents(
                                            deflection_histogram,
                                            # dt=0.1, 
                                            recordStep=recordStep, 
                                            # dt=sim.cfg.recordStep, 
                                            # dt=sim.cfg.dt, 
                                            grouped_mean_voltages = None,
                                            savefig_name='debug_spikeHist_mean.png',
                                            savefig_dpi=500,
            )
    plt.show()
    import sys
    sys.exit()

    baseline_deflections_list   = []
    ON_deflections_list         = []
    OFF_deflections_list        = []

    # Create a bar plot
    plotData=False
    if plotData:
        plt.figure(figsize=(15, 5))
        plt.xlabel('Sample Index')
        plt.ylabel('Height')
        plt.title('Bar Plot of Sampled Heights with Double Exponential Decay')

    deflection_data_dict={'pop1':{}}
    n_dists=10
    for i in range(n_dists):
        # Parameters
        mean_height = 0.75 + 0.5*np.random.uniform()
        std_dev = 0.05
        num_uniform_samples = 50
        num_modified_samples = 150
        initial_offset  = 5 * mean_height
        initial_offset2 = 3 * mean_height
        # Generate the first 2000 samples from the uniform distribution
        uniform_samples = np.random.uniform(mean_height - std_dev, mean_height + std_dev, size=num_uniform_samples)
        # Generate the decay using a double exponential function
        t = np.arange(num_modified_samples)
        decay = initial_offset * np.exp(-t / (num_modified_samples / 5)) * np.exp(-t / (num_modified_samples / 7.5))
        # Generate the next 8000 samples from the uniform distribution with the double exponential decay
        modified_samples  = np.random.uniform(mean_height - std_dev, mean_height + std_dev, size=num_modified_samples) + decay
        # Generate the decay using a double exponential function
        decay2 = initial_offset2 * np.exp(-t / (num_modified_samples / 5)) * np.exp(-t / (num_modified_samples / 10))
        # Generate the next 8000 samples from the uniform distribution with the double exponential decay
        modified_samples2  = np.random.uniform(mean_height - std_dev, mean_height + std_dev, size=num_modified_samples) + decay2
        # Combine the samples
        samples = np.concatenate((uniform_samples, modified_samples, modified_samples2))
        # plt.subplot(n_dists,1,i+1)
        if plotData:plt.bar(range(len(samples)), samples, width=1.0, alpha=0.1)

        baseline_deflections_list.append(   list(uniform_samples))
        ON_deflections_list.append(         list(modified_samples))
        OFF_deflections_list.append(        list(modified_samples2))

        deflection_data_dict['pop1'].update({   str(i):{
                                                    'baseline': { 'bins':[], 'data':uniform_samples },
                                                    'ON':       { 'bins':[], 'data':modified_samples, },
                                                    'OFF':      { 'bins':[], 'data':modified_samples2,}
                                                    },
                                                })
    if plotData:plt.show()

    AnalyzeDeflection.analyzeEvents(deflection_histogram=deflection_data_dict, recordStep=0.025, savefig_name='spikeHist_mean_debug_fig.png')
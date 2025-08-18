################################################################################################################################################

import os
import pickle
import netpyne
# from netpyne import sim
import json

import scipy.stats as stats

import sys
from matplotlib import pyplot as plt
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

import matplotlib.patches as patches
import numpy as np

########################################################################################################################################################################################################
# --- Setting the path for loading the datasets
data_path = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data'
# batch_sim_group = 'c_0003'
# batch_sim_label = 'c_0003_000_sleep_wake'
# sim_output_folder = 'sim_output'

# batch_sim_group = 'paper_000'
# batch_sim_group = 'paper_001'
batch_sim_group = 'paper_002'
batch_sim_label = batch_sim_group+'_'+'barreloid_batch_fig02_sleep_wake'
sim_output_folder = 'sim_output'

simulation_data_path = data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/'+sim_output_folder+'/'

########################################################################################################################################################################################################
# --- Path to save analysis figures
save_figures_path = data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/analysis_figures/'
if not os.path.exists(save_figures_path): os.makedirs(save_figures_path)

########################################################################################################################################################################################################
def significance_symbol(p_value):
    """
    Returns a significance symbol based on the p-value.
    Parameters: p_value (float): The p-value from a statistical test.
    Returns:    str: A string representing significance level ('ns', '*', '**', '***', '****').
    """
    if p_value >= 0.05:
        return 'ns'
    elif p_value < 0.0001:
        return '****'
    elif p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    else:
        return '*'
################################################################################################################################################

popColors_dict={}
for pop in ['VPM__pop', 'TRN__pop', 'TRN_ring__pop', 
            'MLe__pop', 
            'L6A_activated__pop',
            'CTvirtual_uniform__pop',
            'CTvirtual_activated__pop',
            ]:
    if   'VPM' in           pop:    popColors_dict.update({pop:'g'})
    elif 'TRN__pop' in      pop:    popColors_dict.update({pop:'b'})
    elif 'ring' in          pop:    popColors_dict.update({pop:'mediumvioletred'})
    elif 'CTvirtual_' in    pop:    popColors_dict.update({pop:'r'})
    elif 'L6A' in           pop:    popColors_dict.update({pop:'r'})
    elif 'L6A_activated' in pop:    popColors_dict.update({pop:'r'})
    elif 'MLe' in           pop:    popColors_dict.update({pop:'k'})
    else:                           popColors_dict.update({pop:'gray'})
    
pop_color_maps={
    'VPM__pop':'summer',
    'TRN__pop':'winter',
                }

################################################################################################################################################
# --- Figure properties
plt.rcParams.update({'font.size': 15})

sim_index_mapping   = [-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0]
spindle_freq_range  = [8,16] # Hz
y_lims              = {'VPM__pop':[-1,51],'TRN__pop':[-1,151]}
y_text_offset       = {'VPM__pop':1,'TRN__pop':4}
limit_frequency     = 50 # Hz

bar_y_lims          = {'VPM__pop':[-1,160],'TRN__pop':[-1,410]}
bar_y_ticks         = {'VPM__pop':[0,75,150],'TRN__pop':[0,200,400]}

plot_pops=['TRN__pop','VPM__pop']

################################################################################################################################################
#####   Figure 02: Sleep vs Wake PSD data   #####
# 0: Tuned Sleep network                            (selected)
# 1: Tuned Awake network - 10 Hz                    (selected)
## Removed from simulations
# 2: Tuned Awake network - 20 Hz - No deflection
# 3: Tuned Awake network - 20 Hz - With 180 deg deflection
# 4: Default BBP params

sleep_vs_wake_select_sims=[0,1]

# --- Data path
sleep_vs_wake_spike_psd_data_path   = simulation_data_path+'spikePSD_custom/'
sleep_vs_wake_data_folder = sleep_vs_wake_spike_psd_data_path+'data/'

sleep_vs_wake_list_files = os.listdir(sleep_vs_wake_data_folder)
sleep_vs_wake_list_files_PSD = [file for file in sleep_vs_wake_list_files if '_ratePSD_values.json' in file]
sleep_vs_wake_list_files_PSD.sort()

store_sleep_vs_wake_PSD_dict={}
for sleep_vs_wake_data_file_PSD_ind, sleep_vs_wake_data_file_PSD in enumerate(sleep_vs_wake_list_files_PSD):
    filename_json = sleep_vs_wake_data_folder+sleep_vs_wake_data_file_PSD
    with open(filename_json, 'r') as file: sleep_vs_wake_PSD_dict = json.load(file)
    store_sleep_vs_wake_PSD_dict.update({'sim_'+str(sleep_vs_wake_data_file_PSD_ind):sleep_vs_wake_PSD_dict})

# --- Figure 02: PSD traces
plt.figure(figsize=[5,5])
# plt.figure(figsize=[10,10])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)

    square = patches.Rectangle((spindle_freq_range[0], y_lims[plot_pop][0]), spindle_freq_range[1]-spindle_freq_range[0], y_lims[plot_pop][1], edgecolor='k', facecolor='w', alpha=0.1)
    ax=plt.gca()
    ax.add_patch(square)

    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_sleep_vs_wake_PSD_dict.keys()):
        if sleep_vs_wake_select_sims is not None:
            if sim_index not in sleep_vs_wake_select_sims: continue

        if (sim_index == 0):
            line_alpha = 0.75
            line_width = 4
            # line_width = 3
            line_color = 'mediumturquoise'
            # line_color = cmap(0.2)
        elif (sim_index == 1):
            line_alpha = 0.75
            line_width = 2
            # line_width = 3
            line_color = 'mediumpurple'
            # line_color = cmap(0.8)
        else:
            line_alpha = 0.75
            line_width = 1
            line_color = 'grey'
        
        plt.plot(store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['freqs'],store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['signal'],
                                                                          color=line_color,
                                                                        #   color=cmap(sim_index/len(store_sleep_vs_wake_PSD_dict.keys())),
                                                                          linewidth=line_width, 
                                                                        #   linewidth=2+(0.15*sim_index), 
                                                                          alpha=line_alpha)
    plt.ylim(y_lims[plot_pop])
    plt.ylabel('Power')
    
    if plot_pop_ind==len(plot_pops)-1:
        plt.xticks([0,12,24,36,48])
        plt.xlabel('Frequency (Hz)')
    else:
        plt.xticks([])
        plt.xlabel('')
    
    plt.xlim([0,limit_frequency])
    
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()

# plt.savefig(sleep_vs_wake_spike_psd_data_path+'spikePSD_signal.png',dpi=200)
plt.savefig(save_figures_path+'spikePSD_signal.png',dpi=200)

sleep_vs_wake_sim_index_mapping   = [0,1,2,3,4]
sleep_vs_wake_y_text_offset       = {'VPM__pop':1,'TRN__pop':6}
# --- PSD bar plot of max power value and frequency of max power
plt.figure(figsize=[3,5])
# plt.figure(figsize=[3,10])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_sleep_vs_wake_PSD_dict.keys()):
        if sleep_vs_wake_select_sims is not None:
            if sim_index not in sleep_vs_wake_select_sims: continue 

        slice_data  = [(freq,store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
        max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
        
        full_data = [(freq,store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_sleep_vs_wake_PSD_dict[sim_name][plot_pop]['freqs']) if freq<=limit_frequency]
        max_power_signal        = max(full_data, key=lambda x: x[1])

        # Extract frequency and power values
        slice_freqs, slice_power = zip(*slice_data)
        # Compute area under the curve using the trapezoidal rule
        power_in_band = np.trapz(slice_power, slice_freqs)

        bar_value = power_in_band
        bar_string = ''
        bar_string_color = 'k'
        font_weight='bold'

        if (sim_index == 0):
            bar_color = 'mediumturquoise'
            # bar_color = cmap(0.2)
            bar_alpha = 1.0
            edge_color = 'k'
        elif (sim_index == 1):
            bar_color = 'mediumpurple'
            # bar_color = cmap(0.8)
            bar_alpha = 1.0
            edge_color = 'k'
        else:
            bar_color = 'grey'
            bar_alpha = 0.5
            edge_color = None

        plt.bar(sleep_vs_wake_sim_index_mapping[sim_index],bar_value,width=0.84,color=bar_color, edgecolor = edge_color, alpha=bar_alpha)
        # plt.text(sleep_vs_wake_sim_index_mapping[sim_index]-0.25, bar_value+sleep_vs_wake_y_text_offset[plot_pop], bar_string, color=bar_string_color, rotation=0, size=15, weight=font_weight, alpha = bar_alpha)

        if plot_pop_ind==len(plot_pops)-1:
            plt.xticks([0,1],['Sleep','Wake'],rotation=45)
            plt.xlabel('')
        else:
            plt.xticks([0,1],['',''])
            plt.xticks([])
            plt.xlabel('')
        plt.xlim([-1,4])

    plt.ylim(bar_y_lims[plot_pop])
    plt.yticks(bar_y_ticks[plot_pop])
    plt.ylabel('Power 8-16 Hz')
    
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()

# plt.savefig(sleep_vs_wake_spike_psd_data_path+'spikePSD_peakFrequency.png',dpi=200)
plt.savefig(save_figures_path+'spikePSD_peakFrequency.png',dpi=200)

################################################################################################################################################
#####   LFP PSD data   #####
# --- Data path
sleep_vs_wake_lfp_psd_data_path   = simulation_data_path+'lfpPSD/'
sleep_vs_wake_data_folder_lfpPSD = sleep_vs_wake_lfp_psd_data_path+'data/'

sleep_vs_wake_list_files_lfpData = os.listdir(sleep_vs_wake_data_folder_lfpPSD)
sleep_vs_wake_list_files_lfpPSD = [file for file in sleep_vs_wake_list_files_lfpData if '_lfpPSD_values.json' in file]
sleep_vs_wake_list_files_lfpPSD.sort()

store_sleep_vs_wake_lfpPSD_dict={}
for sleep_vs_wake_data_file_lfpPSD_ind, sleep_vs_wake_data_file_lfpPSD in enumerate(sleep_vs_wake_list_files_lfpPSD):
    filename_json = sleep_vs_wake_data_folder_lfpPSD+sleep_vs_wake_data_file_lfpPSD
    with open(filename_json, 'r') as file: sleep_vs_wake_lfpPSD_dict = json.load(file)
    store_sleep_vs_wake_lfpPSD_dict.update({'sim_'+str(sleep_vs_wake_data_file_lfpPSD_ind):sleep_vs_wake_lfpPSD_dict})

################################################################################################################################################

lfpPSD_y_lims              = {'VPM__pop':[-0.1e-5,3.6e-5],'TRN__pop':[-0.1e-5,3.6e-5]}
lfpPSD_y_lims_bar          = {'VPM__pop':[-0.1e-5,12.0e-5],'TRN__pop':[-0.1e-5,12.0e-5]}

# --- LFP traces
plt.figure(figsize=[5,5])
# plt.figure(figsize=[10,10])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)

    square = patches.Rectangle((spindle_freq_range[0], lfpPSD_y_lims[plot_pop][0]), spindle_freq_range[1]-spindle_freq_range[0], lfpPSD_y_lims[plot_pop][1], edgecolor='k', facecolor='w', alpha=0.1)
    ax=plt.gca()
    ax.add_patch(square)

    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_sleep_vs_wake_lfpPSD_dict.keys()):
        if sleep_vs_wake_select_sims is not None:
            if sim_index not in sleep_vs_wake_select_sims: continue 
        
        if (sim_index == 0):
            line_alpha = 0.75
            line_width = 4
            line_color = 'mediumturquoise'
            # line_color = cmap(0.2)
        elif (sim_index == 1):
            line_alpha = 0.75
            line_width = 2
            line_color = 'mediumpurple'
            # line_color = cmap(0.8)
        else:
            line_alpha = 0.75
            line_width = 1
            line_color = 'lightgrey'

        plt.plot(store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['freqs'],store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['signal'],
                                                                          color=line_color,
                                                                        #   color=cmap(sim_index/len(store_sleep_vs_wake_lfpPSD_dict.keys())),
                                                                          linewidth=line_width, 
                                                                        #   linewidth=2+(0.15*sim_index), 
                                                                          alpha=line_alpha)
    plt.ylim(lfpPSD_y_lims[plot_pop])
    plt.ylabel('Power')
    
    if plot_pop_ind==len(plot_pops)-1:
        plt.xticks([0,12,24,36,48])
        plt.xlabel('Frequency (Hz)')
    else:
        plt.xticks([])
        plt.xlabel('')
    
    plt.yticks([0,1.5e-5,3.0e-5])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
    plt.xlim([0,limit_frequency])
    
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()

# plt.savefig(sleep_vs_wake_lfp_psd_data_path+'lfpPSD_signal.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_signal.png',dpi=200)

# --- PSD bar plot of max power value and frequency of max power
plt.figure(figsize=[3,5])
# plt.figure(figsize=[3,10])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_sleep_vs_wake_lfpPSD_dict.keys()): 
        if sleep_vs_wake_select_sims is not None:
            if sim_index not in sleep_vs_wake_select_sims: continue 
        
        slice_data  = [(freq,store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
        max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
        
        full_data = [(freq,store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_sleep_vs_wake_lfpPSD_dict[sim_name][plot_pop]['freqs']) if freq<=limit_frequency]
        max_power_signal        = max(full_data, key=lambda x: x[1])


        # Extract frequency and power values
        slice_freqs, slice_power = zip(*slice_data)
        # Compute area under the curve using the trapezoidal rule
        power_in_band = np.trapz(slice_power, slice_freqs)

        bar_value = power_in_band
        bar_string = ''
        bar_string_color = 'k'
        font_weight='bold'

        if (sim_index == 0):
            bar_color = 'mediumturquoise'
            # bar_color = cmap(0.2)
            bar_alpha = 1.0
            edge_color = 'k'
        elif (sim_index == 1):
            bar_color = 'mediumpurple'
            # bar_color = cmap(0.8)
            bar_alpha = 1.0
            edge_color = 'k'
        else:
            bar_color = 'lightgrey'
            bar_alpha = 0.5
            edge_color = None


        # if max_power_spindleRange[0]==max_power_signal[0]:  
        #     if (sim_index == 2):    bar_color = cmap(0.8)
        #     else:                   bar_color = cmap(0.2)
        #     # bar_color = cmap(sim_index/len(store_sleep_vs_wake_lfpPSD_dict.keys()))
        #     bar_value = max_power_spindleRange[1]
        #     bar_string = max_power_spindleRange[0]
        #     bar_string_color = 'k'
        #     edge_color = bar_color
        # else:                                               
        #     bar_color = 'w'
        #     bar_value = max_power_signal[1]
        #     bar_string = ''
        #     bar_string_color = 'r'
        #     edge_color = 'k'

        # if (sim_index == 1) or (sim_index == 2): 
        #     edge_color='k'
        #     font_weight='bold'
        #     bar_alpha = 1.0
        # else: 
        #     font_weight='normal'
        #     bar_alpha = 0.5

        plt.bar(sleep_vs_wake_sim_index_mapping[sim_index],bar_value,width=0.84,color=bar_color, edgecolor = edge_color, alpha=bar_alpha)
        # plt.text(sleep_vs_wake_sim_index_mapping[sim_index]-0.25, bar_value+0.175e-5, bar_string, color=bar_string_color, rotation=0, size=15, weight=font_weight, alpha = bar_alpha)

        if plot_pop_ind==len(plot_pops)-1:
            plt.xticks([0,1],['Sleep','Wake'],rotation=45)
            plt.xlabel('')
        else:
            plt.xticks([0,1],['',''])
            plt.xticks([])
            plt.xlabel('')
        plt.xlim([-1,4])

    plt.ylim(lfpPSD_y_lims_bar[plot_pop])
    plt.ylabel('Power 8-16 Hz')
    plt.yticks([0,5.0e-5,10.0e-5])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()

# plt.savefig(sleep_vs_wake_lfp_psd_data_path+'lfpPSD_peakFrequency.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_peakFrequency.png',dpi=200)

################################################################################################################################################
#####   LFP PSD windowed data   #####
# --- Data path
lfp_psd_windowed_data_path   = simulation_data_path+'lfpPSD_windows/'
data_folder_lfpPSD_windowed = lfp_psd_windowed_data_path+'data/'

list_files_lfpData = os.listdir(data_folder_lfpPSD_windowed)
list_files_lfpPSD_windowed = [file for file in list_files_lfpData if '_lfpPSD_values_windows.json' in file]
list_files_lfpPSD_windowed.sort()

store_lfpPSD_windowed_dict={}
for data_file_lfpPSD_windowed_ind, data_file_lfpPSD_windowed in enumerate(list_files_lfpPSD_windowed):
    filename_json = data_folder_lfpPSD_windowed+data_file_lfpPSD_windowed
    with open(filename_json, 'r') as file: lfpPSD_dict = json.load(file)
    store_lfpPSD_windowed_dict.update({'sim_'+str(data_file_lfpPSD_windowed_ind):lfpPSD_dict})

################################################################################################################################################

lfpPSD_window_y_lims              = {'VPM__pop':[-0.1e-5,1.5e-4],'TRN__pop':[-0.1e-5,1.5e-4]}
lfpPSD_window_y_lims_bar          = {'VPM__pop':[-0.1e-5,1.5e-4],'TRN__pop':[-0.1e-5,1.5e-4]}


store_window_power_overTime_dict={} # power in a 8-16 Hz band (area under the curve)
store_window_power_allFreqs_dict={} # power in a 8-16 Hz band (raw values)
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    auc_store_window_powers = []
    auc_store_window_powers_Thresh=[]

    store_window_power_overTime_dict.update({plot_pop:{}})
    store_window_power_allFreqs_dict.update({plot_pop:{}})

    for sim_index, sim_name in enumerate(store_lfpPSD_windowed_dict.keys()): 
        
        # if (sim_index!=2) and (sim_index!=4) and (sim_index!=6):continue
        
        store_window_power_overTime = []
        store_window_time_overTime=[]
        
        store_window_power_overTime_allFreqs = []
        store_window_freqs_overTime_allFreqs = []
        store_window_time_overTime_allFreqs=[]

        for time_window in store_lfpPSD_windowed_dict[sim_name][plot_pop].keys():
            
            # Separating power data in the 8-16 Hz frequency band
            slice_data  = [(freq,store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
            max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
            
            # Extract frequency and power values
            slice_freqs, slice_power = zip(*slice_data)
            # Compute area under the curve using the trapezoidal rule
            power_in_band = np.trapz(slice_power, slice_freqs)

            store_window_power_overTime.append(power_in_band)
            store_window_time_overTime.append((int(time_window.split('__')[1])+int(time_window.split('__')[0]))/2)


            # Separating power data in the 0-50 Hz frequency band (full spectrum)
            full_data = [(freq,store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['freqs']) if freq<=limit_frequency]
            max_power_signal        = max(full_data, key=lambda x: x[1])

            # Extract frequency and power values
            slice_allfreqs, slice_allpower = zip(*full_data)

            store_window_power_overTime_allFreqs.append(slice_allpower)
            # store_window_freqs_overTime_allFreqs.append(slice_allfreqs)
            store_window_time_overTime_allFreqs.append((int(time_window.split('__')[1])+int(time_window.split('__')[0]))/2)

        store_window_power_overTime_dict[plot_pop].update({sim_name:{'t':store_window_time_overTime, 'power': store_window_power_overTime}})
        
        store_window_power_allFreqs_dict[plot_pop].update({sim_name:{'t':store_window_time_overTime_allFreqs, 'powers': store_window_power_overTime_allFreqs, 'freqs': slice_allfreqs}})


# --- Plot PSDs of windows with median and quartile vals
plt.figure(figsize=(6,10))
# plt.suptitle('Filtered LFP envelope Mean power')
for plot_pop_ind, plot_pop in enumerate(store_window_power_allFreqs_dict.keys()):
    plt.subplot(2,1,plot_pop_ind+1)

    square = patches.Rectangle((spindle_freq_range[0], 0.02e-5), spindle_freq_range[1]-spindle_freq_range[0], 4.98e-5, edgecolor='none', facecolor='k', alpha=0.05)
    ax=plt.gca()
    ax.add_patch(square)

    for sim_index, sim_name in enumerate(store_window_power_allFreqs_dict[plot_pop].keys()): 

        store_powers={freq:[] for freq in store_window_power_allFreqs_dict[plot_pop][sim_name]['freqs']}

        for freq_powers in store_window_power_allFreqs_dict[plot_pop][sim_name]['powers']:
            for freq_ind, freq in enumerate(store_window_power_allFreqs_dict[plot_pop][sim_name]['freqs']):
                store_powers[freq].append(freq_powers[freq_ind])


        if (sim_index == 0):    
            color = 'mediumturquoise'
            marker='.'
        elif (sim_index == 1):    
            color = 'mediumpurple'
            marker='.'


        medians = [np.median(    store_powers[freq])     for freq in store_powers.keys()]
        q25     = [np.percentile(store_powers[freq], 25) for freq in store_powers.keys()]
        q75     = [np.percentile(store_powers[freq], 75) for freq in store_powers.keys()]
        conds   = list(store_powers.keys())

        plt.plot(conds, medians, marker='none', color=color,   label='Median')  # Line for medians
        plt.fill_between(conds,  q25, q75,   color=color,   alpha=0.4, label='IQR (25th-75th)')  # Shaded IQR
        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
    
    if plot_pop_ind==len(plot_pops)-1:
        plt.xticks([0,12,24,36,48])
        plt.xlabel('Frequency (Hz)')
    else:
        plt.xticks([])
        plt.xlabel('')
    plt.ylabel(f'Power (mV$^2$)')
    plt.ylim([0,5e-5])
    plt.yticks([0,2.0e-5,4.0e-5])
    # plt.ylim([0,2.0e-5])
    # plt.yticks([0,1.0e-5,2.0e-5])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.savefig(lfp_psd_windowed_data_path+'lfpPSD_windowed_powerTraces.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_windowed_powerTraces.png',dpi=200)

for pop_ind,plot_pop in enumerate(store_window_power_overTime_dict.keys()):
    print('\n\n Comparing sim differences in wake vs sleep')
    print(f"Plot pop: {plot_pop}")
    
    test_altenative = 'greater'
    test_altenative = 'less'
    # test_altenative = 'two-sided'

    for sim_combination in [['sim_0','sim_1','two-sided'],
                            ]:
        print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}, {sim_combination[2]}]")

        groupA = store_window_power_overTime_dict[plot_pop][sim_combination[0]]['power']  # Replace with the first group of interest
        groupB = store_window_power_overTime_dict[plot_pop][sim_combination[1]]['power']  # Replace with the second group of interest

        # U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
        U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
        # if p_mw<0.05:   significance = '*'
        # else:           significance = 'ns'
        significance = significance_symbol(p_mw)
        print(f"\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")

plt.figure(figsize=[5,10])
for plot_pop_ind, plot_pop in enumerate(store_window_power_overTime_dict.keys()):
    plt.subplot(2,1,plot_pop_ind+1)

    for sim_index, sim_name in enumerate(store_window_power_overTime_dict[plot_pop].keys()): 

        if (sim_index == 0) or (sim_index == 1):
            if (sim_index == 0):    
                medianline_color = 'w'
                box_color = ['mediumturquoise']
            elif (sim_index == 1):    
                medianline_color = 'w'
                box_color = ['mediumpurple']
            alpha = 1.0
        else:
            alpha = 1.0

        ax = plt.gca()
        boxes = ax.boxplot(store_window_power_overTime_dict[plot_pop][sim_name]['power'], 
                        #    labels=[state_label], 
                        #    notch=True,
                            # medianprops=dict(color="k", alpha=0.7),
                            medianprops=dict(color=medianline_color, alpha=1.0), 
                            patch_artist=True, 
                            positions=np.array([sim_index]),
                            widths=0.15,
                            )
        # Set the face color of each box
        for box, color in zip(boxes['boxes'], box_color):
            box.set_facecolor(color)
            box.set_alpha(alpha)

    # plt.plot([0,1],[1.8e-4,1.8e-4],'-k') # Stats significance line

    plt.xlim([-1,2])
    plt.xticks([0,1],['Sleep','Wake'],rotation=45)
    plt.ylim([0,3.5e-4])
    plt.yticks([0,1.0e-4,2.0e-4,3.0e-4])
    # plt.ylim([0,2.0e-4])
    # plt.yticks([0,1.0e-4,2.0e-4])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.ylabel(f'Band power [8-16 Hz] (mV$^2$)')
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
# plt.savefig(lfp_psd_windowed_data_path+'lfpPSD_windowed_powerInBand_8_16_Hz_distribution.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_windowed_powerInBand_8_16_Hz_distribution.png',dpi=200)

################################################################################################################################################
#####   Figure 02: LFP Traces data   #####
# --- Data path
# sleep_vs_wake_lfp_traces_data_path  = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401a3/sim_output/lfpTraces/'
sleep_vs_wake_lfp_traces_data_path  = simulation_data_path+'lfpTraces/'
sleep_vs_wake_data_folder_lfpTraces = sleep_vs_wake_lfp_traces_data_path+'data/'

sleep_vs_wake_list_files_lfpTraces = os.listdir(sleep_vs_wake_data_folder_lfpTraces)
sleep_vs_wake_list_files_lfpTraces = [file for file in sleep_vs_wake_list_files_lfpTraces if '_LFP_data.json' in file]
sleep_vs_wake_list_files_lfpTraces.sort()

store_sleep_vs_wake_lfpTraces_dict={}
for sleep_vs_wake_data_file_lfpTraces_ind, sleep_vs_wake_data_file_lfpTraces in enumerate(sleep_vs_wake_list_files_lfpTraces):
    filename_json = sleep_vs_wake_data_folder_lfpTraces+sleep_vs_wake_data_file_lfpTraces
    with open(filename_json, 'r') as file: sleep_vs_wake_lfpTraces_dict = json.load(file)
    store_sleep_vs_wake_lfpTraces_dict.update({'sim_'+str(sleep_vs_wake_data_file_lfpTraces_ind):sleep_vs_wake_lfpTraces_dict})

################################################################################################################################################

# store_sleep_vs_wake_lfpTraces_dict[<simulation>][<recording_frequency>][dict_keys(['data', 't', 'names'])][<electrode_inxed>: 0='all'; 1=0; 2=1]

# recording_frequency = '0'   # [1, 400] Hz
recording_frequency = '1'   # [1,  50] Hz

select_sims=['sim_0','sim_1']
sim_cmap_color=[0,0.2,0.8,0,0]
sim_color=['mediumturquoise','mediumpurple','lightgrey','lightgrey','lightgrey']
sim_linewidth =[1.5,1.5,1,1,1]
# sim_linewidth =[1,2,2,1,1]
# select_sims=None

for sim_index, sim_name in enumerate(store_sleep_vs_wake_lfpTraces_dict.keys()):
    
    if select_sims is not None:
        if sim_name not in select_sims: continue

    plt.figure(figsize=[20,5])

    # plt.suptitle(sim_name)
    # TRN
    cmap = matplotlib.colormaps[pop_color_maps['TRN__pop']]
    line_color = sim_color[sim_index]
    # line_color = cmap(sim_cmap_color[sim_index])
    line_width = sim_linewidth[sim_index]
    plt.subplot(2,1,1)
    plt.plot(store_sleep_vs_wake_lfpTraces_dict[sim_name][recording_frequency]['t'],store_sleep_vs_wake_lfpTraces_dict[sim_name][recording_frequency]['data'][1],color=line_color, linewidth=line_width)
    plt.ylim([-0.02,0.02])
    plt.yticks([])
    # plt.yticks([-0.02,0,0.02])

    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # plt.xlim([1000,8000])
    plt.xlim([1500,16000])


    plt.xticks([])
    plt.xlabel('')

    # VPM
    print(sim_index,'VPM__pop',pop_color_maps['VPM__pop'],cmap,line_color)
    cmap = matplotlib.colormaps[pop_color_maps['VPM__pop']]
    line_color = sim_color[sim_index]
    # line_color = cmap(sim_cmap_color[sim_index])
    line_width = sim_linewidth[sim_index]
    plt.subplot(2,1,2)
    plt.plot(store_sleep_vs_wake_lfpTraces_dict[sim_name][recording_frequency]['t'],store_sleep_vs_wake_lfpTraces_dict[sim_name][recording_frequency]['data'][2],color=line_color, linewidth=line_width)
    plt.ylim([-0.02,0.02])
    plt.yticks([])
    # plt.yticks([-0.02,0,0.02])

    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    if sim_name == 'sim_1': 
        plt.plot([1500,2000],[-0.015,-0.015],'k',linewidth=3)
        plt.plot([1500,1500],[ 0.005,-0.015],'k',linewidth=3)
        plt.text(1320, -0.020, '500 ms', color='k', rotation=0, size=15, weight=font_weight, alpha = 1.0)
        plt.text(1125, -0.013, '0.02 mV', color='k', rotation=90,size=15, weight=font_weight, alpha = 1.0)
    
    plt.xticks([])
    # plt.xlim([1000,8000])
    plt.xlim([1500,16000])
    # plt.xlabel('Time (ms)')
    # plt.savefig(sleep_vs_wake_lfp_traces_data_path+'lfpTraces_'+sim_name+'.png',dpi=200)
    plt.savefig(save_figures_path+'lfpTraces_'+sim_name+'.png',dpi=200)

################################################################################################################################################
#####   Figure 02: Mean Voltage Traces data   #####
# --- Data path
# sleep_vs_wake_meanVoltage_data_path  = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401a3/sim_output/mean_voltage/'
sleep_vs_wake_meanVoltage_data_path  = simulation_data_path+'mean_voltage/'
sleep_vs_wake_data_folder_meanVoltage = sleep_vs_wake_meanVoltage_data_path+'data/'

sleep_vs_wake_list_files_meanVoltage = os.listdir(sleep_vs_wake_data_folder_meanVoltage)
sleep_vs_wake_list_files_meanVoltage = [file for file in sleep_vs_wake_list_files_meanVoltage if '_meanVoltage__data.json' in file]
sleep_vs_wake_list_files_meanVoltage.sort()

store_sleep_vs_wake_meanVoltage_dict={}
for sleep_vs_wake_data_file_meanVoltage_ind, sleep_vs_wake_data_file_meanVoltage in enumerate(sleep_vs_wake_list_files_meanVoltage):
    filename_json = sleep_vs_wake_data_folder_meanVoltage+sleep_vs_wake_data_file_meanVoltage
    with open(filename_json, 'r') as file: sleep_vs_wake_meanVoltage_dict = json.load(file)
    store_sleep_vs_wake_meanVoltage_dict.update({'sim_'+str(sleep_vs_wake_data_file_meanVoltage_ind):sleep_vs_wake_meanVoltage_dict})

plt.show()
sys.exit()
################################################################################################################################################
import pywt
# def compute_psd_morlet(time_series, dt, freqs):
#     """
#     Computes the power spectral density (PSD) of a time series using the Morlet wavelet.

#     Parameters:
#     - time_series: 1D numpy array of voltage values (mean voltage trace)
#     - dt: Time step (sampling interval) in ms
#     - freqs: 1D numpy array of frequencies for the wavelet transform (in Hz)

#     Returns:
#     - power: 2D numpy array of wavelet power (rows: frequencies, cols: time points)
#     """
#     scales = 1 / (2 * np.pi * freqs * dt * 1e-3)  # Convert Hz to scale
#     cwt_matrix, _ = pywt.cwt(time_series, scales, 'cmor', sampling_period=dt)
#     power = np.abs(cwt_matrix) ** 2  # Compute power

#     return power

def compute_psd_morlet(time_series, dt, freqs, wavelet='cmor1.5-1.0'):
    """
    Computes the power spectral density (PSD) of a time series using the Morlet wavelet.

    Parameters:
    - time_series: 1D numpy array of voltage values (mean voltage trace)
    - dt: Time step (sampling interval) in ms
    - freqs: 1D numpy array of frequencies for the wavelet transform (in Hz)
    - wavelet: Morlet wavelet name with adjustable bandwidth and center frequency

    Returns:
    - power: 2D numpy array of wavelet power (rows: frequencies, cols: time points)
    """
    scales = 1 / (2 * np.pi * freqs * dt * 1e-3)  # Convert Hz to scale
    cwt_matrix, _ = pywt.cwt(time_series, scales, wavelet, sampling_period=dt)
    
    power = np.abs(cwt_matrix) ** 2  # Compute power
    power /= np.max(power)  # Normalize power for better visualization

    return power

# recording_frequency = '0'   # [1, 400] Hz
recording_frequency = '1'   # [1,  50] Hz

select_sims=['sim_0','sim_1']
# select_sims=None

dt = 0.1

binSize         = 1, 
minFreq         = 1, 
maxFreq         = 100, 
transformMethod = 'morlet', 
stepFreq        = 1, 
NFFT            = 256, 
noverlap        = 128, 
smooth          = 0, 
# plt.psd(x, NFFT=256, Fs=1/(0.1e-3), Fc=None, detrend=None, window=None, noverlap=128, pad_to=None, sides=None, scale_by_freq=None, return_line=None, *, data=None, **kwargs)

freqs = np.linspace(1, 50, 50)  # Frequencies from 1 to 100 Hz
skip_time_ms = [1500,500] # ms
    
for sim_index, sim_name in enumerate(store_sleep_vs_wake_meanVoltage_dict.keys()):
    
    if select_sims is not None:
        if sim_name not in select_sims: continue

    t_window        = [0,16000]
    # t_window        = [0,8000]
    # t_window        = [1500,8000]

    t_vec_TRN       = np.array(store_sleep_vs_wake_meanVoltage_dict[sim_name]['TRN__pop']['t'][int(t_window[0]/dt):int(t_window[1]/dt)])
    mean_v_vec_TRN  = np.array(store_sleep_vs_wake_meanVoltage_dict[sim_name]['TRN__pop']['mean_v'][int(t_window[0]/dt):int(t_window[1]/dt)])

    t_vec_VPM       = np.array(store_sleep_vs_wake_meanVoltage_dict[sim_name]['VPM__pop']['t'][int(t_window[0]/dt):int(t_window[1]/dt)])
    mean_v_vec_VPM  = np.array(store_sleep_vs_wake_meanVoltage_dict[sim_name]['VPM__pop']['mean_v'][int(t_window[0]/dt):int(t_window[1]/dt)])

    plt.figure(figsize=[20,5])

    # plt.suptitle(sim_name)

    # TRN
    plt.subplot(2,2,1)
    # psd = compute_psd_morlet(mean_v_vec_TRN, dt, freqs, wavelet='cmor1.5-1.0')
    # plt.imshow(psd, aspect='auto', extent=[t_vec_TRN[0], t_vec_TRN[-1], freqs[0], freqs[-1]], origin='lower', cmap='jet')
    # plt.colorbar(label="Power")
    
    # Extract time step (dt in ms)
    dt = np.mean(np.diff(t_vec_TRN))  # Should be ~0.1 ms
    T_sec = (t_vec_TRN[-1] - t_vec_TRN[0]) / 1000  # Convert total time to seconds
    
    # Define frequency range (log-spaced for better resolution)
    freqs = np.logspace(np.log10(1), np.log10(50), 100)  # 1 to 50 Hz, 100 steps
    # freqs = np.linspace(1, 50, 50)  # Frequencies from 1 to 100 Hz
    
    # Compute scales for Morlet wavelet
    scales = 1 / (2 * np.pi * freqs * dt * 1e-3)  # Convert Hz to scale
    
    # Compute Continuous Wavelet Transform (CWT)
    cwt_matrix, _ = pywt.cwt(mean_v_vec_TRN, scales, 'cmor1.5-1.0', sampling_period=dt)
    
    # Compute power (square of magnitude)
    power = np.abs(cwt_matrix) ** 2
    power /= np.max(power)  # Normalize for visualization
    
    # Remove first and last 100 ms from the time vector
    t_start, t_end = t_vec_TRN[0] + skip_time_ms[0], t_vec_TRN[-1] - skip_time_ms[1]
    
    # Find the indices that correspond to the cropped time range
    valid_indices = np.where((t_vec_TRN >= t_start) & (t_vec_TRN <= t_end))[0]
    
    # Crop power and time vector
    t_vec_cropped = t_vec_TRN[valid_indices]
    power_cropped = power[:, valid_indices]  # Keep all frequencies, crop only time
    
    # Plotting
    plt.imshow(power_cropped, aspect='auto', extent=[t_vec_cropped[0], t_vec_cropped[-1], freqs[0], freqs[-1]],
               origin='lower', cmap='jet', interpolation='bilinear')
    plt.colorbar(label="Normalized Power")
    plt.xlabel("Time (ms)")
    plt.ylabel("Frequency (Hz)")
    plt.yscale("log")  # Log scale for frequency

    # plt.psd(store_sleep_vs_wake_meanVoltage_dict[sim_name]['TRN__pop']['mean_v'], NFFT=256, Fs=1/(0.1e-3), Fc=12, detrend=None, window=None, noverlap=128, pad_to=None, sides=None, scale_by_freq=None)
    # plt.plot(store_sleep_vs_wake_meanVoltage_dict[sim_name][recording_frequency]['t'],store_sleep_vs_wake_meanVoltage_dict[sim_name][recording_frequency]['data'][1],color='blue')
    # plt.ylim([-0.02,0.02])
    # plt.yticks([])
    # # plt.yticks([-0.02,0,0.02])

    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # plt.xlim([1000,8000])
    plt.xlim([1500,16000])

    plt.xticks([])
    plt.xlabel('')

    plt.subplot(2,2,2)
    plt.plot(t_vec_TRN,mean_v_vec_TRN)


    # VPM
    plt.subplot(2,2,3)
    # psd = compute_psd_morlet(mean_v_vec_VPM, dt, freqs, wavelet='cmor1.5-1.0')
    # plt.imshow(psd, aspect='auto', extent=[t_vec_VPM[0], t_vec_VPM[-1], 
    #                                        freqs[0], freqs[-1]], origin='lower', cmap='jet')
    # plt.colorbar(label="Power")


    # Extract time step (dt in ms)
    dt = np.mean(np.diff(t_vec_VPM))  # Should be ~0.1 ms
    T_sec = (t_vec_VPM[-1] - t_vec_VPM[0]) / 1000  # Convert total time to seconds
        
    # Define frequency range (log-spaced for better resolution)
    freqs = np.logspace(np.log10(1), np.log10(50), 100)  # 1 to 50 Hz, 100 steps
    # freqs = np.linspace(1, 50, 50)  # Frequencies from 1 to 100 Hz

    # Compute scales for Morlet wavelet
    scales = 1 / (2 * np.pi * freqs * dt * 1e-3)  # Convert Hz to scale
    
    # Compute Continuous Wavelet Transform (CWT)
    cwt_matrix, _ = pywt.cwt(mean_v_vec_VPM, scales, 'cmor1.5-1.0', sampling_period=dt)
    
    # Compute power (square of magnitude)
    power = np.abs(cwt_matrix) ** 2
    power /= np.max(power)  # Normalize for visualization
    
    # Remove first and last 100 ms from the time vector
    t_start, t_end = t_vec_VPM[0] + skip_time_ms[0], t_vec_VPM[-1] - skip_time_ms[1]
    
    # Find the indices that correspond to the cropped time range
    valid_indices = np.where((t_vec_VPM >= t_start) & (t_vec_VPM <= t_end))[0]
    
    # Crop power and time vector
    t_vec_cropped = t_vec_VPM[valid_indices]
    power_cropped = power[:, valid_indices]  # Keep all frequencies, crop only time
    
    # Plotting
    plt.imshow(power_cropped, aspect='auto', extent=[t_vec_cropped[0], t_vec_cropped[-1], freqs[0], freqs[-1]],
               origin='lower', cmap='jet', interpolation='bilinear')
    plt.colorbar(label="Normalized Power")
    plt.xlabel("Time (ms)")
    plt.ylabel("Frequency (Hz)")
    plt.yscale("log")  # Log scale for frequency

    # plt.psd(store_sleep_vs_wake_meanVoltage_dict[sim_name]['VPM__pop']['mean_v'], NFFT=256, Fs=1/(0.1e-3), Fc=12, detrend=None, window=None, noverlap=128, pad_to=None, sides=None, scale_by_freq=None)
    # plt.plot(store_sleep_vs_wake_meanVoltage_dict[sim_name][recording_frequency]['t'],store_sleep_vs_wake_meanVoltage_dict[sim_name][recording_frequency]['data'][2],color='green')
    # plt.ylim([-0.02,0.02])
    # plt.yticks([])
    # # plt.yticks([-0.02,0,0.02])

    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)



    plt.subplot(2,2,4)
    plt.plot(t_vec_VPM,mean_v_vec_VPM)


    # if sim_name == 'sim_1': 
    #     plt.plot([1250,1750],[-0.015,-0.015],'k',linewidth=3)
    #     plt.plot([1250,1250],[ 0.005,-0.015],'k',linewidth=3)
    #     plt.text(1320, -0.013, '500 ms', color='k', rotation=0, size=15, weight=font_weight, alpha = 1.0)
    #     plt.text(1125, -0.013, '0.02 mV', color='k', rotation=90,size=15, weight=font_weight, alpha = 1.0)
    
    # plt.xticks([])
    # plt.xlim([1000,8000])
    # plt.xlabel('Time (ms)')
    # plt.savefig(sleep_vs_wake_meanVoltage_data_path+'meanVoltage_'+sim_name+'.png',dpi=200)
    plt.savefig(save_figures_path+'meanVoltage_'+sim_name+'.png',dpi=200)


################################################################################################################################################
plt.show()

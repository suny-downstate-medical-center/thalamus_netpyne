################################################################################################################################################

import os
import pickle
import netpyne
# from netpyne import sim
import json

import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp

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

import math

########################################################################################################################################################################################################
# --- Setting the path for loading the datasets
data_path = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data'
# batch_sim_group = 'paper_sims'
# batch_sim_label = 'paperSim_Fig07_sleep'
# sim_output_folder = 'sim_output'

batch_sim_group = 'paper_002'
batch_sim_label = batch_sim_group+'_'+'barreloid_batch_fig06_sleep'

# # batch_sim_group = 'paper_00X_3'
# batch_sim_group = 'paper_00X_3_1'
# # batch_sim_group = 'paper_00X_3_2'
# batch_sim_label = batch_sim_group+'_'+'test_TRN_RMP'


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

import matplotlib.pyplot as plt

def draw_vertical_bracket(ax, x, y, height, is_left=True, line_width=1, color='black'):
    """
    Draws a square bracket on the given axes.

    Args:
        ax: The matplotlib axes to draw on.
        x: The x-coordinate of the bracket's base.
        y: The y-coordinate of the bracket's base.
        height: The height of the bracket.
        is_left: If True, draws a left bracket; otherwise, draws a right bracket.
        line_width: The line width of the bracket.
        color: The color of the bracket.
    """
    if is_left:
        ax.plot([x, x], [y, y + height], lw=line_width, color=color)
        ax.plot([x, x + height / 5], [y, y], lw=line_width, color=color)
        ax.plot([x, x + height / 5], [y + height, y + height], lw=line_width, color=color)
    else:
        ax.plot([x, x], [y, y + height], lw=line_width, color=color)
        ax.plot([x, x - height / 5], [y, y], lw=line_width, color=color)
        ax.plot([x, x - height / 5], [y + height, y + height], lw=line_width, color=color)

# # Example usage
# fig, ax = plt.subplots()

# # Set plot limits and remove ticks
# ax.set_xlim(0, 10)
# ax.set_ylim(0, 10)
# ax.set_xticks([])
# ax.set_yticks([])

# # Draw left and right square brackets
# draw_square_bracket(ax, 2, 2, 3, is_left=True, line_width=2, color='red')
# draw_square_bracket(ax, 8, 2, 3, is_left=False, line_width=2, color='blue')

# plt.show()

import matplotlib.pyplot as plt

def draw_horizontal_bracket(ax, x, y, width, v_height, is_top=True, line_width=1, color='black', alpha=1.0):
    """
    Draws a horizontal square bracket with specified vertical line height.

    Args:
        ax: The matplotlib axes to draw on.
        x: The x-coordinate of the bracket's starting point.
        y: The y-coordinate of the bracket's base.
        width: The width of the bracket (horizontal line length).
        v_height: The height of the vertical lines.
        is_top: If True, draws a top bracket; otherwise, draws a bottom bracket.
        line_width: The line width of the bracket.
        color: The color of the bracket.
    """
    if is_top:
        ax.plot([x, x + width],             [y, y],             lw=line_width, color=color, alpha=alpha)  # Top horizontal line
        ax.plot([x, x],                     [y, y - v_height],  lw=line_width, color=color, alpha=alpha)  # Left vertical line
        ax.plot([x + width, x + width],     [y, y - v_height],  lw=line_width, color=color, alpha=alpha)  # Right vertical line
    else:
        ax.plot([x, x + width],             [y, y],             lw=line_width, color=color, alpha=alpha)  # Bottom horizontal line
        ax.plot([x, x],                     [y, y + v_height],  lw=line_width, color=color, alpha=alpha)  # Left vertical line
        ax.plot([x + width, x + width],     [y, y + v_height],  lw=line_width, color=color, alpha=alpha)  # Right vertical line

# # Example usage
# fig, ax = plt.subplots()

# # Set plot limits and remove ticks
# ax.set_xlim(0, 10)
# ax.set_ylim(0, 10)
# ax.set_xticks([])
# ax.set_yticks([])

# # Draw top and bottom square brackets with different vertical heights
# draw_horizontal_bracket(ax, 2, 8, 4, 1, is_top=True, line_width=2, color='red')
# draw_horizontal_bracket(ax, 2, 2, 4, 2, is_top=False, line_width=2, color='blue')

# plt.show()

################################################################################################################################################
#####   Spike PSD data   #####
# --- Data path
spike_psd_data_path   = simulation_data_path+'spikePSD_custom/'
data_folder = spike_psd_data_path+'data/'

list_files = os.listdir(data_folder)
list_files_PSD = [file for file in list_files if '_ratePSD_values.json' in file]
list_files_PSD.sort()

store_PSD_dict={}
for data_file_PSD_ind, data_file_PSD in enumerate(list_files_PSD):
    filename_json = data_folder+data_file_PSD
    with open(filename_json, 'r') as file: PSD_dict = json.load(file)
    store_PSD_dict.update({'sim_'+str(data_file_PSD_ind):PSD_dict})

################################################################################################################################################
# --- Figure properties
plt.rcParams.update({'font.size': 15})

sim_index_mapping   = [-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0]
spindle_freq_range  = [8,16] # Hz
# y_lims              = {'VPM__pop':[-1,26],'TRN__pop':[-1,91]}
y_lims              = {'VPM__pop':[-0.1,2.6],'TRN__pop':[-0.1,11]}
y_text_offset       = {'VPM__pop':1,'TRN__pop':4}
limit_frequency     = 50 # Hz

bar_y_lims          = {'VPM__pop':[-1,11],'TRN__pop':[-1,41]}
bar_y_ticks         = {'VPM__pop':[0,5,10],'TRN__pop':[0,20,40]}

plot_pops=['TRN__pop','VPM__pop']

################################################################################################################################################

select_sims = ['sim_2','sim_6']
# --- PSD traces
plt.figure(figsize=[5,5])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)

    square = patches.Rectangle((spindle_freq_range[0], y_lims[plot_pop][0]), spindle_freq_range[1]-spindle_freq_range[0], y_lims[plot_pop][1], edgecolor='k', facecolor='w', alpha=0.1)
    ax=plt.gca()
    ax.add_patch(square)

    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_PSD_dict.keys()):  
        if select_sims is not None:
            if sim_name not in select_sims: continue
       
        if (sim_index == 2) or (sim_index == 6):
            if (sim_index == 2):    line_color = 'k'
            else:                   line_color = 'grey'
            line_alpha = 1.0
            line_width = 4
        else:
            line_alpha = 0.5
            line_width = 0.5
            line_color = 'k'
        
        plt.plot(store_PSD_dict[sim_name][plot_pop]['freqs'],store_PSD_dict[sim_name][plot_pop]['signal'],
                                                                          color=line_color,
                                                                        #   color=cmap(sim_index/len(store_PSD_dict.keys())),
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

# plt.savefig(spike_psd_data_path+'spikePSD_signal.png',dpi=200)
plt.savefig(save_figures_path+'spikePSD_signal.png',dpi=200)

# --- PSD bar plot of max power value and frequency of max power
plt.figure(figsize=[3,5])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_PSD_dict.keys()): 
        if select_sims is not None:
            if sim_name not in select_sims: continue

        slice_data  = [(freq,store_PSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_PSD_dict[sim_name][plot_pop]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
        max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
        
        full_data = [(freq,store_PSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_PSD_dict[sim_name][plot_pop]['freqs']) if freq<=limit_frequency]
        max_power_signal        = max(full_data, key=lambda x: x[1])


        # Extract frequency and power values
        slice_freqs, slice_power = zip(*slice_data)
        # Compute area under the curve using the trapezoidal rule
        power_in_band = np.trapz(slice_power, slice_freqs)

        bar_value = power_in_band
        bar_string = ''
        bar_string_color = 'k'
        font_weight='bold'
        
        if (sim_index == 2) or (sim_index == 6):
            if (sim_index == 2):    bar_color = 'k'
            else:                   bar_color = 'grey'
            bar_alpha = 1.0
            edge_color = 'k'
        else:
            bar_color = 'w'
            bar_alpha = 0.25
            edge_color = 'k'

        plt.bar(sim_index_mapping[sim_index],bar_value,width=0.21,color=bar_color, edgecolor = edge_color, alpha=bar_alpha)
        # plt.text(sim_index_mapping[sim_index]-0.075, bar_value+y_text_offset[plot_pop], bar_string, color=bar_string_color, rotation=0, size=15, weight=font_weight, alpha = bar_alpha)

    plt.ylim(bar_y_lims[plot_pop])
    plt.yticks(bar_y_ticks[plot_pop])
    plt.ylabel('Power 8-16 Hz')
    
    if plot_pop_ind==len(plot_pops)-1:
        # plt.xticks([-1,0,1])
        plt.xticks([-0.5,0.5])
        # plt.xlabel('$P_{T}$')
        plt.xlabel('CSI')
    else:
        plt.xticks([-0.5,0.5],['',''])
        plt.xlabel('')
    plt.xlim([-1.1,1.1])
    
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    plt.tight_layout()

# plt.savefig(spike_psd_data_path+'spikePSD_peakFrequency.png',dpi=200)
plt.savefig(save_figures_path+'spikePSD_peakFrequency.png',dpi=200)









################################################################################################################################################
#####   LFP PSD data   #####
# --- Data path
# lfp_psd_data_path   = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401dd3/sim_output/lfpPSD/'
# lfp_psd_data_path   = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401dd4/sim_output/lfpPSD/'
# lfp_psd_data_path   = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401dd4_2dot5Hz/sim_output/lfpPSD/'
# lfp_psd_data_path   = '/Users/joao/Research/Models/BBP/thalamus_netpyne/paper_analysis/bWgt_9401/bWgt_9401dd4_1dot0Hz/sim_output/lfpPSD/'
lfp_psd_data_path   = simulation_data_path+'lfpPSD/'
data_folder_lfpPSD = lfp_psd_data_path+'data/'

list_files_lfpData = os.listdir(data_folder_lfpPSD)
list_files_lfpPSD = [file for file in list_files_lfpData if '_lfpPSD_values.json' in file]
list_files_lfpPSD.sort()

store_lfpPSD_dict={}
for data_file_lfpPSD_ind, data_file_lfpPSD in enumerate(list_files_lfpPSD):
    filename_json = data_folder_lfpPSD+data_file_lfpPSD
    with open(filename_json, 'r') as file: lfpPSD_dict = json.load(file)
    store_lfpPSD_dict.update({'sim_'+str(data_file_lfpPSD_ind):lfpPSD_dict})

################################################################################################################################################

lfpPSD_y_lims              = {'VPM__pop':[-0.01e-5,0.5e-5],'TRN__pop':[-0.01e-5,0.5e-5]}
lfpPSD_y_lims_bar          = {'VPM__pop':[-0.01e-5,2.6e-5],'TRN__pop':[-0.01e-5,2.6e-5]}

# --- LFP traces
plt.figure(figsize=[5,5])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)

    square = patches.Rectangle((spindle_freq_range[0], lfpPSD_y_lims[plot_pop][0]), spindle_freq_range[1]-spindle_freq_range[0], lfpPSD_y_lims[plot_pop][1], edgecolor='k', facecolor='w', alpha=0.1)
    ax=plt.gca()
    ax.add_patch(square)

    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_lfpPSD_dict.keys()):
        if (sim_index == 2) or (sim_index == 6):
            if (sim_index == 2):    line_color = 'k'
            else:                   line_color = 'grey'
            line_alpha = 1.0
            line_width = 4
        else:
            line_alpha = 0.5
            line_width = 0.5
            line_color = 'k'
        
        plt.plot(store_lfpPSD_dict[sim_name][plot_pop]['freqs'],store_lfpPSD_dict[sim_name][plot_pop]['signal'],
                                                                          color=line_color,
                                                                        #   color=cmap(sim_index/len(store_lfpPSD_dict.keys())),
                                                                          linewidth=line_width, 
                                                                        #   linewidth=2+(0.15*sim_index), 
                                                                          alpha=line_alpha)
    plt.ylim(lfpPSD_y_lims[plot_pop])
    # plt.yticks([0,1.25e-5,2.5e-5])
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

# plt.savefig(lfp_psd_data_path+'lfpPSD_signal.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_signal.png',dpi=200)

# --- PSD bar plot of max power value and frequency of max power
plt.figure(figsize=[3,5])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,plot_pop_ind+1)
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    for sim_index, sim_name in enumerate(store_lfpPSD_dict.keys()): 

        slice_data  = [(freq,store_lfpPSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_dict[sim_name][plot_pop]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
        max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
        
        full_data = [(freq,store_lfpPSD_dict[sim_name][plot_pop]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_dict[sim_name][plot_pop]['freqs']) if freq<=limit_frequency]
        max_power_signal        = max(full_data, key=lambda x: x[1])


        # Extract frequency and power values
        slice_freqs, slice_power = zip(*slice_data)
        # Compute area under the curve using the trapezoidal rule
        power_in_band = np.trapz(slice_power, slice_freqs)

        bar_value = power_in_band
        bar_string = ''
        bar_string_color = 'k'
        font_weight='bold'
        
        if (sim_index == 2) or (sim_index == 6):
            if (sim_index == 2):    bar_color = 'k'
            else:                   bar_color = 'grey'
            bar_alpha = 1.0
            edge_color = 'k'
        else:
            bar_color = 'w'
            bar_alpha = 0.25
            edge_color = 'k'


        plt.bar(sim_index_mapping[sim_index],bar_value,width=0.21,color=bar_color, edgecolor = edge_color, alpha=bar_alpha)

    plt.ylim(lfpPSD_y_lims_bar[plot_pop])
    plt.ylabel('Power 8-16 Hz')
    # plt.yticks([0,5.0e-5,10.0e-5])
    
    if plot_pop_ind==len(plot_pops)-1:
        plt.xticks([-1,0,1])
        # plt.xlabel('$P_{T}$')
        plt.xlabel('CSI')
    else:
        plt.xticks([])
        plt.xlabel('')
    
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()

# plt.savefig(lfp_psd_data_path+'lfpPSD_peakFrequency.png',dpi=200)
plt.savefig(save_figures_path+'lfpPSD_peakFrequency.png',dpi=200)

################################################################################################################################################
#####   LFP Traces data   #####
# --- Data path
lfp_traces_data_path   = simulation_data_path+'lfpTraces/'
data_folder_lfpTraces = lfp_traces_data_path+'data/'

list_files_lfpTraces = os.listdir(data_folder_lfpTraces)
list_files_lfpTraces = [file for file in list_files_lfpTraces if '_LFP_data.json' in file]
list_files_lfpTraces.sort()

store_lfpTraces_dict={}
for data_file_lfpTraces_ind, data_file_lfpTraces in enumerate(list_files_lfpTraces):
    filename_json = data_folder_lfpTraces+data_file_lfpTraces
    with open(filename_json, 'r') as file: lfpTraces_dict = json.load(file)
    store_lfpTraces_dict.update({'sim_'+str(data_file_lfpTraces_ind):lfpTraces_dict})

################################################################################################################################################

from scipy.signal import butter, filtfilt
import numpy as np

def bandpass_filter(data, fs, lowcut, highcut, order=4):
    """
    Apply a bandpass Butterworth filter to the LFP signal.
    
    Parameters:
        data (array-like): LFP signal.
        fs (float): Sampling frequency in Hz.
        lowcut (float): Low cutoff frequency in Hz.
        highcut (float): High cutoff frequency in Hz.
        order (int): Order of the Butterworth filter.

    Returns:
        np.array: Filtered signal.
    """
    # Check for NaNs or Infs
    data = np.array(data)
    if np.any(np.isnan(data)) or np.any(np.isinf(data)):
        print("Warning: LFP signal contains NaNs or Infs. Replacing them.")
        data = np.nan_to_num(data)

    # Ensure valid frequency range
    nyquist = 0.5 * fs
    lowcut = max(0.1, lowcut)  # Avoid 0 Hz issues
    highcut = min(0.9 * nyquist, highcut)  # Avoid Nyquist limit

    low = lowcut / nyquist
    high = highcut / nyquist

    if not (0 < low < 1) or not (0 < high < 1) or low >= high:
        raise ValueError(f"Invalid cutoff frequencies: low={low}, high={high}. Check fs and cutoffs.")

    # Ensure data length is sufficient
    if len(data) < 3 * order:
        raise ValueError(f"LFP signal is too short for filtering. Minimum required: {3 * order} samples.")

    # Design Butterworth filter
    b, a = butter(3, [low, high], btype='band')

    # Apply zero-phase filtering
    filtered_data = filtfilt(b, a, data, padtype='odd')

    return filtered_data

def moving_sum_avg(signal, fs, window_ms=500, overlap_ms=50, start_time=1500):
    """
    Computes the moving sum and moving average of a signal.
    
    Parameters:
        signal (array-like): The input signal.
        fs (int): Sampling frequency in Hz.
        window_ms (int): Window size in milliseconds. Default is 500 ms.
        overlap_ms (int): Overlap size in milliseconds. Default is 50 ms.

    Returns:
        tuple: (moving_sums, moving_averages, time_centers)
    """
    signal = np.array(signal)
    
    # Convert window and overlap from ms to samples
    window_size = int(window_ms * fs / 1000)  # 500 ms → 5000 samples
    step_size = window_size - int(overlap_ms * fs / 1000)  # 50 ms overlap → step of 4500 samples

    # Define indices for moving window
    num_windows = (len(signal) - window_size) // step_size + 1  # Number of windows
    moving_sums = np.zeros(num_windows)
    moving_averages = np.zeros(num_windows)
    time_centers = np.zeros(num_windows)

    for i in range(num_windows):
        start = i * step_size
        end = start + window_size
        segment = signal[start:end]

        moving_sums[i] = np.sum(segment)
        moving_averages[i] = np.mean(segment)
        time_centers[i] = start_time + ((start + end) * 1000) / (2 * fs)  # Convert to mili seconds

    return moving_sums, moving_averages, time_centers

from scipy.signal import hilbert

def compute_envelope(signal):
    """
    Computes the envelope of a signal using the Hilbert transform.

    Parameters:
        signal (array-like): Input signal.

    Returns:
        array: The envelope of the signal.
    """
    analytic_signal = hilbert(signal)  # Compute analytic signal
    envelope = np.abs(analytic_signal)  # Envelope is the magnitude

    return envelope

################################################################################################################################################

# store_lfpTraces_dict[<simulation>][<recording_frequency>][dict_keys(['data', 't', 'names'])][<electrode_inxed>: 0='all'; 1=0; 2=1]

# recording_frequency = '0'   # [1, 400] Hz
recording_frequency = '1'   # [1,  50] Hz

# select_sims=['sim_0','sim_2','sim_4','sim_6','sim_8']
# select_sims=['sim_2','sim_6']
select_sims=None

lowcut = 8  # Lower bound of frequency band (e.g., theta band 4-8 Hz)
highcut = 16  # Upper bound of frequency band
fs = 1000/0.1
window_ms=50
overlap_ms=5

sim_color=['lightgrey','lightgrey','k','lightgrey','lightgrey','lightgrey','grey','lightgrey','lightgrey']
sim_linewidth =1.5

skip_up_time = 0
ct_up_times=[   [ 2500, 3000],  [ 4500, 5000],  [ 6500, 7000],  [ 8500, 9000],  [10500,11000],  [12500,13000],  [14500,15000],]
import copy
ct_up_times=[   [up_start+skip_up_time,up_end] for [up_start,up_end] in ct_up_times]
ct_up_times_draw_window = copy.deepcopy(ct_up_times) # draws the window without considering the skip time during LFP data slicing
# ct_up_times=[]
# skip_down_time = 200 # ms
# ct_down_times=[ [ 3000, 4500],  [ 5000, 6500],  [ 7000, 8500],  [ 9000,10500],  [11000,12500],  [13000,14500],  [15000,16000]]
skip_down_time = 0 # ms
ct_down_times=[ [ 3500, 4000],  [ 5500, 6000],  [ 7500, 8000],  [ 9500,10000],  [11500,12000],  [13500,14000],  [15500,16000]] # testing using only 500 ms of activity in the middle of the down-state
ct_down_times = [[donw_start+skip_down_time,down_end-skip_down_time] for [donw_start,down_end] in ct_down_times]
ct_down_times_draw_window = copy.deepcopy(ct_down_times) # draws the window without considering the skip time during LFP data slicing
# ct_down_times=[]
mean_envelope_dict={}
store_lfp_rms_dict={}
store_lfp_mean_dict={}
pop_electrodes={'TRN__pop':1,'VPM__pop':2}

for sim_index, sim_name in enumerate(store_lfpTraces_dict.keys()):
    
    if select_sims is not None:
        if sim_name not in select_sims: continue

    mean_envelope_dict.update({sim_name:{}})
    store_lfp_rms_dict.update({sim_name:{}})
    store_lfp_mean_dict.update({sim_name:{}})

    for plot_pop in ['TRN__pop', 'VPM__pop']:
        store_lfp_rms_dict[sim_name].update({plot_pop:{}})
        store_lfp_mean_dict[sim_name].update({plot_pop:{}})

        plt.figure(figsize=[25,3])
        # plt.title(sim_name+'__'+plot_pop)
        # plt.subplot(2,1,1)
        # UP times square patches
        if len(ct_up_times)>0:
            for ct_up_time in ct_up_times_draw_window:        
                square = patches.Rectangle((ct_up_time[0], -0.02), ct_up_time[1]-ct_up_time[0], 0.02-(-0.02), edgecolor=None, facecolor='r', alpha=0.05)
                ax=plt.gca()
                ax.add_patch(square)
        # DOWN times square patches
        if len(ct_down_times_draw_window)>0:
            for ct_down_time in ct_down_times_draw_window:        
                square = patches.Rectangle((ct_down_time[0], -0.02), ct_down_time[1]-ct_down_time[0], 0.02-(-0.02), edgecolor=None, facecolor='mediumpurple', alpha=0.05)
                ax=plt.gca()
                ax.add_patch(square)

        # TRN
        cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
        if   'TRN' in plot_pop: line_color = 'b'
        elif 'VPM' in plot_pop: line_color = 'g'
        else:                   line_color = 'k'
        # line_color = sim_color[sim_index]
        time_signal = store_lfpTraces_dict[sim_name][recording_frequency]['t']
        lfp_signal  = store_lfpTraces_dict[sim_name][recording_frequency]['data'][pop_electrodes[plot_pop]] # Replace with your LFP signal list
        plt.plot(time_signal,lfp_signal,color=line_color, linewidth=1.0,alpha = 0.75)
        # plt.plot(time_signal,lfp_signal,color=line_color, linewidth=sim_linewidth)
        plt.ylim([-0.02,0.02])
        plt.yticks([])

        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        plt.plot([1420,1420],[-0.01,-0.0075],'k')
        plt.plot([1420,1920],[-0.01,-0.01],'k')
        plt.text(1380,-0.014,'500 ms')
        plt.text(1150,-0.014,'2.5 µV',rotation=90)

        plt.xlim([1000,16000])
        plt.xticks([])
        plt.xlabel('')

        # Example usage    
        filtered_lfp = bandpass_filter(lfp_signal, fs, lowcut, highcut)

        # Compute envelope
        lfp_envelope = compute_envelope(filtered_lfp)

        # plt.savefig(lfp_traces_data_path+'lfpTraces_'+sim_name+'_'+plot_pop+'.png',dpi=200)
        plt.savefig(save_figures_path+'lfpTraces_'+sim_name+'_'+plot_pop+'.png',dpi=200)

        #######################################################################################################################################

        # slicing LFP indexes based on UP and DOWN times
        if len(ct_up_times)>0:      ct_up_indexes   = [[i for i, t in enumerate(time_signal) if start <= t <= end] for start, end in ct_up_times]
        if len(ct_down_times)>0:    ct_down_indexes = [[i for i, t in enumerate(time_signal) if start <= t <= end] for start, end in ct_down_times]

        # slicing LFP data based on UP and DOWN times
        time_signal_up   = [[time_signal[i] for i in indexes] for indexes in ct_up_indexes]
        lfp_signal_up    = [[lfp_signal[i]  for i in indexes] for indexes in ct_up_indexes]
        time_signal_down = [[time_signal[i] for i in indexes] for indexes in ct_down_indexes]
        lfp_signal_down  = [[lfp_signal[i]  for i in indexes] for indexes in ct_down_indexes]

        lfp_filtered_up   = [[filtered_lfp[i] for i in indexes] for indexes in ct_up_indexes]
        lfp_filtered_down = [[filtered_lfp[i] for i in indexes] for indexes in ct_down_indexes]

        lfp_envelope_up   = [[lfp_envelope[i] for i in indexes] for indexes in ct_up_indexes]
        lfp_envelope_down = [[lfp_envelope[i] for i in indexes] for indexes in ct_down_indexes]

        plt.figure(figsize=[25,3])
        # UP times square patches
        if len(ct_up_times_draw_window)>0:
            for ct_up_time in ct_up_times_draw_window:        
                square = patches.Rectangle((ct_up_time[0], -0.02), ct_up_time[1]-ct_up_time[0], 0.02-(-0.02), edgecolor=None, facecolor='r', alpha=0.05)
                ax=plt.gca()
                ax.add_patch(square)
        # DOWN times square patches
        if len(ct_down_times_draw_window)>0:
            for ct_down_time in ct_down_times_draw_window:        
                square = patches.Rectangle((ct_down_time[0], -0.02), ct_down_time[1]-ct_down_time[0], 0.02-(-0.02), edgecolor=None, facecolor='mediumpurple', alpha=0.05)
                ax=plt.gca()
                ax.add_patch(square)

        # TRN
        cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
        if   'TRN' in plot_pop: line_color = 'b'
        elif 'VPM' in plot_pop: line_color = 'g'
        else:                   line_color = 'k'
        # line_color = sim_color[sim_index]
        time_signal = store_lfpTraces_dict[sim_name][recording_frequency]['t']
        lfp_signal  = store_lfpTraces_dict[sim_name][recording_frequency]['data'][pop_electrodes[plot_pop]] # Replace with your LFP signal list
        
        plt.plot(time_signal,filtered_lfp,color=line_color, linewidth=1.0,alpha = 0.75)
        # plt.plot(time_signal,lfp_signal,color=line_color, linewidth=1.0,alpha = 0.5)

        # # plot colored traces on top of LFP signal
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_up[idx],    lfp_signal_up[idx],      '--',color='firebrick', linewidth=sim_linewidth)
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_down[idx],  lfp_signal_down[idx],    '--',color='firebrick', linewidth=sim_linewidth)

        # # plot colored traces on top of LFP signal
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_up[idx],    lfp_filtered_up[idx],      '--',color='firebrick', linewidth=sim_linewidth)
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_down[idx],  lfp_filtered_down[idx],    '--',color='firebrick', linewidth=sim_linewidth)

        # # plot colored traces on top of LFP signal
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_up[idx],    lfp_envelope_up[idx],      '-',color='firebrick', linewidth=sim_linewidth)
        # for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_down[idx],  lfp_envelope_down[idx],    '-',color='lightcoral', linewidth=sim_linewidth)
        # plot colored traces on top of LFP signal
        for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_up[idx],    lfp_envelope_up[idx],      '-',color='r',              linewidth=1.5)
        for idx,tvals in enumerate(time_signal_up): plt.plot(time_signal_down[idx],  lfp_envelope_down[idx],    '-',color='mediumpurple',   linewidth=1.5)


        store_lfp_mean_up_envelope=[]
        store_lfp_mean_down_envelope=[]

        for idx,tvals in enumerate(time_signal_up): 
            # mean LFP - UP
            mean_lfp_timerange   = [time_signal_up[idx][0],time_signal_up[idx][-1]]
            store_lfp_mean_up_envelope.append(np.mean(lfp_envelope_up[idx]))
            
            # # plt.plot(rms_lfp_timerange, [rms_lfp,rms_lfp],                   '--',  color='firebrick', linewidth=0.5)
            # # plt.plot(rms_lfp_timerange, [rms_lfp_filtered,rms_lfp_filtered], '-',   color='firebrick',  linewidth=1.0)
            # # plt.plot(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='k',  linewidth=1.0)
            # plt.fill_between(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='k',  linewidth=1.0, alpha=0.2, edgecolor='none')
            plt.fill_between(mean_lfp_timerange, [np.mean(lfp_envelope_up[idx]),np.mean(lfp_envelope_up[idx])], '-',   color='k',  linewidth=1.0, alpha=0.3, edgecolor='none')

        for idx,tvals in enumerate(time_signal_down): 
            # mean LFP - DOWN
            mean_lfp_timerange   = [time_signal_down[idx][0],time_signal_down[idx][-1]]
            store_lfp_mean_down_envelope.append(np.mean(lfp_envelope_down[idx]))
            
            # # plt.plot(rms_lfp_timerange, [rms_lfp,rms_lfp],                   '--',  color='lightcoral', linewidth=0.5)
            # # plt.plot(rms_lfp_timerange, [rms_lfp_filtered,rms_lfp_filtered], '-',   color='lightcoral',  linewidth=1.0)
            # # plt.plot(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='grey',  linewidth=1.0)
            # plt.fill_between(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='grey',  linewidth=1.0, alpha=0.2, edgecolor='none')
            plt.fill_between(mean_lfp_timerange, [np.mean(lfp_envelope_down[idx]),np.mean(lfp_envelope_down[idx])], '-',   color='k',  linewidth=1.0, alpha=0.3, edgecolor='none')

        # plot LFP RMS value
        store_lfp_rms_up=[]
        store_lfp_rms_down=[]

        store_lfp_rms_up_raw=[]
        store_lfp_rms_up_filtered=[]
        store_lfp_rms_up_envelope=[]
        store_lfp_rms_up_ratio=[]
        store_lfp_rms_up_ratio2=[]
        
        store_lfp_rms_down_raw=[]
        store_lfp_rms_down_filtered=[]
        store_lfp_rms_down_envelope=[]
        store_lfp_rms_down_ratio=[]
        store_lfp_rms_down_ratio2=[]

        for idx,tvals in enumerate(time_signal_up): 
            # LFP RMS
            rms_lfp             = math.sqrt(sum(x**2 for x in lfp_signal_up[idx])   / len(lfp_signal_up[idx]))
            rms_lfp_filtered    = math.sqrt(sum(x**2 for x in lfp_filtered_up[idx]) / len(lfp_filtered_up[idx]))
            rms_lfp_envelope    = math.sqrt(sum(x**2 for x in lfp_envelope_up[idx]) / len(lfp_envelope_up[idx]))
            rms_lfp_timerange   = [time_signal_up[idx][0],time_signal_up[idx][-1]]
            
            store_lfp_rms_up.append([rms_lfp_filtered,rms_lfp])
            store_lfp_rms_up_raw.append(rms_lfp)
            store_lfp_rms_up_filtered.append(rms_lfp_filtered)
            store_lfp_rms_up_ratio.append(rms_lfp_filtered/rms_lfp)
            store_lfp_rms_up_ratio2.append(rms_lfp_filtered/(rms_lfp-rms_lfp_filtered))
            #
            store_lfp_rms_up_envelope.append(rms_lfp_envelope)
            
            # plt.plot(rms_lfp_timerange, [rms_lfp,rms_lfp],                   '--',  color='firebrick', linewidth=0.5)
            # plt.plot(rms_lfp_timerange, [rms_lfp_filtered,rms_lfp_filtered], '-',   color='firebrick',  linewidth=1.0)
            # plt.plot(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='k',  linewidth=1.0)
            # plt.fill_between(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='k',  linewidth=1.0, alpha=0.2, edgecolor='none')
        
        for idx,tvals in enumerate(time_signal_down): 
            rms_lfp             = math.sqrt(sum(x**2 for x in lfp_signal_down[idx])   / len(lfp_signal_down[idx]))
            rms_lfp_filtered    = math.sqrt(sum(x**2 for x in lfp_filtered_down[idx]) / len(lfp_filtered_down[idx]))
            rms_lfp_envelope    = math.sqrt(sum(x**2 for x in lfp_envelope_down[idx]) / len(lfp_envelope_down[idx]))
            rms_lfp_timerange   = [time_signal_down[idx][0],time_signal_down[idx][-1]]
            
            store_lfp_rms_down.append([rms_lfp_filtered,rms_lfp])
            store_lfp_rms_down_raw.append(rms_lfp)
            store_lfp_rms_down_filtered.append(rms_lfp_filtered)
            store_lfp_rms_down_ratio.append(rms_lfp_filtered/rms_lfp)
            store_lfp_rms_down_ratio2.append(rms_lfp_filtered/(rms_lfp-rms_lfp_filtered))
            #
            store_lfp_rms_down_envelope.append(rms_lfp_envelope)
            
            # plt.plot(rms_lfp_timerange, [rms_lfp,rms_lfp],                   '--',  color='lightcoral', linewidth=0.5)
            # plt.plot(rms_lfp_timerange, [rms_lfp_filtered,rms_lfp_filtered], '-',   color='lightcoral',  linewidth=1.0)
            # plt.plot(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='grey',  linewidth=1.0)
            # plt.fill_between(rms_lfp_timerange, [rms_lfp_envelope,rms_lfp_envelope], '-',   color='grey',  linewidth=1.0, alpha=0.2, edgecolor='none')

        plt.ylim([-0.01,0.01])
        plt.yticks([])

        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        plt.plot([1420,1420],[-0.005,-0.0025],'k')
        plt.plot([1420,1920],[-0.005,-0.005],'k')
        plt.text(1380,-0.007,'500 ms')
        plt.text(1150,-0.007,'2.5 µV',rotation=90)

        plt.xlim([1000,16000])
        plt.xticks([])
        plt.xlabel('')

        # plt.savefig(lfp_traces_data_path+'lfpTraces_'+sim_name+'_'+plot_pop+'__upDown.png',dpi=200)
        plt.savefig(save_figures_path+'lfpTraces_'+sim_name+'_'+plot_pop+'__upDown.png',dpi=200)

        mean_envelope_up   = np.mean([item for sublist in lfp_envelope_up   for item in sublist])
        mean_envelope_down = np.mean([item for sublist in lfp_envelope_down for item in sublist])

        mean_envelope_dict[sim_name].update({plot_pop:{'mean_env_up':mean_envelope_up,'mean_env_down':mean_envelope_down}})

        # mean with stats
        std_envelope_up   = np.std([item for sublist in lfp_envelope_up   for item in sublist])
        std_envelope_down = np.std([item for sublist in lfp_envelope_down for item in sublist])
        mean_envelope_dict[sim_name][plot_pop].update({'std_env_up':std_envelope_up,'std_env_down':std_envelope_down})

        # plt.figure(figsize=(5,5))
        # plt.bar([1,2],[mean_envelope_down,mean_envelope_up],width=0.6,color=['b','r'])
        # plt.xticks([1,2],['Down', 'Up'])
        # # plt.savefig(lfp_traces_data_path+'lfpTraces_'+sim_name+'_'+plot_pop+'__upDown_bar.png',dpi=200)
        # plt.savefig(save_figures_path+'lfpTraces_'+sim_name+'_'+plot_pop+'__upDown_bar.png',dpi=200)
    
        store_lfp_rms_dict[sim_name][plot_pop].update({
            'up':{
                'vals':     store_lfp_rms_up,
                'raw':      store_lfp_rms_up_raw,
                'filtered': store_lfp_rms_up_filtered,
                'envelope': store_lfp_rms_up_envelope,
                'ratio':    store_lfp_rms_up_ratio,
                'ratio2':   store_lfp_rms_up_ratio2,
                },
            'down':{
                'vals':     store_lfp_rms_down,
                'raw':      store_lfp_rms_down_raw,
                'filtered': store_lfp_rms_down_filtered,
                'envelope': store_lfp_rms_down_envelope,
                'ratio':    store_lfp_rms_down_ratio,
                'ratio2':   store_lfp_rms_down_ratio2,
                },
                })
        store_lfp_mean_dict[sim_name][plot_pop].update({
            'up':{
                'envelope': store_lfp_mean_up_envelope,
                },
            'down':{
                'envelope': store_lfp_mean_down_envelope,
                },
                })

plt.figure(figsize=(13,9))
plt.suptitle('Filtered LFP envelope RMS power')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    for sim_index, sim_name in enumerate(store_lfp_rms_dict.keys()):    
        # if sim_index==2 or sim_index==6:
        #     alpha = 0.8
        #     down_color = 'mediumpurple'
        #     up_color = 'r'
        #     down_edgecolor=None
        #     up_edgecolor=None
        # else:
        #     alpha = 0.4
        #     down_color = 'w'
        #     up_color = 'w'
        #     down_edgecolor='b'
        #     up_edgecolor='r'
        
        if sim_index==4:
            alpha = 0.8
            down_color = 'purple'
            up_color = 'darkred'
            down_edgecolor=None
            up_edgecolor=None
            medianline_color='w'
        else:
            alpha = 0.5
            down_color = 'mediumpurple'
            up_color = 'r'
            down_edgecolor=None
            up_edgecolor=None
            medianline_color='k'
            
        # plt.bar(     sim_index-0.1,np.mean(store_lfp_rms_dict[sim_name][plot_pop]['down']['envelope']), width=0.2,                                                                color=down_color, alpha=alpha, edgecolor=down_edgecolor)
        # plt.errorbar(sim_index-0.1,np.mean(store_lfp_rms_dict[sim_name][plot_pop]['down']['envelope']), yerr=np.std(store_lfp_rms_dict[sim_name][plot_pop]['down']['envelope']),  color='k', capsize=3)
        
        # plt.bar(     sim_index+0.1,np.mean(store_lfp_rms_dict[sim_name][plot_pop]['up']['envelope']),   width=0.2,                                                                color=up_color,   alpha=alpha, edgecolor=up_edgecolor)
        # plt.errorbar(sim_index+0.1,np.mean(store_lfp_rms_dict[sim_name][plot_pop]['up']['envelope']),   yerr=np.std(store_lfp_rms_dict[sim_name][plot_pop]['up']['envelope']),    color='k', capsize=3)

        # plt.boxplot(store_lfp_rms_dict[sim_name][plot_pop]['down']['envelope'], positions=np.array([sim_index-0.3]))
        # plt.boxplot(store_lfp_rms_dict[sim_name][plot_pop]['up']['envelope'],   positions=np.array([sim_index+0.3]))
        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.1
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.1
                colors = [down_color]
                state_label='D'

            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_rms_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label],
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([sim_index+box_shift]))
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(alpha)
            
            plt.ylim([0,9e-3])
            plt.yticks([0,2e-3,4e-3,6e-3])
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')


# for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
#     print('\n\n Comparing sim differences to uniform connectivity')
#     print(f"Plot pop: {plot_pop}")

#     for state in ['up','down']:
#         print(f' - state: {state}')
#         for sim_combination in [['sim_4','sim_0'],
#                                 ['sim_4','sim_1'],
#                                 ['sim_4','sim_2'],
#                                 ['sim_4','sim_3'],
#                                 ['sim_4','sim_5'],
#                                 ['sim_4','sim_6'],
#                                 ['sim_4','sim_7'],
#                                 ['sim_4','sim_8'],
#                                 ]:
#             print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}]")

                
#             groupA = store_lfp_rms_dict[sim_combination[0]][plot_pop][state]['envelope']  # Replace with the first group of interest
#             groupB = store_lfp_rms_dict[sim_combination[1]][plot_pop][state]['envelope']  # Replace with the second group of interest

#             U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
#             # if p_mw<0.05:   significance = '*'
#             # else:           significance = 'ns'
#             significance = significance_symbol(p_mw)
#             print(f"\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")

#     for sim_combination in [['sim_2','sim_6'],
#                             ]:
#         print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}]")
#         for state in ['up','down']:
                
#             groupA = store_lfp_rms_dict[sim_combination[0]][plot_pop][state]['envelope']  # Replace with the first group of interest
#             groupB = store_lfp_rms_dict[sim_combination[1]][plot_pop][state]['envelope']  # Replace with the second group of interest

#             U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
#             # if p_mw<0.05:   significance = '*'
#             # else:           significance = 'ns'
#             significance = significance_symbol(p_mw)
#             print(f"\t\tstate: {state} \t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")


#     print(f"Comparing states: up vs down")
#     for sim_name in ['sim_0','sim_1','sim_2','sim_3','sim_4','sim_5','sim_6','sim_7','sim_8']:
#         print(f"\tsim name: [{sim_name}]")
#         groupA = store_lfp_rms_dict[sim_name][plot_pop]['down']['envelope']  # Replace with the first group of interest
#         groupB = store_lfp_rms_dict[sim_name][plot_pop]['up']['envelope']  # Replace with the second group of interest

#         U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
#         # if p_mw<0.05:   significance = '*'
#         # else:           significance = 'ns'
#         significance = significance_symbol(p_mw)
#         print(f" Plot pop: {plot_pop} - sim name: [{sim_name}] \t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")


select_sims = ['sim_2', 'sim_6']
plt.figure(figsize=(5,8))
# plt.suptitle('Filtered LFP envelope RMS power - sim_2 and sim_6')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    sim_shift=0
    plt.title(plot_pop.split('__')[0])
    for sim_index, sim_name in enumerate(store_lfp_rms_dict.keys()):    

        if select_sims is not None:
            if sim_name not in select_sims:continue

        if sim_index==2 or sim_index==6:    
            alpha = 0.8
            down_color = 'mediumpurple'
            up_color = 'r'
            down_edgecolor=None
            up_edgecolor=None
        else:                               
            alpha = 0.4
            down_color = 'w'
            up_color = 'w'
            down_edgecolor='b'
            up_edgecolor='r'

        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.15
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.15
                colors = [down_color]
                state_label='D'

            if select_sims is not None:
                x_tick = sim_shift+box_shift
            else:
                x_tick = sim_index+box_shift


            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_rms_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label], 
                            #    notch=True,
                               medianprops=dict(color="k", alpha=0.7),
                               patch_artist=True, 
                               positions=np.array([x_tick]),
                               widths=0.15,
                               )
            
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(0.5)
        
        sim_shift+=1

        plt.plot([0.15,1.15],[0.0060,0.0060],'k')
        plt.text(0.7,0.0061,'*')

        plt.plot([-0.15,0.15],[0.0052,0.0052],'k')
        plt.text(0.0,0.0053,'*')

        plt.plot([0.85,1.15],[0.0052,0.0052],'k')
        plt.text(1.0,0.0053,'*')
        
        plt.ylim([0,6.5e-3])
        plt.yticks([0,2e-3,4e-3,6e-3])
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.ylabel('Envelope RMS Power')
        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        if select_sims is not None: plt.xticks([0,1],[-0.5,0.5])
        else:                       plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')
plt.tight_layout()
# plt.savefig(lfp_traces_data_path+'lfpEnvelopeRMS_1.png',dpi=200)
plt.savefig(save_figures_path+'lfpEnvelopeRMS_1.png',dpi=200)
# plt.show()

select_sims = ['sim_2', 'sim_4', 'sim_6']
plt.figure(figsize=(5,9))
plt.suptitle('Filtered LFP envelope RMS power - sim_2, sim_4 and sim_6')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    sim_shift=0
    # plt.title(plot_pop)
    for sim_index, sim_name in enumerate(store_lfp_rms_dict.keys()):    

        if select_sims is not None:
            if sim_name not in select_sims:continue

        alpha = 0.8
        down_color = 'mediumpurple'
        up_color = 'r'
        down_edgecolor=None
        up_edgecolor=None

        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.2
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.2
                colors = [down_color]
                state_label='D'

            if select_sims is not None:
                x_tick = sim_shift+box_shift
            else:
                x_tick = sim_index+box_shift


            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_rms_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label], 
                            #    notch=True,
                               medianprops=dict(color="k", alpha=0.7),
                               patch_artist=True, 
                               positions=np.array([x_tick]),
                               widths=0.15,
                               )
            
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(0.5)
        
        sim_shift+=1

        if plot_pop == 'TRN__pop':
            # sim_2 - sim_4
            plt.plot([0.2,2.2],[0.0058,0.0058],'k')
            plt.text(1.2,0.0059,'*')

        if plot_pop == 'VPM__pop':
            # sim_2 - sim_4
            plt.plot([0.2,2.2],[0.0058,0.0058],'k')
            plt.text(1.2,0.0059,'*')

            # sim_4 - sim_6
            plt.plot([1.2,2.2],[0.0052,0.0052],'k')
            plt.text(1.7,0.0053,'*')

        plt.ylim([0,6.5e-3])

        if select_sims is not None: plt.xticks([0,1,2],[-0.5,0,0.5])
        else:                       plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')
plt.tight_layout()
# plt.savefig(lfp_traces_data_path+'lfpEnvelopeRMS_2.png',dpi=200)
plt.savefig(save_figures_path+'lfpEnvelopeRMS_2.png',dpi=200)
# plt.show()

#######################################################################################################################################

# Filtered LFP envelope Mean power
plt.figure(figsize=(13,9))
plt.suptitle('Filtered LFP envelope Mean power')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    for sim_index, sim_name in enumerate(store_lfp_mean_dict.keys()):    
        if sim_index==4:
            alpha = 0.8
            down_color = 'purple'
            up_color = 'darkred'
            down_edgecolor=None
            up_edgecolor=None
            medianline_color='w'
        else:
            alpha = 0.5
            down_color = 'mediumpurple'
            up_color = 'r'
            down_edgecolor=None
            up_edgecolor=None
            medianline_color='k'
            
        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.1
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.1
                colors = [down_color]
                state_label='D'

            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_mean_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label],
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([sim_index+box_shift]))
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(alpha)
            
            plt.ylim([0,9e-3])
            plt.yticks([0,2e-3,4e-3,6e-3])
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax1=plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            
    plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
# plt.savefig(lfp_traces_data_path+'lfpEnvelopeMean_0.png',dpi=200)
plt.savefig(save_figures_path+'lfpEnvelopeMean_0.png',dpi=200)

# Filtered LFP envelope Mean power
plt.figure(figsize=(3.5,11))
plt.suptitle('Filtered LFP envelope Mean power')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    store_xticks=[]
    for sim_index, sim_name in enumerate(['sim_4']):    
        alpha = 0.5
        down_color = 'mediumpurple'
        up_color = 'r'
        down_edgecolor=None
        up_edgecolor=None
            
        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.1
                colors = [up_color]
                state_label='U'
                medianline_color='r'
            else:               
                box_shift = -0.1
                colors = [down_color]
                state_label='D'
                medianline_color='mediumpurple'

            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_mean_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label],
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([sim_index+box_shift]))
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                # box.set_facecolor(color)
                box.set_facecolor('w')
                box.set_alpha(alpha)
            
            store_xticks.append(sim_index+box_shift)
            
            plt.ylim([0,6e-3])
            plt.yticks([0,2e-3,4e-3,6e-3])
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax1=plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)

        # plt.xticks([0],['0.0'])
        plt.xticks(store_xticks,['up','down'])
plt.xlabel('$P_{T}$=0.0')
# plt.savefig(lfp_traces_data_path+'lfpEnvelopeMean_1.png',dpi=200)
plt.savefig(save_figures_path+'lfpEnvelopeMean_1.png',dpi=200)



for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    print('\n\n Comparing sim differences to uniform connectivity')
    print(f"Plot pop: {plot_pop}")
    
    test_altenative = 'greater'
    test_altenative = 'less'
    # test_altenative = 'two-sided'

    # for sim_combination in [['sim_4','sim_0','less'],
    #                         ['sim_4','sim_1','less'],
    #                         ['sim_4','sim_2','less'],
    #                         ['sim_4','sim_3','less'],
    #                         ['sim_4','sim_5','greater'],
    #                         ['sim_4','sim_6','greater'],
    #                         ['sim_4','sim_7','greater'],
    #                         ['sim_4','sim_8','greater'],
    #                         ]:

    for sim_combination in [['sim_4','sim_0','two-sided'],
                            ['sim_4','sim_1','two-sided'],
                            ['sim_4','sim_2','two-sided'],
                            ['sim_4','sim_3','two-sided'],
                            ['sim_4','sim_5','two-sided'],
                            ['sim_4','sim_6','two-sided'],
                            ['sim_4','sim_7','two-sided'],
                            ['sim_4','sim_8','two-sided'],
                            ]:
        print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}, {sim_combination[2]}]")
        for state in ['up','down']:
                
            groupA = store_lfp_mean_dict[sim_combination[0]][plot_pop][state]['envelope']  # Replace with the first group of interest
            groupB = store_lfp_mean_dict[sim_combination[1]][plot_pop][state]['envelope']  # Replace with the second group of interest

            # U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
            U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
            # if p_mw<0.05:   significance = '*'
            # else:           significance = 'ns'
            significance = significance_symbol(p_mw)
            print(f"\t\tstate: {state} \t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")

    print('\n')
    for sim_combination in [
                            ['sim_3','sim_5','two-sided'],
                            ['sim_2','sim_6','two-sided'],
                            ['sim_1','sim_7','two-sided'],
                            ['sim_0','sim_8','two-sided'],
                            ]:
        print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}]")
        for state in ['up','down']:
                
            groupA = store_lfp_mean_dict[sim_combination[0]][plot_pop][state]['envelope']  # Replace with the first group of interest
            groupB = store_lfp_mean_dict[sim_combination[1]][plot_pop][state]['envelope']  # Replace with the second group of interest

            U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
            # if p_mw<0.05:   significance = '*'
            # else:           significance = 'ns'
            significance = significance_symbol(p_mw)
            print(f"\t\tstate: {state} \t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")


    print(f"Comparing states: up vs down")
    for sim_name in ['sim_0','sim_1','sim_2','sim_3','sim_4','sim_5','sim_6','sim_7','sim_8']:
        print(f"\tsim name: [{sim_name}]")
        groupA = store_lfp_mean_dict[sim_name][plot_pop]['down']['envelope']  # Replace with the first group of interest
        groupB = store_lfp_mean_dict[sim_name][plot_pop]['up']['envelope']  # Replace with the second group of interest

        U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative='two-sided')
        # if p_mw<0.05:   significance = '*'
        # else:           significance = 'ns'
        significance = significance_symbol(p_mw)
        print(f" Plot pop: {plot_pop} - sim name: [{sim_name}] \t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}")



plt.figure(figsize=(6,9))
plt.suptitle('Filtered LFP envelope Mean power')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    for state in ['up','down']:
        if state=='up': color='r'
        else:           color='mediumpurple'

        medians = [np.median(    store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'])     for sim_name in store_lfp_mean_dict.keys()]
        q25     = [np.percentile(store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'], 25) for sim_name in store_lfp_mean_dict.keys()]
        q75     = [np.percentile(store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'], 75) for sim_name in store_lfp_mean_dict.keys()]
        conds   = list(store_lfp_mean_dict.keys())
        
        plt.plot(conds, medians, marker='o', color=color,   label='Median')  # Line for medians
        plt.fill_between(conds,  q25, q75,   color=color,   alpha=0.3, label='IQR (25th-75th)')  # Shaded IQR

    plt.title(plot_pop)
    plt.ylim([0,6e-3])
    plt.yticks([0,2e-3,4e-3,6e-3])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    # plt.xticks([])
    plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')

# --- Main analysis figure for paper
plt.figure(figsize=(10,11))
# plt.suptitle('Filtered LFP envelope Mean power')
for pop_ind,plot_pop in enumerate(['TRN__pop', 'VPM__pop']):
    plt.subplot(2,1,pop_ind+1)
    for state in ['up','down']:

        if state == 'up':   
            box_shift = 0.1
            colors = ['none']
            state_label='U'
            medianline_color='darkred'
        else:               
            box_shift = -0.1
            colors = ['none']
            state_label='D'
            medianline_color='purple'

        if state=='up': 
            color_='r'
            box_shift = 0.1
        else:           
            color_='mediumpurple'
            box_shift = -0.1

        medians = [np.median(    store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'])     for sim_name in store_lfp_mean_dict.keys()]
        q25     = [np.percentile(store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'], 25) for sim_name in store_lfp_mean_dict.keys()]
        q75     = [np.percentile(store_lfp_mean_dict[sim_name][plot_pop][state]['envelope'], 75) for sim_name in store_lfp_mean_dict.keys()]
        conds   = [x_point+box_shift for x_point in range(len(store_lfp_mean_dict.keys()))]
        
        plt.plot(conds, medians, marker='.', color=color_,   label='Median')  # Line for medians
        plt.fill_between(conds,  q25, q75,   color=color_,   alpha=0.3, label='IQR (25th-75th)')  # Shaded IQR

        for sim_index, sim_name in enumerate(store_lfp_mean_dict.keys()):
            ax = plt.gca()
            boxes = ax.boxplot([store_lfp_mean_dict[sim_name][plot_pop][state]['envelope']], 
                            #    labels=[state_label],
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([sim_index+box_shift]))
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(alpha)

    # plt.title(plot_pop.split('__')[0])
    plt.ylim([0,10e-3])
    plt.yticks([0,2.5e-3,5e-3,7.5e-3])
    # plt.yticks([0,2e-3,4e-3,6e-3])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    # plt.xticks([])
    plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('LFP envelope mean voltage (μV)')
# plt.savefig(lfp_traces_data_path+'lfpEnvelopeMean_2_linePlot.png',dpi=200)
plt.savefig(save_figures_path+'lfpEnvelopeMean_2_linePlot.png',dpi=200)
# plt.show()



# plt.show()
# sys.exit()

# plt.figure(figsize=(10,5))
# for sim_index, sim_name in enumerate(mean_envelope_dict.keys()):
#     for pop_ind,plot_pop in enumerate(mean_envelope_dict[sim_name].keys()):
#         plt.subplot(2,1,pop_ind+1)
#         plt.title(plot_pop)
#         plt.bar([sim_index-0.2,sim_index+0.2],[mean_envelope_dict[sim_name][plot_pop]['mean_env_down'],mean_envelope_dict[sim_name][plot_pop]['mean_env_up']],width=0.4,color=['b','r'])
#         plt.errorbar(sim_index-0.2,mean_envelope_dict[sim_name][plot_pop]['mean_env_down'],
#                      yerr=mean_envelope_dict[sim_name][plot_pop]['std_env_down'],
#                      color='k')
#         plt.errorbar(sim_index+0.2,mean_envelope_dict[sim_name][plot_pop]['mean_env_up'],
#                      yerr=mean_envelope_dict[sim_name][plot_pop]['std_env_up'],
#                      color='k')
#         plt.xticks([],[])
# plt.xticks(list(range(len(mean_envelope_dict.keys()))),list(mean_envelope_dict.keys()))
# plt.tight_layout()
# # plt.savefig(lfp_traces_data_path+'lfpTraces__upDown_bar.png',dpi=200)
# plt.savefig(save_figures_path+'lfpTraces__upDown_bar.png',dpi=200)


#######################################################################################################################################
spectrogram_windows_data_path   = simulation_data_path+'spectrogram/'
data_folder_spectrogram_windows = spectrogram_windows_data_path+'data_windows/'

list_files_spectrogram_windows = os.listdir(data_folder_spectrogram_windows)
# list_files_spectrogram_windows = [file for file in list_files_spectrogram_windows if '_spectrogram_values_windows_upStates.pkl' in file]
list_files_spectrogram_windows = [file for file in list_files_spectrogram_windows if '_spectrogram_values_windows.pkl' in file]
list_files_spectrogram_windows.sort()

store_spectrogram_windows_dict={}
for data_file_spectrogram_windows_ind, data_file_spectrogram_windows in enumerate(list_files_spectrogram_windows):
    filename_pkl = data_folder_spectrogram_windows+data_file_spectrogram_windows
    with open(filename_pkl, 'rb') as file:spectrogram_windows_dict = pickle.load(file)
    store_spectrogram_windows_dict.update({'sim_'+str(data_file_spectrogram_windows_ind):spectrogram_windows_dict})

lowcut = 8  # Lower bound of frequency band (e.g., theta band 4-8 Hz)
highcut = 16  # Upper bound of frequency band

#######################################################################################################################################

# select_sims = ['sim_2']
select_sims = ['sim_4']
select_sims = None

store_rms_dict={}
store_mean_dict={}

skip_time_up = 0
sampling_rate=0.1

store_band_power_dict={}
for sim_index, sim_name in enumerate(store_spectrogram_windows_dict.keys()):
    if select_sims is not None:
        if sim_name not in select_sims:continue

    store_rms_dict.update({sim_name:{}})
    store_mean_dict.update({sim_name:{}})
    store_band_power_dict.update({sim_name:{}})

    for plot_pop in store_spectrogram_windows_dict[sim_name].keys():
        store_rms_dict[sim_name].update({plot_pop:{'up':[],'down':[]}})
        store_mean_dict[sim_name].update({plot_pop:{'up':[],'down':[]}})
        store_band_power_dict[sim_name].update({plot_pop:{}})
        for time_window in store_spectrogram_windows_dict[sim_name][plot_pop].keys():

            time_range = store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['extent'][0][0:2] # 1st and 2nd vals
            freq_range = store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['extent'][0][2:4] # 3rd and 4th vals

            # time
            num_vals_x_axis = len(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][0])
            # frequencies
            num_vals_y_axis = len(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0])

            time_array = np.arange(time_range[0],time_range[1],  ((time_range[1]-time_range[0])/(num_vals_x_axis)))
            freq_array = np.arange(freq_range[0],freq_range[1]+1,((freq_range[1]-freq_range[0])/(num_vals_y_axis-1)))

            store_freqs=[]
            store_freq_bands = []
            for freq_index,freq_val in enumerate(freq_array):
                if (freq_val>=lowcut) and (freq_val<=highcut):
                    store_freqs.append(freq_val)
                    store_freq_bands.append(list(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][freq_index]))

            # transposes the list of lists to be organized from list of power by frequency (where each list is a frequency) to list of powers for each timestep (where each list is a timestep)
            store_freq_bands_transpose = list(map(list, zip(*store_freq_bands)))

            if 'u_' in time_window: 
                skip_samples = int(skip_time_up/sampling_rate)
            elif 'd_' in time_window: 
                skip_samples = 0
            else:
                skip_samples = 0
            # print('skipping samples for time window ', time_window, ' | samples: ', skip_samples)

            store_freq_bands_transpose_sum = [sum(band_powers) for band_powers_ind,band_powers in enumerate(store_freq_bands_transpose) if band_powers_ind>skip_samples]

            store_freq_bands_transpose_sum_mean = np.mean(store_freq_bands_transpose_sum)

            store_freq_bands_transpose_sum_rms = math.sqrt(sum(x**2 for x in store_freq_bands_transpose_sum) / len(store_freq_bands_transpose_sum))
            
            # plt.figure()
            # plt.plot(list(range(len(store_freq_bands_transpose_sum))),store_freq_bands_transpose_sum,'-k')
            # plt.plot([0,len(store_freq_bands_transpose)],[store_freq_bands_transpose_sum_rms,store_freq_bands_transpose_sum_rms],'--',color='firebrick',)
            # plt.plot([0,len(store_freq_bands_transpose)],[store_freq_bands_transpose_sum_mean,store_freq_bands_transpose_sum_mean],'--',color='lightcoral',)
            # plt.ylim([0,1.5e-4])
            # plt.title('Event name - '+ time_window)

            if   'u_' in time_window:   store_rms_dict[sim_name][plot_pop]['up'].append(  store_freq_bands_transpose_sum_rms)
            elif 'd_' in time_window:   store_rms_dict[sim_name][plot_pop]['down'].append(store_freq_bands_transpose_sum_rms)
            else:                       continue

            if   'u_' in time_window:   store_mean_dict[sim_name][plot_pop]['up'].append(  store_freq_bands_transpose_sum_mean)
            elif 'd_' in time_window:   store_mean_dict[sim_name][plot_pop]['down'].append(store_freq_bands_transpose_sum_mean)
            else:                       continue

cuttof_freq = 50
select_spectrograms=['sim_2','sim_4','sim_6']
for sim_index, sim_name in enumerate(store_spectrogram_windows_dict.keys()):
    if select_spectrograms is not None:
        if sim_name not in select_spectrograms:continue
    plt.figure(figsize=(22,9))
    plt.suptitle(f'Full Spectrogram - {sim_name}')
    for pop_ind,plot_pop in enumerate(['TRN__pop','VPM__pop']):

        time_windows_up   = [time_window for time_window in store_spectrogram_windows_dict[sim_name][plot_pop].keys() if 'u_' in time_window]
        time_windows_down = [time_window for time_window in store_spectrogram_windows_dict[sim_name][plot_pop].keys() if 'd_' in time_window]

        time_windows_up_samples = [(int((int(time_window_up.split('__')[1])-1500)/0.1),int((int(time_window_up.split('__')[2])-1500)/0.1)) for time_window_up in time_windows_up]

        plt.subplot(2,1,pop_ind+1)
        plt.title(plot_pop)
        dataset = store_spectrogram_windows_dict[sim_name][plot_pop]['a__1500__16000']['electrodes']['spectrogram']['morlet'][0][0:cuttof_freq]
        # Convert list of lists to a NumPy array
        data_array = np.array(dataset)

        # Display the data using imshow
        # plt.imshow(data_array, aspect='auto', cmap='viridis')
        plt.imshow(data_array, aspect='auto', cmap='viridis', vmax=4e-5)
        plt.colorbar(label='Power')
        plt.title(plot_pop.split('__')[0])
        plt.ylabel('Frequency (Hz)')
        plt.gca().invert_yaxis()
        plt.xticks([])

        plt.yticks(np.arange(9,cuttof_freq,10),[str(y_val+1) for y_val in np.arange(9,cuttof_freq,10)])
        
        for time_windows_up_sample in time_windows_up_samples:
            # square = patches.Rectangle((time_windows_up_sample[0], lowcut-1), time_windows_up_sample[1]-time_windows_up_sample[0], highcut-lowcut-1, edgecolor='w', facecolor='none', alpha=0.2)
            # ax=plt.gca()
            # ax.add_patch(square)

            ax=plt.gca()
            draw_horizontal_bracket(ax, time_windows_up_sample[0], 16, time_windows_up_sample[1]-time_windows_up_sample[0], v_height=0.5, is_top=True,  line_width=1.0, color='r', alpha=1.0)
            draw_horizontal_bracket(ax, time_windows_up_sample[0],  8, time_windows_up_sample[1]-time_windows_up_sample[0], v_height=0.5, is_top=False, line_width=1.0, color='r', alpha=1.0)
            # draw_horizontal_bracket(ax, time_windows_up_sample[0], 16, time_windows_up_sample[1]-time_windows_up_sample[0], v_height=0.5, is_top=True,  line_width=0.35, color='r', alpha=1.0)
            # draw_horizontal_bracket(ax, time_windows_up_sample[0],  8, time_windows_up_sample[1]-time_windows_up_sample[0], v_height=0.5, is_top=False, line_width=0.35, color='r', alpha=1.0)
        
    plt.xticks(np.arange(10000,len(dataset[0])+1,20000),[int((t_val*0.1)+1500) for t_val in np.arange(10000,len(dataset[0])+1,20000)])
    plt.xlabel('Time (ms)')
    # plt.savefig(spectrogram_windows_data_path+'Spectrogram_'+sim_name+'_'+plot_pop+'.png',dpi=200)
    plt.savefig(save_figures_path+'Spectrogram_'+sim_name+'_'+plot_pop+'.png',dpi=200)

plt.show()
sys.exit()






# plot_pops=['TRN__pop','VPM__pop']
# plt.figure(figsize=(9,9))
# plt.suptitle('RMS power in band')
# for pop_ind,plot_pop in enumerate(plot_pops):
#     plt.subplot(2,1,pop_ind+1)
#     for sim_index,sim_name in enumerate(store_rms_dict.keys()):
#         if sim_index==2 or sim_index==6:    
#             alpha = 0.8
#             down_color = 'mediumpurple'
#             up_color = 'r'
#             down_edgecolor=None
#             up_edgecolor=None
#         else:                               
#             alpha = 0.4
#             down_color = 'w'
#             up_color = 'w'
#             down_edgecolor='b'
#             up_edgecolor='r'
        
#         plt.bar(     sim_index-0.2,np.mean(store_rms_dict[sim_name][plot_pop]['down']), width=0.4,                                                 color=down_color, alpha=alpha, edgecolor=down_edgecolor)
#         plt.errorbar(sim_index-0.2,np.mean(store_rms_dict[sim_name][plot_pop]['down']), yerr=np.std(store_rms_dict[sim_name][plot_pop]['down']),   color='k', capsize=3)
#         plt.bar(     sim_index+0.2,np.mean(store_rms_dict[sim_name][plot_pop]['up']),   width=0.4,                                                 color=up_color,   alpha=alpha, edgecolor=up_edgecolor)
#         plt.errorbar(sim_index+0.2,np.mean(store_rms_dict[sim_name][plot_pop]['up']),   yerr=np.std(store_rms_dict[sim_name][plot_pop]['up']),     color='k', capsize=3)

plt.figure(figsize=(9,9))
plt.suptitle('RMS power in band - Boxplot')
for pop_ind,plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,pop_ind+1)
    for sim_index,sim_name in enumerate(store_rms_dict.keys()):

        alpha = 0.8
        down_color = 'mediumpurple'
        up_color = 'r'
        down_edgecolor=None
        up_edgecolor=None

        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.2
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.2
                colors = [down_color]
                state_label='D'

            if select_sims is not None:
                x_tick = sim_shift+box_shift
            else:
                x_tick = sim_index+box_shift


            ax = plt.gca()
            boxes = ax.boxplot([store_rms_dict[sim_name][plot_pop][state]], 
                            #    labels=[state_label], 
                            #    notch=True,
                               medianprops=dict(color="k", alpha=0.7),
                               patch_artist=True, 
                               positions=np.array([x_tick]),
                               widths=0.15,
                               )
            
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(0.5)
        plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')


# plot_pops=['TRN__pop','VPM__pop']
# plt.figure(figsize=(9,9))
# plt.suptitle('Mean power in band')
# for pop_ind,plot_pop in enumerate(plot_pops):
#     plt.subplot(2,1,pop_ind+1)
#     for sim_index,sim_name in enumerate(store_mean_dict.keys()):
#         plt.bar(     sim_index-0.2,np.mean(store_mean_dict[sim_name][plot_pop]['down']), width=0.4,                                               color='b')
#         plt.errorbar(sim_index-0.2,np.mean(store_mean_dict[sim_name][plot_pop]['down']), yerr=np.std(store_mean_dict[sim_name][plot_pop]['down']), color='k')
#         plt.bar(     sim_index+0.2,np.mean(store_mean_dict[sim_name][plot_pop]['up']),   width=0.4,                                               color='r')
#         plt.errorbar(sim_index+0.2,np.mean(store_mean_dict[sim_name][plot_pop]['up']),   yerr=np.std(store_mean_dict[sim_name][plot_pop]['up']),   color='k')

plt.figure(figsize=(9,9))
plt.suptitle('Mean power in band - Boxplot')
for pop_ind,plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,pop_ind+1)
    for sim_index,sim_name in enumerate(store_mean_dict.keys()):

        alpha = 0.8
        down_color = 'mediumpurple'
        up_color = 'r'
        down_edgecolor=None
        up_edgecolor=None

        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.2
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.2
                colors = [down_color]
                state_label='D'

            if select_sims is not None:
                x_tick = sim_shift+box_shift
            else:
                x_tick = sim_index+box_shift


            ax = plt.gca()
            boxes = ax.boxplot([store_mean_dict[sim_name][plot_pop][state]], 
                            #    labels=[state_label], 
                            #    notch=True,
                               medianprops=dict(color="k", alpha=0.7),
                               patch_artist=True, 
                               positions=np.array([x_tick]),
                               widths=0.15,
                               )
            
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(0.5)
        plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')

# plt.show()
# sys.exit()
        
        
        
    
    


#######################################################################################################################################

store_band_power_dict={}
store_band_stats={}
store_freq_bands_transpose_dict={}
for sim_index, sim_name in enumerate(store_spectrogram_windows_dict.keys()):
    store_band_power_dict.update({sim_name:{}})
    store_band_stats.update({sim_name:{}})
    store_freq_bands_transpose_dict.update({sim_name:{}})
    for plot_pop in store_spectrogram_windows_dict[sim_name].keys():
        store_band_power_dict[sim_name].update({plot_pop:{}})
        store_band_stats[sim_name].update({plot_pop:{}})
        store_freq_bands_transpose_dict[sim_name].update({plot_pop:{}})
        for time_window in store_spectrogram_windows_dict[sim_name][plot_pop].keys():

            time_range = store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['extent'][0][0:2] # 1st and 2nd vals
            freq_range = store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['extent'][0][2:4] # 3rd and 4th vals

            # time
            num_vals_x_axis = len(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][0])
            # frequencies
            num_vals_y_axis = len(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0])

            time_array = np.arange(time_range[0],time_range[1],  ((time_range[1]-time_range[0])/(num_vals_x_axis)))
            freq_array = np.arange(freq_range[0],freq_range[1]+1,((freq_range[1]-freq_range[0])/(num_vals_y_axis-1)))

            store_freqs=[]
            store_freq_bands = []
            for freq_index,freq_val in enumerate(freq_array):
                if (freq_val>=lowcut) and (freq_val<=highcut):
                    store_freqs.append(freq_val)
                    store_freq_bands.append(list(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][freq_index]))

            # transposes the list of lists to be organized from list of power by frequency (where each list is a frequency) to list of powers for each timestep (where each list is a timestep)
            store_freq_bands_transpose = list(map(list, zip(*store_freq_bands)))

            # stores the values to process all together
            store_freq_bands_transpose_dict[sim_name][plot_pop].update({time_window:store_freq_bands_transpose})

            ############
            # area under the curve for each list of powers in the frequency range
            power_in_band_byTimestep=[]
            for t_step in range(len(store_freq_bands_transpose)):
                power_in_band_byTimestep.append(np.trapz(store_freq_bands_transpose[t_step], store_freqs))

            # mean area under the curve
            mean_power_in_band = np.mean(power_in_band_byTimestep)            
            std_power_in_band  = np.std(power_in_band_byTimestep)            

            store_band_stats[sim_name][plot_pop].update({time_window:{'mean':mean_power_in_band,'std':std_power_in_band}})


            sum_band_power=0
            for freq_index,freq_val in enumerate(freq_array):
                if (freq_val>=lowcut) and (freq_val<=highcut):
                    # print('freq_index:', freq_index, ' | freq_val: ', freq_val)
                    sum_freq_band=sum(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][freq_index])
                    sum_band_power+=sum_freq_band

            # sum_band_power = sum([sum(store_spectrogram_windows_dict[sim_name][plot_pop][time_window]['electrodes']['spectrogram']['morlet'][0][freq_index]) for freq_index,freq_val in enumerate(freq_array) if (freq_val>=lowcut) and (freq_val<=highcut)])
            
            store_band_power_dict[sim_name][plot_pop].update({time_window:sum_band_power})

# store_freq_bands_transpose_dict_pooled={}
# for sim_index, sim_name in enumerate(store_freq_bands_transpose_dict.keys()):
#     store_freq_bands_transpose_dict_pooled.update({sim_name:{}})
#     for plot_pop in store_freq_bands_transpose_dict[sim_name].keys():
#         store_freq_bands_transpose_dict_pooled[sim_name].update({plot_pop:{'up':[],'down':[]}})
#         for time_window in store_freq_bands_transpose_dict[sim_name][plot_pop].keys():
#             if   'u_' in time_window: store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['up'].extend(store_freq_bands_transpose_dict[sim_name][plot_pop][time_window])
#             elif 'd_' in time_window: store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['down'].extend(store_freq_bands_transpose_dict[sim_name][plot_pop][time_window])
#             else: continue

# summing the power on the frequnecy band for each timestep
store_freq_bands_transpose_dict_sum={}
for sim_index, sim_name in enumerate(store_freq_bands_transpose_dict.keys()):
    store_freq_bands_transpose_dict_sum.update({sim_name:{}})
    for plot_pop in store_freq_bands_transpose_dict[sim_name].keys():
        store_freq_bands_transpose_dict_sum[sim_name].update({plot_pop:{'up':[],'down':[]}})
        for time_window in store_freq_bands_transpose_dict[sim_name][plot_pop].keys():
            store_freq_bands_transpose_dict_sum[sim_name][plot_pop][time_window] = [sum(pwr) for pwr in store_freq_bands_transpose_dict[sim_name][plot_pop][time_window]]

store_freq_bands_transpose_dict_pooled={}
for sim_index, sim_name in enumerate(store_freq_bands_transpose_dict_sum.keys()):
    store_freq_bands_transpose_dict_pooled.update({sim_name:{}})
    for plot_pop in store_freq_bands_transpose_dict_sum[sim_name].keys():
        store_freq_bands_transpose_dict_pooled[sim_name].update({plot_pop:{'up':[],'down':[]}})
        for time_window in store_freq_bands_transpose_dict_sum[sim_name][plot_pop].keys():
            if   'u_' in time_window: store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['up'].extend(store_freq_bands_transpose_dict_sum[sim_name][plot_pop][time_window])
            elif 'd_' in time_window: store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['down'].extend(store_freq_bands_transpose_dict_sum[sim_name][plot_pop][time_window])
            else: continue

plt.figure(figsize=(9,9))
plt.suptitle('Power in band - Boxplot - Pooled raw data')
for pop_ind,plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,pop_ind+1)
    sim_shift=0
    for sim_index,sim_name in enumerate(store_freq_bands_transpose_dict_pooled.keys()):

        alpha = 0.8
        down_color = 'mediumpurple'
        up_color = 'r'
        down_edgecolor=None
        up_edgecolor=None

        for state in ['up','down']:    
            if state == 'up':   
                box_shift = 0.2
                colors = [up_color]
                state_label='U'
            else:               
                box_shift = -0.2
                colors = [down_color]
                state_label='D'

            if select_sims is not None:
                x_tick = sim_shift+box_shift
            else:
                x_tick = sim_index+box_shift


            ax = plt.gca()
            boxes = ax.boxplot([store_freq_bands_transpose_dict_pooled[sim_name][plot_pop][state]], 
                            #    labels=[state_label], 
                            #    notch=True,
                               medianprops=dict(color="k", alpha=0.7),
                               patch_artist=True, 
                               positions=np.array([x_tick]),
                               widths=0.15,
                               )
            sim_shift+=1
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], colors):
                box.set_facecolor(color)
                box.set_alpha(0.5)
        plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')


plt.figure(figsize=(9,9))
plt.suptitle('Power in band - Boxplot - Pooled raw data')
plot_counter=0
for pop_ind,plot_pop in enumerate(plot_pops):
    # plt.subplot(2,1,pop_ind+1)
    for sim_index,sim_name in enumerate(store_freq_bands_transpose_dict_pooled.keys()):
        plot_counter+=1
        plt.subplot(2,9,plot_counter)
        plt.hist(store_freq_bands_transpose_dict_pooled[sim_name][plot_pop][state],bins=50)

# sum freq bands

store_freq_bands_transpose_dict_stats={}
for sim_index, sim_name in enumerate(store_freq_bands_transpose_dict_pooled.keys()):
    store_freq_bands_transpose_dict_stats.update({sim_name:{}})
    
    for plot_pop in store_freq_bands_transpose_dict_pooled[sim_name].keys():
        store_freq_bands_transpose_dict_stats[sim_name].update({plot_pop:{
            'up':{
                'mean':np.mean(store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['up']),
                'std': np.std( store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['up']),
                },
            'down':{
                'mean':np.mean(store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['down']),
                'std': np.std( store_freq_bands_transpose_dict_pooled[sim_name][plot_pop]['down']),
                },
        }})

plot_pops=['TRN__pop','VPM__pop']
plt.figure(figsize=(9,9))
for pop_ind,plot_pop in enumerate(plot_pops):
    plt.subplot(2,1,pop_ind+1)
    for sim_index, sim_name in enumerate(store_freq_bands_transpose_dict_stats.keys()):
        plt.bar(sim_index-0.2,store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['down']['mean'],width=0.4,color='b')
        plt.errorbar(sim_index-0.2,store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['down']['mean'],yerr=store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['down']['std'],color='k')
        plt.bar(sim_index+0.2,store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['up']['mean'],width=0.4,color='r')
        plt.errorbar(sim_index+0.2,store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['up']['mean'],yerr=store_freq_bands_transpose_dict_stats[sim_name][plot_pop]['up']['std'],color='k')

''''
In [202]: store_freq_bands_transpose_dict_stats
Out[202]: 
    {
        'sim_0': {  'VPM__pop': {'mean': 6.1654298137068675e-06,    'std': 7.083697651202317e-06},
                    'TRN__pop': {'mean': 4.608061360688848e-06,     'std': 4.999927279941659e-06}},
        'sim_1': {  'VPM__pop': {'mean': 4.920601457813169e-06,     'std': 5.348115760907324e-06},
                    'TRN__pop': {'mean': 4.267929904468301e-06,     'std': 4.792009001969788e-06}},
        'sim_2': {  'VPM__pop': {'mean': 4.538391127011015e-06,     'std': 5.023375532412966e-06},
                    'TRN__pop': {'mean': 4.490865150649334e-06,     'std': 4.8580793806475926e-06}},
        'sim_3': {  'VPM__pop': {'mean': 2.447494035758693e-06,     'std': 2.6410350058762637e-06},
                    'TRN__pop': {'mean': 2.936900800036073e-06,     'std': 3.062411376176449e-06}},
        'sim_4': {  'VPM__pop': {'mean': 3.220274418727559e-06,     'std': 3.957469702136875e-06},
                    'TRN__pop': {'mean': 3.2350284002051895e-06,    'std': 3.936983823640933e-06}},
        'sim_5': {  'VPM__pop': {'mean': 1.1002532202149505e-06,    'std': 1.1169628382412444e-06},
                    'TRN__pop': {'mean': 2.5963558491703634e-06,    'std': 2.658850074039185e-06}},
        'sim_6': {  'VPM__pop': {'mean': 9.391285276407625e-07,     'std': 9.402138688418937e-07},
                    'TRN__pop': {'mean': 2.5719798668192042e-06,    'std': 2.771534934889074e-06}},
        'sim_7': {  'VPM__pop': {'mean': 5.256289153851991e-07,     'std': 5.700584040273944e-07},
                    'TRN__pop': {'mean': 1.775007561129801e-06,     'std': 2.0204609128819862e-06}},
        'sim_8': {  'VPM__pop': {'mean': 6.469712643595799e-07,     'std': 7.42818979855675e-07},
                    'TRN__pop': {'mean': 1.1550436514381913e-06,    'std': 1.3297207570974143e-06}
                    }
    }
    {
        'sim_0': {  'VPM__pop': {   'up': {     'mean': 6.1654298137068675e-06,     'std': 7.083697651202317e-06},
                                    'down': {   'mean': 6.418015592369788e-07,      'std': 6.708071668356058e-07}},
                    'TRN__pop': {   'up': {     'mean': 4.608061360688848e-06,      'std': 4.999927279941659e-06},
                                    'down': {   'mean': 3.6556157514894656e-07,     'std': 3.742782006557369e-07}}},
        'sim_1': {  'VPM__pop': {   'up': {     'mean': 4.920601457813169e-06,      'std': 5.348115760907324e-06},
                                    'down': {   'mean': 7.182037364746858e-07,      'std': 8.63136848306883e-07}},
                    'TRN__pop': {   'up': {     'mean': 4.267929904468301e-06,      'std': 4.792009001969788e-06},
                                    'down': {   'mean': 3.3150735948434846e-07,     'std': 3.6443945557695605e-07}}},
        'sim_2': {  'VPM__pop': {   'up': {     'mean': 4.538391127011015e-06,      'std': 5.023375532412966e-06},
                                    'down': {   'mean': 6.424014650079329e-07,      'std': 9.783043079179346e-07}},
                    'TRN__pop': {   'up': {     'mean': 4.490865150649334e-06,      'std': 4.8580793806475926e-06},
                                    'down': {   'mean': 2.591914649456824e-07,      'std': 3.132578234210637e-07}}},
        'sim_3': {  'VPM__pop': {   'up': {     'mean': 2.447494035758693e-06,      'std': 2.6410350058762637e-06},
                                    'down': {   'mean': 6.179228959121295e-07,      'std': 7.398465494851129e-07}},
                    'TRN__pop': {   'up': {     'mean': 2.936900800036073e-06,      'std': 3.062411376176449e-06},
                                    'down': {   'mean': 1.7997819979001596e-07,     'std': 1.7501412685932592e-07}}}}
    
    
    '''
        


# store_band_power_stats_dict={}
# for sim_index, sim_name in enumerate(store_spectrogram_windows_dict.keys()):
#     store_band_power_stats_dict.update({sim_name:{}})
#     for plot_pop in store_spectrogram_windows_dict[sim_name].keys():
#         window_powers = list(store_band_power_dict[sim_name][plot_pop].values())
#         store_band_power_stats_dict[sim_name].update({plot_pop:{'mean':np.mean(window_powers),'std':np.std(window_powers)}})

#######################################################################################################################################

# sys.exit()


    # plt.figure(figsize=[15,8])

    # # VPM
    # cmap = matplotlib.colormaps[pop_color_maps['VPM__pop']]
    # # line_color = cmap(sim_index/len(store_PSD_dict.keys()))
    # line_color = sim_color[sim_index]
    # plt.subplot(2,1,1)
    # time_signal = store_lfpTraces_dict[sim_name][recording_frequency]['t']
    # lfp_signal  = store_lfpTraces_dict[sim_name][recording_frequency]['data'][2] # Replace with your LFP signal list
    # plt.plot(time_signal,lfp_signal,color=line_color, linewidth=sim_linewidth)
    # # plt.plot(store_lfpTraces_dict[sim_name][recording_frequency]['t'],store_lfpTraces_dict[sim_name][recording_frequency]['data'][2],color='green')
    # plt.ylim([-0.02,0.02])
    # plt.yticks([])
    # # plt.yticks([-0.02,0,0.02])

    # ax1=plt.gca()
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_visible(False)
    # ax1.spines['left'].set_visible(False)

    # if sim_name == 'sim_6': 
    #     plt.plot([1250,1750],[-0.015,-0.015],'k',linewidth=3)
    #     plt.plot([1250,1250],[ 0.005,-0.015],'k',linewidth=3)
    #     plt.text(1320, -0.013, '500 ms', color='k', rotation=0, size=15, weight=font_weight, alpha = 1.0)
    #     plt.text(1125, -0.013, '0.02 mV', color='k', rotation=90,size=15, weight=font_weight, alpha = 1.0)
    
    # plt.xticks([])
    # plt.xlim([1000,16000])
    # # plt.xlabel('Time (ms)')

    # # --- 

    # # Example usage    
    # filtered_lfp = bandpass_filter(lfp_signal, fs, lowcut, highcut)

    # plt.subplot(2,1,2)
    # plt.plot(time_signal,filtered_lfp,color=line_color, linewidth=sim_linewidth)
    # plt.ylim([-0.02,0.02])
    # plt.yticks([])
    # # plt.yticks([-0.02,0,0.02])

    # # Compute envelope
    # lfp_envelope = compute_envelope(filtered_lfp)
    # plt.plot(time_signal,lfp_envelope,'--',color='firebrick', linewidth=0.5)

    # ax1=plt.gca()
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_visible(False)
    # ax1.spines['left'].set_visible(False)

    # plt.xlim([1000,16000])

    # plt.xticks([])
    # plt.xlabel('')

    # # # --- 

    # # filtered_lfp_abs = [abs(pwr) for pwr in filtered_lfp]

    # # moving_sums, moving_averages, time_centers = moving_sum_avg(filtered_lfp_abs, fs, window_ms, overlap_ms, start_time=1500)
    # # plt.subplot(6,1,6)
    # # plt.plot([1500,16000], [1.0, 1.0], '--',color='firebrick', linewidth=0.5)
    # # plt.plot(time_centers,moving_sums,color=line_color, linewidth=sim_linewidth)

    # # # plt.ylim([-0.02,0.02])
    # # plt.yticks([])
    # # # plt.yticks([-0.02,0,0.02])

    # # ax1=plt.gca()
    # # ax1.spines['top'].set_visible(False)
    # # ax1.spines['right'].set_visible(False)
    # # ax1.spines['bottom'].set_visible(False)
    # # ax1.spines['left'].set_visible(False)

    # # plt.xlim([1000,16000])
    # # plt.xticks([])
    # # plt.xlabel('')

    # # --- 

    # plt.savefig(lfp_traces_data_path+'lfpTraces_'+sim_name+'_VPM__pop.png',dpi=200)
    # # plt.savefig(lfp_traces_data_path+'lfpTraces_'+sim_name+'.png',dpi=200)


# plt.figure(figsize=(5,5))
# for sim_index, sim_name in enumerate(mean_envelope_dict.keys()):
#     for pop_ind,plot_pop in enumerate(['TRN__pop','VPM__pop']):
#         plt.subplot(2,1,pop_ind+1)
#         plt.bar([sim_index-0.2,sim_index+0.2],[mean_envelope_dict[sim_name]['mean_env_down'],mean_envelope_dict[sim_name]['mean_env_up']],width=0.4,color=['b','r'])

# plt.savefig(lfp_traces_data_path+'lfpTraces__TRN__pop__upDown_bar.png',dpi=200)


sys.exit()


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

# --- PSD bar plot of max power value and frequency of max power
store_window_power_overTime_dict={}
plt.figure(figsize=[15,5])
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(4,1,plot_pop_ind+1)
    cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
    auc_store_window_powers = []
    auc_store_window_powers_Thresh=[]

    store_window_power_overTime_dict.update({plot_pop:{}})

    for sim_index, sim_name in enumerate(store_lfpPSD_windowed_dict.keys()): 
        
        # if (sim_index!=2) and (sim_index!=4) and (sim_index!=6):continue
        
        store_window_power_overTime = []
        store_window_time_overTime=[]
        for time_window in store_lfpPSD_windowed_dict[sim_name][plot_pop].keys():

            slice_data  = [(freq,store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['freqs']) if freq>=spindle_freq_range[0] and freq<=spindle_freq_range[1]]
            max_power_spindleRange  = max(slice_data, key=lambda x: x[1])
            
            full_data = [(freq,store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['signal'][freq_ind]) for freq_ind,freq in enumerate(store_lfpPSD_windowed_dict[sim_name][plot_pop][time_window]['freqs']) if freq<=limit_frequency]
            max_power_signal        = max(full_data, key=lambda x: x[1])

            # Extract frequency and power values
            slice_freqs, slice_power = zip(*slice_data)
            # Compute area under the curve using the trapezoidal rule
            power_in_band = np.trapz(slice_power, slice_freqs)

            store_window_power_overTime.append(power_in_band)
            store_window_time_overTime.append((int(time_window.split('__')[1])+int(time_window.split('__')[0]))/2)

            pwrThresh=5e-5
            store_window_power_overTime_aboveThresh_inds = [pwr_ind for pwr_ind,pwr in enumerate(store_window_power_overTime) if pwr>=pwrThresh]
            store_window_power_overTime_aboveThresh = [store_window_power_overTime[pwr_ind] for pwr_ind in store_window_power_overTime_aboveThresh_inds]
            store_window_time_overTime_aboveThresh = [store_window_time_overTime[pwr_ind] for pwr_ind in store_window_power_overTime_aboveThresh_inds]

        auc_store_window_powers.append(np.trapz(store_window_power_overTime, store_window_time_overTime))
        
        auc_store_window_powers_Thresh.append(np.trapz(store_window_power_overTime_aboveThresh, store_window_time_overTime_aboveThresh))

        if (sim_index == 2) or (sim_index == 6):
            if (sim_index == 2):    line_color = 'k'
            else:                   line_color = 'grey'
            line_alpha = 1.0
            line_width = 4
        else:
            line_alpha = 0.5
            line_width = 2
            line_color = 'k'
        
        plt.plot(store_window_time_overTime,store_window_power_overTime,
                                            # color=line_color,
                                            color=cmap(sim_index/len(store_lfpPSD_windowed_dict.keys())),
                                            linewidth=line_width, 
                                            alpha=line_alpha)
        
        store_window_power_overTime_dict[plot_pop].update({sim_name:{'t':store_window_time_overTime, 'power': store_window_power_overTime}})

    plt.ylim(lfpPSD_window_y_lims_bar[plot_pop])
    plt.ylabel('Power 8-16 Hz')
    plt.yticks([0,5.0e-5,10.0e-5,15.0e-5,20.0e-5])
    
    plt.subplot(4,1,plot_pop_ind+3)
    # plt.bar(list(store_lfpPSD_windowed_dict.keys()), auc_store_window_powers)
    plt.bar(list(store_lfpPSD_windowed_dict.keys()), auc_store_window_powers_Thresh)

# n_bins = 50
# for plot_pop_ind, plot_pop in enumerate(plot_pops):
#     plt.figure(figsize=[15,5])
#     cmap = matplotlib.colormaps[pop_color_maps[plot_pop]]
#     plot_count=0
#     for sim_index, sim_name in enumerate(store_lfpPSD_windowed_dict.keys()): 
#         if (sim_index!=0) and (sim_index!=2) and (sim_index!=4) and (sim_index!=6) and (sim_index!=8):continue
#         plt.subplot(1,5,plot_count+1)
#         # plt.subplot(1,len(store_lfpPSD_windowed_dict.keys()),plot_count+1)
#         plt.hist(store_window_power_overTime_dict[plot_pop][sim_name]['power'],n_bins,color='w',edgecolor=cmap(sim_index/len(store_lfpPSD_windowed_dict.keys())),histtype='step')
#         plot_count+=1
#         plt.ylim(0,50)
#         plt.xlim(0,3.1e-4)
#     # plt.tight_layout()

# plt.show()






#         bar_value = power_in_band
#         bar_string = ''
#         bar_string_color = 'k'
#         font_weight='bold'
        
#         if (sim_index == 2) or (sim_index == 6):
#             if (sim_index == 2):    bar_color = 'k'
#             else:                   bar_color = 'grey'
#             bar_alpha = 1.0
#             edge_color = 'k'
#         else:
#             bar_color = 'w'
#             bar_alpha = 0.25
#             edge_color = 'k'


#         plt.bar(sim_index_mapping[sim_index],bar_value,width=0.21,color=bar_color, edgecolor = edge_color, alpha=bar_alpha)

#     plt.ylim(lfpPSD_y_lims_bar[plot_pop])
#     plt.ylabel('Power 8-16 Hz')
#     plt.yticks([0,5.0e-5,10.0e-5])
    
#     if plot_pop_ind==len(plot_pops)-1:
#         plt.xticks([-1,0,1])
#         # plt.xlabel('$P_{T}$')
#         plt.xlabel('CSI')
#     else:
#         plt.xticks([])
#         plt.xlabel('')
    
#     ax1=plt.gca()
#     ax1.spines['top'].set_visible(False)
#     ax1.spines['right'].set_visible(False)
#     ax1.spines['bottom'].set_visible(False)
#     ax1.spines['left'].set_visible(False)

#     plt.tight_layout()

# plt.savefig(lfp_psd_data_path+'lfpPSD_peakFrequency.png',dpi=200)



# ###############################################################################################################################################
plt.show()



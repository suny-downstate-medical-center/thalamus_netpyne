'''
script replacing AnalyzeAngularTuningStats_CT.py
'''

import os
from cProfile import label
import json
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import sys
from matplotlib import pyplot as plt
from sympy import plot
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
    
# ########################################################################################################################################################################################################

# data_path = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data'
# sim_output_folder = 'sim_output'
# _addVersion='_v04'

# # batch_sim_group = 'c_0003'
# # batch_sim_group = 'paper_000' # closed-loop CT feedback by default
# batch_sim_group = 'paper_001' # uniform CT feedback by default

# # batch_sim_label = batch_sim_group+'_'+'barreloid_batch_fig05_sleep_wake'

# if batch_sim_group == 'paper_000':
#     # paper_000
#     batch_sim_labels = [
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),

#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_0'      , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_10'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_25'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_50'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_100'    , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),

#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#                 ]
# else:
#     # paper_001
#     batch_sim_labels = [
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_0_uniform_ct_modulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),

#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_0'      , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_10'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_25'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_50'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
#                 (batch_sim_group+'_barreloid_batch_fig07_1_mixed_ct_modulation_100'    , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),

#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
#                 (batch_sim_group+'_barreloid_batch_fig07_2_closed_ct_modulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
#                 ]

# simulation_data_paths = [data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/sim_output/' for (batch_sim_label,connTemplate) in batch_sim_labels]

# conn_templates_set = [  'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform', 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50', 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform']

# ########################################################################################################################################################################################################

# # --- Path to save analysis figures
# save_figures_path = data_path+'/'+batch_sim_group+'/'+batch_sim_group+'_barreloid_batch_fig07_analysis_figures/'
# if not os.path.exists(save_figures_path): os.makedirs(save_figures_path)

########################################################################################################################################################################################################

data_path = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data'
sim_output_folder = 'sim_output'
_addVersion='_v04'

batch_sim_group = 'paper_002'

# analyze_figure = 'e'  # closed-loop CT feedback by default
analyze_figure = 'f'  # uniform CT feedback by default

if analyze_figure == 'e':
    # closed-loop CT feedback by default
    batch_sim_labels = [
                (batch_sim_group+'_barreloid_batch_fig07e_0_uniform_ct_modulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07e_0_uniform_ct_modulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07e_0_uniform_ct_modulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07e_0_uniform_ct_modulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07e_0_uniform_ct_modulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),

                (batch_sim_group+'_barreloid_batch_fig07e_1_mixed_ct_modulation_0'      , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07e_1_mixed_ct_modulation_10'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07e_1_mixed_ct_modulation_25'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07e_1_mixed_ct_modulation_50'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07e_1_mixed_ct_modulation_100'    , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),

                (batch_sim_group+'_barreloid_batch_fig07e_2_closed_ct_modulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
                (batch_sim_group+'_barreloid_batch_fig07e_2_closed_ct_modulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
                (batch_sim_group+'_barreloid_batch_fig07e_2_closed_ct_modulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
                (batch_sim_group+'_barreloid_batch_fig07e_2_closed_ct_modulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
                (batch_sim_group+'_barreloid_batch_fig07e_2_closed_ct_modulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
                ]
else:
    # uniform CT feedback by default
    batch_sim_labels = [
                (batch_sim_group+'_barreloid_batch_fig07f_0_uniform_ct_modulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_0_uniform_ct_modulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_0_uniform_ct_modulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_0_uniform_ct_modulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_0_uniform_ct_modulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),

                (batch_sim_group+'_barreloid_batch_fig07f_1_mixed_ct_modulation_0'      , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07f_1_mixed_ct_modulation_10'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07f_1_mixed_ct_modulation_25'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07f_1_mixed_ct_modulation_50'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),
                (batch_sim_group+'_barreloid_batch_fig07f_1_mixed_ct_modulation_100'    , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50'),

                (batch_sim_group+'_barreloid_batch_fig07f_2_closed_ct_modulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_2_closed_ct_modulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_2_closed_ct_modulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_2_closed_ct_modulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
                (batch_sim_group+'_barreloid_batch_fig07f_2_closed_ct_modulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform'),
                ]

simulation_data_paths = [data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/sim_output/' for (batch_sim_label,connTemplate) in batch_sim_labels]

conn_templates_set = [  'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform', 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50', 'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform']

########################################################################################################################################################################################################

# --- Path to save analysis figures
save_figures_path = data_path+'/'+batch_sim_group+'/'+batch_sim_group+'_barreloid_batch_fig07'+analyze_figure+'_analysis_figures/'
if not os.path.exists(save_figures_path): os.makedirs(save_figures_path)

########################################################################################################################################################################################################

# # simFolder='c_0003'
# # _addVersion='_v04'
# # 100% uniform
# # conn_template = 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'
# simNames = [
#             ('c_0003_300_awake_0_uniform_CTmodulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),]

# # # 50% closed / 50% uniform
# # conn_template = 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'
# # simNames = [
# #             ('c_0003_300_awake_1_mixed5050_CTmodulation_0'      , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
# #             ('c_0003_300_awake_1_mixed5050_CTmodulation_10'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
# #             ('c_0003_300_awake_1_mixed5050_CTmodulation_25'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
# #             ('c_0003_300_awake_1_mixed5050_CTmodulation_50'     , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
# #             ('c_0003_300_awake_1_mixed5050_CTmodulation_100'    , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),]

# # # 100% closed
# # conn_template = 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'
# # simNames = [
# #             ('c_0003_300_awake_2_closed_CTmodulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
# #             ('c_0003_300_awake_2_closed_CTmodulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
# #             ('c_0003_300_awake_2_closed_CTmodulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
# #             ('c_0003_300_awake_2_closed_CTmodulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
# #             ('c_0003_300_awake_2_closed_CTmodulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),]

# simNames = [
#             ('c_0003_300_awake_0_uniform_CTmodulation_0'    , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_10'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_25'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_50'   , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),
#             ('c_0003_300_awake_0_uniform_CTmodulation_100'  , 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform'),

#             ('c_0003_300_awake_1_mixed5050_CTmodulation_0'  , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#             ('c_0003_300_awake_1_mixed5050_CTmodulation_10' , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#             ('c_0003_300_awake_1_mixed5050_CTmodulation_25' , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#             ('c_0003_300_awake_1_mixed5050_CTmodulation_50' , 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),
#             ('c_0003_300_awake_1_mixed5050_CTmodulation_100', 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50'),

#             ('c_0003_300_awake_2_closed_CTmodulation_0'     , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#             ('c_0003_300_awake_2_closed_CTmodulation_10'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#             ('c_0003_300_awake_2_closed_CTmodulation_25'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#             ('c_0003_300_awake_2_closed_CTmodulation_50'    , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#             ('c_0003_300_awake_2_closed_CTmodulation_100'   , 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed'),
#             ]

# conn_templates_set = [  'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform', 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50', 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed']

########################################################################################################################################################################################################

boxplot_colors          = ['lightgray','dimgray','k']
angular_tuning_windows  = ['20ms']
ct_feedbacks            = ['ct_uniform']
dataset_type            = 'normalized'

angulart_tuning_pooled_data_dict={}
mean_angulart_tuning_pooled_data_dict={}

angulart_tuning_pooled_data_dict_absolute={}
mean_angulart_tuning_pooled_data_dict_absolute={}

for angular_tuning_window in angular_tuning_windows:
    angulart_tuning_pooled_data_dict.update({angular_tuning_window:{}})
    mean_angulart_tuning_pooled_data_dict.update({angular_tuning_window:{}})
    angulart_tuning_pooled_data_dict_absolute.update({angular_tuning_window:{}})
    mean_angulart_tuning_pooled_data_dict_absolute.update({angular_tuning_window:{}})

    for ct_feedback in ct_feedbacks:
        angulart_tuning_pooled_data_dict[angular_tuning_window].update({ct_feedback:{}})
        mean_angulart_tuning_pooled_data_dict[angular_tuning_window].update({ct_feedback:{}})
        angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].update({ct_feedback:{}})
        mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].update({ct_feedback:{}})
        # for conn_template in conn_templates:
        for source_folder_ind, source_folder in enumerate(simulation_data_paths):
            try: conn_template = batch_sim_labels[source_folder_ind][1]+_addVersion
            except: print('Adding conn version failed')

            ############################################################################################################################################
            if angular_tuning_window == '20ms':
                file_flag = 'angularTuning3_spikes_ON'
                dataFolder      = source_folder+'angular_tuning_ON/'
            ############################################################################################################################################
            elif angular_tuning_window == '250ms':
                file_flag = 'angularTuning3_spikes_all'
                dataFolder      = source_folder+'angular_tuning_all/'
            ############################################################################################################################################
            # normalized data
            angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].update({batch_sim_labels[source_folder_ind][0]:{}})
            mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].update({batch_sim_labels[source_folder_ind][0]:{}})
            filename = dataFolder+'angular_tuning_Statistics__'+dataset_type+'__'+angular_tuning_window+'__'+ct_feedback+'__conn_'+conn_template+'.json'
            with open(filename, 'r') as file: ang_tun_data = json.load(file)

            angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]]=ang_tun_data

            for pop in ang_tun_data:
                mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]].update({pop:{}})
                max_mean = max([np.mean(ang_tun_data[pop][angle]) for angle in ang_tun_data[pop].keys()])
                print(batch_sim_labels[source_folder_ind][0], pop, max_mean)
                for angle in ang_tun_data[pop].keys():
                    mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][pop].update({angle:{'mean':np.mean(ang_tun_data[pop][angle])/max_mean,'std':np.std(ang_tun_data[pop][angle])/max_mean}})

            # absolute data
            angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].update({batch_sim_labels[source_folder_ind][0]:{}})
            mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].update({batch_sim_labels[source_folder_ind][0]:{}})
            filename_absolute = dataFolder+'angular_tuning_Statistics__absolute__'+angular_tuning_window+'__'+ct_feedback+'__conn_'+conn_template+'.json'
            with open(filename_absolute, 'r') as file: ang_tun_data_absolute = json.load(file)

            angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]]=ang_tun_data_absolute

            for pop in ang_tun_data_absolute:
                mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]].update({pop:{}})
                max_mean = max([np.mean(ang_tun_data_absolute[pop][angle]) for angle in ang_tun_data_absolute[pop].keys()])
                print(batch_sim_labels[source_folder_ind][0], pop, max_mean)
                for angle in ang_tun_data_absolute[pop].keys():
                    mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][pop].update({
                            angle:{
                                'mean':np.mean(ang_tun_data_absolute[pop][angle])/max_mean,
                                'std':np.std(ang_tun_data_absolute[pop][angle])/max_mean,

                                'mean_notNormalized':np.mean(ang_tun_data_absolute[pop][angle]),
                                'std_notNormalized':np.std(ang_tun_data_absolute[pop][angle]),
                                }
                            }
                        )

########################################################################################################################################################################################################

angular_tuning_hartings={
                        'MLe__pop':{
                            '270':  1.0,
                            '315':  0.6704982467108909,
                            '0':    0.33882794754053996,
                            '45':   0.22636951177424108,
                            '90':   0.19764487038033318,
                            '135':  0.17639018632728395,
                            '180':  0.3105937087413645,
                            '225':  0.610603008575053,},
                        'VPM__pop':{
                            '270':  1.0,
                            '315':  367.0387667424478/589.0280369402608,
                            '0':    294.51401847013045/589.0280369402608,
                            '45':   241.4544744541332/589.0280369402608,
                            '90':   231.9652440864924/589.0280369402608,
                            '135':  270.7449775143726/589.0280369402608,
                            '180':  332.24549916261424/589.0280369402608,
                            '225':  394.3636706579752/589.0280369402608,},
                        'TRN__pop':{
                            '270':  1.0,
                            '315':  465.9616573671514/589.0280369402608,
                            '0':    394.82515964684/589.0280369402608,
                            '45':   340.35391028959623/589.0280369402608,
                            '90':   340.52806291144435/589.0280369402608,
                            '135':  371.63102179991336/589.0280369402608,
                            '180':  406.7796123474516/589.0280369402608,
                            '225':  478.96760294949587/589.0280369402608,},
    }
angular_tuning_hartings_meanVals={
                        'MLe__pop':{
                            '270':  1.0,
                            '315':  (0.6704982467108909+0.610603008575053)/2,
                            '0':    (0.33882794754053996+0.3105937087413645)/2,
                            '45':   (0.22636951177424108+0.17639018632728395)/2,
                            '90':   0.19764487038033318,
                            '135':  (0.22636951177424108+0.17639018632728395)/2,
                            '180':  (0.33882794754053996+0.3105937087413645)/2,
                            '225':  (0.6704982467108909+0.610603008575053)/2,},
                        'VPM__pop':{
                            '270':  1.0,
                            '315':  ((367.0387667424478/589.0280369402608)+(394.3636706579752/589.0280369402608))/2,
                            '0':    ((294.51401847013045/589.0280369402608)+(332.24549916261424/589.0280369402608))/2,
                            '45':   ((241.4544744541332/589.0280369402608)+(270.7449775143726/589.0280369402608))/2,
                            '90':   231.9652440864924/589.0280369402608,
                            '135':  ((241.4544744541332/589.0280369402608)+(270.7449775143726/589.0280369402608))/2,
                            '180':  ((294.51401847013045/589.0280369402608)+(332.24549916261424/589.0280369402608))/2,
                            '225':  ((367.0387667424478/589.0280369402608)+(394.3636706579752/589.0280369402608))/2,},
                        'TRN__pop':{
                            '270':  1.0,
                            '315':  ((465.9616573671514/589.0280369402608)+(478.96760294949587/589.0280369402608))/2,
                            '0':    ((394.82515964684/589.0280369402608)+(406.7796123474516/589.0280369402608))/2,
                            '45':   ((340.35391028959623/589.0280369402608)+(371.63102179991336/589.0280369402608))/2,
                            '90':   340.52806291144435/589.0280369402608,
                            '135':  ((340.35391028959623/589.0280369402608)+(371.63102179991336/589.0280369402608))/2,
                            '180':  ((394.82515964684/589.0280369402608)+(406.7796123474516/589.0280369402608))/2,
                            '225':  ((465.9616573671514/589.0280369402608)+(478.96760294949587/589.0280369402608))/2,},
    }

angles_sequence=['135', '180', '225', '270', '315', '0', '45', '90']

pop_colors={
    'MLe__pop':'gist_gray',
    # 'VPM__pop':'gist_yarg',
    # 'TRN__pop':'gist_yarg',
    'VPM__pop':'YlGn',
    'TRN__pop':'Blues',
    # 'VPM__pop':'summer',
    # 'TRN__pop':'winter',
    'CTvirtual_uniform__pop':'Reds',
    'other': 'Purples',
            }

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

ct_feedback = 'ct_uniform'
plot_angular_tuning_windows  = angular_tuning_windows
print('\t - CT feedback: ', ct_feedback)

plot_pops = ['TRN__pop', 'VPM__pop']

########################################################################################################################################################################################################

# --- Calculating linear metrics to estimate the deviance between the model and experimental angular tuning
remove_270deg=True # removing 270 deg from the calculation, since it is the normalization reference and will always return error=0, altering the results
store_RMSE_dict={}
store_MAE_dict={}
for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
    store_RMSE_dict.update({plot_pop:{}})
    store_MAE_dict.update({plot_pop:{}})
    # print(plot_pop)
    for angular_tuning_window in plot_angular_tuning_windows:
        store_RMSE_dict[plot_pop].update({angular_tuning_window:{}})
        store_MAE_dict[plot_pop].update({angular_tuning_window:{}})
        # print(angular_tuning_window)
        # for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
        for source_folder_ind, source_folder in enumerate(simulation_data_paths):
            store_RMSE_dict[plot_pop][angular_tuning_window].update({batch_sim_labels[source_folder_ind][0]:{}})
            store_MAE_dict[plot_pop][angular_tuning_window].update({batch_sim_labels[source_folder_ind][0]:{}})
            # print(batch_sim_labels[source_folder_ind][0])
            for angle_ind, angle in enumerate(angles_sequence):
                if remove_270deg:
                    if angle=='270':continue
                # angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][plot_pop][angle]
                # angular_tuning_hartings_meanVals[plot_pop][angle]
                # 1) RMSE 
                rmse = np.sqrt(np.mean((np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][plot_pop][angle]) - np.array(angular_tuning_hartings_meanVals[plot_pop][angle]))**2)) 
                # 2) MAE (manual version)
                mae = np.mean(np.abs(np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][plot_pop][angle]) - np.array(angular_tuning_hartings_meanVals[plot_pop][angle])))
                
                # print(angle, '\tRMSE',rmse, '\tMAE',mae)
                
                store_RMSE_dict[plot_pop][angular_tuning_window][batch_sim_labels[source_folder_ind][0]].update({angle:rmse})
                store_MAE_dict[plot_pop][angular_tuning_window][batch_sim_labels[source_folder_ind][0]].update({angle:mae})

# --- RMSE - Currently being used
simName_CT_0=batch_sim_labels[0][0]

for pop_ind,plot_pop in enumerate(plot_pops):

    if   plot_pop=='TRN__pop': test_reference='two-sided'
    elif plot_pop=='VPM__pop': test_reference='two-sided'
    else:                      test_reference='two-sided'
    # if   plot_pop=='TRN__pop': test_reference='greater'
    # elif plot_pop=='VPM__pop': test_reference='less'
    # else:                      test_reference='greater'

    paiwise_test_reference = 'two-sided'

    print(f"\nPlot pop: {plot_pop}")
    for angular_tuning_window in plot_angular_tuning_windows:
        for sim_combination in [
                                [simName_CT_0,batch_sim_labels[1][0],  test_reference],
                                [simName_CT_0,batch_sim_labels[2][0],  test_reference],
                                [simName_CT_0,batch_sim_labels[3][0],  test_reference],
                                [simName_CT_0,batch_sim_labels[4][0],  test_reference],
                                ]:
            print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}, {sim_combination[2]}]")

            groupA = list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_combination[0]].values())  # Replace with the first group of interest
            groupB = list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_combination[1]].values())  # Replace with the second group of interest
            
            U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
            # if p_mw<0.05:   significance = '*'
            # else:           significance = 'ns'
            significance = significance_symbol(p_mw)
            print(f"\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")


list_batch_sim_labels=[batch_sim_labels[source_folder_ind][0] for source_folder_ind, source_folder in enumerate(simulation_data_paths)]
y_lim   = {'TRN__pop':[0,0.71],     'VPM__pop':[0,0.71]}
y_ticks = {'TRN__pop':[0,0.25,0.5], 'VPM__pop':[0,0.25,0.5]}

fill_colors={'VPM__pop':'teal','TRN__pop':'mediumslateblue'}
angular_tuning_window='20ms'

plt.figure(figsize=(10,10))
# plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'font.size': 25,'mathtext.default': 'regular'})
for plot_pop_ind, plot_pop in enumerate(['TRN__pop','VPM__pop']):
    plt.subplot(2,1,plot_pop_ind+1)
    for conn_template_ind, conn_template in enumerate(conn_templates_set):
        
        median  = [np.median(    list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()))     for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        q25     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 25) for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        q75     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 75) for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        conds   = [sim_name for sim_name in store_RMSE_dict[plot_pop][angular_tuning_window].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        
        
        x_shift = (conn_template_ind*0.2)-0.15
        fill_x_pos = [x_val+x_shift for x_val in list(range(len(q25)))]

        plt.fill_between(fill_x_pos,  q25, q75,   color=fill_colors[plot_pop],   alpha=0.3+(0.3*conn_template_ind), edgecolor='k', linewidth=1)
        # plt.fill_between(fill_x_pos,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.2+(0.4*conn_template_ind), edgecolor='k', linewidth=1)
        
        x_val=-1
        for sim_name_ind, sim_name in enumerate(list_batch_sim_labels):
            
            if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) not in sim_name: continue
            else: x_val+=1
            
            x_pos = x_val+x_shift

            medianline_color = 'w'
            # alpha = 0.2+(0.3*conn_template_ind)

            alpha           = 0.6+(0.2*conn_template_ind)
            cmap            = mpl.colormaps[pop_colors[plot_pop]]
            cmap_mapping    = 0.2+(0.2*(sim_name_ind%5))                

            ax = plt.gca()
            boxes = ax.boxplot(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 
                            #    labels=[state_label], 
                            #    notch=True,
                                # medianprops=dict(color="k", alpha=0.7),
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                boxprops=dict(edgecolor='k', alpha=1.0),
                                patch_artist=True, 
                                positions=np.array([x_pos]),
                                widths=0.15,
                                )
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                # box.set_facecolor('k')
                box.set_alpha(1.0)
                box.set_facecolor(cmap(cmap_mapping))
                # box.set_facecolor(boxplot_colors[conn_template_ind])

    if   'TRN' in plot_pop: plt.ylim([0,0.66])
    elif 'VPM' in plot_pop: plt.ylim([0,0.21])
    else:                   plt.ylim([0,1.0])
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)        
    plt.xticks([0,1,2,3,4],[0,10,25,50,100])
    plt.xlabel('')

    plt.ylabel('Angular tuning RMSE')

plt.xticks([0,1,2,3,4],[0,10,25,50,100])
plt.xlabel('CT modulation (%)')
plt.tight_layout()
plt.savefig(save_figures_path+'CT_RMSE_byPop'+'_'+dataset_type+'.png',dpi=200)




plt.figure(figsize=(15,10))
# plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'font.size': 25,'mathtext.default': 'regular'})
plot_ind=0
for conn_template_ind, conn_template in enumerate(conn_templates_set):
    plt.subplot(1,len(conn_templates_set),plot_ind+1)
    plot_ind+=1
    
    for plot_pop_ind, plot_pop in enumerate(['TRN__pop','VPM__pop']):
        median  = [np.median(    list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()))     for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        q25     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 25) for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        q75     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 75) for sim_name in list_batch_sim_labels if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        conds   = [sim_name for sim_name in store_RMSE_dict[plot_pop][angular_tuning_window].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
        
        
        x_shift = (conn_template_ind*0.2)-0.15
        fill_x_pos = [x_val+x_shift for x_val in list(range(len(q25)))]

        plt.fill_between(fill_x_pos,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.3, edgecolor='none', linewidth=1)
        # plt.fill_between(fill_x_pos,  q25, q75,   color=fill_colors[plot_pop],   alpha=0.3+(0.3*conn_template_ind), edgecolor='k', linewidth=1)
        # plt.fill_between(fill_x_pos,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.2+(0.4*conn_template_ind), edgecolor='k', linewidth=1)
        
        x_val=-1
        for sim_name_ind, sim_name in enumerate(list_batch_sim_labels):
            
            if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) not in sim_name: continue
            else: x_val+=1
            
            x_pos = x_val+x_shift

            medianline_color = 'w'
            # alpha = 0.2+(0.3*conn_template_ind)

            alpha           = 0.6+(0.2*conn_template_ind)
            cmap            = mpl.colormaps[pop_colors[plot_pop]]
            cmap_mapping    = 0.2+(0.2*(sim_name_ind%5))                

            ax = plt.gca()
            boxes = ax.boxplot(list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values()), 
                            #    labels=[state_label], 
                            #    notch=True,
                                # medianprops=dict(color="k", alpha=0.7),
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                boxprops=dict(edgecolor='k', alpha=1.0),
                                patch_artist=True, 
                                positions=np.array([x_pos]),
                                widths=0.15,
                                )
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                # box.set_facecolor('k')
                box.set_alpha(1.0)
                box.set_facecolor(cmap(cmap_mapping))
                # box.set_facecolor(boxplot_colors[conn_template_ind])

        if   'TRN' in plot_pop: plt.ylim([0,0.66])
        elif 'VPM' in plot_pop: plt.ylim([0,0.21])
        else:                   plt.ylim([0,1.0])
        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)        
        plt.xticks([0,1,2,3,4],[0,10,25,50,100])
        plt.xlabel('')

        if conn_template_ind==0:plt.yticks([0.0,0.2,0.4,0.6])
        else:                   plt.yticks([])
        plt.ylim([0.0,0.65])

        if conn_template_ind==0:plt.ylabel('Angular tuning RMSE')
        if conn_template_ind==1:plt.xlabel('CT modulation (%)')
plt.xticks([0,1,2,3,4],[0,10,25,50,100])
plt.tight_layout()
plt.savefig(save_figures_path+'CT_RMSE2_byPop'+'_'+dataset_type+'.png',dpi=200)






print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t This test shows how CT modulation changes the RMSE of the angular tuning - see test explanation in code')
print('------------------------------------------------------------------------------------------------------------------------------------------------')

for plot_pop_ind, plot_pop in enumerate(plot_pops):
    print('\n\nplot pop: ',plot_pop)
    for angular_tuning_window in angulart_tuning_pooled_data_dict_absolute.keys():
        for ct_feedback in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].keys():
            for conn_template_ind, conn_template in enumerate(conn_templates_set):
                print('\t\tconn template',conn_template)

                test_data=[list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values())      for sim_name in store_RMSE_dict[plot_pop][angular_tuning_window].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
                
                for data_ind, data_combination in enumerate(test_data):
                    if data_ind == 0: continue

                    # print('\t\t\t\n',test_data[0],'\t',data_combination)

                    sim_combination = [test_data[0],data_combination,'two-sided']

                    U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                    
                    # if p_mw<0.05:   significance = '*'
                    # else:           significance = 'ns'
                    significance = significance_symbol(p_mw)
                    print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

                    '''
                    - This test shows how CT modulation changes the RMSE of the angular tuning
                    - The main result is: 
                        no significant change in RMSE for the same connectivity with increasing CT modulation
                        shows that CT modulation affects the error uniformly, even if it increases tuning or not (which will be shown later)
                    ___________________________________________________________________________________

                    Figure mapping to statistical tests:
                    
                        x   x   x   x   x       x = uniform
                        y   y   y   y   y       y = mixed
                        z   z   z   z   z       z = closed
                        0   10  25  50  100
                                CT mod
                    ___________________________________________________________________________________
                    
                    Test:
                        x_0 vs [x_10, x_25, x_50, x_100] - [for all x, y and z]
                        - varying CT modulation for the same angle/same connectivity

                    '''




########################################################################################################################################################################################################
angular_tuning_window='20ms'
angular_tuning_processed={}
for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
    print('\t\t - Population plotted: ', plot_pop)
    angular_tuning_processed.update({plot_pop:{}})
    try:    pop_color_map=pop_colors[plot_pop]
    except: pop_color_map=pop_colors['other']
    for angular_tuning_window in plot_angular_tuning_windows:

        angular_tuning_processed[plot_pop].update({angular_tuning_window:{}})

        print('\n - Angular tuning window: ', angular_tuning_window)
        # for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
        for sim_name_ind, sim_name in enumerate(list_batch_sim_labels):
            
            angular_tuning_processed[plot_pop][angular_tuning_window].update({sim_name:{}})

            print('\t\t\t - Sim name: ', sim_name)

            plot_means = []
            plot_stds  = []
            plot_x_vals = []
            plot_x_labels = []
            plot_medians=[]
            plot_q25s=[]
            plot_q75s=[]
            # for angle in angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys():
            for angle_ind, angle in enumerate(angles_sequence):
                plot_data = angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][sim_name][plot_pop][angle]

                plot_means.append(np.mean(plot_data))
                plot_stds.append(np.std(plot_data))
                plot_x_vals.append(angle_ind)
                plot_x_labels.append(int(angle))
                # adding non-parametric stats
                plot_medians.append(np.median(plot_data))
                plot_q25s.append(np.percentile(plot_data,25))
                plot_q75s.append(np.percentile(plot_data,75))

            plot_x_vals_array = np.array(plot_x_vals)            
            plot_means_array  = np.array(plot_means)/max(plot_means)
            plot_stds_array   = np.array(plot_stds)/max(plot_means)
            # adding non-parametric stats
            plot_medians_array= np.array(plot_medians)/max(plot_medians)
            plot_q25s         = np.array(plot_q25s)/max(plot_medians)
            plot_q75s         = np.array(plot_q75s)/max(plot_medians)

            angular_tuning_processed[plot_pop][angular_tuning_window][sim_name].update({'inds':plot_x_vals_array,'angles':angles_sequence,
                                                                                             'mean':plot_means_array,'std':plot_stds_array,
                                                                                             'median':plot_medians_array, 'q25':plot_q25s,'q75':plot_q75s,
                                                                                             })

########################################################################################################################################################################################################

def calculate_polygon_area(coords):
    """
    Calculate the area of a polygon given its (x, y) coordinates.

    Parameters:
        coords (list of tuples): List of (x, y) coordinates of the polygon vertices.
                                 The coordinates should be ordered sequentially around the polygon
                                 (either clockwise or counter-clockwise).
                                 The polygon can be open or closed.

    Returns:
        float: The absolute area of the polygon.
    """
    if len(coords) < 3:
        raise ValueError("A polygon must have at least 3 vertices.")
    
    # Ensure the polygon is closed by appending the first point to the end if needed
    if coords[0] != coords[-1]:
        coords.append(coords[0])

    # Apply the Shoelace formula
    n = len(coords)
    area = 0
    for i in range(n - 1):  # Iterate through vertices
        x1, y1 = coords[i]
        x2, y2 = coords[i + 1]
        area += x1 * y2 - x2 * y1

    return abs(area) / 2

import math
angular_tuning_window='20ms'
store_area = {}
store_area_std_up = {}
store_area_std_down = {}
# for conn_template in conn_templates_ordered:
for sim_name in list_batch_sim_labels:
    store_area.update({sim_name:{}})
    store_area_std_up.update({sim_name:{}})
    store_area_std_down.update({sim_name:{}})
    for plot_pop in angular_tuning_processed.keys():
        coords =[
            (   angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind] * math.cos(np.deg2rad(int(angle))),
                angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind] * math.sin(np.deg2rad(int(angle)))
            ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['angles'])
        ]
        polygon_area = calculate_polygon_area(coords)
        store_area[sim_name].update({plot_pop:polygon_area})

        # area of std values (up and down boundaries)
        coords_std_up =[
            (   (angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind]+angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['std'][angle_ind]) * math.cos(np.deg2rad(int(angle))),
                (angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind]+angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['std'][angle_ind]) * math.sin(np.deg2rad(int(angle)))
            ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['angles'])
        ]
        polygon_area_std_up = calculate_polygon_area(coords_std_up)
        store_area_std_up[sim_name].update({plot_pop:polygon_area_std_up})
        coords_std_down =[
            (   (angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind]-angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['std'][angle_ind]) * math.cos(np.deg2rad(int(angle))),
                (angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['mean'][angle_ind]-angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['std'][angle_ind]) * math.sin(np.deg2rad(int(angle)))
            ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][sim_name]['angles'])
        ]
        polygon_area_std_down = calculate_polygon_area(coords_std_down)
        store_area_std_down[sim_name].update({plot_pop:polygon_area_std_down})

area_ratio = {}
# for conn_template in conn_templates_ordered:
for sim_name in list_batch_sim_labels:
    area_ratio.update({sim_name:{}})
    for plot_pop in angular_tuning_processed.keys():
        if plot_pop=='MLe__pop':continue
        area_ratio[sim_name].update({plot_pop:store_area[sim_name][plot_pop]/store_area[sim_name]['MLe__pop']})
    # # ratio between TRN/VPM areas
    # area_ratio[sim_name].update({'TRN-VPM':store_area[sim_name]['TRN__pop']/store_area[sim_name]['VPM__pop']})

area_ratio_std_up = {}
# for conn_template in conn_templates_ordered:
for sim_name in list_batch_sim_labels:
    area_ratio_std_up.update({sim_name:{}})
    for plot_pop in angular_tuning_processed.keys():
        if plot_pop=='MLe__pop':continue
        area_ratio_std_up[sim_name].update({plot_pop:store_area_std_up[sim_name][plot_pop]/store_area[sim_name]['MLe__pop']}) # normalize to respective mean value of MLe area
area_ratio_std_down = {}
# for conn_template in conn_templates_ordered:
for sim_name in list_batch_sim_labels:
    area_ratio_std_down.update({sim_name:{}})
    for plot_pop in angular_tuning_processed.keys():
        if plot_pop=='MLe__pop':continue
        area_ratio_std_down[sim_name].update({plot_pop:store_area_std_down[sim_name][plot_pop]/store_area[sim_name]['MLe__pop']}) # normalize to respective mean value of MLe area


store_control_area={}
for plot_pop in angular_tuning_hartings.keys():
        coords =[
            (   angular_tuning_hartings[plot_pop][angle] * math.cos(np.deg2rad(int(angle))),
                angular_tuning_hartings[plot_pop][angle] * math.sin(np.deg2rad(int(angle)))
            ) for angle_ind, angle in enumerate(angular_tuning_hartings[plot_pop].keys())
        ]
        polygon_area = calculate_polygon_area(coords)
        store_control_area.update({plot_pop:polygon_area})
control_area_ratio={}
for plot_pop in store_control_area.keys():
    if plot_pop=='MLe__pop':continue
    control_area_ratio.update({plot_pop:store_control_area[plot_pop]/store_control_area['MLe__pop']})
# control_area_ratio.update({'TRN-VPM':store_control_area['TRN__pop']/store_control_area['VPM__pop']})






# Bar plot
plt.figure(figsize=(13,10))
plt.rcParams.update({'font.size': 25,'mathtext.default': 'regular'})
# plt.title('Area - Model vs Data')
plot_control=[1]
plot_control.extend(list(control_area_ratio.values()))

# Hartings Control reference lines
bar_colors=['k','g','b']
for line_ind,line_heigth in enumerate(plot_control):
    if line_ind==0:continue # skip plotting line for MLe
    plt.plot([-1.5,8.5],[line_heigth,line_heigth],'--k',linewidth=2,alpha=0.25,color=bar_colors[line_ind])

# Hartings Control bars
plt.bar([-2-0.375,-2-0.125,-2+0.125],plot_control,width=0.25,
        # color=['k','g','b'],
        edgecolor=['k','g','b'],
        linewidth=2,
        color=['w','w','w'],
        label=['MLe','VPM/MLe','TRN/MLe'])

# Reference bar for model
plt.bar([-0.375],[1],width=0.25,color=['k'],edgecolor='k', linewidth=1)

for conn_template_ind, conn_template in enumerate(conn_templates_set):
    for sim_name_ind, sim_name in enumerate(list_batch_sim_labels):

        # alpha = 0.6+(0.2*conn_template_ind)
        alpha = 1.0

        cmap_TRN            = mpl.colormaps[pop_colors['TRN__pop']]
        cmap_mapping_TRN    = 0.2+(0.2*(sim_name_ind%5))
        cmap_VPM            = mpl.colormaps[pop_colors['VPM__pop']]
        cmap_mapping_VPM    = 0.2+(0.2*(sim_name_ind%5))

        
        if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name:
            
            # x_tick_main = (2*(sim_name_ind%(len(list_batch_sim_labels)/len(conn_templates_set))))+2 # doubles the step, and adds 2 to shift everyone right
            # x_bar_width = 0.5
            # x_tick_shift = ((conn_template_ind*x_bar_width)/2)-x_bar_width-0.875

            # x_slide_bars_VPM = len(conn_templates_set)*x_bar_width/2
            # x_slide_bars_TRN = 0

            x_tick_main = 1.5*conn_template_ind
            # (3*(conn_template_ind%(len(conn_templates_set)/len(list_batch_sim_labels))))+5 # doubles the step, and adds 2 to shift everyone right
            x_bar_width = 0.5
            x_tick_shift = ((sim_name_ind*x_bar_width)/2)-x_bar_width-0.875

            # x_slide_bars_VPM = 0
            x_slide_bars_VPM = 1.5
            x_slide_bars_TRN = 2.75

            xtick_VPM = x_tick_main+x_tick_shift+x_slide_bars_VPM
            xtick_TRN = x_tick_main+x_tick_shift+x_slide_bars_TRN

            plt.bar([xtick_VPM,xtick_TRN],area_ratio[sim_name].values(),width=0.25,color=[cmap_VPM(cmap_mapping_VPM),cmap_TRN(cmap_mapping_TRN)], alpha=alpha, edgecolor='k', linewidth=1)

            plotNormalizedBarErrors=False # removed because it doesn't scale linearly - upper bounds result in a much bigger difference in the areas than the lower bounds - must use a linear metric instead
            if plotNormalizedBarErrors:
                plt.errorbar([sim_name_ind-0.125,sim_name_ind+0.125],list(area_ratio[sim_name].values()), yerr=list(area_ratio_std_up[sim_name].values()), 
                            lolims=True, 
                            uplims=False,
                            fmt='.', capsize=5,color='grey',
                            )
                plt.errorbar([sim_name_ind-0.125,sim_name_ind+0.125],list(area_ratio[sim_name].values()), yerr=list(area_ratio_std_down[sim_name].values()), 
                            uplims=True, 
                            lolims=False,
                            fmt='.', capsize=5,color='orange',
                            )
plt.ylim([-0.1,4.1])
plt.yticks([0,1,2,3])
plt.xlim([-3,11])
plt.xticks([-2.125, -0.375, 1.25, 4.0, 6.75], ['Exp.', 'MLe', '0.0', '0.5', '1.0'])
plt.xlabel('CSI')
plt.ylabel('Area ratio - norm. to MLe')
ax1 = plt.gca()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.tight_layout()
plt.savefig(save_figures_path+'CT_barPlot'+'_'+dataset_type+'.png',dpi=200)



########################################################################################################################################################################################################

# --- Plotting firing - absolute angular tuning
# Create a recursive generator function to yield all values in a nested dictionary
def flatten(d):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from flatten(v)
        else:
            yield v

try:    max_val = max(max(flatten(angulart_tuning_pooled_data_dict_absolute)))+200
except: max_val = 2000

# list_angles         = ['135', '180', '225', '270', '315', '0', '45', '90']
# list_angles_labels  = ['-3π/4', '-π/2', '-π/4', 'Target', 'π/4', 'π/2', '3π/4', 'π']
list_angles         = ['270', '90']
list_angles_labels  = ['270 (Target)', '90 (Opposite)']

angular_tuning_window = '20ms'
ct_feedback = 'ct_uniform'
# plt.figure(figsize=(15,10))
# plt.rcParams.update({'font.size': 15})
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.figure(figsize=(20,10))
    plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})

    for conn_template_ind, conn_template in enumerate(conn_templates_set):
        subplot_index = (conn_template_ind+1)
        plt.subplot(1,len(conn_templates_set),subplot_index)
        for ang_ind, ang in enumerate(list_angles):
        
            median  = [np.median(    list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]))      for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
            q25     = [np.percentile(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]), 25)  for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
            q75     = [np.percentile(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]), 75)  for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
            conds   = [sim_name for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]

            x_shift = (ang_ind*0.2)-0.10
            # x_shift = (conn_template_ind*0.2)-0.15
            fill_x_pos = [x_val+x_shift for x_val in list(range(len(q25)))]

            plt.fill_between(fill_x_pos,  q25, q75,   color='k',   alpha=1-(0.3+(0.4*ang_ind)), edgecolor='none', linewidth=1)
            # plt.fill_between(fill_x_pos,  q25, q75,   color=popColors_dict[plot_pop],   alpha=1-(0.3+(0.4*ang_ind)), edgecolor='k', linewidth=1)
            
            # plt.plot(fill_x_pos, median, color='k', alpha=0.2+(0.2*conn_template_ind), linewidth = 3)
            # plt.plot(fill_x_pos, median, color='k', alpha=0.2+(0.2*conn_template_ind), linewidth = 3, marker = 'o')

            x_val=-1
            for sim_name_ind, sim_name in enumerate(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()):
                
                if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) not in sim_name: continue
                else: x_val+=1
                
                x_pos = x_val+x_shift

                medianline_color = 'w'
                alpha = 0.2+(0.3*conn_template_ind)

                cmap = mpl.colormaps[pop_colors[plot_pop]]
                cmap_mapping = 0.2+(0.2*(sim_name_ind%5))
                # c=cmap(cmap_mapping)

                ax = plt.gca()
                boxes = ax.boxplot(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]), 
                                #    labels=[state_label], 
                                #    notch=True,
                                    # medianprops=dict(color="k", alpha=0.7),
                                    medianprops=dict(color=medianline_color, alpha=1.0), 
                                    boxprops=dict(edgecolor='k', alpha=1.0),
                                    patch_artist=True, 
                                    positions=np.array([x_pos]),
                                    widths=0.15,
                                    )
                # Set the face color of each box
                for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                    # box.set_facecolor('k')
                    box.set_alpha(1.0)
                    box.set_facecolor(cmap(cmap_mapping))
                    # box.set_facecolor(boxplot_colors[conn_template_ind])

        if   'TRN' in plot_pop: 
            plt.ylim([0,625])
            plt.yticks([0,200,400,600])
        elif 'VPM' in plot_pop: 
            plt.ylim([0,1550])
            plt.yticks([0,500,1000,1500])
        else:                   
            plt.ylim([0,1750])
        if conn_template_ind==0: plt.ylabel('Number of spikes')
        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)        
        if conn_template_ind!=0:  plt.yticks([])
        plt.xticks([0,1,2,3,4],[0,10,25,50,100])
        plt.xlabel('CT modulation (%)')
    plt.tight_layout()
    plt.savefig(save_figures_path+'CT_spikingActivity_'+str(plot_pop_ind)+'_'+plot_pop+'_'+dataset_type+'.png',dpi=200)

########################################################################################################################################################################################################


print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t Table 4iii and S2iii - This test shows how CT modulation changes the firing for each angle - see test explanation in code')
print('------------------------------------------------------------------------------------------------------------------------------------------------')
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    print('\n\nplot pop: ',plot_pop)
    for ang_ind, ang in enumerate(list_angles):
        print('\n\tangle: ',ang)
        for angular_tuning_window in angulart_tuning_pooled_data_dict_absolute.keys():
            for ct_feedback in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].keys():
                for conn_template_ind, conn_template in enumerate(conn_templates_set):
                    print('\t\tconn template',conn_template)

                    test_data=[angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]      for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) in sim_name]
                    
                    for data_ind, data_combination in enumerate(test_data):
                        if data_ind == 0: continue

                        print('\t\t\t',test_data[0],'\t',data_combination)

                        sim_combination = [test_data[0],data_combination,'two-sided']

                        U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                        
                        # if p_mw<0.05:   significance = '*'
                        # else:           significance = 'ns'
                        significance = significance_symbol(p_mw)
                        print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

                        '''
                        - This test shows how CT modulation changes the firing for each angle
                        - The main result is: 
                            for 270 deg, the change is apparent for all networks
                            at 90 deg the full closed loop network shows no difference, but the others do. 
                                This can be interpreted as if the activity in the target direction increases, but in the opposite direction it doesn't
                                Allows us to claim that the network is becoming more tuned to the angle
                        
                        ___________________________________________________________________________________

                        Figure mapping to statistical tests:
                                                Angles
                        135     180     225     270     315     0       45      90
                        xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx       x = uniform
                        yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy       y = mixed
                        zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz       z = closed
                        0-100   0-100   0-100   0-100   0-100   0-100   0-100   0-100
                                                CT mod
                        ___________________________________________________________________________________
                        
                        Test:
                            x_0 vs [x_10, x_25, x_50, x_100] - [for all x, y and z]
                            - varying CT modulation for the same angle/same connectivity

                        '''

print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t Table 4i and S2i - This test shows how CT modulation changes the firing rate for different angles in the closed loop and mixed connectivity - see test explanation in code')
print('------------------------------------------------------------------------------------------------------------------------------------------------')
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    print('\n\nplot pop:',plot_pop)
    # for ang_ind, ang in enumerate(list_angles):
    #     print('\n\tangle: ',ang)
    for angular_tuning_window in angulart_tuning_pooled_data_dict_absolute.keys():
        for ct_feedback in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].keys():
            for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys():
                print('\n\t\tsim name: ',sim_name)
                
                # for conn_template_ind, conn_template in enumerate(conn_templates_set):
                #     print('\t\t\tconn template',conn_template)
                #     # for ang_ind, ang in enumerate(list_angles):

                test_data_ang=[angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]      for ang_ind, ang in enumerate(list_angles)]
                # test_data_ang=[angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]      for conn_template_ind, conn_template in enumerate(conn_templates_set)]
                
                for data_ind, data_combination in enumerate(test_data_ang):
                    if data_ind == 0: continue
                    test_ang = list_angles[data_ind]
                    # test_ang = list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys())[data_ind]
                    print('\t\t\t',test_ang, ' vs ',list_angles[0] ,'\t',test_data_ang[0],'\t',data_combination)

                    sim_combination = [test_data_ang[0],data_combination,'two-sided']

                    U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                    
                    # if p_mw<0.05:   significance = '*'
                    # else:           significance = 'ns'
                    significance = significance_symbol(p_mw)
                    print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

                    '''
                    - This test shows how CT modulation changes the firing rate for different angles in the closed loop and mixed connectivity
                    - The main result is: 
                        for closed loop and mixed, the data is significantly different across angles, for the same CT modulation
                        TRN data for uniform connectivity is not significantly different across all angles, for all CT modulations
                            shows that closed loop and mixed conns are modulated by cortex, but uniform connections are not
                        results are less obivious in VPM, since it is less sensitive to modulation because it has fixed MLe inputs
                                            
                    ___________________________________________________________________________________

                    Figure mapping to statistical tests:
                                            Angles
                    135     180     225     270     315     0       45      90
                    xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx       x = uniform
                    yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy       y = mixed
                    zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz       z = closed
                    0-100   0-100   0-100   0-100   0-100   0-100   0-100   0-100
                                            CT mod
                    ___________________________________________________________________________________
                    
                    Test:
                    x_270 vs [x_135 ... x_90] - [for all x, y and z]
                    - varying angle for the same modulation/same connectivity

                    '''

print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t Table 4ii and S2ii - This test shows how different connectivity shapes the response at different angles, for a given CT modulation - see test explanation in code')
print('------------------------------------------------------------------------------------------------------------------------------------------------')
conn_templates_set_reduced=['uniform','mixed','closed']
modulation_values=['0','10','25','50','100']
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    print('\n\nplot_pop: ',plot_pop)
    for ang_ind, ang in enumerate(list_angles):
        print('\n\tangle: ',ang)
        for angular_tuning_window in angulart_tuning_pooled_data_dict_absolute.keys():
            for ct_feedback in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].keys():
                for modulation_value in modulation_values:
                    print('\n\t\tCT modulation: ',modulation_value)
                    test_data_modulation=[angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][ang]      for sim_name in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys() if int(sim_name.split('_')[-1]) == int(modulation_value)]
                        
                    for data_ind, data_combination in enumerate(test_data_modulation):
                        if data_ind == 0: continue
                        print('\t\t\t',conn_templates_set_reduced[0],' vs ',conn_templates_set_reduced[data_ind], '\t',test_data_modulation[0],'\t',data_combination)

                        sim_combination = [test_data_modulation[0],data_combination,'two-sided']

                        U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                        
                        # if p_mw<0.05:   significance = '*'
                        # else:           significance = 'ns'
                        significance = significance_symbol(p_mw)
                        print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

                        '''
                        - This test shows how different connectivity shapes the response at different angles, for a given CT modulation
                        - The main result is: 
                            network response overlaps at lateral angles, but are significantly opposite at target and opposite angle
                            TRN data showcases that difference more explicitly, but VPM also shows the same effect
                        
                            ___________________________________________________________________________________

                            Figure mapping to statistical tests:
                                                    Angles
                            135     180     225     270     315     0       45      90
                            xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx   xxxxx       x = uniform
                            yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy   yyyyy       y = mixed
                            zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz   zzzzz       z = closed
                            0-100   0-100   0-100   0-100   0-100   0-100   0-100   0-100
                                                    CT mod
                            ___________________________________________________________________________________
                            
                            Test:
                            x_uniform vs [x_mixed, x_closed] - [for all x, y and z]
                            - varying connectivity for the same modulation/same angle

                            '''

########################################################################################################################################################################################################
# --- Plotting angular tuning maps for the Normalized dataset for visualization in the paper figure

# conn_templates_set = [  'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform', 'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50', 'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed']
# angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][batch_sim_labels[source_folder_ind][0]][plot_pop][angle]

ct_feedback = 'ct_uniform'
angular_tuning_window='20ms'
plot_pops=['TRN__pop','VPM__pop']
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    for conn_template_ind, conn_template in enumerate(conn_templates_set):
        plt.figure(figsize=(25,25))
        plt.rcParams.update({'font.size': 80})
        plt.subplot(1,1,1,projection='polar')
        line_width=1
        for sim_name_ind, sim_name in enumerate(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].keys()):
            if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) not in sim_name: continue

            line_width+=1

            # keys
            angles = [int(angle) for angle in mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys()]
            mean_values = [mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][sim_name][plot_pop][angle]['mean'] for angle in mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys()]

            angles.append(angles[0])
            mean_values.append(mean_values[0])
            
            angles_=np.array(angles)
            mean_values_=np.array(mean_values)

            alpha = 0.2+(0.3*conn_template_ind)
            cmap = mpl.colormaps[pop_colors[plot_pop]]
            cmap_mapping = 0.2+(0.2*(sim_name_ind%5))

            plt.plot(np.deg2rad(angles_),mean_values_,c=cmap(cmap_mapping), alpha=1.0,linewidth=5)
            # plt.plot(np.deg2rad(angles_),mean_values_,c=popColors_dict[plot_pop], alpha=alpha,linewidth=line_width)

        plt.yticks([0.5, 1.0])
        plt.ylim([0,1.5])
        plt.xticks(np.deg2rad([0,90,180,270]))
        # plt.grid(False)

        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine

        plt.savefig(save_figures_path+'CT_angularTuningPlots_'+str(plot_pop_ind)+'_'+plot_pop+'_'+str(conn_template_ind)+'_'+conn_template+'_'+dataset_type+'.png',dpi=200)

# -- Absolute angular tuning maps
y_lim_absolute      ={'VPM__pop':[0,2000],'TRN__pop':[0,400]}
y_lim_absolute_ticks={'VPM__pop':[1000,2000],'TRN__pop':[200,400]}
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    for conn_template_ind, conn_template in enumerate(conn_templates_set):
        plt.figure(figsize=(25,25))
        plt.rcParams.update({'font.size': 80})
        plt.subplot(1,1,1,projection='polar')
        line_width=1
        for sim_name_ind, sim_name in enumerate(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()):
            if 'barreloid_batch_fig07'+analyze_figure+'_'+str(conn_template_ind) not in sim_name: continue

            line_width+=1

            # keys
            angles = [int(angle) for angle in mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys()]
            mean_values = [mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop][angle]['mean_notNormalized'] for angle in mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][sim_name][plot_pop].keys()]

            angles.append(angles[0])
            mean_values.append(mean_values[0])
            
            angles_=np.array(angles)
            mean_values_=np.array(mean_values)

            alpha = 0.2+(0.3*conn_template_ind)
            cmap = mpl.colormaps[pop_colors[plot_pop]]
            cmap_mapping = 0.2+(0.2*(sim_name_ind%5))

            plt.plot(np.deg2rad(angles_),mean_values_,c=cmap(cmap_mapping), alpha=1.0,linewidth=5)
            # plt.plot(np.deg2rad(angles_),mean_values_,c=popColors_dict[plot_pop], alpha=alpha,linewidth=line_width)

        plt.ylim(y_lim_absolute[plot_pop])
        plt.yticks(y_lim_absolute_ticks[plot_pop])
        plt.xticks(np.deg2rad([0,90,180,270]))
        # plt.grid(False)

        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine
        plt.savefig(save_figures_path+'CT_angularTuningPlots_absolute_'+str(plot_pop_ind)+'_'+plot_pop+'_'+str(conn_template_ind)+'_'+conn_template+'_absolute.png',dpi=200)





########################################################################################################################################################################################################

plt.show()

########################################################################################################################################################################################################


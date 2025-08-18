'''
script replacing AnalyzeAngularTuningStats.py
'''

import os
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

from cProfile import label
import json
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

########################################################################################################################################################################################################
def significance_symbol(p_value):
    """
    Returns a significance symbol based on the p-value.
    Parameters: p_value (float): The p-value from a statistical test.
    Returns:    str: A string representing significance level ('ns', '*', '**', '***', '****').
    """
    if (p_value >= 0.05) or np.isnan(p_value):
        return 'ns'
    elif p_value < 0.0001:
        return '****'
    elif p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    else:
        return '*'

########################################################################################################################################################################################################

data_path = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data'
# batch_sim_group = 'paper_sims'
# batch_sim_label = 'paperSim_Fig05_awake'
# sim_output_folder = 'sim_output'

# batch_sim_group = 'paper_000'
# batch_sim_group = 'paper_001'
batch_sim_group = 'paper_002'
batch_sim_label = batch_sim_group+'_'+'barreloid_batch_fig05_awake'
sim_output_folder = 'sim_output'

simulation_data_path = data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/'+sim_output_folder+'/'

########################################################################################################################################################################################################

# --- Path to save analysis figures
save_figures_path = data_path+'/'+batch_sim_group+'/'+batch_sim_label+'/analysis_figures/'
if not os.path.exists(save_figures_path): os.makedirs(save_figures_path)

########################################################################################################################################################################################################

_addVersion='_v04'
# conn_templates=[ 
#                     'VPM_ff_closed__TRN_fb_open__L6A_fb_closed',
#                     'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25',
#                     'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50',
#                     'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75',
#                     'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
#                     'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75',
#                     'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
#                     'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25',
#                     'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',
#                     ]

conn_templates=[    
                    'VPM_ff_closed__TRN_fb_open__L6A_fb_uniform',
                    'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_75_25',
                    'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_50_50',
                    'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_25_75',
                    'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                    'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_25_75',
                    'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',
                    'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_75_25',
                    'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',
                    ]


try: conn_templates = [c+_addVersion for c in conn_templates]
except: print('Adding conn version failed')

#######################################################################

angular_tuning_windows  = ['20ms']
ct_feedbacks            = ['ct_uniform']

# dataset_type = 'absolute' # do not use it here, just 'normalized'
dataset_type = 'normalized'

angulart_tuning_pooled_data_dict={}
mean_angulart_tuning_pooled_data_dict={}

angulart_tuning_pooled_data_dict_absolute={}
mean_angulart_tuning_pooled_data_dict_absolute={}

for angular_tuning_window in angular_tuning_windows:
    ############################################################################################################################################
    if angular_tuning_window == '20ms':
        file_flag = 'angularTuning3_spikes_ON'
        dataFolder      = simulation_data_path+'angular_tuning_ON/'
    ############################################################################################################################################
    elif angular_tuning_window == '250ms':
        file_flag = 'angularTuning3_spikes_all'
        dataFolder      = simulation_data_path+'angular_tuning_all/'
    ############################################################################################################################################
    angulart_tuning_pooled_data_dict.update({angular_tuning_window:{}})
    mean_angulart_tuning_pooled_data_dict.update({angular_tuning_window:{}})
    angulart_tuning_pooled_data_dict_absolute.update({angular_tuning_window:{}})
    mean_angulart_tuning_pooled_data_dict_absolute.update({angular_tuning_window:{}})

    for ct_feedback in ct_feedbacks:
        angulart_tuning_pooled_data_dict[angular_tuning_window].update({ct_feedback:{}})
        mean_angulart_tuning_pooled_data_dict[angular_tuning_window].update({ct_feedback:{}})
        angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].update({ct_feedback:{}})
        mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].update({ct_feedback:{}})
        for conn_template in conn_templates:
            angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].update({conn_template:{}})
            mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].update({conn_template:{}})

            filename = dataFolder+'angular_tuning_Statistics__'+dataset_type+'__'+angular_tuning_window+'__'+ct_feedback+'__conn_'+conn_template+'.json'
            with open(filename, 'r') as file: ang_tun_data = json.load(file)

            angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template]=ang_tun_data

            for pop in ang_tun_data:
                mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template].update({pop:{}})
                max_mean = max([np.mean(ang_tun_data[pop][angle]) for angle in ang_tun_data[pop].keys()])
                print(conn_template, pop, max_mean)
                for angle in ang_tun_data[pop].keys():
                    mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][pop].update({angle:{'mean':np.mean(ang_tun_data[pop][angle])/max_mean,'std':np.std(ang_tun_data[pop][angle])/max_mean}})
            
            # absolute data
            angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].update({conn_template:{}})
            mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].update({conn_template:{}})
            filename_absolute = dataFolder+'angular_tuning_Statistics__absolute__'+angular_tuning_window+'__'+ct_feedback+'__conn_'+conn_template+'.json'
            with open(filename_absolute, 'r') as file: ang_tun_data_absolute = json.load(file)

            angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template]=ang_tun_data_absolute

            for pop in ang_tun_data_absolute:
                mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template].update({pop:{}})
                max_mean = max([np.mean(ang_tun_data_absolute[pop][angle]) for angle in ang_tun_data_absolute[pop].keys()])
                print(conn_template, pop, max_mean)
                for angle in ang_tun_data_absolute[pop].keys():
                    mean_angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][pop].update({
                            angle:{
                                'mean':np.mean(ang_tun_data_absolute[pop][angle])/max_mean,
                                'std':np.std(ang_tun_data_absolute[pop][angle])/max_mean,

                                'mean_notNormalized':np.mean(ang_tun_data_absolute[pop][angle]),
                                'std_notNormalized':np.std(ang_tun_data_absolute[pop][angle]),
                                }
                            }
                        )

# mean_angulart_tuning_pooled_data_dict={}

# for angular_tuning_window in mean_angulart_tuning_pooled_data_dict.keys():
#     for ct_feedback in mean_angulart_tuning_pooled_data_dict[angular_tuning_window].keys():
#         for conn_template in mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback].keys():

# mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][test_pop]

# VPM vs MLe and TRN vs MLe
for angular_tuning_window in angular_tuning_windows:
    print('\n - Angular tuning window: ', angular_tuning_window)
    for ct_feedback in ct_feedbacks:
        print('\t - CT feedback: ', ct_feedback)
        for conn_template in conn_templates:
            print('\t\t - Conn template: ', conn_template)
            for test_pop in ['VPM__pop','TRN__pop']:
                print('\t\t\t - Population tested: MLe__pop vs ', test_pop)
                for angle in angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template]['MLe__pop'].keys():

                    data1 = angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template]['MLe__pop'][angle]
                    data2 = angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][test_pop][angle]
                    if data1==data2:
                        print('\t\t\t', angle, '\tvalues are the same')
                    else:
                        statistics = stats.ttest_ind(data1,data2)
                        if      statistics[1]>0.05:                                 significance = 'ns'
                        elif    (statistics[1]<=0.05)  and (statistics[1]>0.005):   significance = '*'
                        elif    (statistics[1]<=0.005) and (statistics[1]>0.0005):  significance = '**'
                        elif    (statistics[1]<=0.0005):                            significance = '***'


                        print('\t\t\t',angle, '\t',statistics[1], '\t',significance)

####################################################################################################################################################################################
# Plot open traces for angular tuning, grouped by pop

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
    'VPM__pop':'summer',
    'TRN__pop':'winter',
    # 'MLe__pop':'Greys',
    # 'VPM__pop':'Greens',
    # 'TRN__pop':'Blues',
    'CTvirtual_uniform__pop':'Reds',
    'other': 'Purples',
            }
# pop_colors={
#     'MLe__pop':'k',
#     'VPM__pop':'g',
#     'TRN__pop':'b',
#     'CTvirtual_uniform__pop':'r',
#     'other': 'cyan',
#             }

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
# plot_angular_tuning_windows  = ['20ms']
print('\t - CT feedback: ', ct_feedback)

# conn_templates_ordered = [ 'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',

#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_10_90',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_20_80',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_30_70',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_40_60',
                            
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_60_40',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_70_30',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_80_20',
#                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_90_10',

#                             'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',
#                             ]; _addVersion='_v03'
# try: conn_templates_ordered = [c+_addVersion for c in conn_templates_ordered]
# except: print('Adding conn version failed')

conn_templates_ordered = conn_templates

########################################################################################################################
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
        for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
            store_RMSE_dict[plot_pop][angular_tuning_window].update({conn_template:{}})
            store_MAE_dict[plot_pop][angular_tuning_window].update({conn_template:{}})
            # print(conn_template)
            for angle_ind, angle in enumerate(angles_sequence):
                if remove_270deg:
                    if angle=='270':continue
                # angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle]
                # angular_tuning_hartings_meanVals[plot_pop][angle]
                # 1) RMSE 
                rmse = np.sqrt(np.mean((np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle]) - np.array(angular_tuning_hartings_meanVals[plot_pop][angle]))**2)) 
                # 2) MAE (manual version)
                mae = np.mean(np.abs(np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle]) - np.array(angular_tuning_hartings_meanVals[plot_pop][angle])))
                
                # print(angle, '\tRMSE',rmse, '\tMAE',mae)
                
                store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].update({angle:rmse})
                store_MAE_dict[plot_pop][angular_tuning_window][conn_template].update({angle:mae})

# plt.figure(figsize=(10,10))
# for plot_pop_ind, plot_pop in enumerate(['VPM__pop','TRN__pop']):
#     for angular_tuning_window in plot_angular_tuning_windows:
#         for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
#             plt.bar(     conn_template_ind+((plot_pop_ind*0.2)-0.1), np.mean(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values())),width=0.2,color=popColors_dict[plot_pop])
#             plt.errorbar(conn_template_ind+((plot_pop_ind*0.2)-0.1), np.mean(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values())), np.std(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values())), color='grey')


# plt.figure(figsize=(10,10))
# for plot_pop_ind, plot_pop in enumerate(['VPM__pop','TRN__pop']):
#     for angular_tuning_window in plot_angular_tuning_windows:
#         for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
#             plt.bar(     conn_template_ind+((plot_pop_ind*0.2)-0.1), np.mean(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values())),width=0.2,color=popColors_dict[plot_pop])
#             plt.errorbar(conn_template_ind+((plot_pop_ind*0.2)-0.1), np.mean(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values())), np.std(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values())), color='grey')


# # --- MAE
# conn_template_uniform='VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform_v04'
# for pop_ind,plot_pop in enumerate(['VPM__pop','TRN__pop']):

#     if   plot_pop=='TRN__pop': test_reference='two-sided'
#     elif plot_pop=='VPM__pop': test_reference='two-sided'
#     else:                      test_reference='two-sided'
#     # if   plot_pop=='TRN__pop': test_reference='greater'
#     # elif plot_pop=='VPM__pop': test_reference='less'
#     # else:                      test_reference='greater'

#     paiwise_test_reference = 'two-sided'

#     print(f"\nPlot pop: {plot_pop}")
#     for angular_tuning_window in plot_angular_tuning_windows:
#         for sim_combination in [
#                                 [conn_template_uniform,'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed_v04',                   test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25_v04',       test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50_v04',       test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75_v04',       test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75_v04',   test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50_v04',   test_reference],
#                                 [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25_v04',   test_reference],
#                                 [conn_template_uniform,'VPM_ff_closed__TRN_fb_open__L6A_fb_closed_v04',                     test_reference],

#                                 ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed_v04',             'VPM_ff_closed__TRN_fb_open__L6A_fb_closed_v04',                    paiwise_test_reference],
#                                 ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25_v04',  paiwise_test_reference],
#                                 ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50_v04',  paiwise_test_reference],
#                                 ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75_v04',  paiwise_test_reference],
#                                 ]:
#             print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}, {sim_combination[2]}]")

#             groupA = list(store_MAE_dict[plot_pop][angular_tuning_window][sim_combination[0]].values())  # Replace with the first group of interest
#             groupB = list(store_MAE_dict[plot_pop][angular_tuning_window][sim_combination[1]].values())  # Replace with the second group of interest
            
#             U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
#             if p_mw<0.05:   significance = '*'
#             else:           significance = 'ns'
#             print(f"\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

# plt.figure(figsize=(7,10))
# plt.rcParams.update({'font.size': 15})
# for plot_pop_ind, plot_pop in enumerate(['VPM__pop','TRN__pop']):
#     for angular_tuning_window in plot_angular_tuning_windows:

#         q25     = [np.percentile(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 25) for conn_template in conn_templates_ordered]
#         q75     = [np.percentile(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 75) for conn_template in conn_templates_ordered]
#         conds   = conn_templates_ordered
#         plt.fill_between(conds,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.3)

#         for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
#             medianline_color = 'w'
#             alpha = 0.8

#             ax = plt.gca()
#             boxes = ax.boxplot(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 
#                             #    labels=[state_label], 
#                             #    notch=True,
#                                 # medianprops=dict(color="k", alpha=0.7),
#                                 medianprops=dict(color=medianline_color, alpha=1.0), 
#                                 patch_artist=True, 
#                                 positions=np.array([conn_template_ind+((plot_pop_ind*0.2)-0.1)]),
#                                 widths=0.15,
#                                 )
#             # Set the face color of each box
#             for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
#                 box.set_facecolor(color)
#                 box.set_alpha(alpha)
#             plt.xticks([])

# plt.ylim([0,0.51])
# plt.yticks([0,0.25,0.5])
# ax1=plt.gca()
# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['bottom'].set_visible(False)
# ax1.spines['left'].set_visible(False)
# plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# # plt.xlabel('$P_{T}$')
# plt.xlabel('CSI')
# plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_MAE_all'+'_'+dataset_type+'.png',dpi=200)

# y_lim   = {'TRN__pop':[0,0.51],'VPM__pop':[0,0.26]}
# y_ticks = {'TRN__pop':[0,0.25,0.5],'VPM__pop':[0,0.1,0.2]}
# plt.figure(figsize=(7,10))
# plt.rcParams.update({'font.size': 15})
# for plot_pop_ind, plot_pop in enumerate(['TRN__pop','VPM__pop']):
#     plt.subplot(2,1,plot_pop_ind+1)
#     for angular_tuning_window in plot_angular_tuning_windows:

#         q25     = [np.percentile(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 25) for conn_template in conn_templates_ordered]
#         q75     = [np.percentile(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 75) for conn_template in conn_templates_ordered]
#         conds   = conn_templates_ordered
#         plt.fill_between(conds,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.3)

#         for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
#             medianline_color = 'w'
#             alpha = 0.8

#             ax = plt.gca()
#             boxes = ax.boxplot(list(store_MAE_dict[plot_pop][angular_tuning_window][conn_template].values()), 
#                             #    labels=[state_label], 
#                             #    notch=True,
#                                 # medianprops=dict(color="k", alpha=0.7),
#                                 medianprops=dict(color=medianline_color, alpha=1.0), 
#                                 patch_artist=True, 
#                                 positions=np.array([conn_template_ind]),
#                                 widths=0.15,
#                                 )
#             # Set the face color of each box
#             for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
#                 box.set_facecolor(color)
#                 box.set_alpha(alpha)
#             plt.xticks([])

#     plt.ylim(y_lim[plot_pop])
#     plt.yticks(y_ticks[plot_pop])
#     # plt.ylim([0,0.51])
#     # plt.yticks([0,0.25,0.5])
#     ax1=plt.gca()
#     ax1.spines['top'].set_visible(False)
#     ax1.spines['right'].set_visible(False)
#     ax1.spines['bottom'].set_visible(False)
#     ax1.spines['left'].set_visible(False)
# plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# # plt.xlabel('$P_{T}$')
# plt.xlabel('CSI')
# plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_MAE_byPop'+'_'+dataset_type+'.png',dpi=200)


# --- RMSE - Currently being used

print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t Paper - Figure 5G - Comparing RMSE statistics === ATTENTION, THESE RESULTS ARE REVERSED - START CLOSED- AND PROGRESSES TO OPEN-LOOP')
print('------------------------------------------------------------------------------------------------------------------------------------------------')

conn_template_uniform='VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform_v04'
for pop_ind,plot_pop in enumerate(['TRN__pop','VPM__pop']):

    if   plot_pop=='TRN__pop': test_reference='two-sided'
    elif plot_pop=='VPM__pop': test_reference='two-sided'
    else:                      test_reference='two-sided'
    # if   plot_pop=='TRN__pop': test_reference='greater'
    # elif plot_pop=='VPM__pop': test_reference='less'
    # else:                      test_reference='greater'

    paiwise_test_reference = 'two-sided'

    print(f"\nPlot pop: {plot_pop}")
    for angular_tuning_window in plot_angular_tuning_windows:
        # for sim_combination in [
        #                         [conn_template_uniform,'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed_v04',                   test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25_v04',       test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50_v04',       test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75_v04',       test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75_v04',   test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50_v04',   test_reference],
        #                         [conn_template_uniform,'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25_v04',   test_reference],
        #                         [conn_template_uniform,'VPM_ff_closed__TRN_fb_open__L6A_fb_closed_v04',                     test_reference],

        #                         ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed_v04',             'VPM_ff_closed__TRN_fb_open__L6A_fb_closed_v04',                    paiwise_test_reference],
        #                         ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25_v04',  paiwise_test_reference],
        #                         ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50_v04',  paiwise_test_reference],
        #                         ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75_v04', 'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75_v04',  paiwise_test_reference],
        #                         ]:
        for sim_combination in [
                                [conn_template_uniform,   conn_templates[8],   test_reference],
                                [conn_template_uniform,   conn_templates[7],   test_reference],
                                [conn_template_uniform,   conn_templates[6],   test_reference],
                                [conn_template_uniform,   conn_templates[5],   test_reference],
                                [conn_template_uniform,   conn_templates[3],   test_reference],
                                [conn_template_uniform,   conn_templates[2],   test_reference],
                                [conn_template_uniform,   conn_templates[1],   test_reference],
                                [conn_template_uniform,   conn_templates[0],   test_reference],

                                [conn_templates[8], conn_templates[0],  paiwise_test_reference],
                                [conn_templates[7], conn_templates[1],  paiwise_test_reference],
                                [conn_templates[6], conn_templates[2],  paiwise_test_reference],
                                [conn_templates[5], conn_templates[3],  paiwise_test_reference],
                                ]:
            print(f"\tsim combination: [{sim_combination[0]}, {sim_combination[1]}, {sim_combination[2]}]")

            groupA = list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_combination[0]].values())  # Replace with the first group of interest
            groupB = list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_combination[1]].values())  # Replace with the second group of interest
            
            U, p_mw = stats.mannwhitneyu(groupA, groupB, alternative=sim_combination[2])
            # if p_mw<0.05:   significance = '*'
            # else:           significance = 'ns'
            significance = significance_symbol(p_mw)
            print(f"\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

plt.figure(figsize=(7,10))
plt.rcParams.update({'font.size': 15})
for plot_pop_ind, plot_pop in enumerate(['VPM__pop','TRN__pop']):
    for angular_tuning_window in plot_angular_tuning_windows:

        q25     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 25) for conn_template in conn_templates_ordered]
        q75     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 75) for conn_template in conn_templates_ordered]
        conds   = conn_templates_ordered
        plt.fill_between(conds,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.3)

        for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
            medianline_color = 'w'
            alpha = 0.8

            ax = plt.gca()
            boxes = ax.boxplot(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 
                            #    labels=[state_label], 
                            #    notch=True,
                                # medianprops=dict(color="k", alpha=0.7),
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([conn_template_ind+((plot_pop_ind*0.2)-0.1)]),
                                widths=0.15,
                                )
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                box.set_facecolor(color)
                box.set_alpha(alpha)
            plt.xticks([])

plt.ylim([0,0.51])
plt.yticks([0,0.25,0.5])
ax1=plt.gca()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
# plt.xlabel('$P_{T}$')
plt.xlabel('CSI')
plt.ylabel('RMSE')
plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_RMSE_all'+'_'+dataset_type+'.png',dpi=200)

y_lim   = {'TRN__pop':[0,0.61],'VPM__pop':[0,0.33]}
y_ticks = {'TRN__pop':[0,0.25,0.5],'VPM__pop':[0,0.15,0.30]}
plt.figure(figsize=(10,10))
plt.rcParams.update({'font.size': 20})
for plot_pop_ind, plot_pop in enumerate(['TRN__pop','VPM__pop']):
    plt.subplot(2,1,plot_pop_ind+1)
    for angular_tuning_window in plot_angular_tuning_windows:

        q25     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 25) for conn_template in conn_templates_ordered]
        q75     = [np.percentile(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 75) for conn_template in conn_templates_ordered]
        conds   = conn_templates_ordered
        plt.fill_between(conds,  q25, q75,   color=popColors_dict[plot_pop],   alpha=0.3, edgecolor='none', linewidth=1)

        for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
            medianline_color = 'w'
            alpha = 0.8

            ax = plt.gca()
            boxes = ax.boxplot(list(store_RMSE_dict[plot_pop][angular_tuning_window][conn_template].values()), 
                            #    labels=[state_label], 
                            #    notch=True,
                                # medianprops=dict(color="k", alpha=0.7),
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                patch_artist=True, 
                                positions=np.array([conn_template_ind]),
                                widths=0.15,
                                )
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                box.set_facecolor(color)
                box.set_alpha(alpha)
            plt.xticks([])

    plt.ylim(y_lim[plot_pop])
    plt.yticks(y_ticks[plot_pop])
    # plt.ylim([0,0.51])
    # plt.yticks([0,0.25,0.5])
    ax1=plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('RMSE')
plt.tight_layout()
plt.savefig(save_figures_path+'_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_RMSE_byPop'+'_'+dataset_type+'.png',dpi=200)
# plt.show()

########################################################################################################################

angular_tuning_processed={}

for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
    print('\t\t - Population plotted: ', plot_pop)

    angular_tuning_processed.update({plot_pop:{}})

    try:    pop_color_map=pop_colors[plot_pop]
    except: pop_color_map=pop_colors['other']
    # try:    pop_color=pop_colors[plot_pop]
    # except: pop_color=pop_colors['other']

    legends=[]
    legends_Rc=[]

    for angular_tuning_window in plot_angular_tuning_windows:

        angular_tuning_processed[plot_pop].update({angular_tuning_window:{}})

        print('\n - Angular tuning window: ', angular_tuning_window)
        plt.figure(figsize=(15,15))
        plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
        for conn_template_ind, conn_template in enumerate(conn_templates_ordered):
            
            angular_tuning_processed[plot_pop][angular_tuning_window].update({conn_template:{}})

            print('\t\t\t - Conn template: ', conn_template)

            if '_tr_' in conn_template:
                c_fraction=conn_template.split('_tr_')[1].split('_')[0]
                u_fraction=conn_template.split('_tr_')[1].split('_')[1]
                Rcu = c_fraction+':'+u_fraction
                Rc  = str(int(c_fraction)/(int(c_fraction)+int(u_fraction)))
            elif 'uniform' in conn_template:
                Rcu = '100:0'
                Rc  = '0.0'
            elif 'closed' in conn_template:
                Rcu = '0:100'
                Rc  = '1.0'
            else:
                Rcu = 'unknown'
                Rc  = 'unknown'
            
            legends.append(Rcu)
            legends_Rc.append(Rc)

            plot_means = []
            plot_stds  = []
            plot_x_vals = []
            plot_x_labels = []
            plot_medians=[]
            plot_q25s=[]
            plot_q75s=[]
            # for angle in angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop].keys():
            for angle_ind, angle in enumerate(angles_sequence):
                plot_data = angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle]

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

            angular_tuning_processed[plot_pop][angular_tuning_window][conn_template].update({'inds':plot_x_vals_array,'angles':angles_sequence,
                                                                                             'mean':plot_means_array,'std':plot_stds_array,
                                                                                             'median':plot_medians_array, 'q25':plot_q25s,'q75':plot_q75s,
                                                                                             })

            cmap = mpl.colormaps[pop_color_map]

            # Take colors at regular intervals spanning the colormap.
            pop_color = cmap((7+conn_template_ind)/(len(conn_templates_ordered)+10))

            # Plot the mean line
            # linewidth = 3
            linewidth = 2+(conn_template_ind*0.25)
            
            plt.plot(plot_x_vals_array, plot_means_array, color=pop_color, linewidth = linewidth)
            
            plt.errorbar(plot_x_vals, plot_means_array, plot_stds_array, color=pop_color, alpha=0.1, capsize=10, elinewidth=3)
            # # Fill the area between the mean ± standard deviation
            # plt.fill_between(plot_x_vals, plot_means_array - plot_stds_array, plot_means_array + plot_stds_array, color=pop_color, alpha=0.1)
            
        # Reference values from Hartings (MLe, VPM and TRN)
        plot_x_ref_vals = []
        plot_y_ref_vals = []
        for angle_ind, angle in enumerate(angles_sequence):
            plot_x_ref_vals.append(angle_ind)
            plot_y_ref_vals.append(angular_tuning_hartings[plot_pop][angle])
        plt.plot(plot_x_ref_vals,plot_y_ref_vals, color='grey', linewidth=10)
        legends.append('control')
        legends_Rc.append('control')
        
        # Mean reference values from Hartings (MLe, VPM and TRN)
        plot_x_ref_vals = []
        plot_y_ref_vals_mean = []
        for angle_ind, angle in enumerate(angles_sequence):
            plot_x_ref_vals.append(angle_ind)
            plot_y_ref_vals_mean.append(angular_tuning_hartings_meanVals[plot_pop][angle])
        plt.plot(plot_x_ref_vals,plot_y_ref_vals_mean, color='red', linewidth=2)
        legends.append('control_mean')
        legends_Rc.append('control_mean')
        
        # plt.legend(legends,loc='upper right')
        plt.legend(legends_Rc,loc='upper right')

        plt.ylim([-0.1,1.6])
        plt.yticks([-0,0.5,1.0,1.5])

        plt.xticks(plot_x_vals_array, ['','-π/2','','0','','π/2','','π'])
        # plt.xticks(plot_x_vals_array, plot_x_labels)
        plt.xlabel('Angle (deg)')
        plt.ylabel('Angular tuning')

        # plt.title('Plot with Shaded Standard Deviation')
        
        ax1 = plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_test_plot_linear_angular_tuning_'+plot_pop+'_'+dataset_type+'.png')
        except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_test_plot_linear_angular_tuning_'+plot_pop+'_'+dataset_type+'.png')

# import sys; sys.exit()

# ####################################################################################################################################################################################
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
# angular_tuning_window='20ms'
for angular_tuning_window in plot_angular_tuning_windows:
    store_area = {}
    store_area_std_up = {}
    store_area_std_down = {}
    for conn_template in conn_templates_ordered:
        store_area.update({conn_template:{}})
        store_area_std_up.update({conn_template:{}})
        store_area_std_down.update({conn_template:{}})
        for plot_pop in angular_tuning_processed.keys():
            coords =[
                (   angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind] * math.cos(np.deg2rad(int(angle))),
                    angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind] * math.sin(np.deg2rad(int(angle)))
                ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles'])
            ]
            polygon_area = calculate_polygon_area(coords)
            store_area[conn_template].update({plot_pop:polygon_area})

            # area of std values (up and down boundaries)
            coords_std_up =[
                (   (angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]+angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std'][angle_ind]) * math.cos(np.deg2rad(int(angle))),
                    (angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]+angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std'][angle_ind]) * math.sin(np.deg2rad(int(angle)))
                ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles'])
            ]
            polygon_area_std_up = calculate_polygon_area(coords_std_up)
            store_area_std_up[conn_template].update({plot_pop:polygon_area_std_up})
            coords_std_down =[
                (   (angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]-angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std'][angle_ind]) * math.cos(np.deg2rad(int(angle))),
                    (angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]-angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std'][angle_ind]) * math.sin(np.deg2rad(int(angle)))
                ) for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles'])
            ]
            polygon_area_std_down = calculate_polygon_area(coords_std_down)
            store_area_std_down[conn_template].update({plot_pop:polygon_area_std_down})

    area_ratio = {}
    for conn_template in conn_templates_ordered:
        area_ratio.update({conn_template:{}})
        for plot_pop in angular_tuning_processed.keys():
            if plot_pop=='MLe__pop':continue
            area_ratio[conn_template].update({plot_pop:store_area[conn_template][plot_pop]/store_area[conn_template]['MLe__pop']})
        # # ratio between TRN/VPM areas
        # area_ratio[conn_template].update({'TRN-VPM':store_area[conn_template]['TRN__pop']/store_area[conn_template]['VPM__pop']})

    area_ratio_std_up = {}
    for conn_template in conn_templates_ordered:
        area_ratio_std_up.update({conn_template:{}})
        for plot_pop in angular_tuning_processed.keys():
            if plot_pop=='MLe__pop':continue
            area_ratio_std_up[conn_template].update({plot_pop:store_area_std_up[conn_template][plot_pop]/store_area[conn_template]['MLe__pop']}) # normalize to respective mean value of MLe area
    area_ratio_std_down = {}
    for conn_template in conn_templates_ordered:
        area_ratio_std_down.update({conn_template:{}})
        for plot_pop in angular_tuning_processed.keys():
            if plot_pop=='MLe__pop':continue
            area_ratio_std_down[conn_template].update({plot_pop:store_area_std_down[conn_template][plot_pop]/store_area[conn_template]['MLe__pop']}) # normalize to respective mean value of MLe area


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
    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 20,'mathtext.default': 'regular'})
    # plt.title('Area - Model vs Data')
    plot_control=[1]
    plot_control.extend(list(control_area_ratio.values()))

    # Hartings Control reference lines
    bar_colors=['k','g','b']
    for line_ind,line_heigth in enumerate(plot_control):
        if line_ind==0:continue # skip plotting line for MLe
        plt.plot([-1.5,8],[line_heigth,line_heigth],'--k',linewidth=2,alpha=0.25,color=bar_colors[line_ind])

    # Hartings Control bars
    plt.bar([-2-0.375,-2-0.125,-2+0.125],plot_control,width=0.25,
            # color=['k','g','b'],
            edgecolor=['k','g','b'],
            linewidth=2,
            color=['w','w','w'],
            label=['MLe','VPM/MLe','TRN/MLe'])
    
    # Reference bar for model
    plt.bar([-0.375],[1],width=0.25,color=['k'])
    
    for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
        plt.bar([conn_template_ind-0.125,conn_template_ind+0.125],area_ratio[conn_template].values(),width=0.25,color=['g','b'])

        plotNormalizedBarErrors=False # removed because it doesn't scale linearly - upper bounds result in a much bigger difference in the areas than the lower bounds - must use a linear metric instead
        if plotNormalizedBarErrors:
            plt.errorbar([conn_template_ind-0.125,conn_template_ind+0.125],list(area_ratio[conn_template].values()), yerr=list(area_ratio_std_up[conn_template].values()), 
                        lolims=True, 
                        uplims=False,
                        fmt='.', capsize=5,color='grey',
                        )
            plt.errorbar([conn_template_ind-0.125,conn_template_ind+0.125],list(area_ratio[conn_template].values()), yerr=list(area_ratio_std_down[conn_template].values()), 
                        uplims=True, 
                        lolims=False,
                        fmt='.', capsize=5,color='orange',
                        )
    
    # plt.legend()
    plt.xlim([-3,len(conn_templates_ordered)])
    # plt.xticks([-2, 0, 5, 10, 15, 20], ['Control', -1.0, -0.5, 0.0, 0.5, 1.0])
    plt.xticks([-2, 0, 2, 4, 6, 8], ['Exper.', -1.0, -0.5, 0.0, 0.5, 1.0])
    # plt.xticks([-2, 0, 5, 10], ['Control', 0, 0.5, 1.0])
    plt.ylim([-0.1,4.1])
    plt.yticks([0,1,2,3,4])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('Area ratio - norm. to MLe')
    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    plt.tight_layout()
    try:    plt.savefig(save_figures_path+'_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'.png',dpi=200)
    except: plt.savefig('_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'.png')

    # Diff bar plot
    plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
    # plt.title('Area Ratio - Model vs Data')
    plt.bar([-2-0.125,-2+0.125],[0,0],width=0.25,color=['g','b'],label=['diff VPM/MLe','diff TRN/MLe'])
    for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
        area_ratio_array = np.array(list(area_ratio[conn_template].values()))
        control_area_ratio_array=np.array(list(control_area_ratio.values()))
        plot_diff = area_ratio_array-control_area_ratio_array
        plt.bar([conn_template_ind-0.125,conn_template_ind+0.125],plot_diff,width=0.25,color=['g','b'])
    # plt.legend()
    plt.xlim([-3,len(conn_templates_ordered)])
    plt.ylim([-1.1,2])
    plt.yticks([-1.0,0,1.0,2.0])
    # plt.xticks([-2, 0, 5, 10, 15, 20], ['', -1.0, -0.5, 0.0, 0.5, 1.0])
    # plt.xticks([-2, 0, 2, 4, 6, 8], ['Control', -1.0, -0.5, 0.0, 0.5, 1.0])
    plt.xticks([-2, 0, 2, 4, 6, 8], ['Exper.', -1.0, -0.5, 0.0, 0.5, 1.0])
    # plt.xticks([-2, 0, 5, 10], ['', 0, 0.5, 1.0])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('Area ratio - norm. to MLe')
    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_diff.png')
    except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_diff.png')





    # Absolute Difference Line plot
    plot_colors=['g','b']

    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
    # plt.title('Area Ratio - Model vs Data')
    plt.bar([-2-0.125,-2+0.125],[0,0],width=0.25,color=['g','b'],label=['abs diff VPM/MLe','abs diff TRN/MLe'])

    plt.plot([-1,9],[0,0],'--k',linewidth=2, alpha=0.25)

    for pop_diff_i, pop_diff in enumerate(['diff VPM/MLe','diff TRN/MLe']):
        pop_diff_values=[]
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            pop_diff_values.append(abs(np.array(list(area_ratio[conn_template].values()))[pop_diff_i]-np.array(list(control_area_ratio.values()))[pop_diff_i]))
        plt.plot(pop_diff_values,color=plot_colors[pop_diff_i],linewidth=2)

    for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
        area_ratio_array = np.array(list(area_ratio[conn_template].values()))
        control_area_ratio_array=np.array(list(control_area_ratio.values()))
        plot_diff = abs(area_ratio_array-control_area_ratio_array)
        # if (conn_template_ind==2) or (conn_template_ind==6):    markerfacecolor = 'lightgrey'
        # else:                                                   markerfacecolor = 'w'
        markerfacecolor = 'w'
        for i in range(len(plot_diff)):
            plt.plot(conn_template_ind,plot_diff[i],marker='d',markerfacecolor=markerfacecolor,markeredgecolor=plot_colors[i],markeredgewidth=2,markersize=10)
        
    # plt.legend()
    plt.xlim([-3,len(conn_templates_ordered)])
    plt.ylim([-0.1,4.1])
    plt.yticks([0.0,1.0,2.0])
    # plt.xticks([0, 5, 10, 15, 20], [-1.0, -0.5, 0.0, 0.5, 1.0])
    plt.xticks([0, 2, 4, 6, 8], [-1.0, -0.5, 0.0, 0.5, 1.0])
    # plt.xticks([0, 5, 10], [0, 0.5, 1.0])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('Difference in area ratio')
    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    plt.tight_layout()
    # plt.grid(True)
    try:    plt.savefig(save_figures_path+'_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_diff_Line.png',dpi=200)
    except: plt.savefig('_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_diff_Line.png')




    # Pair-wise distance and Mean Pair-wise distance plots
    plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
    # plt.title('Mean Pairwise Distance - Model vs Data')
    plt.bar([-2-0.125,-2+0.125,-2+0.375],[0,0,0],width=0.25,color=['g','b','purple'],label=['MLe','VPM','TRN'])
    for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
        mean_vals_diff={}
        for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
            mean_vals_diff.update({plot_pop:[]})
            for angle_ind, angle in enumerate(angles_sequence):
                mean_vals_diff[plot_pop].append(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]-angular_tuning_hartings[plot_pop][angle])

        plot_mean_vals_diff=[list(mean_vals_diff[plot_pop]) for plot_pop in mean_vals_diff.keys()]

        plt.bar([conn_template_ind-0.375,conn_template_ind-0.125,conn_template_ind+0.125],np.mean(plot_mean_vals_diff,1),width=0.25,color=['k','g','b'])
    plt.xlim([-3,len(conn_templates_ordered)])
    plt.ylim([-0.26,0.26])
    plt.yticks([-0.2,0,0.2])
    # plt.xticks([-2, 0, 5, 10], ['Control', 0, 0.5, 1.0])
    # plt.xticks([-2, 0, 2, 4], ['Control', 0.0, 0.5, 1.0])
    plt.xticks([-2, 0, 2, 4], ['Exper.', 0.0, 0.5, 1.0])
    # plt.xlabel('$P_{T}$')
    plt.xlabel('CSI')
    plt.ylabel('Area ratio - norm. to MLe')
    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_meanVals_diff.png')
    except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_test_bar_plot_angTun_v2_'+'_'+dataset_type+'_meanVals_diff.png')


    print('Plotting normalized and non-normalized angular tuning')
    # --- Keep debugging angular tuning figure
    plt.figure(figsize=(20,25))
    plt.rcParams.update({'font.size': 80})
    plt.subplot(1,1,1,projection='polar')
    for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
        # keys
        polar_list_keys=[int(ang) for ang in angular_tuning_hartings[plot_pop]]
        polar_list_keys.extend([polar_list_keys[0]])
        # 
        control_values = list(angular_tuning_hartings[plot_pop].values())
        control_values.extend([control_values[0]])
        try:    plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=0.5,linewidth=10)
        except: plt.plot(np.deg2rad(polar_list_keys),control_values)
    plt.suptitle('Angular tuning plot - Control')
    plt.yticks([0.5, 1.0])
    plt.ylim([0,1.5])
    plt.xticks(np.deg2rad([0,90,180,270]))
    # plt.grid(False)

    ax1 = plt.gca()
    ax1.spines['polar'].set_visible(False)  # Hide the circular spine
    try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_test_angularTuning_Hartings_'+dataset_type+'.png')
    except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_test_angularTuning_Hartings_'+dataset_type+'.png')

    print('Plotting normalized and non-normalized angular tuning')
    # --- Keep debugging angular tuning figure
    plt.figure(figsize=(20,25))
    plt.rcParams.update({'font.size': 80})
    plt.subplot(1,1,1,projection='polar')
    for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
        if   'MLe' in plot_pop: plot_marker='^'
        elif 'VPM' in plot_pop: plot_marker='s'
        elif 'TRN' in plot_pop: plot_marker='D'
        else:                   plot_marker='o'

        # keys
        polar_list_keys=[int(ang) for ang in angular_tuning_hartings[plot_pop]]
        polar_list_keys.extend([polar_list_keys[0]])
        # 
        control_values = list(angular_tuning_hartings[plot_pop].values())
        control_values.extend([control_values[0]])
        try:    
            plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=0.5,linewidth=10)
            plt.plot(np.deg2rad(polar_list_keys),control_values,plot_marker,c=popColors_dict[plot_pop], alpha=0.8,markersize=25)
        except: 
            plt.plot(np.deg2rad(polar_list_keys),control_values)
            plt.plot(np.deg2rad(polar_list_keys),control_values,plot_marker,markersize=25)
    plt.yticks([0.5, 1.0])
    plt.ylim([0,1.5])
    plt.xticks(np.deg2rad([0,90,180,270]))
    # plt.grid(False)

    ax1 = plt.gca()
    ax1.spines['polar'].set_visible(False)  # Hide the circular spine
    try:    plt.savefig(save_figures_path+'_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_angularTuning_Hartings_'+dataset_type+'_paperFig.png')
    except: plt.savefig('_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_test_angularTuning_Hartings_'+dataset_type+'_paperFig.png')


    for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
        print('\n\nplot pop: ',plot_pop)
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            print('\t\tconn template',conn_template)
            for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles']):
                mean_val = angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'][angle_ind]
                std_val  = angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std'][angle_ind]
                ang_tun_val = angular_tuning_hartings[plot_pop][angle]
                if mean_val==ang_tun_val:
                    print('\t\t\t', angle, '\tvalues are the same')
                else:
                    statistic,pval = stats.ttest_1samp(mean_val, ang_tun_val)
                    symbol = significance_symbol(pval)
                    print('\t\t\t', angle, '\t\t\t\t\t\t\t test val: ', ang_tun_val,'\t | mean:', mean_val, ' - std: ', std_val, '\t')
                    print(f"\t\t\t|\t{symbol} \t| ttest_1samp statistic: {statistic}, p-value: {pval}\n")
        print('\n')

    # One-sample t-test
    for plot_pop in ['MLe__pop']:
        print('\n\nplot pop: ',plot_pop)
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            if conn_template_ind>0:continue # all the same, so only run once
            print('\t\tconn template',conn_template)
            for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles']):
                if (angle==270) or (angle=='270'): 
                    print(f"\t\t\t|\t{angle} \t|\t{'ns'} \t| same value\n")
                    continue
                model_outputs = np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle])
                target_mean = angular_tuning_hartings[plot_pop][angle]
                # 1) One-sample t-test
                # H0: The mean of model_outputs is equal to target_mean
                # HA: The mean of model_outputs differs from target_mean
                t_stat, p_val_t = stats.ttest_1samp(model_outputs, popmean=target_mean)
                t_symbol = significance_symbol(p_val_t)
                print(f"\t\t\t|\t{angle} \t|\t{t_symbol} \t| One-sample t-test: {t_stat}, p-value: {p_val_t}\n")
    print('\n')
    # PAPER STATS - Wilcoxon signed-rank test
    store_ns=0; store_s=0
    store_ns_exp=0; store_s_exp=0
    for plot_pop in ['MLe__pop']:
        print('\n\nplot pop: ',plot_pop)
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            if conn_template_ind>0:continue # all the same, so only run once
            print('\t\tconn template',conn_template)
            for angle_ind, angle in enumerate(angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['angles']):
                if (angle==270) or (angle=='270'): 
                    print(f"\t\t\t|\t{angle} \t|\t{'ns'} \t| same value\n")
                    continue
                model_outputs = np.array(angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][plot_pop][angle])
                target_mean = angular_tuning_hartings[plot_pop][angle]
                # 2) Wilcoxon signed-rank test
                # H0: The median of model_outputs is equal to target_mean
                # HA: The median of model_outputs differs from target_mean
                # We compare each data point against the target_mean by subtracting
                # and passing that to wilcoxon().
                w_stat, p_val_w = stats.wilcoxon(model_outputs - target_mean)
                w_symbol = significance_symbol(p_val_w)
                print(f"\t\t\t|\t{angle} \t|\t{w_symbol} \t| Wilcoxon signed-rank test: {w_stat}, p-value: {p_val_w}\n")

                
                if w_symbol=='ns':
                    store_ns+=np.mean(model_outputs)
                    store_ns_exp+=target_mean
                else:
                    store_s+=np.mean(model_outputs)
                    store_s_exp+=target_mean
        print('\n')
        print('ns spiking |\t goal: ', (store_ns_exp*100)/(store_ns_exp+store_s_exp),'%\t|\t observed: ', (store_ns*100)/(store_ns+store_s),'%\t|\t difference: ', abs(((store_ns*100)/(store_ns+store_s))-((store_ns_exp*100)/(store_ns_exp+store_s_exp))),'%')
        print('s  spiking |\t goal: ', (store_s_exp*100)/(store_ns_exp+store_s_exp), '%\t|\t observed: ', (store_s*100)/(store_ns+store_s), '%\t|\t difference: ', abs(((store_s*100)/(store_ns+store_s))-((store_s_exp*100)/(store_ns_exp+store_s_exp))),'%')
        
        '''
        plot pop:  MLe__pop
			|	135 	|	** 	| Wilcoxon signed-rank test: 0.0, p-value: 0.0078125
			|	180 	|	** 	| Wilcoxon signed-rank test: 0.0, p-value: 0.0078125
			|	225 	|	ns 	| Wilcoxon signed-rank test: 6.0, p-value: 0.109375
			|	270 	|	ns 	| same value
			|	315 	|	ns 	| Wilcoxon signed-rank test: 8.0, p-value: 0.1953125
			|	0 	    |	ns 	| Wilcoxon signed-rank test: 15.0, p-value: 0.7421875
			|	45 	    |	ns 	| Wilcoxon signed-rank test: 12.0, p-value: 0.4609375
			|	90 	    |	* 	| Wilcoxon signed-rank test: 2.0, p-value: 0.0234375

            not significant vs significant spiking - showing that the area that don't match the statistics only accounts for 3.5% of spikes, so can be ignored
            ns spiking |	 goal:  72.94949101285448 %	    |	 observed:  69.47646191126958 %	|	 difference:  3.473029101584899 %
            s  spiking |	 goal:  27.050508987145523 %	|	 observed:  30.52353808873042 %	|	 difference:  3.4730291015848955 %
        '''

    # MLe stats figure
    for plot_pop in ['MLe__pop']:
        try:    pop_color_map=pop_colors[plot_pop]
        except: pop_color_map=pop_colors['other']
        
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            
            if (plot_pop=='MLe__pop') and conn_template_ind>=1: continue # plot MLe only once, since it doesn't change across templates

            plt.figure(figsize=(20,25))
            plt.rcParams.update({'font.size': 80})
            plt.subplot(1,1,1,projection='polar')

            # keys
            polar_list_keys=[int(ang) for ang in angular_tuning_hartings[plot_pop]]
            polar_list_keys.extend([polar_list_keys[0]])
            # 
            control_values = list(angular_tuning_hartings[plot_pop].values())
            control_values.extend([control_values[0]])
            # plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=0.5,linewidth=10)
            plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=1.0,linewidth=3)
            plt.plot(np.deg2rad(polar_list_keys),control_values,'^',markersize=25,c=popColors_dict[plot_pop], alpha=1.0,linewidth=3)

            cmap = mpl.colormaps[pop_color_map]
            pop_color = cmap((7+conn_template_ind)/(len(conn_templates_ordered)+10))
            
            # plot_inds=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['inds']
            # plot_inds=np.append(plot_inds, plot_inds[0])
            plot_inds=np.deg2rad([135, 180, 225, 270, 315, 0, 45, 90, 135])
            
            plot_ang_tun=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean']
            plot_ang_tun=np.append(plot_ang_tun,plot_ang_tun[0])
            
            plot_ang_tun_hi=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'] + angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std']
            plot_ang_tun_hi=np.append(plot_ang_tun_hi,plot_ang_tun_hi[0])

            plot_ang_tun_lo=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'] - angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std']
            plot_ang_tun_lo=np.append(plot_ang_tun_lo,plot_ang_tun_lo[0])

            plt.fill_between(plot_inds, plot_ang_tun_hi, plot_ang_tun_lo, 
                            color='magenta', alpha=0.2)
            plt.plot(plot_inds,plot_ang_tun,linewidth=4,color='magenta', alpha=0.7)

            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            # plt.grid(False)

            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine

            try:    plt.savefig(save_figures_path+'_PaperFig_'+batch_sim_label+'_'+angular_tuning_window+'_angularTuning_individual_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')
            except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_angularTuning_individual_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')
        

    # Individual figures with Angular tuning measured and control
    for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
        try:    pop_color_map=pop_colors[plot_pop]
        except: pop_color_map=pop_colors['other']
        
        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            
            if (plot_pop=='MLe__pop') and conn_template_ind>=1: continue # plot MLe only once, since it doesn't change across templates

            plt.figure(figsize=(20,25))
            plt.rcParams.update({'font.size': 80})
            plt.subplot(1,1,1,projection='polar')

            # keys
            polar_list_keys=[int(ang) for ang in angular_tuning_hartings[plot_pop]]
            polar_list_keys.extend([polar_list_keys[0]])
            # 
            control_values = list(angular_tuning_hartings[plot_pop].values())
            control_values.extend([control_values[0]])
            plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=0.5,linewidth=10)

            cmap = mpl.colormaps[pop_color_map]
            pop_color = cmap((7+conn_template_ind)/(len(conn_templates_ordered)+10))
            
            # plot_inds=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['inds']
            # plot_inds=np.append(plot_inds, plot_inds[0])
            plot_inds=np.deg2rad([135, 180, 225, 270, 315, 0, 45, 90, 135])
            
            plot_ang_tun=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean']
            plot_ang_tun=np.append(plot_ang_tun,plot_ang_tun[0])
            
            plot_ang_tun_hi=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'] + angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std']
            plot_ang_tun_hi=np.append(plot_ang_tun_hi,plot_ang_tun_hi[0])

            plot_ang_tun_lo=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean'] - angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['std']
            plot_ang_tun_lo=np.append(plot_ang_tun_lo,plot_ang_tun_lo[0])

            plt.fill_between(plot_inds, plot_ang_tun_hi, plot_ang_tun_lo, 
                            color=popColors_dict[plot_pop], alpha=0.1)
            plt.plot(plot_inds,plot_ang_tun,linewidth=5,color=popColors_dict[plot_pop])

            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            # plt.grid(False)

            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine

            try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_angularTuning_individual_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')
            except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_angularTuning_individual_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')
        

    # Angular tuning grouped by pop
    for plot_pop in ['MLe__pop','VPM__pop','TRN__pop']:
        
        plt.figure(figsize=(20,25))
        plt.rcParams.update({'font.size': 80})
        plt.subplot(1,1,1,projection='polar')

        # keys
        polar_list_keys=[int(ang) for ang in angular_tuning_hartings[plot_pop]]
        polar_list_keys.extend([polar_list_keys[0]])
        # 
        control_values = list(angular_tuning_hartings[plot_pop].values())
        control_values.extend([control_values[0]])
        plt.plot(np.deg2rad(polar_list_keys),control_values,c='grey', alpha=0.25,linewidth=10)
        # plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[plot_pop], alpha=0.5,linewidth=10)

        for conn_template_ind,conn_template in enumerate(conn_templates_ordered):
            
            if (plot_pop=='MLe__pop') and conn_template_ind>=1: continue # plot MLe only once, since it doesn't change across templates


            try:    pop_color_map=pop_colors[plot_pop]
            except: pop_color_map=pop_colors['other']


            cmap = mpl.colormaps[pop_color_map]
            pop_color = cmap((7+conn_template_ind)/(len(conn_templates_ordered)+10))
            
            # plot_inds=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['inds']
            # plot_inds=np.append(plot_inds, plot_inds[0])
            plot_inds=np.deg2rad([135, 180, 225, 270, 315, 0, 45, 90, 135])
            
            plot_ang_tun=angular_tuning_processed[plot_pop][angular_tuning_window][conn_template]['mean']
            plot_ang_tun=np.append(plot_ang_tun,plot_ang_tun[0])
            
            plt.plot(plot_inds,plot_ang_tun,linewidth=2,color=pop_color)
            # plt.fill_between(plot_inds, plot_ang_tun_hi, plot_ang_tun_lo, 
            #                  color=popColors_dict[plot_pop], alpha=0.1)

        plt.yticks([0.5, 1.0])
        plt.ylim([0,1.5])
        plt.xticks(np.deg2rad([0,90,180,270]))
        # plt.grid(False)

        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine

        try:    plt.savefig(save_figures_path+batch_sim_label+'_'+angular_tuning_window+'_angularTuning_grouped_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')
        except: plt.savefig(batch_sim_label+'_'+angular_tuning_window+'_angularTuning_grouped_'+plot_pop+'__'+str(conn_template_ind)+'_'+dataset_type+'.png')


# import sys; sys.exit()



####################################################################################################################################################################################

plot_pops = ['TRN__pop', 'VPM__pop']
ct_feedback = 'ct_uniform'

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

plt.figure(figsize=(10,10))
# plt.rcParams.update({'font.size': 20,'mathtext.default': 'regular'})
plt.rcParams.update({'font.size': 20,'mathtext.default': 'regular'})
for plot_pop_ind, plot_pop in enumerate(plot_pops):
    plt.subplot(len(plot_pops),1,plot_pop_ind+1)
    for ang_ind, ang in enumerate(list_angles):
        median  = [np.median(    list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop][ang]))      for conn_template in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()]
        q25     = [np.percentile(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop][ang]), 25)  for conn_template in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()]
        q75     = [np.percentile(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop][ang]), 75)  for conn_template in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()]
        conds   = list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys())

        x_shift = (ang_ind*0.2)-0.10
        fill_x_pos = [x_val+x_shift for x_val in list(range(len(q25)))]
        plt.fill_between(fill_x_pos,  q25, q75,   color='k',   alpha=1-(0.3+(0.4*ang_ind)), edgecolor='none', linewidth=1)
        x_val=-1
        for conn_template_ind, conn_template in enumerate(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()):
            x_val+=1
            x_pos = x_val+x_shift

            medianline_color = 'w'
            # alpha = 0.2+(0.3*conn_template_ind)
            alpha = 0.8

            cmap = mpl.colormaps[pop_colors[plot_pop]]
            cmap_mapping = 0.2+(0.2*(conn_template_ind%5))
            # c=cmap(cmap_mapping)
            ax = plt.gca()
            boxes = ax.boxplot(list(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop][ang]), 
                                medianprops=dict(color=medianline_color, alpha=1.0), 
                                boxprops=dict(edgecolor='k', alpha=1.0),
                                patch_artist=True, 
                                positions=np.array([x_pos]),
                                widths=0.15,
                                )
            # Set the face color of each box
            for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
                box.set_facecolor(color)
                box.set_alpha(alpha)
            # # Set the face color of each box
            # for box, color in zip(boxes['boxes'], [popColors_dict[plot_pop]]):
            #     # box.set_facecolor('k')
            #     box.set_alpha(1.0)
            #     box.set_facecolor(cmap(cmap_mapping))
            #     # box.set_facecolor(boxplot_colors[conn_template_ind])

        if   'TRN' in plot_pop: 
            plt.ylim([0,330])
            plt.yticks([0,75,150,225,300])
        elif 'VPM' in plot_pop: 
            plt.ylim([0,1100])
            plt.yticks([0,250,500,750,1000])

        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)        
        plt.xticks(list(range(9)),[(-1) + i*((1) - (-1))/(9 - 1) for i in range(9)])
        # plt.xlabel('$P_{T}$')
        plt.xlabel('CSI')
        plt.ylabel('Number of spikes')
        # if plot_pop_ind==0: plt.ylabel('Number of spikes')

plt.tight_layout()
plt.savefig(save_figures_path+'_PaperFig_'+'Awake_spikingActivity_distribution.png',dpi=200)

print('\n\n------------------------------------------------------------------------------------------------------------------------------------------------')
print('\t Paper - Figure 5H - Comparing spiking in the awake condition for each connectivity schema - 270 vs 90 deg - This test shows that topoogy increases tuning')
print('------------------------------------------------------------------------------------------------------------------------------------------------')

for plot_pop_ind, plot_pop in enumerate(plot_pops):
    print('\n\nplot pop: ',plot_pop)
    for angular_tuning_window in angulart_tuning_pooled_data_dict_absolute.keys():
        for ct_feedback in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window].keys():
            print('\tct_feedback: ', ct_feedback)
            for conn_template in angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys():
                print('\t\tconn_template: ', conn_template)
                
                sim_combination = [angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop]['270'],
                                   angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback][conn_template][plot_pop]['90'],
                                   'two-sided']
                U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                significance = significance_symbol(p_mw)
                print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")
                
                
                # reference_data_index = 4
                # for data_ind, data_combination in enumerate(test_data):
                #     if data_ind == reference_data_index: continue

                #     # print('\t\t\t\n',test_data[reference_data_index],'\t',data_combination)

                #     sim_combination = [test_data[reference_data_index],data_combination,'two-sided']

                #     U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                    
                #     # if p_mw<0.05:   significance = '*'
                #     # else:           significance = 'ns'
                #     significance = significance_symbol(p_mw)
                #     print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")



            # for conn_template_ind, conn_template in enumerate(angulart_tuning_pooled_data_dict_absolute[angular_tuning_window][ct_feedback].keys()):
            # # for conn_template_ind, conn_template in enumerate(conn_templates_set):
            #     print('\t\tconn template',conn_template)

            #     test_data=[list(store_RMSE_dict[plot_pop][angular_tuning_window][sim_name].values())      for sim_name in store_RMSE_dict[plot_pop][angular_tuning_window].keys() if 'awake_'+str(conn_template_ind) in sim_name]
                
            #     for data_ind, data_combination in enumerate(test_data):
            #         if data_ind == 0: continue

            #         print('\t\t\t\n',test_data[0],'\t',data_combination)

            #         sim_combination = [test_data[0],data_combination,'two-sided']

            #         U, p_mw = stats.mannwhitneyu(sim_combination[0], sim_combination[1], alternative=sim_combination[2])
                    
            #         # if p_mw<0.05:   significance = '*'
            #         # else:           significance = 'ns'
            #         significance = significance_symbol(p_mw)
            #         print(f"\t\t\t|\t{significance} \t| Mann-Whitney U test statistic: {U}, p-value: {p_mw}\n")

            #         '''
            #         - This test shows how CT modulation changes the RMSE of the angular tuning
            #         - The main result is: 
            #             no significant change in RMSE for the same connectivity with increasing CT modulation
            #             shows that CT modulation affects the error uniformly, even if it increases tuning or not (which will be shown later)
            #         ___________________________________________________________________________________

            #         Figure mapping to statistical tests:
                    
            #             x   x   x   x   x       x = uniform
            #             y   y   y   y   y       y = mixed
            #             z   z   z   z   z       z = closed
            #             0   10  25  50  100
            #                     CT mod
            #         ___________________________________________________________________________________
                    
            #         Test:
            #             x_0 vs [x_10, x_25, x_50, x_100] - [for all x, y and z]
            #             - varying CT modulation for the same angle/same connectivity

            #         '''
























































import sys; sys.exit()


####################################################################################################################################################################################
from scipy.spatial import ConvexHull

# Function to calculate metrics
def calculate_metrics(values, angles):
    # Convert polar to Cartesian coordinates
    x = values * np.cos(angles)
    y = values * np.sin(angles)
    points = np.column_stack((x, y))
    
    # Convex Hull for accurate area and perimeter
    hull = ConvexHull(points)
    area = hull.volume  # Area of the polygon
    perimeter = hull.area  # Perimeter of the polygon
    
    # Compactness metric
    compactness = (4 * np.pi * area) / (perimeter ** 2) if perimeter > 0 else 0
    
    return {
        "area": area,
        "perimeter": perimeter,
        "compactness": compactness,
        "points": points,
    }

# Calculate differences
def calculate_differences(values1, values2):
    radial_diff = np.abs(values1 - values2)
    angular_diff = np.arccos(
        np.cos(angles[:-1]) * np.cos(angles[:-1]) + np.sin(angles[:-1]) * np.sin(angles[:-1])
    )
    return {
        "radial_diff_mean": np.mean(radial_diff),
        "radial_diff_std": np.std(radial_diff),
        "angular_diff_mean": np.mean(angular_diff),
        "angular_diff_std": np.std(angular_diff),
    }



shape_stats = {}
mean_angulart_tuning_pooled_data_dict['20ms']['ct_uniform']['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_10_90_v03']['VPM__pop']
for angular_tuning_window in angular_tuning_windows:
    shape_stats.update({angular_tuning_window:{}})
    print('\n - Angular tuning window: ', angular_tuning_window)
    for ct_feedback in ct_feedbacks:
        shape_stats[angular_tuning_window].update({ct_feedback:{}})
        print('\t - CT feedback: ', ct_feedback)
        for conn_template in conn_templates:
            shape_stats[angular_tuning_window][ct_feedback].update({conn_template:{}})
            print('\t\t - Conn template: ', conn_template)
            for test_pop in ['VPM__pop','TRN__pop']:
                
                angles = [int(angle) for angle in mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][test_pop].keys()]
                mean_values = [mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][test_pop][angle]['mean'] for angle in mean_angulart_tuning_pooled_data_dict[angular_tuning_window][ct_feedback][conn_template][test_pop].keys()]

                angles.append(angles[0])
                mean_values.append(mean_values[0])
                
                angles_=np.array(angles)
                mean_values_=np.array(mean_values)

                shape_metrics=calculate_metrics(mean_values_,angles_)
                shape_stats[angular_tuning_window][ct_feedback][conn_template].update({test_pop:shape_metrics})



shape_stats_hartings = {}
for pop in angular_tuning_hartings.keys():
        
    angles_hartings = [int(angle) for angle in angular_tuning_hartings[pop].keys()]
    mean_values_hartings = [angular_tuning_hartings[pop][angle] for angle in angular_tuning_hartings[pop].keys()]

    angles_hartings.append(angles_hartings[0])
    mean_values_hartings.append(mean_values_hartings[0])


    angles_hartings_=np.array(angles_hartings)
    mean_values_hartings_=np.array(mean_values_hartings)

    shape_metrics_hartings=calculate_metrics(mean_values_hartings_,angles_hartings_)
    shape_stats_hartings.update({pop:shape_metrics_hartings})


# compare to Hartings data
compare_to_hartings={}
for pop in shape_stats_hartings.keys():
    compare_to_hartings.update({pop:{}})
    for angular_tuning_window in angular_tuning_windows:
        compare_to_hartings[pop].update({angular_tuning_window:{}})
        print('\n - Angular tuning window: ', angular_tuning_window)
        for ct_feedback in ct_feedbacks:
            compare_to_hartings[pop][angular_tuning_window].update({ct_feedback:{}})
            print('\t - CT feedback: ', ct_feedback)
            for conn_template in conn_templates:
                print('\t\t - Conn template: ', conn_template)

                pop_differences = calculate_differences(shape_stats_hartings[pop], shape_stats[angular_tuning_window][ct_feedback][conn_template][pop])

                compare_to_hartings[pop][angular_tuning_window][ct_feedback].update({conn_template:pop_differences})




                '''Continue updating chat gpt code here'''

                # Compute shape differences
                differences = calculate_differences(values_schema1[:-1], values_schema2[:-1])

                # Print results
                print("Schema 1 Metrics:", metrics_schema1)
                print("Schema 2 Metrics:", metrics_schema2)
                print("Differences:", differences)

                # Plotting the polygons
                fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})
                ax.plot(angles, values_schema1, label='Schema 1', color='blue', linewidth=2)
                ax.fill(angles, values_schema1, color='blue', alpha=0.3)
                ax.plot(angles, values_schema2, label='Schema 2', color='red', linewidth=2)
                ax.fill(angles, values_schema2, color='red', alpha=0.3)

                # Customizations
                ax.set_theta_offset(np.pi / 2)  # Start at top
                ax.set_theta_direction(-1)     # Clockwise
                ax.set_ylim(0, 1)              # Normalize radius
                ax.set_xticks(angles[:-1])     # Show angle markers
                ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1))

                plt.title("Polygon Comparison in Polar Coordinates")
plt.show()

###################################################################################################################################################################################

'''
Output:

# Used for stats in SfN 2024

- Angular tuning window:  20ms
     - CT feedback:  ct_uniform
         - Conn template:  uniform
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.7521647670540255      ns
             0      0.8160781034438156      ns
             45      0.6461295673045595      ns
             90      0.33567025476613344      ns
             135      0.5211147244744442      ns
             180      0.7440833309073283      ns
             225      0.8407216342733626      ns
             - Population tested: MLe__pop vs  TRN__pop
<ipython-input-37-ce7293c8b83d>:17: RuntimeWarning: Precision loss occurred in moment calculation due to catastrophic cancellation. This occurs when the data are nearly identical. Results may be unreliable.
  statistics = stats.ttest_ind(data1,data2)
             270      0.003519942901928932      **
             315      0.04713623710269477      *
             0      4.113968290590229e-07      ***
             45      3.093227441229925e-09      ***
             90      2.1459879621579646e-08      ***
             135      4.890368944248449e-07      ***
             180      3.2475849053072085e-06      ***
             225      0.004517099970446306      **
         - Conn template:  open
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.8409149922165972      ns
             0      0.6707005592939309      ns
             45      0.9039393865753051      ns
             90      0.6269059289848079      ns
             135      0.9543490442830542      ns
             180      0.8883366897667537      ns
             225      0.7565389880203184      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.0807360045427442      ns
             315      0.32287926960417046      ns
             0      0.6909197324824439      ns
             45      0.7045685693227325      ns
             90      0.8643933834537688      ns
             135      0.8292541979703243      ns
             180      0.7902790669413556      ns
             225      0.4901406798321766      ns
         - Conn template:  closed
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.9754254502666918      ns
             0      0.5919158962766029      ns
             45      0.5776536252303559      ns
             90      0.06813736413042203      ns
             135      0.33811383653364235      ns
             180      0.7850241196330295      ns
             225      0.6919794747889988      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.3342819433946589      ns
             315      0.08954726643851267      ns
             0      0.17993112602572528      ns
             45      0.3283147519632821      ns
             90      0.10617175918271851      ns
             135      0.20824131459945489      ns
             180      0.1626404304671703      ns
             225      0.058632679331827885      ns
         - Conn template:  mixed
             - Population tested: MLe__pop vs  VPM__pop
             270      0.33428194339465434      ns
             315      0.5255560672755276      ns
             0      0.2512161186950974      ns
             45      0.22641319319117156      ns
             90      0.2461629280751862      ns
             135      0.47364948017220865      ns
             180      0.7104002716633144      ns
             225      0.9958568913761077      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.06990595013781015      ns
             315      0.09079596992950438      ns
             0      0.1480519436350724      ns
             45      0.3398973426048145      ns
             90      0.11939101896951501      ns
             135      0.6460006436736326      ns
             180      0.37667221635933      ns
             225      0.1578127988236212      ns
     - CT feedback:  ct_active_fb
         - Conn template:  uniform
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.9159692017935898      ns
             0      0.70801099108505      ns
             45      0.5568583980828119      ns
             90      0.09099982572225208      ns
             135      0.3055851887624087      ns
             180      0.9586380779124999      ns
             225      0.6078614384518157      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.006310914661733997      *
             315      0.0077096042680869226      *
             0      1.886219141348052e-08      ***
             45      1.732267726982455e-09      ***
             90      4.659752596515513e-10      ***
             135      1.8485649680479874e-09      ***
             180      3.7517296215144927e-07      ***
             225      0.0012713425508146224      **
         - Conn template:  open
             - Population tested: MLe__pop vs  VPM__pop
             270      0.3342819433946581      ns
             315      0.9309540222305736      ns
             0      0.7502920912747777      ns
             45      0.6634030493885488      ns
             90      0.5458559918699506      ns
             135      0.577382609834701      ns
             180      0.8898129562339198      ns
             225      0.7080594982102526      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.05980867328195896      ns
             315      0.3965648239324808      ns
             0      0.5182417721541434      ns
             45      0.8861892746440201      ns
             90      0.5049956323911747      ns
             135      0.8174736695735035      ns
             180      0.74082183863775      ns
             225      0.2878049770073335      ns
         - Conn template:  closed
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.6801840189461942      ns
             0      0.7253903375973777      ns
             45      0.5710410965705911      ns
             90      0.20408177350098236      ns
             135      0.20776085438936825      ns
             180      0.39370491831650933      ns
             225      0.6218620436748303      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.1495075777270052      ns
             315      0.32900006031726525      ns
             0      0.4304480029207195      ns
             45      0.968388099088163      ns
             90      0.301707092856313      ns
             135      0.49738038825167075      ns
             180      0.17092726778897238      ns
             225      0.12043871004659591      ns
         - Conn template:  mixed
             - Population tested: MLe__pop vs  VPM__pop
             270      0.3342819433946561      ns
             315      0.4189098699939853      ns
             0      0.4113505687164942      ns
             45      0.23088227354059157      ns
             90      0.23615742254862038      ns
             135      0.28264953070975435      ns
             180      0.5714889148323132      ns
             225      0.9863454550930799      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.20865804075248506      ns
             315      0.08886752361027642      ns
             0      0.3626172621505077      ns
             45      0.8093006000035265      ns
             90      0.4149821729216967      ns
             135      0.36655681330203116      ns
             180      0.21105085873384777      ns
             225      0.46780638897528715      ns

 - Angular tuning window:  250ms
     - CT feedback:  ct_uniform
         - Conn template:  uniform
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.012622791363871391      *
             0      0.016270904676565688      *
             45      0.05251095504576225      ns
             90      0.9231855944630367      ns
             135      0.05887130045556171      ns
             180      0.019793809072945232      *
             225      0.062233255224648475      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      4.523788342626824e-05      ***
             315      0.15617214436186322      ns
             0      1.9563200286997452e-05      ***
             45      6.066912331781352e-08      ***
             90      1.5920461770643077e-08      ***
             135      1.5000844204300226e-06      ***
             180      1.9081794114168628e-05      ***
             225      0.6733925930060165      ns
         - Conn template:  open
             - Population tested: MLe__pop vs  VPM__pop
             270      0.1645387887431976      ns
             315      0.5425191899392674      ns
             0      0.8478508766181552      ns
             45      0.21356280007761433      ns
             90      0.12909798977059053      ns
             135      0.2510760956481604      ns
             180      0.700662608850043      ns
             225      0.48453975156745044      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.11377626701669141      ns
             315      0.3616824082203851      ns
             0      0.3231349644606416      ns
             45      9.451524403872792e-05      ***
             90      3.1246857731962244e-05      ***
             135      0.0009682641849897999      **
             180      0.35764168026753307      ns
             225      0.4345966311006293      ns
         - Conn template:  closed
             - Population tested: MLe__pop vs  VPM__pop
             270      0.0379515531175728      *
             315      0.17504210003728374      ns
             0      0.020828709586738024      *
             45      0.010347020245842815      *
             90      0.02911435013971513      *
             135      0.021988151482412382      *
             180      0.02845645370396743      *
             225      0.10236542267889964      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.33428194339465733      ns
             315      0.5668663693468178      ns
             0      0.7722587430611713      ns
             45      0.04701670297461397      *
             90      0.3640710866549921      ns
             135      0.1234955558723532      ns
             180      0.34713911761578986      ns
             225      0.666071723775566      ns
         - Conn template:  mixed
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.03879530138172978      *
             0      0.007696019315064937      *
             45      0.12718754477180796      ns
             90      0.491817691436363      ns
             135      0.09513213053005522      ns
             180      0.034829894557695446      *
             225      0.06269358070799357      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.3342819433946572      ns
             315      0.6640379296977164      ns
             0      0.0851874308302789      ns
             45      0.00019610938979277437      ***
             90      0.0005132889740781802      **
             135      9.796326015249533e-05      ***
             180      0.032382117991461117      *
             225      0.8128806606713808      ns
     - CT feedback:  ct_active_fb
         - Conn template:  uniform
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.001302103474073535      **
             0      4.3958582184454936e-05      ***
             45      2.2359483616847987e-06      ***
             90      2.925211151593771e-06      ***
             135      5.593565142623789e-06      ***
             180      1.2369832370814e-05      ***
             225      0.004169884669974172      **
             - Population tested: MLe__pop vs  TRN__pop
             270      2.294644411752482e-05      ***
             315      0.018265527994845092      *
             0      1.1388893319951115e-07      ***
             45      2.942132152330028e-09      ***
             90      1.232488720769387e-10      ***
             135      1.670372217510306e-09      ***
             180      1.7723348982494484e-06      ***
             225      0.07187436414605468      ns
         - Conn template:  open
             - Population tested: MLe__pop vs  VPM__pop
             270      0.07480403213000102      ns
             315      0.42693851648903947      ns
             0      0.28275600553802066      ns
             45      0.4280831513633977      ns
             90      0.49942401332819797      ns
             135      0.4703868772420895      ns
             180      0.36980400311243034      ns
             225      0.569778240450557      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.052128583667100645      ns
             315      0.16934248746116795      ns
             0      0.05066368065438297      ns
             45      0.296520589616069      ns
             90      0.518904000744726      ns
             135      0.40883204238283455      ns
             180      0.13592937382348133      ns
             225      0.2858731886238332      ns
         - Conn template:  closed
             - Population tested: MLe__pop vs  VPM__pop
             270      0.0021728834090289674      **
             315      0.2221817331817116      ns
             0      0.024360988377415946      *
             45      0.015295228334090753      *
             90      0.029565194150196438      *
             135      0.013181845090746117      *
             180      0.025235951222648986      *
             225      0.18543604392609048      ns
             - Population tested: MLe__pop vs  TRN__pop
             270      0.2020255416083232      ns
             315      0.02934464928216502      *
             0      0.0015056246686059856      **
             45      0.000749128766645593      **
             90      0.000694993205377088      **
             135      0.0010125642894514937      **
             180      0.002466319208073631      **
             225      0.019127092351960053      *
         - Conn template:  mixed
             - Population tested: MLe__pop vs  VPM__pop
             270     values are the same
             315      0.0004077845182760703      ***
             0      8.477844505587531e-05      ***
             45      0.0005931388692210781      **
             90      0.0007205936168932276      **
             135      0.00035835778559547584      ***
             180      0.00011112408265378173      ***
             225      0.0005615978187869388      **
             - Population tested: MLe__pop vs  TRN__pop
             270      0.05192705329708703      ns
             315      0.10196914746030047      ns
             0      0.032257032352717485      *
             45      0.18387950481852952      ns
             90      0.20392443828593101      ns
             135      0.13618196879423528      ns
             180      0.036038614723105857      *
             225      0.13803915437627692      ns



'''


# angular_tuning = {  'ind':[  0,      1,                  2,                  3,                   4],
#                     'mean':[  1.0,    0.6982148208637337, 0.5676691729323308, 0.44418423989406364, 0.4282115869017633]}



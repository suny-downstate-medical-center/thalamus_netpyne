import json
import numpy as np

load_files=[
     ('../conn/conn_validation/connectivity_data__VPM_ff_closed__TRN_fb_closed__L6A_fb_closed_v04.json',        'closed'),
     ('../conn/conn_validation/connectivity_data__VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform_v04.json',   'uniform'),
     ('../conn/conn_validation/connectivity_data__VPM_ff_closed__TRN_fb_open__L6A_fb_closed_v04.json',      'open'),
]
convergence_stats_dict={'bbp':      {
                            'VPM': {
                                'MLe': {'mean': 143.4382153702406,      'std': 48.286062396690696},  
                                'TRN': {'mean': 520.8336684539768,      'std': 446.4196876699123},  
                                'L6A': {'mean': 1044.720254720143,      'std': 452.7204081938519}}, 
                            'TRN': {
                                'VPM': {'mean': 188.48339784070075,     'std': 149.87060656327023}, 
                                'TRN': {'mean': 71.76451415766958,      'std': 59.63718275874621},
                                'L6A': {'mean': 458.2288861689106,      'std': 162.177222938063},}, 
                                }}

fontsize=60

for (load_file,load_file_key) in load_files:
    # Open a JSON file and load its contents into a Python dictionary
    with open(load_file, 'r') as file: model_data = json.load(file)
    # --- Reads the saved dictionary and modifies the keys to rename CTvirtual_uniform to L6A, and merge TRN and TRN_ring into a single pop
    output = {}
    for target_pop, targets in model_data['cell_convergence_mean_std'].items():
        if target_pop not in ['VPM','TRN']:continue
        post_pop = 'TRN' if target_pop in ['TRN', 'TRN_ring'] else target_pop
        if post_pop not in output: output[post_pop] = {}

        for source_pop, value in targets.items():
            if source_pop in ['TRN', 'TRN_ring']:     pre_pop = 'TRN' 
            elif source_pop == 'VPM':                 pre_pop = 'VPM' 
            elif source_pop == 'MLe':                 pre_pop = 'MLe' 
            elif source_pop == 'CTvirtual_uniform':   pre_pop = 'L6A'
            else:continue
            if pre_pop not in output[post_pop]: output[post_pop][pre_pop] = {'mean':0,'std':0}
            output[post_pop][pre_pop]['mean']+=value[0]
            output[post_pop][pre_pop]['std']+=value[1]
    # --- Stores the values from different network templates
    convergence_stats_dict.update({load_file_key:output})


import matplotlib.pyplot as plt
networks = list(convergence_stats_dict.keys())
offset = 0.15
target_pops = list(convergence_stats_dict[networks[0]].keys())
colors_dict={
                'bbp':      'cornflowerblue',
                'closed':   'orange',
                'uniform':  'peru',
                'open':     'saddlebrown',
                }

plt.figure(figsize=(15,15))
plt.rcParams.update({'font.size': fontsize,'mathtext.default': 'regular'})
# plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
for target_pop_ind,target_pop in enumerate(['TRN','VPM']):
# for target_pop_ind,target_pop in enumerate(target_pops):

    source_pops=list(convergence_stats_dict[networks[0]][target_pop].keys())
    # replaces (VPM by TC) and (L6A by CT)
    source_pops_name=[]
    for s_pop_name in source_pops:
        if      'VPM' in s_pop_name:    source_pops_name.append('TC')
        elif    'L6A' in s_pop_name:    source_pops_name.append('CT')
        else:                           source_pops_name.append(s_pop_name)


    print(target_pop)
    print(source_pops,'\n')

    plt.subplot(len(target_pops),1,target_pop_ind+1)
    for source_pop_ind, source_pop in enumerate(source_pops):
        for network_ind, network in enumerate(networks):
            plt.bar(source_pop_ind+(network_ind*offset),convergence_stats_dict[network][target_pop][source_pop]['mean'], width=0.14, color=colors_dict[network])
            plt.errorbar(source_pop_ind+(network_ind*offset),
                         convergence_stats_dict[network][target_pop][source_pop]['mean'],
                         yerr=convergence_stats_dict[network][target_pop][source_pop]['std'],
                        fmt='', linestyle='',
                        capsize=7,elinewidth=4,
                        color='k')
    # plt.title(target_pop)
    plt.box(False)
    plt.xticks([i+2*offset for i in range(3)], source_pops_name)
    plt.yticks([0, 750, 1500])
    plt.ylim([0,1700])
    plt.ylabel('Convergence')
    if target_pop_ind+1==len(target_pops): plt.xlabel('Source pop')
    plt.tight_layout()
    plt.savefig('../sim_figs/fig01_combined_bar_plot.png')
    # plt.savefig('figs_conn_validation/combined_bar_plot_stats2.png')

pathway_ratios_dict={}
for network_template in convergence_stats_dict.keys():
    ratios_dict = {
        'L6A/MLe ->VPM':convergence_stats_dict[network_template]['VPM']['L6A']['mean']/convergence_stats_dict[network_template]['VPM']['MLe']['mean'],
        'TRN/MLe ->VPM':convergence_stats_dict[network_template]['VPM']['TRN']['mean']/convergence_stats_dict[network_template]['VPM']['MLe']['mean'],
        'L6A-> VPM/TRN':convergence_stats_dict[network_template]['VPM']['L6A']['mean']/convergence_stats_dict[network_template]['TRN']['L6A']['mean'],    
        'L6A/VPM ->TRN':convergence_stats_dict[network_template]['TRN']['L6A']['mean']/convergence_stats_dict[network_template]['TRN']['VPM']['mean'],
        'TRN/VPM ->TRN':convergence_stats_dict[network_template]['TRN']['TRN']['mean']/convergence_stats_dict[network_template]['TRN']['VPM']['mean'],
    }
    pathway_ratios_dict.update({network_template:ratios_dict})

pathway_ratios= list(pathway_ratios_dict[list(pathway_ratios_dict.keys())[0]].keys())
pathway_ratios_labels = ['i','ii','iii','iv','v']

plt.figure(figsize=(15,15))
plt.rcParams.update({'font.size': fontsize,'mathtext.default': 'regular'})
# plt.rcParams.update({'font.size': 30,'mathtext.default': 'regular'})
for pathway_ratio_ind,pathway_ratio in enumerate(pathway_ratios):
    for network_ind, network in enumerate(pathway_ratios_dict.keys()):
            # plotting separately to create legend with a single instance of each network name on label
            if pathway_ratio_ind == 0:  plt.bar(pathway_ratio_ind+(network_ind*offset), pathway_ratios_dict[network][pathway_ratio], width=0.14, color=colors_dict[network],label=network)
            else:                       plt.bar(pathway_ratio_ind+(network_ind*offset), pathway_ratios_dict[network][pathway_ratio], width=0.14, color=colors_dict[network])
plt.box(False)
# plt.xticks([i+2*offset for i in range(len(pathway_ratios))], pathway_ratios, rotation=10)
plt.xticks([i+2*offset for i in range(len(pathway_ratios))], pathway_ratios_labels)
plt.yticks([0, 2.5, 5, 7.5])
plt.ylabel('Convergence ratio')
legend = plt.legend()
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
if network_ind+1==len(pathway_ratios_dict.keys()): plt.xlabel('Pathways')
plt.tight_layout()
plt.savefig('../sim_figs/fig01_combined_pathway_ratios.png')
# plt.savefig('figs_conn_validation/combined_bar_plot_pathway_ratios2.png')


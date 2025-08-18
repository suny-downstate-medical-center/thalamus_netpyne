'''
script replacing PlotAngularTuning.py
'''

import sys
import os

try:
    # Add project root to sys.path
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if project_root not in sys.path: sys.path.append(project_root)
    from sim.Build_Net import BuildNetwork
except: 
    print('Failed to import network template from /sim/Build_Net.py. Setting network properties manually')
    diam_default = {
        'VPM': [450.0, 550.0],
        'TRN': [450.0, 550.0],
        'MLe': [472.5, 527.5],
        'TRN_ring': [430.0, 570.0],
        'L6A': [325.0, 675.0],
        'L6A_activated': [325.0, 675.0],
        'L6A_suppressed': [325.0, 675.0],
        'L6A_silent': [325.0, 675.0],
        'L6A_sparse': [325.0, 675.0],
        'CTvirtual': [325.0, 675.0],
        'CTvirtual_activated': [325.0, 675.0],
        'CTvirtual_suppressed': [325.0, 675.0],
        'CTvirtual_sparse': [325.0, 675.0],
        'CTvirtual_silent': [325.0, 675.0],
        'CTvirtual_uniform': [325.0, 675.0],
        'BaselineDriver': [400.0, 600.0]}
    height_default = {
        'VPM': [3300, 4000],
        'TRN': [2550, 2800],
        'MLe': [4500, 5700],
        'TRN_ring': [2550, 2800],
        'L6A': [890, 1090],
        'L6A_activated': [890, 1090],
        'L6A_suppressed': [890, 1090],
        'L6A_silent': [890, 1090],
        'L6A_sparse': [890, 1090],
        'CTvirtual': [890, 1090],
        'CTvirtual_activated': [890, 1090],
        'CTvirtual_suppressed': [890, 1090],
        'CTvirtual_sparse': [890, 1090],
        'CTvirtual_silent': [890, 1090],
        'CTvirtual_uniform': [890, 1090],
        'BaselineDriver': [2500, 2550]}

    pass

import json
import numpy as np
from matplotlib import pyplot as plt
# from ..sim.Build_Net import BuildNetwork
from src.ProcessAngles import AngleProcessor
angle_processor = AngleProcessor()

class PlotAngularTuning():
    def plotAngularTuning(angular_tuning_dict, 
                          barreloid_network_cell_properties='../network_template/barreloid_network_cell_properties.json', 
                          center_point=500,
                          select_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop', 'L6A_activated__pop'],
                          save_fig_name = None,
                          select_deflection=['all'],
                          deflection_angle = 0,
                          ):

        # --- Load saved network template - to make plotting compatible in parallel mode (MPI)
        with open(barreloid_network_cell_properties, 'r') as file: allCells = json.load(file)

        # angular_tuning_dict = {'spkTimes':spkTimes_byDeflection, 'spkInds':spkInds_byDeflection, 'spkGids':spkGids_byDeflection}

        spkTimes_byDeflection = angular_tuning_dict['spkTimes']
        spkGids_byDeflection  = angular_tuning_dict['spkGids']

        angular_tuning_plot_data_dict = {}

        # --- Plot angular tuning
        for deflection_flag in spkTimes_byDeflection.keys():
            if deflection_flag not in select_deflection: 
                # print('Skipping deflection ', deflection_flag)
                continue

            # dict to store plotted data
            angular_tuning_plot_data_dict.update({deflection_flag:{}})

            count_spikes_grouped        = PlotAngularTuning.GroupCountSpikes(allCells, spkGids=spkGids_byDeflection[deflection_flag], select_pops= select_pops, center_point = center_point)
            # normalized_spiking_grouped  = PlotAngularTuning.NormalizeSpiking(count_spikes_grouped, select_pops)
            normalized_spiking_grouped  = PlotAngularTuning.NormalizeSpiking_v2(count_spikes_grouped, select_pops, deflection_angle) # Normalizing based on the deflection angle, instead of the max value for a given pop
            PlotAngularTuning.plotFigure(deflection_flag, normalized_spiking_grouped, count_spikes_grouped, save_fig_name)

            for cell_pop in normalized_spiking_grouped.keys():
                angular_tuning_plot_data_dict[deflection_flag].update({cell_pop:{
                                                                                    'normalized':{
                                                                                        'keys':     list(normalized_spiking_grouped[cell_pop].keys()),
                                                                                        'values':   list(normalized_spiking_grouped[cell_pop].values()),
                                                                                    },
                                                                                    'absolute':{
                                                                                        'keys':     list(count_spikes_grouped[cell_pop].keys()),
                                                                                        'values':   list(count_spikes_grouped[cell_pop].values()),
                                                                                    },
                                                                                }})
        return angular_tuning_plot_data_dict

    def GroupCountSpikes(allCells, spkGids, select_pops, center_point):

        # creates a dictionary with spike times for each cell
        spkGids_int = [int(spkGid_) for spkGid_ in spkGids]
        from collections import Counter
        count_spks = Counter(spkGids_int)

        # dict with all gids for the selected pops and their spike counts
        store_count_spks_ON={}
        store_spikes_grouped_ON={}
        count_spikes_grouped={}
        for cell_ind in allCells.keys():
            cell_gid = allCells[str(cell_ind)]['gid']
            cell_pop = allCells[str(cell_ind)]['tags']['pop']
            if cell_pop not in select_pops:continue
            if cell_pop not in store_count_spks_ON.keys():      store_count_spks_ON.update({cell_pop:{}})
            if int(cell_gid) in count_spks.keys():              store_count_spks_ON[cell_pop].update({int(cell_gid):count_spks[int(cell_gid)]})
            else:                                               store_count_spks_ON[cell_pop].update({int(cell_gid):0})
            #
            if cell_pop not in store_spikes_grouped_ON.keys(): 
                angles_dict={angle:[] for angle in angle_processor.reference_angles_}
                store_spikes_grouped_ON.update({cell_pop:angles_dict}) 
                angles_dict_counter={angle:0 for angle in angle_processor.reference_angles_}
                count_spikes_grouped.update({cell_pop:angles_dict_counter}) 
            #
            if ('MLe' in cell_pop) or ('VPM' in cell_pop) or ('TRN' in cell_pop):
                # diam,height         = BuildNetwork.getNetworkSize(sim.cfg.center_point)
                try:
                    diam,height         = BuildNetwork.getNetworkSize(center_point)
                except:
                    diam    = diam_default
                    height  = height_default

                # --- debug boundaries --- seem inverted
                cell_upper_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_lower_boundary = height[cell_pop.split('__pop')[0]][0]
                # # --- debug boundaries --- seem inverted
                # cell_upper_boundary = height[cell_pop.split('__pop')[0]][0]
                # cell_lower_boundary = height[cell_pop.split('__pop')[0]][1]
                cell_position   = allCells[str(cell_ind)]['tags']['y']
            else: # needs to be tested, because L6 is deactivated for now
                cell_upper_boundary=360
                cell_lower_boundary=0
                cell_position_x   = allCells[str(cell_ind)]['tags']['x']
                cell_position_z   = allCells[str(cell_ind)]['tags']['z']
                cell_position=(np.remainder((((np.arctan2(cell_position_x-500,cell_position_z-500))*(180/np.pi))+360),360)/360)
            cell_position_relative = (cell_position-cell_lower_boundary)/(cell_upper_boundary-cell_lower_boundary)
            matching_angle = angle_processor.get_closest_angle(cell_position_relative)

            count_spikes_grouped[cell_pop][matching_angle]+=store_count_spks_ON[cell_pop][cell_gid]

        return count_spikes_grouped

    def NormalizeSpiking(count_spikes_grouped, select_pops):
        print(' - Normalizing spiking activity - ')
        # --- calculates the normalized values
        normalized_spiking_grouped  = {cell_pop:{} for cell_pop in count_spikes_grouped.keys()}
        pop_spiking                 = {cell_pop:0  for cell_pop in count_spikes_grouped.keys()}
        for cell_pop in count_spikes_grouped.keys():
            if cell_pop in select_pops:
                # Counts the number of spikes for a given pop
                pop_spiking[cell_pop]=sum(count_spikes_grouped[cell_pop].values())
                for matching_angle in count_spikes_grouped[cell_pop].keys():
                    # print(matching_angle)
                    if max(count_spikes_grouped[cell_pop].values())<=0:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:0})
                    else:normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped[cell_pop][matching_angle]/max(count_spikes_grouped[cell_pop].values())})
        return normalized_spiking_grouped

    def NormalizeSpiking_v2(count_spikes_grouped, select_pops, deflection_angle):
        print(' - Normalizing spiking activity - ')
        # --- calculates the normalized values
        normalized_spiking_grouped  = {cell_pop:{} for cell_pop in count_spikes_grouped.keys()}
        pop_spiking                 = {cell_pop:0  for cell_pop in count_spikes_grouped.keys()}
        for cell_pop in count_spikes_grouped.keys():
            if cell_pop in select_pops:
                # Counts the number of spikes for a given pop
                pop_spiking[cell_pop]=sum(count_spikes_grouped[cell_pop].values())
                for matching_angle in count_spikes_grouped[cell_pop].keys():
                    # print(cell_pop,matching_angle)
                    # print('value - matching angle: ', count_spikes_grouped[cell_pop][matching_angle])
                    # print('value - deflect  angle: ', count_spikes_grouped[cell_pop][deflection_angle])
                    # print('value - max val  angle: ', max(count_spikes_grouped[cell_pop].values()))
                    if max(count_spikes_grouped[cell_pop].values())<=0:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:0})
                    else:
                        normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped[cell_pop][matching_angle]/count_spikes_grouped[cell_pop][deflection_angle]})
                        # normalized_spiking_grouped[cell_pop].update({matching_angle:count_spikes_grouped[cell_pop][matching_angle]/max(count_spikes_grouped[cell_pop].values())})
        return normalized_spiking_grouped
            

    def plotFigure(deflection_flag, normalized_spiking_grouped, count_spikes_grouped, save_fig_name):

        print(' - Plotting Anuglar Tuning Figure - ', deflection_flag)

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
        
        print('Plotting normalized and non-normalized angular tuning')
        print(deflection_flag)
        # --- Keep debugging angular tuning figure
        plt.figure(figsize=(20,40))
        plt.rcParams.update({'font.size': 80})
        # plt.rcParams.update({'font.size': 40})
        plt.subplot(2,1,1,projection='polar')
        for cell_pop in normalized_spiking_grouped.keys():
            polar_list_values=list(normalized_spiking_grouped[cell_pop].values())
            polar_list_values.extend([polar_list_values[0]])
            polar_list_keys=list(normalized_spiking_grouped[cell_pop].keys())
            polar_list_keys.extend([360])
            # plt.plot(list(normalized_spiking_grouped[cell_pop].keys()),list(normalized_spiking_grouped[cell_pop].values()))
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
        
        plt.yticks([0.5, 1.0])
        plt.ylim([0,2.0])
        plt.xticks(np.deg2rad([0,90,180,270]))
        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine        

        plt.subplot(2,1,2,projection='polar')
        for cell_pop in count_spikes_grouped.keys():
            polar_list_values=list(count_spikes_grouped[cell_pop].values())
            polar_list_values.extend([polar_list_values[0]])
            polar_list_keys=list(count_spikes_grouped[cell_pop].keys())
            polar_list_keys.extend([360])
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values)
        plt.legend(list(normalized_spiking_grouped.keys()),loc='upper right', bbox_to_anchor=(0,0.2))

        plt.suptitle('Angular tuning plot - '+deflection_flag)

        # plt.yticks([0.5, 1.0])
        # plt.ylim([0,1.5])
        plt.xticks(np.deg2rad([0,90,180,270]))
        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine        

        if save_fig_name is not None:
            if '.png' in save_fig_name: save_fig_name=save_fig_name.split('.png')[0]
        
        save_fig_name_=save_fig_name+'__'+deflection_flag+'.png'

        try:    plt.savefig(str(save_fig_name_))
        except: plt.savefig(str('angularTuning3__'+deflection_flag+'.png'))


    def add_rotated_angularTuning_data(angular_tuning_plot_data_dict,deflection_angle,target_rotation_angle=270):
        
        for deflection_flag in angular_tuning_plot_data_dict.keys():
            for cell_pop in angular_tuning_plot_data_dict[deflection_flag].keys():
                keys            = angular_tuning_plot_data_dict[deflection_flag][cell_pop]['normalized']['keys']
                rotation_diff   = int(target_rotation_angle)-int(deflection_angle)
                keys_rotated    = [int(k)+rotation_diff for k in keys]
                
                keys_rotated_=[]
                for kr in keys_rotated:
                    if   kr>= 360:  keys_rotated_.append(kr-360)
                    elif kr<0:      keys_rotated_.append(kr+360)
                    else:           keys_rotated_.append(kr)


                # keys_rotated_   = [kr-360 for kr in keys_rotated if kr>=360]

                angular_tuning_plot_data_dict[deflection_flag][cell_pop]['normalized'].update({'keys_rotated':keys_rotated_})
                angular_tuning_plot_data_dict[deflection_flag][cell_pop]['absolute'].update(  {'keys_rotated':keys_rotated_})

        return angular_tuning_plot_data_dict
    
    def parseTuningValues(store_angular_tuning_plot_data_dict, data_type = 'absolute'):
        # separating the operations carried by combine_rotated_angularTuning_data and combine_rotated_angularTuning_data_absolute
        # makes a dictionary aligned to 270 degrees
        parse_data_dict = {}
        for deflection_angle in store_angular_tuning_plot_data_dict.keys():
            for cell_pop in store_angular_tuning_plot_data_dict[deflection_angle]['all'].keys():
                values          = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop][data_type]['values']
                keys_rotated    = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop][data_type]['keys_rotated']
                
                if cell_pop not in parse_data_dict.keys():parse_data_dict.update({cell_pop:{}})
                parse_data_dict[cell_pop].update({deflection_angle:{'values':values, 'keys_rotated':keys_rotated}})

        deflection_ordered_dict={}
        for cell_pop in parse_data_dict.keys():
            deflection_ordered_dict.update({cell_pop:PlotAngularTuning.sort_keys_and_values(parse_data_dict[cell_pop])})
        
        return deflection_ordered_dict
    
    def calculate_mean_angular_tuning_dict(deflection_ordered_dict):
        mean_angular_tuning_dict={}
        for cell_pop in deflection_ordered_dict.keys():
            tuning_values = [deflection_dict['values']  for deflection_dict in list(deflection_ordered_dict[cell_pop].values())]
            keys_rotated  = list(deflection_ordered_dict[cell_pop].values())[0]['keys_rotated']
            mean_angular_tuning_dict.update({cell_pop:{'keys_rotated':keys_rotated,'mean':np.mean(tuning_values,0),'std':np.std(tuning_values,0)}})
    
        return mean_angular_tuning_dict

    def combine_rotated_angularTuning_data(store_angular_tuning_plot_data_dict):

        parse_data_dict = {}
        for deflection_angle in store_angular_tuning_plot_data_dict.keys():
            for cell_pop in store_angular_tuning_plot_data_dict[deflection_angle]['all'].keys():
                values          = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop]['normalized']['values']
                keys_rotated    = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop]['normalized']['keys_rotated']
                
                if cell_pop not in parse_data_dict.keys():parse_data_dict.update({cell_pop:{}})
                parse_data_dict[cell_pop].update({deflection_angle:{'values':values, 'keys_rotated':keys_rotated}})

        deflection_ordered_dict={}
        for cell_pop in parse_data_dict.keys():
            deflection_ordered_dict.update({cell_pop:PlotAngularTuning.sort_keys_and_values(parse_data_dict[cell_pop])})
            
        mean_angular_tuning_dict={}
        for cell_pop in deflection_ordered_dict.keys():
            tuning_values = [deflection_dict['values']  for deflection_dict in list(deflection_ordered_dict[cell_pop].values())]
            keys_rotated  = list(deflection_ordered_dict[cell_pop].values())[0]['keys_rotated']
            mean_angular_tuning_dict.update({cell_pop:{'keys_rotated':keys_rotated,'mean':np.mean(tuning_values,0),'std':np.std(tuning_values,0)}})
    
        return mean_angular_tuning_dict
    
    def combine_rotated_angularTuning_data_absolute(store_angular_tuning_plot_data_dict):

        parse_data_dict = {}
        for deflection_angle in store_angular_tuning_plot_data_dict.keys():
            for cell_pop in store_angular_tuning_plot_data_dict[deflection_angle]['all'].keys():
                values          = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop]['absolute']['values']
                keys_rotated    = store_angular_tuning_plot_data_dict[deflection_angle]['all'][cell_pop]['absolute']['keys_rotated']
                
                if cell_pop not in parse_data_dict.keys():parse_data_dict.update({cell_pop:{}})
                parse_data_dict[cell_pop].update({deflection_angle:{'values':values, 'keys_rotated':keys_rotated}})

        deflection_ordered_dict={}
        for cell_pop in parse_data_dict.keys():
            deflection_ordered_dict.update({cell_pop:PlotAngularTuning.sort_keys_and_values(parse_data_dict[cell_pop])})
            
        mean_angular_tuning_dict_absolute={}
        for cell_pop in deflection_ordered_dict.keys():
            tuning_values = [deflection_dict['values']  for deflection_dict in list(deflection_ordered_dict[cell_pop].values())]
            keys_rotated  = list(deflection_ordered_dict[cell_pop].values())[0]['keys_rotated']
            mean_angular_tuning_dict_absolute.update({cell_pop:{
                                                                'keys_rotated':keys_rotated,
                                                                'mean':(np.mean(tuning_values,0)/max(np.mean(tuning_values,0))),
                                                                'std':(np.std(tuning_values,0)/max(np.mean(tuning_values,0))),

                                                                'mean_normByMinVal':(np.mean(tuning_values,0)/min(np.mean(tuning_values,0))),
                                                                'std_normByMinVal':(np.std(tuning_values,0)/min(np.mean(tuning_values,0))),

                                                                'mean_notNomalized':(np.mean(tuning_values,0)),
                                                                'std_notNomalized':(np.std(tuning_values,0)),


                                                                }
                                                            }
                                                        )
    
        return mean_angular_tuning_dict_absolute
    
    def plotCombinedAngularTuning(mean_angular_tuning_dict, angular_tuning_hartings, save_fig_name=None):
        
        print(' - Plotting Anuglar Tuning Figure - combined')
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
        
        print('Plotting normalized and non-normalized angular tuning')
        # --- Keep debugging angular tuning figure
        plt.figure(figsize=(20,25))
        plt.rcParams.update({'font.size': 80})
        plt.subplot(1,1,1,projection='polar')
        for cell_pop in mean_angular_tuning_dict.keys():
            # mean
            polar_list_values_mean=list(mean_angular_tuning_dict[cell_pop]['mean'])
            polar_list_values_mean.extend([polar_list_values_mean[0]])
            # std
            polar_list_values_std=list(mean_angular_tuning_dict[cell_pop]['std'])
            polar_list_values_std.extend([polar_list_values_std[0]])
            # keys
            polar_list_keys=list(mean_angular_tuning_dict[cell_pop]['keys_rotated'])
            polar_list_keys.extend([360])

            # # Plot control values from Hartings 2000
            # control_values = list(angular_tuning_hartings[cell_pop].values())
            # control_values.extend([control_values[0]])
            # try:    plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[cell_pop], alpha=0.25,linewidth=10)
            # except: plt.plot(np.deg2rad(polar_list_keys),control_values)

            # Plot angular tuning values from model
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values_mean,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values_mean)

            plt.errorbar(   
                            np.deg2rad(polar_list_keys), polar_list_values_mean, yerr=polar_list_values_std, 
                            fmt='', linestyle='', 
                            capsize=10,elinewidth=5,markeredgewidth=5,
                            color = popColors_dict[cell_pop],
                            )

        # plt.legend(list(mean_angular_tuning_dict.keys()),loc='upper right', bbox_to_anchor=(0,0.2))

        plt.suptitle('Angular tuning plot - combined')
        plt.yticks([0.5, 1.0])
        plt.ylim([0,1.5])
        plt.xticks(np.deg2rad([0,90,180,270]))
        # plt.grid(False)

        ax1 = plt.gca()
        ax1.spines['polar'].set_visible(False)  # Hide the circular spine        
        
        if save_fig_name is not None:
            if '.png' in save_fig_name: save_fig_name=save_fig_name.split('.png')[0]
        
        save_fig_name_=save_fig_name+'__'+'combined'+'.png'

        try:    plt.savefig(str(save_fig_name_),dpi=200)
        except: plt.savefig(str('angularTuning3__'+'combined'+'.png'))

        print('Plotting normalized and non-normalized angular tuning')
        # --- Keep debugging angular tuning figure
        for cell_pop in mean_angular_tuning_dict.keys():
            plt.figure(figsize=(20,25))
            plt.rcParams.update({'font.size': 80})
            plt.subplot(1,1,1,projection='polar')
            # mean
            polar_list_values_mean=list(mean_angular_tuning_dict[cell_pop]['mean'])
            polar_list_values_mean.extend([polar_list_values_mean[0]])
            # std
            polar_list_values_std=list(mean_angular_tuning_dict[cell_pop]['std'])
            polar_list_values_std.extend([polar_list_values_std[0]])
            # keys
            polar_list_keys=list(mean_angular_tuning_dict[cell_pop]['keys_rotated'])
            polar_list_keys.extend([360])

            # Plot control values from Hartings 2000
            control_values = list(angular_tuning_hartings[cell_pop].values())
            control_values.extend([control_values[0]])
            try:    plt.plot(np.deg2rad(polar_list_keys),control_values,c=popColors_dict[cell_pop], alpha=0.25,linewidth=10)
            except: plt.plot(np.deg2rad(polar_list_keys),control_values)

            # Plot angular tuning values from model
            try:    plt.plot(np.deg2rad(polar_list_keys),polar_list_values_mean,c=popColors_dict[cell_pop],linewidth=5)
            except: plt.plot(np.deg2rad(polar_list_keys),polar_list_values_mean)

            plt.errorbar(   
                            np.deg2rad(polar_list_keys), polar_list_values_mean, yerr=polar_list_values_std, 
                            fmt='', linestyle='', 
                            capsize=10,elinewidth=5,markeredgewidth=5,
                            color = popColors_dict[cell_pop],
                            )

            # plt.legend(list(mean_angular_tuning_dict.keys()),loc='upper right', bbox_to_anchor=(0,0.2))

            plt.suptitle('Angular tuning plot - combined')
            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            # plt.grid(False)

            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine        
            
            if save_fig_name is not None:
                if '.png' in save_fig_name: save_fig_name=save_fig_name.split('.png')[0]
            
            save_fig_name_=save_fig_name+'__'+'combined_'+cell_pop+'.png'

            try:    plt.savefig(str(save_fig_name_),dpi=200)
            except: plt.savefig(str('angularTuning3__'+cell_pop+'combined'+'.png'))

    def sort_keys_and_values(data_dict):
        output_dict = {}
        for key, sub_dict in data_dict.items():
            keys_rotated = sub_dict['keys_rotated']
            values = sub_dict['values']
            # Zip keys and values together, sort by the keys, and unzip them
            sorted_pairs = sorted(zip(keys_rotated, values), key=lambda x: x[0])
            sorted_keys, sorted_values = zip(*sorted_pairs)
            # Update the output dictionary with the sorted values
            output_dict[key] = {
                'keys_rotated': list(sorted_keys),
                'values': list(sorted_values)
            }
        return output_dict

if __name__ == "__main__":

    target_angles=[0, 45, 90, 135, 180, 225, 270, 315]
    extra_zeros=''

    angular_tuning_hartings={
                            'MLe__pop':{
                                '0':    0.33882794754053996,
                                '45':   0.22636951177424108,
                                '90':   0.19764487038033318,
                                '135':  0.17639018632728395,
                                '180':  0.3105937087413645,
                                '225':  0.610603008575053,
                                '270':  1.0,
                                '315':  0.6704982467108909,
                                },
                            'VPM__pop':{
                                '0':    294.51401847013045/589.0280369402608,
                                '45':   241.4544744541332/589.0280369402608,
                                '90':   231.9652440864924/589.0280369402608,
                                '135':  270.7449775143726/589.0280369402608,
                                '180':  332.24549916261424/589.0280369402608,
                                '225':  394.3636706579752/589.0280369402608,
                                '270':  1.0,
                                '315':  367.0387667424478/589.0280369402608,
                                },
                            'TRN__pop':{
                                '0':    394.82515964684/589.0280369402608,
                                '45':   340.35391028959623/589.0280369402608,
                                '90':   340.52806291144435/589.0280369402608,
                                '135':  371.63102179991336/589.0280369402608,
                                '180':  406.7796123474516/589.0280369402608,
                                '225':  478.96760294949587/589.0280369402608,
                                '270':  1.0,
                                '315':  465.9616573671514/589.0280369402608,
                                },}
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

    ################################################################################################################################################
    # c_0003_100_awake
    # simFolder='c_0003'; simName='c_0003_100_awake'; _addVersion='_v04'; conn_templates=[ 
    #                                                                                             'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',
    #                                                                                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25',
    #                                                                                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
    #                                                                                             'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75',
    #                                                                                             'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
    #                                                                                             'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75',
    #                                                                                             'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50',
    #                                                                                             'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25',
    #                                                                                             'VPM_ff_closed__TRN_fb_open__L6A_fb_closed',
    #                                                                                             ]  
    simFolder='paper_sims'; simName='paperSim_Fig05_awake'; _addVersion='_v04'; conn_templates=[ 
                                                                                                'VPM_ff_closed__TRN_fb_open__L6A_fb_closed',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75',
                                                                                                'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25',
                                                                                                'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',
                                                                                                ]  

    # #  -  100% uniform
    # simFolder='c_0003'; simName='c_0003_300_awake_0_uniform_CTmodulation_0';     _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    # simFolder='c_0003'; simName='c_0003_300_awake_0_uniform_CTmodulation_10';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    # simFolder='c_0003'; simName='c_0003_300_awake_0_uniform_CTmodulation_25';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    # simFolder='c_0003'; simName='c_0003_300_awake_0_uniform_CTmodulation_50';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    # simFolder='c_0003'; simName='c_0003_300_awake_0_uniform_CTmodulation_100';   _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    # #  -  50% closed / 50% uniform
    # simFolder='c_0003'; simName='c_0003_300_awake_1_mixed5050_CTmodulation_0';  _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    # simFolder='c_0003'; simName='c_0003_300_awake_1_mixed5050_CTmodulation_10'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    # simFolder='c_0003'; simName='c_0003_300_awake_1_mixed5050_CTmodulation_25'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    # simFolder='c_0003'; simName='c_0003_300_awake_1_mixed5050_CTmodulation_50'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    # simFolder='c_0003'; simName='c_0003_300_awake_1_mixed5050_CTmodulation_100';_addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    # #  -  100% closed
    # simFolder='c_0003'; simName='c_0003_300_awake_2_closed_CTmodulation_0';     _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    # simFolder='c_0003'; simName='c_0003_300_awake_2_closed_CTmodulation_10';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    # simFolder='c_0003'; simName='c_0003_300_awake_2_closed_CTmodulation_25';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    # simFolder='c_0003'; simName='c_0003_300_awake_2_closed_CTmodulation_50';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    # simFolder='c_0003'; simName='c_0003_300_awake_2_closed_CTmodulation_100';   _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]

    ################################################################################################################################################
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig05_awake'; _addVersion='_v04'; conn_templates=[ 
                                                                                                'VPM_ff_closed__TRN_fb_open__L6A_fb_closed',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50',
                                                                                                'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75',
                                                                                                'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
                                                                                                'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25',
                                                                                                'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',
                                                                                                ]  
    #  -  100% uniform
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_0';     _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_10';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_25';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_50';    _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_100';   _addVersion='_v04'; conn_templates=['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',]
    #  -  50% closed / 50% uniform
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_0';  _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_10'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_25'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_50'; _addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_100';_addVersion='_v04'; conn_templates=['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',]
    #  -  100% closed
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_2_closed_ct_modulation_0';     _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_2_closed_ct_modulation_10';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_2_closed_ct_modulation_25';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_2_closed_ct_modulation_50';    _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]
    simFolder='paper_000'; simName='paper_000_barreloid_batch_fig07_2_closed_ct_modulation_100';   _addVersion='_v04'; conn_templates=['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',]

    ################################################################################################################################################

    simFolder='paper_000'; _addVersion='_v04'

    sims_dict = {   
                    'paper_000_barreloid_batch_fig05_awake'                         :   ['VPM_ff_closed__TRN_fb_open__L6A_fb_closed',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_75_25',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_50_50',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_remixed_tr_25_75',
                                                                                         'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_25_75',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_75_25',
                                                                                         'VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    'paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_0'     :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_10'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_25'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_50'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_000_barreloid_batch_fig07_0_uniform_ct_modulation_100'   :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_0'       :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    'paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_10'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    'paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_25'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    'paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_50'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    'paper_000_barreloid_batch_fig07_1_mixed_ct_modulation_100'     :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    'paper_000_barreloid_batch_fig07_2_closed_ct_modulation_0'      :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    'paper_000_barreloid_batch_fig07_2_closed_ct_modulation_10'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    'paper_000_barreloid_batch_fig07_2_closed_ct_modulation_25'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    'paper_000_barreloid_batch_fig07_2_closed_ct_modulation_50'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    'paper_000_barreloid_batch_fig07_2_closed_ct_modulation_100'    :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],       
                    }
    ################################################################################################################################################

    simFolder='paper_001'; _addVersion='_v04'

    sims_dict = {   
                    # 'paper_001_barreloid_batch_fig05_awake'                         :   ['VPM_ff_closed__TRN_fb_open__L6A_fb_uniform',
                    #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_75_25',
                    #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_50_50',
                    #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_25_75',
                    #                                                                      'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                    #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_25_75',
                    #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',
                    #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_75_25',
                    #                                                                      'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',], 
                    'paper_001_barreloid_batch_fig07_0_uniform_ct_modulation_0'     :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_0_uniform_ct_modulation_10'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_0_uniform_ct_modulation_25'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_0_uniform_ct_modulation_50'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_0_uniform_ct_modulation_100'   :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_1_mixed_ct_modulation_0'       :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    'paper_001_barreloid_batch_fig07_1_mixed_ct_modulation_10'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    'paper_001_barreloid_batch_fig07_1_mixed_ct_modulation_25'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    'paper_001_barreloid_batch_fig07_1_mixed_ct_modulation_50'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    'paper_001_barreloid_batch_fig07_1_mixed_ct_modulation_100'     :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    'paper_001_barreloid_batch_fig07_2_closed_ct_modulation_0'      :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_2_closed_ct_modulation_10'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_2_closed_ct_modulation_25'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_2_closed_ct_modulation_50'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    'paper_001_barreloid_batch_fig07_2_closed_ct_modulation_100'    :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],   
                    }

    ################################################################################################################################################

    simFolder='paper_002'; _addVersion='_v04'

    sims_dict = {   
                    'paper_002_barreloid_batch_fig05_awake'                         :   ['VPM_ff_closed__TRN_fb_open__L6A_fb_uniform',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_75_25',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_50_50',
                                                                                         'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_25_75',
                                                                                         'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_25_75',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',
                                                                                         'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_75_25',
                                                                                         'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',], 
                    
                    
                    # 'paper_002_barreloid_batch_fig07e_0_uniform_ct_modulation_0'     :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07e_0_uniform_ct_modulation_10'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07e_0_uniform_ct_modulation_25'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07e_0_uniform_ct_modulation_50'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07e_0_uniform_ct_modulation_100'   :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07e_1_mixed_ct_modulation_0'       :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07e_1_mixed_ct_modulation_10'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07e_1_mixed_ct_modulation_25'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07e_1_mixed_ct_modulation_50'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07e_1_mixed_ct_modulation_100'     :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_remixed_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07e_2_closed_ct_modulation_0'      :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    # 'paper_002_barreloid_batch_fig07e_2_closed_ct_modulation_10'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    # 'paper_002_barreloid_batch_fig07e_2_closed_ct_modulation_25'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    # 'paper_002_barreloid_batch_fig07e_2_closed_ct_modulation_50'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],
                    # 'paper_002_barreloid_batch_fig07e_2_closed_ct_modulation_100'    :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_closed',],       
                    

                    # 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_0'     :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_10'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_25'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_50'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_100'   :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_0'       :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_10'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_25'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_50'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_100'     :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
                    # 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_0'      :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_10'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_25'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_50'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
                    # 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_100'    :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],   

                    }


    # sims_dict = {   
    #                 # 'paper_002_barreloid_batch_fig05_awake'                         :   ['VPM_ff_closed__TRN_fb_open__L6A_fb_uniform',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_75_25',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_50_50',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_openRemixed__L6A_fb_uniform_tr_25_75',
    #                 #                                                                      'VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_25_75',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',
    #                 #                                                                      'VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_75_25',
    #                 #                                                                      'VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',], 
    #                 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_0'     :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_10'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_25'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_50'    :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_0_uniform_ct_modulation_100'   :   ['VPM_ff_uniform__TRN_fb_uniform__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_0'       :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
    #                 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_10'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
    #                 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_25'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
    #                 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_50'      :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
    #                 'paper_002_barreloid_batch_fig07f_1_mixed_ct_modulation_100'     :   ['VPM_ff_remixed__TRN_fb_remixed__L6A_fb_uniform_tr_50_50',],
    #                 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_0'      :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_10'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_25'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_50'     :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],
    #                 'paper_002_barreloid_batch_fig07f_2_closed_ct_modulation_100'    :   ['VPM_ff_closed__TRN_fb_closed__L6A_fb_uniform',],   
    #                 }


    ################################################################################################################################################

    for simName,conn_templates in sims_dict.items():

        grouped_params_0='0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0'
        grouped_params_1='1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1'
        grouped_params_2='2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2_2'
        grouped_params_3='3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3'
        grouped_params_4='4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4_4'
        grouped_params_5='5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5_5'
        grouped_params_6='6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6_6'
        grouped_params_7='7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7_7'

        angular_tuning_windows  = ['20ms']
        ct_feedbacks            = ['ct_uniform']

        figure_counter=0
        
        try: conn_templates = [c+_addVersion for c in conn_templates]
        except: print('Adding conn version failed')

        analysis_folder = '/Users/joao/Research/Models/BBP/thalamus_netpyne/sim_data/'
        # source_folder   = analysis_folder+simFolder+'/'+simName+'/'
        source_folder   = analysis_folder+simFolder+'/'+simName+'/sim_output/'

        for angular_tuning_window in angular_tuning_windows:
            ############################################################################################################################################
            if angular_tuning_window == '20ms':
                file_flag = 'angularTuning3_spikes_ON'
                dataFolder      = source_folder+'angular_tuning_ON/'
            ############################################################################################################################################
            elif angular_tuning_window == '250ms':
                file_flag = 'angularTuning3_spikes_all'
                dataFolder      = source_folder+'angular_tuning_all/'
            ############################################################################################################################################
            for ct_feedback in ct_feedbacks:
                for conn_template_ind, conn_template in enumerate(conn_templates):
                    if ct_feedback == 'ct_uniform':
                        sim_identifier = conn_template_ind
                        angular_tuning_files = [
                            simName+'_'+grouped_params_0+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_1+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_2+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_3+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_4+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_5+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_6+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                            simName+'_'+grouped_params_7+'_'+str(sim_identifier)+'_0_0_0_0_0_0_0_0_0_0_'+extra_zeros+file_flag+'.json',
                        ]
                        barreloid_network_cell_properties = '../network_template/barreloid_network_cell_properties__'+conn_template+'.json'

                    ############################################################################################################################################

                    # --- Processing pooled angular tuning

                    deflection_angles = [0, 45, 90, 135, 180, 225, 270, 315]

                    store_angular_tuning_plot_data_dict={}
                    for angular_tuning_file_ind,angular_tuning_file in enumerate(angular_tuning_files):

                        at_filename = dataFolder+angular_tuning_file
                        with open(at_filename, 'r') as openfile: angular_tuning_dict = json.load(openfile)
                            
                        network_name = barreloid_network_cell_properties.split('barreloid_network_cell_properties__')[1].split('.json')[0]

                        save_fig_name = dataFolder+simName+'__'+network_name+'__'+str(angular_tuning_file_ind)

                        angular_tuning_plot_data_dict = PlotAngularTuning.plotAngularTuning(    angular_tuning_dict, 
                                                                                                barreloid_network_cell_properties=barreloid_network_cell_properties, 
                                                                                                center_point=500,
                                                                                                select_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop', 'L6A_activated__pop'],
                                                                                                save_fig_name = save_fig_name,
                                                                                                select_deflection = ['all'], # deflections: ['0', '1', '2', '3', '4', '5', ... , 'all'] - 'all' combines all deflections into a single plot
                                                                                                deflection_angle=deflection_angles[angular_tuning_file_ind]
                                                                                                )
                        
                        angular_tuning_plot_data_dict = PlotAngularTuning.add_rotated_angularTuning_data(angular_tuning_plot_data_dict,deflection_angle=target_angles[angular_tuning_file_ind],target_rotation_angle=270)

                        store_angular_tuning_plot_data_dict.update({target_angles[angular_tuning_file_ind]:angular_tuning_plot_data_dict})

                    # Save plotted data
                    with open(dataFolder+simName+'__'+network_name+'_plot_angular_tuning_data'+'.json', 'w') as fp: json.dump(store_angular_tuning_plot_data_dict, fp, indent=4)
                    
                    ########################################################################################################################
                    # --- Processing normalized values 
                    mean_angular_tuning_dict = PlotAngularTuning.combine_rotated_angularTuning_data(store_angular_tuning_plot_data_dict)
                    # PlotAngularTuning.plotCombinedAngularTuning(mean_angular_tuning_dict, angular_tuning_hartings, save_fig_name=dataFolder+simName+'__'+network_name+'_normalized')
                    PlotAngularTuning.plotCombinedAngularTuning(mean_angular_tuning_dict, angular_tuning_hartings, save_fig_name=dataFolder+'_normalized_'+str(figure_counter)+'__'+simName+'__'+network_name) 

                    deflection_ordered_dict_normalized = PlotAngularTuning.parseTuningValues(store_angular_tuning_plot_data_dict, data_type = 'normalized')
                    grouped_normalized_data = {}
                    for pop_name in deflection_ordered_dict_normalized.keys():
                        grouped_normalized_data.update({pop_name:{'angle':{}}})
                        for deflection_angle in deflection_ordered_dict_normalized[pop_name].keys():
                            for key_angle_ind, key_angle in enumerate(deflection_ordered_dict_normalized[pop_name][deflection_angle]['keys_rotated']):
                                if key_angle not in grouped_normalized_data[pop_name]['angle'].keys(): grouped_normalized_data[pop_name]['angle'].update({key_angle:[]})
                                grouped_normalized_data[pop_name]['angle'][key_angle].append(deflection_ordered_dict_normalized[pop_name][deflection_angle]['values'][key_angle_ind])

                    # stats
                    from scipy.stats import ttest_1samp
                    for pop_name in deflection_ordered_dict_normalized.keys():
                        # print(pop_name)
                        for key_angle in grouped_normalized_data[pop_name]['angle'].keys():
                            # Perform the one-sample t-test
                            t_stat, p_value = ttest_1samp(grouped_normalized_data[pop_name]['angle'][key_angle], angular_tuning_hartings_meanVals[pop_name][str(key_angle)])
                            if      p_value>0.05:                           significance = 'ns'
                            elif    (p_value<=0.05)  and (p_value>0.005):   significance = '*'
                            elif    (p_value<=0.005) and (p_value>0.0005):  significance = '**'
                            elif    (p_value<=0.0005):                      significance = '***'

                            # print('\t\t\t',key_angle, '\t',p_value, '\t',significance, '\t (',t_stat,', ', p_value,')')

                            # print('Input values: ', grouped_normalized_data[pop_name]['angle'][key_angle])
                            # print('Pop mean value: ', angular_tuning_hartings[pop_name][str(key_angle)])

                    # --- Testing processing absolute values 
                    mean_angular_tuning_dict_absolute = PlotAngularTuning.combine_rotated_angularTuning_data_absolute(store_angular_tuning_plot_data_dict)
                    # PlotAngularTuning.plotCombinedAngularTuning(mean_angular_tuning_dict_absolute, angular_tuning_hartings, save_fig_name=dataFolder+simName+'__'+network_name+'_absolute')
                    PlotAngularTuning.plotCombinedAngularTuning(mean_angular_tuning_dict_absolute, angular_tuning_hartings, save_fig_name=dataFolder+'_absolute_'+str(figure_counter)+'__'+simName+'__'+network_name)

                    deflection_ordered_dict_absolute = PlotAngularTuning.parseTuningValues(store_angular_tuning_plot_data_dict, data_type = 'absolute')
                    grouped_absolute_data = {}
                    for pop_name in deflection_ordered_dict_absolute.keys():
                        grouped_absolute_data.update({pop_name:{'angle':{}}})
                        for deflection_angle in deflection_ordered_dict_absolute[pop_name].keys():
                            for key_angle_ind, key_angle in enumerate(deflection_ordered_dict_absolute[pop_name][deflection_angle]['keys_rotated']):
                                if key_angle not in grouped_absolute_data[pop_name]['angle'].keys(): grouped_absolute_data[pop_name]['angle'].update({key_angle:[]})
                                grouped_absolute_data[pop_name]['angle'][key_angle].append(deflection_ordered_dict_absolute[pop_name][deflection_angle]['values'][key_angle_ind])

                    ########################################################################################################################

                    # Extra analysis - processing data to send to the script that pools angular tuning amongst group and run statistical analysis
                    for dataset_type in ['normalized', 'absolute']:

                        store_pooled_data={}

                        angular_tuning_pops = ['MLe__pop', 'VPM__pop', 'TRN__pop']
                        for at_pop in angular_tuning_pops:
                            store_pooled_data.update({at_pop:{}})
                            for deflection_angle in store_angular_tuning_plot_data_dict.keys():
                                for at_angle_ind, at_angle in enumerate(store_angular_tuning_plot_data_dict[deflection_angle]['all'][at_pop][dataset_type]['keys_rotated']):
                                    if at_angle not in  store_pooled_data[at_pop].keys(): store_pooled_data[at_pop].update({at_angle:[]})

                                    store_pooled_data[at_pop][at_angle].append(store_angular_tuning_plot_data_dict[deflection_angle]['all'][at_pop][dataset_type]['values'][at_angle_ind])

                        # Save angular tuning statistics
                        print(f' - Angular tuning window: {angular_tuning_window} | CT feedback: {ct_feedback} | Conn: {conn_template}')
                        with open(dataFolder+'angular_tuning_Statistics__'+dataset_type+'__'+angular_tuning_window+'__'+ct_feedback+'__conn_'+conn_template+'.json', 'w') as fp: json.dump(store_pooled_data, fp, indent=4)

                    figure_counter+=1


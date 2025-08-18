'''
Class to build the stimulation dataset

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

########################################################################################################################
from cProfile import label
import numpy as np
import json
import math
import os
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 40})
import csv
import pandas as pd
from scipy.optimize import curve_fit
########################################################################################################################

class LoadPSTH():
    def load(folderPath, fileNames,plotFigs=False,save_path=''):
        hist_dict={}
        for fileName in fileNames:
            key_name=fileName.split('.')[0]
            hist_dict.update({key_name:{}})
            csvFilePath = folderPath+fileName
        # Initialize an empty list to store the data
            data = []
        # Open the CSV file and read its contents
            with open(csvFilePath, 'r') as file:
                csv_reader = csv.reader(file)
            # Iterate through the rows and append them to the data list
                for row in csv_reader: data.append(row)

            store_x=[]
            store_y=[]
            subtract_first_point=True
            for ind,[x,y] in enumerate(data):
                if key_name!='0stim':   # the first point of the 'stim' signal is not ignored because it lies in the data curve, differently from PrV and VPM histograms
                    if ind==0:continue  # ignores the first point, which is the 'zero reference'
                if subtract_first_point:
                    x0 = float(data[0][0])
                    y0 = float(data[0][1])
                else:
                    x0=0
                    y0=0
                store_x.append(float(x)-x0)
                store_y.append(float(y)-y0)

            hist_dict[key_name].update({'x':store_x,'y':store_y})
            if plotFigs:
                # # Now, 'data' contains the CSV data as a list of lists    
                # Raw data
                plt.figure(figsize=(25,10))
                plt.plot(store_x,store_y,color='red')
                plt.bar(store_x,store_y,color='k')
                plt.ylim([-0.05, 0.4])
                plt.xlabel('Time (ms)')
                plt.ylabel('spike probability - (0.1/ms)')
                plt.savefig(save_path+'1_'+key_name+'.png')

                store_dx=[]
                store_dy=[]
                for i,x in enumerate(store_x):
                    if i>0:dx=store_x[i]-store_x[i-1]
                    else:continue
                    dy=(store_y[i]-store_y[i-1])/dx

                    store_dx.append(dx)
                    store_dy.append(dy)
                hist_dict[key_name].update({'dx':store_dx,'dy':store_dy})

            # First Derivative data
                plt.figure(figsize=(25,20))
                plt.plot(store_x[1::],store_dy,color='k')
                plt.bar(store_x[1::],store_dy)
                plt.ylim([-0.05, 0.4])
                plt.savefig(save_path+'2_'+key_name+'_diff.png')

            # Colors
                if key_name=='PrV':   
                    print('im here')  
                    stim_1=63
                    stim_2=203
                elif key_name=='0stim':   
                    stim_1=38
                    stim_2=122
                else:                   
                    stim_1=63
                    stim_2=191

                print(key_name,stim_1,stim_2)

            # Colored data
                plt.figure(figsize=(25,20))
                for i,x in enumerate(store_x):
                    if i<stim_1:c='r'
                    elif (i>=stim_1) and (i<stim_2):c='b' # 192 and 204
                    else: c='g'
                    plt.bar(x,store_y[i],color=c)
                plt.ylim([-0.05, 0.4])
                plt.savefig(save_path+'3_'+key_name+'_colored.png')

            # Remove baseline
                y_baseline = np.mean(store_y[0:stim_1])
                store_y_noBaseline = [y-y_baseline for y in store_y]
                plt.figure(figsize=(25,20))
                for i,x in enumerate(store_x):
                    if i<stim_1:c='r'
                    elif (i>=stim_1) and (i<stim_2):c='b'
                    else: c='g'
                    plt.bar(x,store_y_noBaseline[i],color=c)
                plt.ylim([-0.05, 0.4])
                plt.savefig(save_path+'4_'+key_name+'_colored_noBaseline.png')
        
        return hist_dict
    
    def formatData(hist_dict,target_nuclei, ON_start, OFF_start, skip_ON, skip_OFF):
        PrV_data={  'baseline'  :{ 'x':hist_dict[target_nuclei]['x'][0:ON_start],                   
                                   'y':hist_dict[target_nuclei]['y'][0:ON_start],},
                    'ON'        :{ 'x':hist_dict[target_nuclei]['x'][ON_start+skip_ON:OFF_start],   
                                   'y':hist_dict[target_nuclei]['y'][ON_start+skip_ON:OFF_start],},
                    'OFF'       :{ 'x':hist_dict[target_nuclei]['x'][OFF_start+skip_OFF::],         
                                   'y':hist_dict[target_nuclei]['y'][OFF_start+skip_OFF::],},
                    'all'       :{ 'x':hist_dict[target_nuclei]['x'],                               
                                   'y':hist_dict[target_nuclei]['y'],},}
        return PrV_data

########################################################################################################################

class FitData():
    def fitPolynomial(x_data,y_data,degree): return np.polyfit(x_data, y_data, degree)
    # def getPolynomialObj(polyFit):           return np.poly1d(polyFit)

    # New methods to improve PSTH data fitting
    @staticmethod
    def fit_double_exp(x_data, y_data, initial_guess=(1, 0.01, 1, 0.01, 0)):
        def double_exp_decay(t, a, b, c, d, f):
            return a * np.exp(-b * t) + c * np.exp(-d * t) + f
        popt, _ = curve_fit(double_exp_decay, x_data, y_data, p0=initial_guess)
        return popt

    @staticmethod
    def fit_triple_exp(x_data, y_data, initial_guess=(1, 0.01, 1, 0.01, 1, 0.01, 0)):
        def triple_exp_decay(t, a, b, c, d, e, f, g):
            return a * np.exp(-b * t) + c * np.exp(-d * t) + e * np.exp(-f * t) + g
        popt, _ = curve_fit(triple_exp_decay, x_data, y_data, p0=initial_guess)
        return popt

    @staticmethod
    def getFittedModel(params, model_type):
        if model_type == 'polynomial':
            return np.poly1d(params)
        elif model_type == 'double_exp':
            def double_exp_decay(t):
                a, b, c, d, f = params
                return a * np.exp(-b * t) + c * np.exp(-d * t) + f
            return double_exp_decay
        elif model_type == 'triple_exp':
            def triple_exp_decay(t):
                a, b, c, d, e, f, g = params
                return a * np.exp(-b * t) + c * np.exp(-d * t) + e * np.exp(-f * t) + g
            return triple_exp_decay
        else:
            raise ValueError("Unsupported model type. Choose 'polynomial', 'double_exp', or 'triple_exp'.")


########################################################################################################################

class SampleData():
    # --- Fit sample from baseline
    def sampleNormal(mean, std, start,stop, dt):
        sampled_dict={'x':[],'y':[]}
        sampled_dict['x'] = list(np.arange(start=start,stop=stop,step=dt))
        sampled_dict['y'] = [np.random.normal(mean,std) for i in range(len(sampled_dict['x']))]
        return sampled_dict
    
    # --- Fit sample from ON and OFF polynomial fits
    def sampleFitFunction(fitFcn, timeSeries, dt,resetTime=False):
        sampled_dict={'x':[],'y':[]}
        for t_sample in np.arange(math.floor(min(timeSeries)),math.ceil(max(timeSeries)),dt):
            sampled_dict['x'].append(t_sample)
            if resetTime: t_value = t_sample-math.floor(min(timeSeries))
            else:         t_value = t_sample
            sampled_dict['y'].append(fitFcn(t_value))
        return sampled_dict
    
    # --- Fit sample from ON and OFF polynomial fits outside the fitting window (e.g.: t>100 for a polynomial thats fit from 0-100 ms)
    #  -  Uses the last probability value on the fitFcn
    def sampleFitFunction_end(fitFcn, end_time, timeSeries, dt,resetTime=False):
        sampled_dict={'x':[],'y':[]}
        for t_sample in np.arange(math.floor(min(timeSeries)),math.ceil(max(timeSeries)),dt):
            sampled_dict['x'].append(t_sample)
            sampled_dict['y'].append(fitFcn(end_time))
        return sampled_dict


    def sampleRaster(data_y_sampled, target_angle, n_cells, prob_compensation, yPoly_angTuning=None):
        if yPoly_angTuning is not None: 
            print('---> Runing angular tuned sampling')
            # yFit_angTuning = FitData.getPolynomialObj(yPoly_angTuning)
            yFit_angTuning = FitData.getFittedModel(params=yPoly_angTuning, model_type='polynomial')
        raster_sampled=[]
        for prob_ind,prob in enumerate(data_y_sampled):
            # --- scaling probability to simulate angular tuning
            if yPoly_angTuning is not None:
                random_var=[]
                for cell_ind in range(n_cells):
                    cell_deg = int(cell_ind*(360/n_cells))
                    x_val = abs(4-(4*abs(cell_deg-target_angle-180)/(180)))
                    y_val = yFit_angTuning(x_val)
                    spike=np.random.binomial(1,prob*y_val*prob_compensation,1)
                    random_var.append(int(spike[0])) # only the 0th element because is a single cell, that will return a 0 or 1
            else:
                # --- no probability scaling
                random_var = np.random.binomial(1,prob,n_cells)

            # --- Stores the raster plot info
            raster_sampled.append(list(random_var))
        return raster_sampled

    # def sampleRaster2(data_y_sampled, deflection_times, dt, target_angle, n_cells, prob_compensation, yPoly_angTuning=None):
    def sampleRaster2(data_y_sampled, deflection_times, dt, target_angle, n_cells, prob_compensation, yPoly_angTuning=None, overshoot=100, max_degree=360):
        '''
        # overshoot: - interval (ms) after release of whisker that should show angular tuning
                     - added to prevent uniform whisker response during OFF stim
        '''
        if yPoly_angTuning is not None: 
            print('---> Runing angular tuned sampling')
            # yFit_angTuning = FitData.getPolynomialObj(yPoly_angTuning)
            yFit_angTuning = FitData.getFittedModel(params=yPoly_angTuning, model_type='polynomial')
        raster_sampled=[]
        for prob_ind,prob in enumerate(data_y_sampled):
            if (prob_ind*dt)%1000==0:print(prob_ind*dt)
            # --- Checks if the sample is within the deflection ranges, otherwise samples uniformly
            # if any(r[0] <= (prob_ind*dt) <= r[1] for r in deflection_times):
            if any(r[0] <= (prob_ind*dt) <= r[1]+overshoot for r in deflection_times):
                # print('-> ', prob_ind*dt, '\tangular tuned sample')
                # --- scaling probability to simulate angular tuning
                if yPoly_angTuning is not None:
                    random_var=[]
                    for cell_ind in range(n_cells):
                        cell_deg = int(cell_ind*(max_degree/n_cells))
                        x_val = abs(4-(4*abs(cell_deg-target_angle-(max_degree/2))/((max_degree/2))))
                        y_val = yFit_angTuning(x_val)
                        spike=np.random.binomial(1,prob*y_val*prob_compensation,1)
                        random_var.append(int(spike[0])) # only the 0th element because is a single cell, that will return a 0 or 1
                else:
                    # --- no probability scaling
                    random_var = np.random.binomial(1,prob,n_cells)
            else:
                # --- no probability scaling
                random_var = np.random.binomial(1,prob,n_cells)

            # --- Stores the raster plot info
            raster_sampled.append(list(random_var))
        return raster_sampled
    
    def countSpikes(n_cells, raster_sampled):
        store_cell_spks=[0 for i in range(n_cells)]
        for time_bin_ind,time_bin in enumerate(raster_sampled):
            for cell_ind,cell_spk in enumerate(time_bin):
                if raster_sampled[time_bin_ind][cell_ind] == 1: store_cell_spks[cell_ind]+=raster_sampled[time_bin_ind][cell_ind]
        return store_cell_spks
    
    def getSpikeTimes(raster_sampled,time_vector,n_cells):
        store_spkts={}
        for cell_gid in range(n_cells):store_spkts.update({cell_gid:[]})
        for sample_ind, sample in enumerate(raster_sampled):
            for cell_gid in range(n_cells):
                if raster_sampled[sample_ind][cell_gid]==1:store_spkts[cell_gid].append(time_vector[sample_ind])
        return store_spkts


    def generateProbabilityArray(deflection_times, baseline_prob, ON_polyFit, OFF_polyFit, resetTime, dt):
        ON_polyFitObj  = FitData.getFittedModel(params=ON_polyFit,  model_type='triple_exp')
        OFF_polyFitObj = FitData.getFittedModel(params=OFF_polyFit, model_type='double_exp')
        # ON_polyFitObj   = FitData.getPolynomialObj(ON_polyFit)
        # OFF_polyFitObj  = FitData.getPolynomialObj(OFF_polyFit)

        prob_array=[]
        events_array=[] # list to record the ON/OFF times
        # sampledProb={}
        print('deflection_times ', deflection_times)
        if deflection_times[0][0]>0: 
            print('0')
            baseline_range = [0,deflection_times[0][0]]
            prob_array.append(SampleData.sampleNormal(mean=baseline_prob[0],std=baseline_prob[1],start=baseline_range[0],stop=baseline_range[1],dt=dt))
            # sampledProb.update({['baseline']:SampleData.sampleNormal(mean=baseline_prob[0],std=baseline_prob[1],dur=baseline_range,dt=dt)})

            # --- Record zeros for baseline
            events_array.append([0 for i in np.arange(start=baseline_range[0],stop=baseline_range[1],step=dt)])

        for ind,deflection_time in enumerate(deflection_times):            
            # --- If interval between ON and OFF is < 200
            if deflection_times[ind][1]-deflection_times[ind][0]<200:
                print('deflection ON')
                prob_array.append(SampleData.sampleFitFunction(ON_polyFitObj , deflection_times[ind], dt, resetTime=resetTime))
            else:
                print('deflection ON - (> 200 ms)')
                prob_array.append(SampleData.sampleFitFunction(ON_polyFitObj , [deflection_times[ind][0],deflection_times[ind][0]+200], dt, resetTime=resetTime))
                prob_array.append(SampleData.sampleNormal(mean=ON_polyFitObj(200),std=baseline_prob[1],start=deflection_times[ind][0]+200,stop=deflection_times[ind][1],dt=dt))
                
            # --- If interval between stims is > 100
            if ind<len(deflection_times)-1:
                if deflection_times[ind+1][0]-deflection_times[ind][1]<100:
                    print('deflection OFF')
                    prob_array.append(SampleData.sampleFitFunction(OFF_polyFitObj, [deflection_times[ind][1],deflection_times[ind+1][0]], dt, resetTime=resetTime))
                else:
                    print('deflection OFF - (new deflection in interval < 100 ms)')
                    prob_array.append(SampleData.sampleFitFunction(OFF_polyFitObj, [deflection_times[ind][1],deflection_times[ind][1]+100], dt, resetTime=resetTime))
                    prob_array.append(SampleData.sampleNormal(mean=OFF_polyFitObj(100),std=baseline_prob[1],start=deflection_times[ind][1]+100,stop=deflection_times[ind+1][0],dt=dt))

            else:
                print('deflection OFF')
                prob_array.append(SampleData.sampleFitFunction(OFF_polyFitObj, [deflection_times[ind][1],deflection_times[ind][1]+100], dt, resetTime=resetTime))
                # prob_array.append(SampleData.sampleNormal(mean=baseline_prob[0],std=baseline_prob[1],start=deflection_times[ind][1]+100,stop=deflection_times[ind+1][0],dt=dt))

        return prob_array
    
    def deflectionEvent(deflection_dict):
        # import sys; sys.exit()
        print(' - Deflection event')
        # --- Local variables
        degree_to_cell = int(deflection_dict['angle_range'][1]/deflection_dict['n_cells'])

        prob_array = SampleData.generateProbabilityArray(   deflection_times    = deflection_dict['deflection_times'],
                                                            baseline_prob       = deflection_dict['baseline_prob'],
                                                            ON_polyFit          = deflection_dict['ON_polyFit'],
                                                            OFF_polyFit         = deflection_dict['OFF_polyFit'],
                                                            resetTime           = deflection_dict['resetTime'],
                                                            dt                  = deflection_dict['dt'])
        
        # --- Removing because it might be slowing down simulations
        # if deflection_dict['plotFigs']: PlotFigs.plotProbArray(prob_array, save_path=deflection_dict['save_figs_folder'])

        prob_hist={'x':[],'y':[]}
        for i in range(len(prob_array)):
            prob_hist['x']+=prob_array[i]['x']
            prob_hist['y']+=prob_array[i]['y']

        print('Rescaling probabilities to match netpyne higher sampling rate - (thus, dataset prob should be divided by %f)'%(int(deflection_dict['paper_sampling']/deflection_dict['dt'])))
        rescale_prob = [prob/(int(deflection_dict['paper_sampling']/deflection_dict['dt'])) for prob in prob_hist['y']]
        prob_hist['y']=rescale_prob
        
        # raster_sampled = SampleData.sampleRaster(   data_y_sampled      = prob_hist['y'],
        #                                             target_angle        = deflection_dict['target_angle'],
        #                                             n_cells             = deflection_dict['n_cells'],
        #                                             prob_compensation   = deflection_dict['prob_compensation'],
        #                                             yPoly_angTuning     = deflection_dict['yPoly_angTuning'])
        
        raster_sampled = SampleData.sampleRaster2(  data_y_sampled      = prob_hist['y'],
                                                    deflection_times    = deflection_dict['deflection_times'],
                                                    dt                  = deflection_dict['dt'],
                                                    target_angle        = deflection_dict['target_angle'],
                                                    n_cells             = deflection_dict['n_cells'],
                                                    prob_compensation   = deflection_dict['prob_compensation'],
                                                    yPoly_angTuning     = deflection_dict['yPoly_angTuning'],
                                                    max_degree          = deflection_dict['angle_range'][1],
                                                    )
        
        sort_spikes=True
        if sort_spikes and deflection_dict['target_angle'] is not None: # skips sorting if target_angle = None
            print('Sorting raster_sampled')

            # # # Creates a list of angle values based on their distance from the target angle (e.g. tgt:180 - list:[180,178,182,176,184,...,0,358])
            # # # sorted_values = sorted(list(range(0, 360, int(360/deflection_dict['n_cells']))), key=lambda x: abs(4-(4*abs(x-deflection_dict['target_angle']-180)/(180))))
            # # sorted_values = sorted(list(range(0, 360, int(360/deflection_dict['n_cells']))), key=lambda x: min(abs(x - deflection_dict['target_angle']), 360 - abs(x - deflection_dict['target_angle'])))
            # sorted_values = sorted(list(range(deflection_dict['angle_range'][0], deflection_dict['angle_range'][1], int(deflection_dict['angle_range'][1]/deflection_dict['n_cells']))), key=lambda x: min(abs(x - deflection_dict['target_angle']), deflection_dict['angle_range'][1] - abs(x - deflection_dict['target_angle'])))

            # --- Testing sorting only a slice of the raster (if distance from target angle <120, else dont sort)
            
            # threshold_distance=45
            # threshold_distance=90 # version that works
            threshold_distance=135

            # Generate the list of angles
            angles = list(range(
                deflection_dict['angle_range'][0],
                deflection_dict['angle_range'][1],
                int(deflection_dict['angle_range'][1] / deflection_dict['n_cells'])
            ))
            # Sort values within the threshold distance
            sorted_within_threshold = sorted(
                [x for x in angles if min(
                    abs(x - deflection_dict['target_angle']),
                    deflection_dict['angle_range'][1] - abs(x - deflection_dict['target_angle'])
                ) <= threshold_distance],
                key=lambda x: min(
                    abs(x - deflection_dict['target_angle']),
                    deflection_dict['angle_range'][1] - abs(x - deflection_dict['target_angle'])
                )
            )
            # # Add the unsorted values (outside the threshold) back to the list
            # unsorted_values = [x for x in angles if x not in sorted_values]
            # Randomize the unsorted values (outside the threshold)
            unsorted_values = [x for x in angles if x not in sorted_within_threshold]
            import random
            random.shuffle(unsorted_values)
            sorted_within_threshold += unsorted_values

            sorted_values=sorted_within_threshold

            # # Uniformly distribute the unsorted values (outside the threshold)
            # unsorted_values = [x for x in angles if x not in sorted_within_threshold]
            # indices = np.linspace(0, len(sorted_within_threshold), len(unsorted_values) + 1, dtype=int)[1:-1]  # Determine insertion indices
            # sorted_values = sorted_within_threshold.copy()

            # for i, val in zip(indices, unsorted_values):
            #     sorted_values.insert(i, val)
            # --- END - Testing sorting only a slice of the raster (if distance from target angle <120, else dont sort)

            store_raster_deflection={}
            store_raster_deflection_index=[]
            store_raster_deflection_list=[]
            for sample_ind,sample_spks in enumerate(raster_sampled):
                if any(r[0] <= (sample_ind*deflection_dict['dt']) <= r[1] for r in deflection_dict['deflection_times']):
                    store_raster_deflection.update({sample_ind:sample_spks})
                    store_raster_deflection_index.append(sample_ind)
                    store_raster_deflection_list.append(sample_spks)
            
            slice_raster_deflection_dict = {}
            for i in range(len(store_raster_deflection_list[0])):
                cell_store_raster_deflection_list=[]
                for j in range(len(store_raster_deflection_list)):
                    cell_store_raster_deflection_list.append(store_raster_deflection_list[j][i])
                # slice_raster_deflection_dict.update({(i*int(360/deflection_dict['n_cells'])):cell_store_raster_deflection_list})
                slice_raster_deflection_dict.update({(i*int(deflection_dict['angle_range'][1]/deflection_dict['n_cells'])):cell_store_raster_deflection_list})

            # Counts the number of spikes for each cell and makes a list of tuples with the format (cell_number:num_spikes)
            # count_slice_raster = [(ind*int(360/deflection_dict['n_cells']),sum(vals)) for ind,vals in enumerate(slice_raster_deflection_dict.values())]
            count_slice_raster = [(ind*int(deflection_dict['angle_range'][1]/deflection_dict['n_cells']),sum(vals)) for ind,vals in enumerate(slice_raster_deflection_dict.values())]

            # Sorts the list of tuples based on the num_spikes - IMPORTANT: with the higher number of spikes first (which is used for resorting the data)
            sort_count_slice_raster = sorted(count_slice_raster, key=lambda x: x[1], reverse=True)

            # Creates a list of tuples that appends sorted_values in front of sort_count_slice_raster, creating a tuple with the (new index,old index, spk count), used for reassigning the spikes
            sort_count_slice_raster_indexed=[]
            for ind,(old_index,spike_count) in enumerate(sort_count_slice_raster):
                sort_count_slice_raster_indexed.append((sorted_values[ind],old_index,spike_count))

            # Reassingns the spikes for each cell based on the distance from the target angle, resulting in a sorted version of slice_raster_dict
            resorted_slice_raster_dict = {new_idx:slice_raster_deflection_dict[old_idx] for (new_idx,old_idx,spike_count) in sort_count_slice_raster_indexed}

            # Formatting the data back into a list of lists for each dt
            resorted_raster_sampled=[]
            for spkind in range(len(resorted_slice_raster_dict[0])):
                cell_indexes = list(resorted_slice_raster_dict.keys())
                cell_indexes.sort()
                resorted_raster_sampled.append([resorted_slice_raster_dict[cell_ind][spkind] for cell_ind in cell_indexes])
            
            # Changing the sorted spike lists into the sampled raster
            for idx, sample_ind in enumerate(store_raster_deflection_index): raster_sampled[sample_ind] = resorted_raster_sampled[idx]
            
            print('Sorting finished!')
            deflection_dict['maximize_target_angle']=False
            print(' --- maximize_target_angle = ', deflection_dict['maximize_target_angle'])

        store_cell_spks = SampleData.countSpikes(   n_cells             = deflection_dict['n_cells'],
                                                    raster_sampled      = raster_sampled)

        ####### This section is skipped in case of sort_spikes=True - no resampling needed
        if deflection_dict['maximize_target_angle']:
            if store_cell_spks[int(deflection_dict['target_angle']/degree_to_cell)+(deflection_dict['target_angle']%degree_to_cell)]>=max(store_cell_spks):
                print('Response peaked at target angle')
            else:
                attempt=0
                while store_cell_spks[int(deflection_dict['target_angle']/degree_to_cell)+(deflection_dict['target_angle']%degree_to_cell)]<max(store_cell_spks):
                    print('Repicking raster until target angle peaks - attempt # ', attempt)
                    # raster_sampled = SampleData.sampleRaster(   data_y_sampled      = prob_hist['y'],
                    #                                             target_angle        = deflection_dict['target_angle'],
                    #                                             n_cells             = deflection_dict['n_cells'],
                    #                                             prob_compensation   = deflection_dict['prob_compensation'],
                    #                                             yPoly_angTuning     = deflection_dict['yPoly_angTuning'])
                    
                    raster_sampled = SampleData.sampleRaster2(  data_y_sampled      = prob_hist['y'],
                                                                deflection_times    = deflection_dict['deflection_times'],
                                                                dt                  = deflection_dict['dt'],
                                                                target_angle        = deflection_dict['target_angle'],
                                                                n_cells             = deflection_dict['n_cells'],
                                                                prob_compensation   = deflection_dict['prob_compensation'],
                                                                yPoly_angTuning     = deflection_dict['yPoly_angTuning'])

                    store_cell_spks = SampleData.countSpikes(   n_cells             = deflection_dict['n_cells'],
                                                                raster_sampled      = raster_sampled)
                    attempt+=1
                print('Picked a raster - spikes at ',deflection_dict['target_angle'],': ',store_cell_spks[int(deflection_dict['target_angle']/degree_to_cell)+(deflection_dict['target_angle']%degree_to_cell)],' | max spikes: ',max(store_cell_spks))
        
        # --- Converting raster into spike times
        store_spkts = SampleData.getSpikeTimes(raster_sampled,time_vector=prob_hist['x'],n_cells=deflection_dict['n_cells'])

        # --- Store raster data in JSON file
        if deflection_dict['save_stim_data']: 
            spikes_dict={   'deflection_dict':  deflection_dict,
                            'target_angle':     deflection_dict['target_angle'],
                            'spike_events':     raster_sampled, 
                            'time_vector':      prob_hist['x'],
                            'spike_prob':       prob_hist['y'],
                            'spike_count':      store_cell_spks,
                            'spkts':            store_spkts,
                            'gids':             [gid for gid in range(deflection_dict['n_cells'])], 
                            'dt':               deflection_dict['dt']
                            }
            # if save_format=='json':
            #     with open(deflection_dict['save_deflection_model_folder']+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us.json', 'w') as fp: json.dump(spikes_dict, fp, indent=4)
            # else:
            #     f_name=deflection_dict['save_deflection_model_folder']+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us.pkl'
            #     with open(f_name, 'wb') as f: pd.to_pickle(spikes_dict, f)

        return spikes_dict

    def LoadJSON(file_path):
        # --- Load the JSON file into a Python dictionary
        with open(file_path, 'r') as fp: loaded_dict = json.load(fp)
        return loaded_dict
    
    def LoadPickle(file_path):
        # --- Load the Pickle file into a Python dictionary
        with open(file_path, 'rb') as f: loaded_dict = pd.read_pickle(f)
        return loaded_dict

    def analyzeDeflectionEvent(spikes_dict):

        # --- You can now access the values in the spikes_dict as needed
        raster_sampled      = spikes_dict['spike_events']
        store_cell_spks     = spikes_dict['spike_count']

        # --- Separates spikes within deflections vs spikes during baseline
        raster_deflection=[]; raster_baseline=[]
        for spike_sample_ind, spike_sample in enumerate(raster_sampled):
            if any(r[0] <= (spike_sample_ind*spikes_dict['dt']) <= r[1] for r in spikes_dict['deflection_dict']['deflection_times']): raster_deflection.append(spike_sample)
            else: raster_baseline.append(spike_sample)

        store_cell_spks_deflection = SampleData.countSpikes(n_cells             = spikes_dict['deflection_dict']['n_cells'],
                                                            raster_sampled      = raster_deflection)
        store_cell_spks_baseline   = SampleData.countSpikes(n_cells             = spikes_dict['deflection_dict']['n_cells'],
                                                            raster_sampled      = raster_baseline)

        deflection_dict = spikes_dict['deflection_dict']
        degree_to_cell = int(360/deflection_dict['n_cells'])

        save_path = spikes_dict['deflection_dict']['save_figs_folder']
    
        print('\t--- Plotting figures for target angle of ', deflection_dict['target_angle'])
        
        if deflection_dict['target_angle'] is not None:
            max_tuning = int(deflection_dict['target_angle']/degree_to_cell)+(deflection_dict['target_angle']%degree_to_cell)
            min_tuning = int((deflection_dict['target_angle']+180)/degree_to_cell)+(deflection_dict['target_angle']%degree_to_cell)
            if min_tuning>=deflection_dict['n_cells']:min_tuning-=deflection_dict['n_cells']
            
            print('Maximum firing for cell: ', max_tuning)
            print('Minimum firing for cell: ', min_tuning)
        else:
            max_tuning = 0
            min_tuning = 0
        ##############################
        # --- Adding the connections based on quadrants
        if deflection_dict['target_angle'] is not None:
            if deflection_dict['target_angle']%2==1:   quadrant_shift=[21,23]
            else:                   quadrant_shift=[22,22]
            quadrant_ = [deflection_dict['target_angle']+angle_shift*45 for angle_shift in np.arange(-16,16)]
        else: 
            quadrant_shift = [23,21]
            quadrant_ =  [-540,-495,-450,-405,-360,-315,-270,-225,-180,-135,-90,-45,0,45,90,135,180,225,270,315,360,405,450,495,540]
        # if deflection_dict['target_angle'] is not None:
        #     if deflection_dict['target_angle']%2==1:   quadrant_shift=[45,45]
        #     else:                   quadrant_shift=[44,46]
        #     quadrant_ = [deflection_dict['target_angle']+angle_shift*90 for angle_shift in np.arange(-16,16)]
        # else: 
        #     quadrant_shift = [44,46]
        #     quadrant_ =  [-540,-450,-360,-270,-180,-90,0,90,180,270,360,450,540]

        quadrant = [q for q in quadrant_ if (q>=0) and (q<360)]
        # --- Range of values to sum angles
        # quadrant_ranges = [[q-quadrant_shift[0],q+quadrant_shift[1]] for q in quadrant]
        quadrant_ranges=[]
        for q in quadrant:
            if q%2==1:  quadrant_ranges.append([q-21,q+23])
            else:       quadrant_ranges.append([q-22,q+22])
        
        # --- List of values to sum angles
        quadrant_values_=[]
        for quadrant_range in quadrant_ranges: quadrant_values_.append([q for q in np.arange(quadrant_range[0],quadrant_range[1],degree_to_cell)])
        # --- Removing values (< 0) and (>= 360)
        quadrant_values=[]
        for values in quadrant_values_:
            values_=[]
            for val in values:
                if val<0:    val+=360
                if val>=360: val-=360
                values_.append(val)
            quadrant_values.append(values_)
        # --- Adding the values in the target ranges
        sum_quadrant_values=[]
        sum_quadrant_values_deflection=[]
        sum_quadrant_values_baseline=[]
        for quad_vals in quadrant_values:
            sum_quadrant            = 0
            sum_quadrant_deflection = 0
            sum_quadrant_baseline   = 0
            for q_val in quad_vals:
                sum_quadrant            += store_cell_spks[int(q_val/degree_to_cell)+(q_val%degree_to_cell)]
                sum_quadrant_deflection += store_cell_spks_deflection[int(q_val/degree_to_cell)+(q_val%degree_to_cell)]
                sum_quadrant_baseline   += store_cell_spks_baseline[int(q_val/degree_to_cell)+(q_val%degree_to_cell)]
                # sum_quadrant+=store_cell_spks[int(q_val/degree_to_cell)]
            sum_quadrant_values.append(sum_quadrant)
            sum_quadrant_values_deflection.append(sum_quadrant_deflection)
            sum_quadrant_values_baseline.append(sum_quadrant_baseline)
        ##############################

        # --- Plotting figures
        plt.figure( figsize=(40,30))
        plt.subplot(3,1,1)
        stim_x  = spikes_dict['time_vector']
        stim_y_ = [0 for i in stim_x]
        for ind, k in enumerate(deflection_dict['deflection_times']):
            t_ON  = deflection_dict['deflection_times'][ind][0]
            t_OFF = deflection_dict['deflection_times'][ind][1]
            for ind,y in enumerate(stim_y_):
                if ind >= int(t_ON/deflection_dict['dt']) and ind < int(t_OFF/deflection_dict['dt']): stim_y_[ind]+=1

        plt.plot(stim_x,stim_y_,color='k')
        plt.ylim(-0.1,1.1)
        plt.yticks([0,1],['OFF','ON'])
        plt.ylabel('Whisker deflection')
        plt.xlim([1850,10000])

        # --- Remove corners
        ax=plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.xticks([])  # --- Remove x-axis ticks
        
        # --- Draws a line in the raster where the firing is maximum (b) and minimum (r)
        plt.subplot(3,1,2)
        # plt.plot([0,int(len(raster_sampled)*deflection_dict['dt'])],[max_tuning,max_tuning],':b', linewidth = 2,alpha=0.5)
        # plt.plot([0,int(len(raster_sampled)*deflection_dict['dt'])],[min_tuning,min_tuning],':r', linewidth = 2,alpha=0.5)
        plt.yticks([0,45,90,135,180], ['0','90','180','270','360'])
        plt.ylabel('Brainstem cell\nAngle (deg)')
        plt.xlim([1850,10000])
        
        # --- Remove corners
        ax=plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.xticks([])  # --- Remove x-axis ticks

        # for prob_ind,prob in enumerate(data_y_sampled):
        for raster_ind,raster in enumerate(raster_sampled):

            random_var_=[]
            for i_r_v,r_v in enumerate(raster_sampled[raster_ind]):
                if r_v==1: random_var_.append(i_r_v)

            # --- Sampled Raster plot
            plt.subplot(3,1,2)
            for val in random_var_: plt.plot(spikes_dict['time_vector'][raster_ind],val,'.k')

        if deflection_dict['paper_sampling']!=deflection_dict['dt']:
            print('paper_sampling: ', deflection_dict['paper_sampling'], '\t', deflection_dict['paper_sampling']/deflection_dict['dt'])
            separateProbPlots=False
            if separateProbPlots:
                # --- Plotting histogram with netpyne sampling rate
                for sample_ind in np.arange(0,len(raster_sampled),1):
                    sampled_raster_chunk = raster_sampled[sample_ind:sample_ind+1]
                    prob_chunk = sum([sum(slice) for slice in sampled_raster_chunk])
                    plt.subplot(6,1,5)
                    plt.bar(spikes_dict['time_vector'][sample_ind],prob_chunk/deflection_dict['n_cells'],color='k')
                    plt.ylabel('Spike probability')
                    plt.xlabel('Time (ms)')

                # --- Plotting histogram with paper sampling rate
                for sample_ind in np.arange(0,len(raster_sampled),int(deflection_dict['paper_sampling']/deflection_dict['dt'])):
                    sampled_raster_chunk = raster_sampled[sample_ind:sample_ind+int(deflection_dict['paper_sampling']/deflection_dict['dt'])]
                    prob_chunk = sum([sum(slice) for slice in sampled_raster_chunk])
                    plt.subplot(6,1,6)
                    plt.bar(spikes_dict['time_vector'][sample_ind],prob_chunk/deflection_dict['n_cells'],color='k')
                    plt.ylabel('Spike probability')
                    plt.xlabel('Time (ms)')
            else:
                # --- Plotting histogram with paper sampling rate
                for sample_ind in np.arange(0,len(raster_sampled),int(deflection_dict['paper_sampling']/deflection_dict['dt'])):
                    sampled_raster_chunk = raster_sampled[sample_ind:sample_ind+int(deflection_dict['paper_sampling']/deflection_dict['dt'])]
                    prob_chunk = sum([sum(slice) for slice in sampled_raster_chunk])
                    plt.subplot(3,1,3)
                    plt.bar(spikes_dict['time_vector'][sample_ind],prob_chunk/deflection_dict['n_cells'],color='k')
                    plt.ylabel('Spike probability')
                    plt.xlabel('Time (ms)')
            
        else:
            # --- Sampled Probability histogram
            for sample_ind in np.arange(0,len(raster_sampled),int(deflection_dict['paper_sampling']/deflection_dict['dt'])):
                sampled_raster_chunk = raster_sampled[sample_ind:sample_ind+int(deflection_dict['paper_sampling']/deflection_dict['dt'])]
                prob_chunk = sum([sum(slice) for slice in sampled_raster_chunk])
                plt.subplot(3,1,3)
                plt.bar(spikes_dict['time_vector'][sample_ind],prob_chunk/deflection_dict['n_cells'],color='k')
                plt.ylabel('Spike probability')
                plt.yticks([0,0.15,0.3])
                plt.xlabel('Time (ms)')

                # --- Remove corners
                ax=plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
        
        plt.xlim([1850,10000])

        plt.savefig(save_path+'sampled_raster/9_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us_'+str(deflection_dict['target_angle'])+deflection_dict['stims_string']+'_degrees_sampled_scatter.png',dpi=300)

        #######################################################################################################################################
        # --- Step polar plot
        cell_degs = [int(cell_ind*(360/deflection_dict['n_cells'])) for cell_ind in range(len(store_cell_spks))]
        cell_rads = [np.radians(int(cell_ind*(360/deflection_dict['n_cells']))) for cell_ind in range(len(store_cell_spks))]

        norm_store_cell_spks            = [cell_spk/max(store_cell_spks)               for cell_spk in store_cell_spks]
        norm_store_cell_spks_deflection = [cell_spk/max(store_cell_spks_deflection)    for cell_spk in store_cell_spks_deflection]
        norm_store_cell_spks_baseline   = [cell_spk/max(store_cell_spks_baseline)      for cell_spk in store_cell_spks_baseline]
        
        plt.figure( figsize=(40,30))
        from scipy.interpolate import interp1d
        plt.subplot(2,2,3,projection='polar')
        # quad_colors = ['cyan','magenta','forestgreen','orange','blue','red','purple','grey','cyan']
        quad_colors = ['blue','darkviolet','mediumvioletred','red','gold','darkkhaki','forestgreen','dodgerblue']
        # quad_colors = ['cyan','magenta','forestgreen','orange']
        for ind_quad_vals,quad_vals in enumerate(quadrant_values):
            lw=5
            if abs(quad_vals[0]-quad_vals[-1])>180: 
                quad_curve_coordinates = [[0,quad_vals[-1]], [1.1,1.1]]
                quad_curve_coordinates[0] = np.deg2rad(quad_curve_coordinates[0])
                x = np.linspace( quad_curve_coordinates[0][0], quad_curve_coordinates[0][1], 500)
                y = interp1d( quad_curve_coordinates[0], quad_curve_coordinates[1])( x)
                plt.plot(x, y, color=quad_colors[ind_quad_vals],alpha=0.5,linewidth=lw)

                quad_curve_coordinates = [[quad_vals[0],360], [1.1,1.1]]

                quad_curve_coordinates[0] = np.deg2rad(quad_curve_coordinates[0])
                x = np.linspace( quad_curve_coordinates[0][0], quad_curve_coordinates[0][1], 500)
                y = interp1d( quad_curve_coordinates[0], quad_curve_coordinates[1])( x)
                plt.plot(x, y, color=quad_colors[ind_quad_vals],alpha=0.5,linewidth=lw)

            else:
                quad_curve_coordinates = [[quad_vals[0],quad_vals[-1]], [1.1,1.1]]            
                quad_curve_coordinates[0] = np.deg2rad(quad_curve_coordinates[0])
                x = np.linspace( quad_curve_coordinates[0][0], quad_curve_coordinates[0][1], 500)
                y = interp1d( quad_curve_coordinates[0], quad_curve_coordinates[1])( x)
                plt.plot(x, y, color=quad_colors[ind_quad_vals],alpha=0.5,linewidth=lw)

            if deflection_dict['addMouseHead']: 
                import matplotlib.image as image 
                # file = save_path+'ims/mouse_image.png'
                # file_coords=[450,700]
                file = save_path+'ims/mouse_image_crop.png'
                file_coords=[2800,2300]
                im = image.imread(file) 
                plt.figimage(im, file_coords[0], file_coords[1], zorder = 0, alpha =.1) 

        for cell_ind,cell_spk in enumerate(store_cell_spks):
            if cell_ind==max_tuning:
                edgecolor='b'
                color='b'

            elif cell_ind==min_tuning:
                edgecolor='r'
                color='r'
            else: continue

            cell_deg = int(cell_ind*(360/deflection_dict['n_cells']))

            # --- Regular contour plot
            plt.subplot(2,1,1)
            plt.bar(cell_deg,cell_spk,color=color,edgecolor=edgecolor,width=0.1)
            plt.xlabel('Angle (deg)')
            plt.ylabel('Number of spikes')
            plt.xticks([0,45,90,135,180,225,270,315,360])
            plt.yticks([0,150,300])

            # --- Remove corners
            ax=plt.gca()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # --- Normalized Polar step plot
            plt.subplot(2,2,3,projection='polar')
            plt.bar(np.radians(cell_deg), cell_spk/max(store_cell_spks), width=0.025, bottom=0.0, color=color, edgecolor=edgecolor, alpha=1)
            ax=plt.gca()
            ax.grid(False)

        # --- Step plots
        plt.subplot(2,1,1)
        plt.step(cell_degs, store_cell_spks,            where='mid',color='k',alpha=1.0,linewidth=5)
        plt.step(cell_degs, store_cell_spks_deflection, where='mid',color='saddlebrown',alpha=0.5,linewidth=5)
        plt.step(cell_degs, store_cell_spks_baseline,   where='mid',color='k',alpha=0.5)

        plt.subplot(2,2,3,projection='polar')
        # plt.step(cell_rads, norm_store_cell_spks,               where='mid',color='k',alpha=0.5)
        plt.step(cell_rads, norm_store_cell_spks_deflection,    where='mid',color='saddlebrown',alpha=0.5,linewidth=4)
        # plt.step(cell_rads, norm_store_cell_spks_baseline,      where='mid',color='k',alpha=0.5)
        
        # --- Normalized cumulative bar plot
        norm_sum_quadrant_values            = [q_val/max(sum_quadrant_values) for q_val in sum_quadrant_values]
        norm_sum_quadrant_values_deflection = [q_val/max(sum_quadrant_values_deflection) for q_val in sum_quadrant_values_deflection]
        norm_sum_quadrant_values_baseline   = [q_val/max(sum_quadrant_values_baseline) for q_val in sum_quadrant_values_baseline]
        plt.subplot(2,2,4)
        # plt.bar(quadrant,norm_sum_quadrant_values,              width=20,color=quad_colors,edgecolor='k')
        plt.bar(quadrant,norm_sum_quadrant_values_deflection,   width=20,color=quad_colors,edgecolor='k')
        # plt.bar(quadrant,norm_sum_quadrant_values_baseline,     width=20,color=quad_colors,edgecolor='k')
        # plt.bar(quadrant,norm_sum_quadrant_values,width=20,color='orange',edgecolor='k')
        
        ## --- Plotting reference angular tuning values from Hartings paper
        # plt.plot(list(np.arange(0,360,45)),[0.19764487038033318, 0.20137985, 0.32471083, 0.64055063, 1.0, 0.64055063, 0.32471083, 0.20137985],'ok')
        
        plt.xticks(quadrant)
        plt.ylim([0,1.1])
        # plt.grid()
        plt.xlabel('Quadrant')
        plt.ylabel('Cumulative Number of spikes')
        # --- Remove corners
        ax=plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        plt.savefig(save_path+'sampled_raster/11_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us_'+str(deflection_dict['target_angle'])+deflection_dict['stims_string']+'_contourPlot.png',dpi=300)

        if deflection_dict['target_angle'] is not None:
            # --- paperFigure - preparing data and colors for plotting angular tuning values
            quadrant_extended = quadrant+[360]
            quadrant_extended_rad = [np.deg2rad(val) for val in quadrant_extended]
            norm_sum_quadrant_values_deflection_extended = norm_sum_quadrant_values_deflection+[norm_sum_quadrant_values_deflection[0]]

            quadrant_extended_colors = []
            for angle in quadrant_extended:
                if   (angle==0) or (angle==360):    quadrant_extended_colors.append('blue')
                elif (angle==45):                   quadrant_extended_colors.append('darkviolet')
                elif (angle==90):                   quadrant_extended_colors.append('mediumvioletred')
                elif (angle==135):                  quadrant_extended_colors.append('red')
                elif (angle==180):                  quadrant_extended_colors.append('gold')
                elif (angle==225):                  quadrant_extended_colors.append('darkkhaki')
                elif (angle==270):                  quadrant_extended_colors.append('forestgreen')
                elif (angle==315):                  quadrant_extended_colors.append('dodgerblue')
                else:                               quadrant_extended_colors.append('grey')

            angular_tuning_shift_angles = deflection_dict['target_angle']-180
            angular_tuning_plot_angles = list([val+angular_tuning_shift_angles for val in np.arange(0,361,45)])
            angular_tuning_plot_angles_positive = []
            for val in angular_tuning_plot_angles:
                if   val<0:     new_val = val+360
                elif val>360:   new_val = val-360
                # elif val==360:  new_val = 0
                else:           new_val = val
                angular_tuning_plot_angles_positive.append(new_val)
            angular_tuning_plot_angles_positive_rad = [np.deg2rad(val) for val in angular_tuning_plot_angles_positive]

            # --- paperFigure - individual angular tuning plots for each deflection dataset

            hartings_mle_angular_tuning = [0.19764487038033318, 0.20137985, 0.32471083, 0.64055063, 1.0, 0.64055063, 0.32471083, 0.20137985, 0.19764487038033318]

            plt.figure(figsize=(25,20))
            plt.subplot(1,1,1,projection='polar')
            plt.rcParams.update({'font.size': 60})  
            plt.plot(quadrant_extended_rad,norm_sum_quadrant_values_deflection_extended, linewidth=5,color='k')        
            for i in range(len(quadrant_extended_rad)): plt.plot(quadrant_extended_rad[i],norm_sum_quadrant_values_deflection_extended[i], 'o', color=quadrant_extended_colors[i], alpha=0.8,markersize=25)
            plt.plot(angular_tuning_plot_angles_positive_rad,hartings_mle_angular_tuning, '-k', alpha=0.5,linewidth=10)
            plt.plot(angular_tuning_plot_angles_positive_rad,hartings_mle_angular_tuning, '^k', alpha=0.8,markersize=25)
            plt.xlim([0,np.deg2rad(360)])
            plt.yticks([0.5, 1.0])
            plt.ylim([0,1.5])
            plt.xticks(np.deg2rad([0,90,180,270]))
            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine
            plt.savefig(save_path+'sampled_raster/12_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us_paperFigure_'+str(deflection_dict['target_angle'])+deflection_dict['stims_string']+'_contourPlot.png',dpi=300)

            # --- paperFigure - angular tuning values only
            plt.figure(figsize=(25,20))
            plt.subplot(1,1,1,projection='polar')
            plt.rcParams.update({'font.size': 60})  
            plt.plot(quadrant_extended_rad,norm_sum_quadrant_values_deflection_extended, linewidth=5,color='k')        
            for i in range(len(quadrant_extended_rad)): plt.plot(quadrant_extended_rad[i],norm_sum_quadrant_values_deflection_extended[i], 'o', color=quadrant_extended_colors[i], alpha=0.8,markersize=25)
            plt.xlim([0,np.deg2rad(360)])
            plt.yticks([])
            plt.ylim([0,1.5])
            plt.xticks([])
            ax1 = plt.gca()
            ax1.spines['polar'].set_visible(False)  # Hide the circular spine
            plt.grid(False)
            plt.savefig(save_path+'sampled_raster/13_samplingBins|'+str(int(deflection_dict['dt']*1000))+'us_paperFigure_'+str(deflection_dict['target_angle'])+deflection_dict['stims_string']+'_contourPlot.png',dpi=300)

            # # stats
            # from scipy.stats import ttest_1samp

            # for at_ind, at_val in enumerate(hartings_mle_angular_tuning):
            #     t_stat, p_value = ttest_1samp(norm_sum_quadrant_values_deflection_extended[at_ind], at_val)
            #     if      p_value>0.05:                           significance = 'ns'
            #     elif    (p_value<=0.05)  and (p_value>0.005):   significance = '*'
            #     elif    (p_value<=0.005) and (p_value>0.0005):  significance = '**'
            #     elif    (p_value<=0.0005):                      significance = '***'
            #     else:                                           significance = '--'
            #     print('\t\t\t',quadrant_extended[at_ind], '\t',p_value, '\t',significance, '\t (',t_stat,', ', p_value,')')


            # # plt.bar(quadrant,norm_sum_quadrant_values,              width=20,color=quad_colors,edgecolor='k')
            # # plt.bar(quadrant,norm_sum_quadrant_values_baseline,     width=20,color=quad_colors,edgecolor='k')
            # # plt.bar(quadrant,norm_sum_quadrant_values,width=20,color='orange',edgecolor='k')
            


            # plt.plot(list(np.arange(0,360,45)),hartings_mle_angular_tuning, alpha=0.5,linewidth=10)
            # plt.plot(list(np.arange(0,360,45)),hartings_mle_angular_tuning, alpha=0.8,markersize=25)
            
            # # plt.plot(list(np.arange(0,360,45)),[0.19764487038033318, 0.20137985, 0.32471083, 0.64055063, 1.0, 0.64055063, 0.32471083, 0.20137985], alpha=0.5,linewidth=10)
            # # plt.plot(list(np.arange(0,360,45)),[0.19764487038033318, 0.20137985, 0.32471083, 0.64055063, 1.0, 0.64055063, 0.32471083, 0.20137985], alpha=0.8,markersize=25)
            
            # plt.xticks(quadrant)
            # plt.ylim([0,1.1])
            # # plt.grid()
            # plt.xlabel('Quadrant')
            # plt.ylabel('Cumulative Number of spikes')
            # # --- Remove corners
            # ax=plt.gca()
            # ax.spines['top'].set_visible(False)
            # ax.spines['right'].set_visible(False)
            # ax.spines['bottom'].set_visible(False)
            # ax.spines['left'].set_visible(False)


########################################################################################################################

class PlotFigs():
    def plotAllTraces(hist_dict, save_path =''):
        plt.figure(figsize=(25,20))
        for key in hist_dict.keys():
            plt.plot(hist_dict[key]['x'],hist_dict[key]['y'],color='k')
        # stim markers
        stim_1=38;stim_2=122
        plt.plot(hist_dict['0stim']['x'][stim_1],hist_dict['0stim']['y'][stim_1],color='r',marker='o')
        plt.plot(hist_dict['0stim']['x'][stim_2],hist_dict['0stim']['y'][stim_2],color='r',marker='o')

        # stim markers
        stim_1=62;stim_2=202
        plt.plot(hist_dict['PrV']['x'][stim_1],hist_dict['PrV']['y'][stim_1],color='cyan',marker='o')
        plt.plot(hist_dict['PrV']['x'][stim_2],hist_dict['PrV']['y'][stim_2],color='cyan',marker='o')

        # stim markers
        stim_1=62;stim_2=190
        plt.plot(hist_dict['VPM']['x'][stim_1],hist_dict['VPM']['y'][stim_1],color='purple',marker='o')
        plt.plot(hist_dict['VPM']['x'][stim_2],hist_dict['VPM']['y'][stim_2],color='purple',marker='o')

        plt.ylabel('spike probability - (0.1/ms)')
        plt.xlabel('time - (ms)')
        plt.ylim([-0.05, 0.4])
        plt.savefig(save_path+'5_allTraces.png',dpi=300)

        for key in hist_dict.keys():
            plt.figure(figsize=(25,20))
            plt.rcParams.update({'font.size': 60})
            if key == '0stim':  plt.plot(hist_dict[key]['x'],hist_dict[key]['y'],color='k')
            else:               plt.fill_between(hist_dict[key]['x'],hist_dict[key]['y'],color='k',step="pre",)
            plt.xlabel('Time (ms)')
            plt.ylabel('Spike probability (0.1/ms)')
            plt.ylim([-0.01, 0.4])
            yticks=[0,0.1,0.2,0.3,0.4]
            plt.yticks(yticks,labels=[10*ytick for ytick in yticks])
            
            plt.xlim([-10,410])
            ON_offset = 106
            OFF_offset = 210
            plt.xticks([ON_offset-100, ON_offset, ON_offset+OFF_offset, ON_offset+OFF_offset+100],labels=[i-ON_offset for i in [ON_offset-100, ON_offset, ON_offset+OFF_offset, ON_offset+OFF_offset+100]])

            ax1 = plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)

            plt.savefig(save_path+'paperFig_5_allTraces_Bar_'+key+'.png',dpi=300)


    def plotSampledDataBar(sampledData, save_path = ''):
        plt.figure( figsize=(25,20))
        plt.bar(sampledData['baseline']['x'], sampledData['baseline']['y'], color='k')
        plt.bar(sampledData['ON']['x'],       sampledData['ON']['y'],       color='green')
        plt.bar(sampledData['OFF']['x'],      sampledData['OFF']['y'],      color='goldenrod')
        plt.savefig(save_path+'8_polyfit_ON_OFF_sampled.png',dpi=300)

    # --- Getting samples for plotting
    def getFitX(xVals,nSamples): return np.linspace(min(xVals), max(xVals), nSamples) # --- sampling over a linespace for plotting below

    def plotFits(PrV_data,yFit_ON,yFit_OFF,yFit_ON_zero,yFit_OFF_zero,nSamples, save_path = ''):
        # --- Getting X samples for plotting
        xFit_ON         = PlotFigs.getFitX(PrV_data['ON']['x'],         nSamples)
        xFit_OFF        = PlotFigs.getFitX(PrV_data['OFF']['x'],        nSamples)
        xFit_ON_zero    = PlotFigs.getFitX(PrV_data['ON_zero']['x'],    nSamples)
        xFit_OFF_zero   = PlotFigs.getFitX(PrV_data['OFF_zero']['x'],   nSamples)

        # --- Plotting data with fit on top
        plt.figure( figsize=(25,20))
        plt.subplot(2,1,1)
        plt.plot(   PrV_data['all']['x'], PrV_data['all']['y'], '-',    linewidth=5,   color= 'k',      alpha=0.1)
        plt.plot(   xFit_ON,              yFit_ON(xFit_ON),     ':',    linewidth=5,   color= 'green',)
        plt.plot(   xFit_OFF,             yFit_OFF(xFit_OFF),   ':',    linewidth=5,   color= 'goldenrod',)
        plt.xlim(-20,440)
        plt.ylim(-0.05,0.4)
        plt.subplot(2,1,2)
        plt.plot(   PrV_data['all']['x'], PrV_data['all']['y'], '-',    linewidth=5,   color= 'k',      alpha=0.1)
        plt.plot(   xFit_ON_zero,      yFit_ON_zero(xFit_ON_zero),   ':',    linewidth=5,   color= 'green',)
        plt.plot(   xFit_OFF_zero,     yFit_OFF_zero(xFit_OFF_zero), ':',    linewidth=5,   color= 'goldenrod',)
        plt.xlim(-20,440)
        plt.ylim(-0.05,0.4)
        plt.savefig(save_path+'7_polyfit_ON_OFF.png',dpi=300)

    def plotFits_zero(PrV_data,yFit_ON_zero,yFit_OFF_zero,nSamples, save_path = ''):
        # --- Getting X samples for plotting
        xFit_ON         = PlotFigs.getFitX(PrV_data['ON']['x'],         nSamples)
        xFit_OFF        = PlotFigs.getFitX(PrV_data['OFF']['x'],        nSamples)
        xFit_ON_zero    = PlotFigs.getFitX(PrV_data['ON_zero']['x'],    nSamples)
        xFit_OFF_zero   = PlotFigs.getFitX(PrV_data['OFF_zero']['x'],   nSamples)

        # --- Plotting data with fit on top
        plt.figure( figsize=(25,20))
        # plt.subplot(2,1,2)
        plt.plot(   PrV_data['all']['x'], PrV_data['all']['y'], '-',    linewidth=5,   color= 'k',      alpha=0.1)
        plt.plot(   xFit_ON,      yFit_ON_zero(xFit_ON_zero),   ':',    linewidth=5,   color= 'green',)
        plt.plot(   xFit_OFF,     yFit_OFF_zero(xFit_OFF_zero), ':',    linewidth=5,   color= 'goldenrod',)
        plt.xlim(-20,440)
        plt.ylim(-0.05,0.4)
        plt.savefig(save_path+'7_polyfit_ON_OFF_zero.png',dpi=300)

        # --- Plotting data with fit on top
        plt.figure( figsize=(25,20))
        # plt.subplot(2,1,2)
        plt.plot(   PrV_data['all']['x'], [y/0.1 for y in PrV_data['all']['y']], '-',    linewidth=5,   color= 'k',      alpha=0.1)
        # plt.plot(   xFit_ON,      yFit_ON_zero(xFit_ON_zero),   ':',    linewidth=5,   color= 'green',)
        # plt.plot(   xFit_OFF,     yFit_OFF_zero(xFit_OFF_zero), ':',    linewidth=5,   color= 'goldenrod',)
        plt.xlim(-20,440)
        # plt.ylim(-0.05,0.4)
        plt.savefig(save_path+'71_polyfit_ON_OFF_zero_paperScale.png',dpi=300)

    def plotFitAndSampled(PrV_data, yFit_ON, yFit_OFF, nSamples, sampledData, save_path = ''):
        # --- Getting X samples for plotting
        xFit_ON         = PlotFigs.getFitX(PrV_data['ON']['x'],         nSamples)
        xFit_OFF        = PlotFigs.getFitX(PrV_data['OFF']['x'],        nSamples)
        # --- Plotting data with fit on top and Sampled Data
        plt.figure( figsize=(25,20))
        plt.subplot(2,1,1)
        plt.plot(   PrV_data['all']['x'],   PrV_data['all']['y'],   '-',    linewidth=5,   color= 'lightgrey',)
        plt.plot(   xFit_ON,                yFit_ON(xFit_ON),       ':',    linewidth=5,   color= 'green',)
        plt.plot(   xFit_OFF,               yFit_OFF(xFit_OFF),     ':',    linewidth=5,   color= 'goldenrod',)
        plt.ylim([-0.05, 0.4])
        plt.ylabel('spike probability - (0.1/ms)')
        plt.subplot(2,1,2)
        plt.bar(sampledData['baseline']['x'],   sampledData['baseline']['y'],   color='k')
        plt.bar(sampledData['ON']['x'],         sampledData['ON']['y'],         color='green')
        plt.bar(sampledData['OFF']['x'],        sampledData['OFF']['y'],        color='goldenrod')
        plt.bar(sampledData['ON_zero']['x'],    sampledData['ON_zero']['y'],    color='limegreen',  alpha=0.25)
        plt.bar(sampledData['OFF_zero']['x'],   sampledData['OFF_zero']['y'],   color='gold',       alpha=0.25)
        plt.ylim([-0.05, 0.4])
        plt.xlabel('Time (ms)')
        plt.ylabel('spike probability - (1/ms)')
        plt.savefig(save_path+'81_polyfit_ON_OFF_sampled.png',dpi=300)

        # --- Plotting data with fit on top and Sampled Data
        plt.figure( figsize=(25,20))
        plt.subplot(2,1,1)
        plt.plot(   PrV_data['all']['x'],   [y/0.1 for y in PrV_data['all']['y']],   '-',    linewidth=5,   color= 'lightgrey',)
        # plt.plot(   xFit_ON,                yFit_ON(xFit_ON),       ':',    linewidth=5,   color= 'green',)
        # plt.plot(   xFit_OFF,               yFit_OFF(xFit_OFF),     ':',    linewidth=5,   color= 'goldenrod',)
        plt.ylim([-0.05, 0.4])
        plt.ylabel('spike probability - (0.1/ms)')
        plt.subplot(2,1,2)
        plt.bar(sampledData['baseline']['x'],   [y/0.1 for y in sampledData['baseline']['y']],  color='k')
        plt.bar(sampledData['ON']['x'],         [y/0.1 for y in sampledData['ON']['y']],        color='green')
        plt.bar(sampledData['OFF']['x'],        [y/0.1 for y in sampledData['OFF']['y']] ,      color='goldenrod')
        plt.bar(sampledData['ON_zero']['x'],    [y/0.1 for y in sampledData['ON_zero']['y']],   color='limegreen',  alpha=0.25)
        plt.bar(sampledData['OFF_zero']['x'],   [y/0.1 for y in sampledData['OFF_zero']['y']],  color='gold',       alpha=0.25)
        plt.ylim([-0.05, 0.4])
        plt.xlabel('Time (ms)')
        plt.ylabel('spike probability - (0.1/ms)')
        plt.savefig(save_path+'812_polyfit_ON_OFF_sampled_paperScale.png',dpi=300)

    def plotFitAndSampled_zero(PrV_data, yFit_ON_zero, yFit_OFF_zero, nSamples, sampledData, save_path = '',zero_times=[0,0],baseline_data=None):
        # --- Getting X samples for plotting
        xFit_ON         = PlotFigs.getFitX(PrV_data['ON']['x'],         nSamples)
        xFit_OFF        = PlotFigs.getFitX(PrV_data['OFF']['x'],        nSamples)
        xFit_ON_zero    = PlotFigs.getFitX(PrV_data['ON_zero']['x'],    nSamples)
        xFit_OFF_zero   = PlotFigs.getFitX(PrV_data['OFF_zero']['x'],   nSamples)

        x_data_ON  = [x+zero_times[0] for x in sampledData['ON_zero']['x']]
        x_data_OFF = [x+zero_times[1] for x in sampledData['OFF_zero']['x']]

        # --- Plotting data with fit on top and Sampled Data
        plt.figure( figsize=(25,20))
        plt.subplot(2,1,1)
        plt.plot(   PrV_data['all']['x'],   PrV_data['all']['y'],   '-',    linewidth=5,   color= 'lightgrey',)
        plt.plot(   xFit_ON,                yFit_ON_zero(xFit_ON_zero),       ':',    linewidth=5,   color= 'green',)
        plt.plot(   xFit_OFF,               yFit_OFF_zero(xFit_OFF_zero),     ':',    linewidth=5,   color= 'goldenrod',)
        plt.ylim([-0.05, 0.4])
        plt.ylabel('spike probability - (0.1/ms)')
        plt.subplot(2,1,2)
        plt.bar(sampledData['baseline']['x'],   sampledData['baseline']['y'],   color='k')
        plt.bar(x_data_ON,                      sampledData['ON_zero']['y'],    color='limegreen',  alpha=0.25)
        plt.bar(x_data_OFF,                     sampledData['OFF_zero']['y'],   color='gold',       alpha=0.25)
        # plt.bar(sampledData['ON']['x'],         sampledData['ON']['y'],         color='green')
        # plt.bar(sampledData['OFF']['x'],        sampledData['OFF']['y'],        color='goldenrod')
        # plt.bar(sampledData['ON_zero']['x'],    sampledData['ON_zero']['y'],    color='limegreen',  alpha=0.25)
        # plt.bar(sampledData['OFF_zero']['x'],   sampledData['OFF_zero']['y'],   color='gold',       alpha=0.25)
        plt.ylim([-0.05, 0.4])
        plt.xlabel('Time (ms)')
        plt.ylabel('spike probability - (1/ms)')
        plt.savefig(save_path+'82_polyfit_ON_OFF_sampled.png',dpi=300)
        
        # --- Plotting data with fit on top and Sampled Data
        plt.figure( figsize=(25,20))
        plt.subplot(2,1,1)
        plt.plot(   PrV_data['all']['x'],   [y/0.1 for y in PrV_data['all']['y']],   '-',    linewidth=5,   color= 'lightgrey',)
        # plt.plot(   xFit_ON,                yFit_ON_zero(xFit_ON_zero),       ':',    linewidth=5,   color= 'green',)
        # plt.plot(   xFit_OFF,               yFit_OFF_zero(xFit_OFF_zero),     ':',    linewidth=5,   color= 'goldenrod',)
        plt.ylim([-0.05, 4])
        plt.ylabel('spike probability - (0.1/ms)')
        plt.subplot(2,1,2)
        plt.bar(sampledData['baseline']['x'],   [y/0.1 for y in sampledData['baseline']['y']],   color='k')
        plt.bar(x_data_ON,                      [y/0.1 for y in sampledData['ON_zero']['y']],    color='limegreen',  alpha=0.25)
        plt.bar(x_data_OFF,                     [y/0.1 for y in sampledData['OFF_zero']['y']],   color='gold',       alpha=0.25)
        # plt.bar(sampledData['ON']['x'],         sampledData['ON']['y'],         color='green')
        # plt.bar(sampledData['OFF']['x'],        sampledData['OFF']['y'],        color='goldenrod')
        # plt.bar(sampledData['ON_zero']['x'],    sampledData['ON_zero']['y'],    color='limegreen',  alpha=0.25)
        # plt.bar(sampledData['OFF_zero']['x'],   sampledData['OFF_zero']['y'],   color='gold',       alpha=0.25)
        plt.ylim([-0.05, 4])
        plt.xlabel('Time (ms)')
        plt.ylabel('spike probability - (0.1/ms)')
        plt.savefig(save_path+'822_polyfit_ON_OFF_sampled.png',dpi=300)



        if baseline_data is not None: 
            baseline_mean = baseline_data[0]
            baseline_std  = baseline_data[1]
            baseline_y_mean     = [baseline_mean for i in range(len(sampledData['baseline']['x']))]
            baseline_y_std_pos  = [baseline_mean+baseline_std for i in range(len(sampledData['baseline']['x']))]
            baseline_y_std_neg  = [baseline_mean-baseline_std for i in range(len(sampledData['baseline']['x']))]





        plt.figure(figsize=(25,20))
        plt.rcParams.update({'font.size': 60})        
        plt.fill_between(PrV_data['all']['x'],   PrV_data['all']['y'],color='k',step="pre", alpha=0.5)
        plt.plot(   xFit_ON,                yFit_ON_zero(xFit_ON_zero),       '-',    linewidth=30,   color= 'darkblue',    alpha = 0.25)
        plt.plot(   xFit_OFF,               yFit_OFF_zero(xFit_OFF_zero),     '-',    linewidth=30,   color= 'darkred',     alpha = 0.25)
        
        plt.plot(sampledData['baseline']['x'], baseline_y_mean, color='darkorange', linewidth=30, alpha = 0.25)
        # y_shift = 0.010
        # plt.errorbar(53, baseline_mean+y_shift, baseline_std, linestyle='None', marker='.', markersize=2, capsize=6, color='darkorange')
        
        plt.xlabel('Time (ms)')
        plt.ylabel('Spike probability (0.1/ms)')
        plt.ylim([-0.01, 0.4])
        yticks=[0,0.1,0.2,0.3,0.4]
        plt.yticks(yticks,labels=[10*ytick for ytick in yticks])

        plt.xlim([-10,410])
        ON_offset = 106
        OFF_offset = 210
        plt.xticks([ON_offset-100, ON_offset, ON_offset+OFF_offset, ON_offset+OFF_offset+100],labels=[i-ON_offset for i in [ON_offset-100, ON_offset, ON_offset+OFF_offset, ON_offset+OFF_offset+100]])

        ax1 = plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        plt.savefig(save_path+'paperFig_82_polyfit_ON_OFF_sampled.png',dpi=300)



        # if key == '0stim':  plt.plot(hist_dict[key]['x'],hist_dict[key]['y'],color='k')
        # else:               plt.fill_between(hist_dict[key]['x'],hist_dict[key]['y'],color='k',step="pre",)
        # plt.xlabel('Time (ms)')
        # plt.ylabel('Spike probability (0.1/ms)')
        # plt.ylim([-0.01, 0.4])
        # yticks=[0,0.1,0.2,0.3,0.4]
        # plt.yticks(yticks,labels=[10*ytick for ytick in yticks])
        
        # plt.xlim([-10,410])
        # ON_offset = 106
        # OFF_offset = 210
        # plt.xticks([ON_offset-100, ON_offset, ON_offset+OFF_offset],labels=[i-ON_offset for i in [ON_offset-100, ON_offset, ON_offset+OFF_offset]])
        
        # ax1 = plt.gca()
        # ax1.spines['top'].set_visible(False)
        # ax1.spines['right'].set_visible(False)
        # ax1.spines['bottom'].set_visible(False)
        # ax1.spines['left'].set_visible(False)

        # plt.savefig(save_path+'5_allTraces_Bar_'+key+'.png',dpi=300)















    def plotAngularTuning(angular_tuning,nSamples,yFit_angTuning, save_path = '', area_params=None):
        xFit_angularTuning = PlotFigs.getFitX(angular_tuning['x'],nSamples)
        # np.linspace(min(angular_tuning['x']), max(angular_tuning['x']), 250)
        
        plt.figure( figsize=(16,16))
        plt.plot(       angular_tuning['x'],   angular_tuning['y'],                 '-',    linewidth=2,   color= 'b',)
        plt.plot(       xFit_angularTuning,    yFit_angTuning(xFit_angularTuning),  ':',    linewidth=5,   color= 'grey',)
        plt.ylim(0,1.1)
        plt.legend(['PrV','polynomial fit'])
        plt.savefig(save_path+'_1_probabilityScaling_polyfit_rescaled.png',dpi=300)
    



        plt.figure(figsize=(25,20))
        plt.rcParams.update({'font.size': 60})    

        # Without angular tuning - used to calculate the are
        plt.plot(       angular_tuning['x'],   [1 for xval in range(len(angular_tuning['x']))],          '-',    linewidth=10,   color= 'grey', alpha=0.5)
        plt.plot(       angular_tuning['x'],   [1 for xval in range(len(angular_tuning['x']))],          'o',    markersize=25,  color= 'grey', alpha=0.9)
        
        # With angular tuning - which has a reduced area
        plt.plot(       angular_tuning['x'],   angular_tuning['y'],                 '-',    linewidth=10,   color= 'k', alpha=0.5)
        plt.plot(       angular_tuning['x'],   angular_tuning['y'],                 '^',    markersize=25,  color= 'k', alpha=0.9)
        
        # With angular tuning + scaling - to compensate for the reduced area
        # plt.plot(       angular_tuning['x'],   angular_tuning['y_scaled'],          '-',    linewidth=10,   color= 'lightseagreen', alpha=0.5)
        plt.plot(       xFit_angularTuning,    yFit_angTuning(xFit_angularTuning),  ':',    linewidth=10,   color= 'lightseagreen', alpha=0.5)
        plt.plot(       angular_tuning['x'],   angular_tuning['y_scaled'],          'p',    markersize=25,  color= 'lightseagreen', alpha=0.9)
        
        plt.ylim(-0.1,3.1)
        plt.yticks([0,1,2,3])
        plt.ylabel('Normalized angular tuning')
        plt.xticks([0,1,2,3,4], ['0','-π/4','-π/2','-3π/4','-π'])
        plt.xlabel('Distance from target angle (rad)')
        # plt.legend(['PrV','polynomial fit'])

        ax1 = plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        plt.savefig(save_path+'paperFigure_1_probabilityRescaling.png',dpi=300)
    


        # With angular tuning + scaling - to compensate for the reduced area
        plt.figure(figsize=(25,20))
        plt.rcParams.update({'font.size': 60})    

        plt.plot(       xFit_angularTuning,    yFit_angTuning(xFit_angularTuning),  ':',    linewidth=20,   color= 'lightseagreen', alpha=0.5)
        plt.plot(       angular_tuning['x'],   angular_tuning['y_scaled'],          'p',    markersize=50,  color= 'lightseagreen', alpha=0.9)
        
        plt.ylim(-0.1,3.1)
        plt.yticks([0,1,2,3])
        plt.ylabel('Normalized angular tuning')
        plt.xticks([0,1,2,3,4], ['0','-π/4','-π/2','-3π/4','-π'])
        plt.xlabel('Distance from target angle (rad)')
        # plt.legend(['PrV','polynomial fit'])

        ax1 = plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        plt.savefig(save_path+'paperFigure_1_probabilityRescaling_rescalingFcnOnly.png',dpi=300)
    
        if area_params is not None:

            plt.figure(figsize=(5,20))
            plt.rcParams.update({'font.size': 60})
            plt.bar([0,1,2],[area_params['auc_ref'],area_params['auc'],area_params['auc_fit']],color=['grey','k','lightseagreen'])

            plt.ylim(-0.1,4.1)
            plt.yticks([0,1,2,3,4])
            plt.ylabel('Area under the curve')
            plt.xticks([0,1,2], ['Unsc','ATun','Comp'],rotation=45)
            plt.xlabel('')
            # plt.legend(['PrV','polynomial fit'])

            ax1 = plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            
            plt.savefig(save_path+'paperFigure_areaUnderCurve_barPlot.png',dpi=300)



    def plotProbArray(prob_array, save_path=''):
        plt.figure( figsize=(30,20))
        for i in range(len(prob_array)): plt.bar(prob_array[i]['x'],prob_array[i]['y'],width=0.75)
        plt.savefig(save_path+'testing.png')






########################################################################################################################

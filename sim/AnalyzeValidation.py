import sys
import os
import json
from matplotlib import pyplot as plt
import NetPyNE_BBP
import numpy as np
from scipy.interpolate import interp1d


class AnalyzeValidation():

    def sort_dict(dictionary): return dict(sorted(dictionary.items(), key=lambda item: item[1], reverse=True))

    def getPearsonCorrelation(files):
        store_corr_coef={}
        for file in files:
            print(file)
            cellNumber = file.split('_')[1]
            with open(folderPath+file, 'r') as openfile: correlation_coefficient_dict = json.load(openfile)
            store_corr_coef.update({cellNumber:correlation_coefficient_dict})
        return store_corr_coef
    
    def find_matching_indices(reference_vector, target_vector):
        # Get the first and last timestamps from reference_vector
        start_time_ref = reference_vector[0]
        end_time_ref = reference_vector[-1]

        # Find the closest timestamps in target_vector to the start and end timestamps of reference_vector
        start_index_target = np.abs(np.array(target_vector) - start_time_ref).argmin()
        end_index_target = np.abs(np.array(target_vector) - end_time_ref).argmin()

        return start_index_target, end_index_target
    
    def resample_vector(vector, original_dt, target_dt):
        original_times = np.arange(0, len(vector) * original_dt, original_dt)
        target_times = np.arange(0, original_times[-1], target_dt)
        resampled_values = np.interp(target_times, original_times, vector)
        return resampled_values.tolist()

    def normalize_dict_values(data, original_min, original_max, new_min=0, new_max=1):
        """
        Normalize the values of dictionaries based on specified ranges.

        Args:
        - data (dict): Dictionary containing key-value pairs to be normalized.
        - original_min (float): Minimum value of the original data range.
        - original_max (float): Maximum value of the original data range.
        - new_min (float): Minimum value of the range for normalization (default: 0).
        - new_max (float): Maximum value of the range for normalization (default: 1).

        Returns:
        - dict: Normalized dictionary with the same keys as the input dictionary.
        """
        # Normalize each value in the dictionary
        normalized_data = {}
        for key, value in data.items():
            normalized_value = ((value - original_min) / (original_max - original_min)) * (new_max - new_min) + new_min
            normalized_data[key] = normalized_value
        
        return normalized_data
    
    def select_keys_with_highest_combined_values(normalized_a, normalized_b, normalized_c, min_b_value):
        # Combine normalized values of datasets "a" and "c"
        combined_values = {key: normalized_a.get(key, 0) + normalized_c.get(key, 0) for key in set(normalized_a) | set(normalized_c)}
        
        # Filter combined dataset to select keys with highest combined values
        highest_combined_keys = [key for key, value in sorted(combined_values.items(), key=lambda x: x[1], reverse=True)]
        
        # Filter dataset "b" to select keys with values above min_b_value
        filtered_b_keys = [key for key, value in normalized_b.items() if value > min_b_value]
        
        # Identify keys that satisfy both conditions
        selected_keys = list(set(highest_combined_keys) & set(filtered_b_keys))
        
        return selected_keys

    def getCummulativeCoeff(store_corr_coef): 
        cummulative_corr_coeff      = [sum(store_corr_coef[cellNumber].values()) for cellNumber in store_corr_coef.keys()]
        cummulative_corr_coeff_dict = {cellNumber:sum(store_corr_coef[cellNumber].values()) for cellNumber in store_corr_coef.keys()}
        cummulative_corr_coeff_dict = AnalyzeValidation.sort_dict(cummulative_corr_coeff_dict)
        return cummulative_corr_coeff,cummulative_corr_coeff_dict
    
    def getMeanCoeff(store_corr_coef): return [sum(store_corr_coef[cellNumber].values())/len(store_corr_coef[cellNumber].values()) for cellNumber in store_corr_coef.keys()]

    def plotPearsonCorrelationHistogram(store_corr_coef, filename):
        plt.figure(figsize=(10,10))
        cummulative_corr_coeff,cummulative_corr_coeff_dict=AnalyzeValidation.getCummulativeCoeff(store_corr_coef)
        mean_corr_coeff=AnalyzeValidation.getMeanCoeff(store_corr_coef)

        plt.hist(cummulative_corr_coeff,20)
        if '.png' not in filename:filename+='.png'
        plt.savefig(filename,dpi=300)
    
    def indices_of_top_percent_values(lst, x):
        if x <= 0 or x > 100: raise ValueError("Percentage must be between 0 and 100.")
        # Calculate the number of elements to consider
        num_elements = len(lst)
        num_to_select = int(num_elements * (x / 100))
        # Use sorted function to get indices of the top x% highest values
        indices = sorted(range(len(lst)), key=lambda i: lst[i], reverse=True)[:num_to_select]
        return indices
    
    def getBestCells(highest_cummulative_corr_coeff,store_corr_coef):
        import numpy as np
        best_cells=[]
        for ind,cellNumber in enumerate(store_corr_coef.keys()):
            if ind in highest_cummulative_corr_coeff: best_cells.append(cellNumber)
        
        best_cells_dict={cellNumber:np.mean(list(store_corr_coef[cellNumber].values())) for cellNumber in store_corr_coef.keys() if cellNumber in best_cells}

        return best_cells,best_cells_dict
    
    def getBBPSpiking(store_corr_coef, target_pop):
        from cell_cfg import cfg
        th_spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=cfg.th_spikes_file, cfg_file=cfg.sonataConfigFile, microcircuit_number=2, showFig=False)
        store_bbp_spiking_dict = {cellNumber:len(th_spk_dict[target_pop][int(cellNumber)]) for cellNumber in store_corr_coef.keys()}
        store_bbp_spiking_dict = AnalyzeValidation.sort_dict(store_bbp_spiking_dict)
        # store_bbp_spiking_dict = dict(sorted(store_bbp_spiking_dict.items(), key=lambda item: item[1], reverse=True))
        return store_bbp_spiking_dict

    def compareSims(cummulative_corr_coeff_dict,store_bbp_spiking_dict,filename = 'fig_compare_simulations.png'):
        x_data = [cummulative_corr_coeff_dict[cellNumber] for cellNumber in cummulative_corr_coeff_dict.keys()]
        y_data = [store_bbp_spiking_dict[cellNumber] for cellNumber in cummulative_corr_coeff_dict.keys()]

        # best_spiking_cells={cellNumber:(store_bbp_spiking_dict[cellNumber],cummulative_corr_coeff_dict[cellNumber]) for cellNumber in store_bbp_spiking_dict.keys() if (cummulative_corr_coeff_dict[cellNumber]>5.9) and store_bbp_spiking_dict[cellNumber]>2}
        
        plt.figure(figsize=(10,10))
        plt.scatter(x_data,y_data,s=10)
        if '.png' not in filename:filename+='.png'
        plt.savefig(filename,dpi=300)
        

    def compareOriginalSpiking(best_cells,store_bbp_spiking_dict):
        best_cells_bbp_spiking = {cellNumber:store_bbp_spiking_dict[cellNumber] for cellNumber in best_cells}
        return best_cells_bbp_spiking
    
    def get_BBP_NetPyNE_firingRate(spk_output_files):
        store_firing_rate={}
        store_correlation_coeff={}
        for file_ind, file in enumerate(spk_output_files):
            cellNumber = file.split('__')[1]
            print(cellNumber, 'remaining: ', len(spk_output_files)-file_ind)
            with open(folder_netpyne_bbp_spiking+file, 'r') as openfile: save_data = json.load(openfile)
            sim_gid = save_data['sim_gid']
            firing_bbp      = len([spkid for spkid in save_data['traces_data']['bbp']['spkt']])
            firing_netpyne  = len([spkid for spkid in save_data['traces_data']['netpyne']['spkt']])
            store_firing_rate.update({cellNumber:(firing_bbp,firing_netpyne)})
        return store_firing_rate

    def characterizeBBPSpiking(spk_output_files, plotFig=True):
        store_cell_IO={}
        store_correlation_coeff={}
        for file_ind, file in enumerate(spk_output_files):
            cellNumber = file.split('__')[1]
            print(cellNumber, 'remaining: ', len(spk_output_files)-file_ind)
            with open(folder_netpyne_bbp_spiking+file, 'r') as openfile: save_data = json.load(openfile)
            sim_gid = save_data['sim_gid']
            cell_input  = len([spkid for spkid in save_data['raster_data']['bbp']['spkid'] if int(spkid)!=int(sim_gid)])
            cell_output = len([spkid for spkid in save_data['raster_data']['bbp']['spkid'] if int(spkid)==int(sim_gid)])
            store_cell_IO.update({cellNumber:(cell_input,cell_output)})
            
            sampling_rate_netpyne = save_data['traces_data']['netpyne']['t'][1]-save_data['traces_data']['netpyne']['t'][0]
            sampling_rate_bbp = 0.1
            # sampling_rate_bbp = save_data['traces_data']['bbp']['t'][1]-save_data['traces_data']['bbp']['t'][0]

            # --- Getting voltage and time vectors for NetPyNE and BBP data
            BBPTime         = save_data['traces_data']['bbp']['t']
            BBPTrace        = save_data['traces_data']['bbp']['v']
            netpyneTime     = save_data['traces_data']['netpyne']['t']
            netpyneTrace    = save_data['traces_data']['netpyne']['v']

            # --- Getting the Start and End indexes of the NetPyNE time vector that match the first and last of the BBP time vector, because sim duration and sampling rate are different
            start_index_target, end_index_target = AnalyzeValidation.find_matching_indices(reference_vector=BBPTime, target_vector=netpyneTime)

            # --- Slicing the NetPyNE time and voltage vectors
            netpyneTrace_slice  = netpyneTrace[start_index_target:end_index_target]
            netpyneTime_slice   = netpyneTime[start_index_target:end_index_target]

            # --- Downsampling the NetPyNE time and voltage vectors
            netpyneTrace_downsampled  = AnalyzeValidation.resample_vector(netpyneTrace_slice, sampling_rate_netpyne, sampling_rate_bbp)
            netpyneTime_downsampled   = AnalyzeValidation.resample_vector(netpyneTime_slice, sampling_rate_netpyne, sampling_rate_bbp)

            # --- Slicing the vectors to make sure both traces have the same number of samples
            if len(netpyneTrace_downsampled)!=len(BBPTrace):
                print('fixing the size of the arrays')
                final_netpyneTrace = netpyneTrace_downsampled[  0:min([len(netpyneTrace_downsampled),len(BBPTrace)])]
                final_BBPTrace     = BBPTrace[                  0:min([len(netpyneTrace_downsampled),len(BBPTrace)])]
                final_netpyneTime  = netpyneTime_downsampled[   0:min([len(netpyneTrace_downsampled),len(BBPTrace)])]
                final_BBPTime      = BBPTrace[                  0:min([len(netpyneTrace_downsampled),len(BBPTrace)])]
            else:
                final_netpyneTrace = netpyneTrace_downsampled
                final_BBPTrace     = BBPTrace
                final_netpyneTime  = netpyneTime_downsampled
                final_BBPTime      = BBPTrace
         
            # Compute the Pearson correlation coefficient
            correlation_coefficient = np.corrcoef(final_netpyneTrace, final_BBPTrace)[0, 1]
            store_correlation_coeff.update({cellNumber:correlation_coefficient})
            store_correlation_coeff=AnalyzeValidation.sort_dict(store_correlation_coeff)

        if plotFig:
            plt.figure(figsize=(10,10))
            x_data = [cell_input for (cell_input,cell_output) in store_cell_IO.values()]
            y_data = [cell_output for (cell_input,cell_output) in store_cell_IO.values()]
            plt.scatter(x_data,y_data,s=10)
            plt.savefig('test_IO.png',dpi=300)

        return store_cell_IO, store_correlation_coeff

########################################################################################################################################################################################

if __name__ == '__main__':
    
    if len(sys.argv) == 1: 
        print ("Plotting selected GIDs for model implementation")
        plotFig=True
        selected_keys = ['37520', '41912', '36219', '36043', '34349', '29892','31315','28925','30252','32049']
        print(selected_keys)

    elif (len(sys.argv) == 2) and (int(sys.argv[1])==0 or int(sys.argv[1])==1):
        if int(sys.argv[1])==0:
            print ("Plotting selected TC GIDs")
            plotFig=True
            post_elected_cells={
                            'good':[37520,41912,36219,36043,34349,],
                            'ok':[37625,37409,],
                            'not_great':[34503,35847,35113,],
                        }
            selected_keys = list(post_elected_cells['good'])

        if int(sys.argv[1])==1: 
            print ("Plotting selected RT GIDs")
            plotFig=True
            post_elected_cells={
                            'good':[29892,31315,28925,30252,32049],
                            'ok':[30872,30007,32411],
                            'not_great':[],
                        }
            selected_keys = list(post_elected_cells['good'])
 
    elif (len(sys.argv) == 2) and (int(sys.argv[1])==2 or int(sys.argv[1])==3):
        if int(sys.argv[1])==2: 
            print('Testing TC gids')
            plotFig=True
            selected_keys = ['37520', '41912', '36219', '36043', '34349']
        if int(sys.argv[1])==3: 
            print('Testing RT gids')
            plotFig=True
            selected_keys = ['29892','31315','28925','30252','32049']
            # selected_keys = ['30872', '29892', '30007', '31315', '28925', '30252', '32411', '32049']
    elif (len(sys.argv) == 2) and (int(sys.argv[1])==4 or int(sys.argv[1])==5):

        if int(sys.argv[1])==4: 
            print ("Running analysis for a single TC cells")
            cellNumbers = [
                            35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847, # sim 01
                            33798, 34368, 36219, 39232, 34389, 34242, 38102, 35703, 38487, 41067, 37463, 38468, 36711, 34932, 38346, 34503, 36248, 41454, 36721, 33741, # sim 02
                            40602, 34274, 41534, 33640, 36881, 34859, 36169, 38276, 37409, 34707, 38440, 41237, 38052, 36302, 33602, 41247, 38036, 39429, 38474, 35824, # sim 03
                            38651, 37968, 40213, 42177, 40168, 40215, 41723, 36655, 38134, 41695, 42422, 42460, 36521, 38775, 35220, 35162, 34349, 36440, 35739, 34954, # sim 04
                            37256, 41168, 39751, 38748, 33967, 35343, 40876, 39755, 36185, 41399, 39299, 38971, 37093, 37917, 37599, 34471, 39745, 39477, 42073, 36043, # sim 05
                            41388, 38169, 34773, 34401, 41379, 37475, 38090, 40659, 37782, 38709, 42405, 41353, 41307, 40641, 37685, 39390, 39239, 35684, 34363, 37548, # sim 06
                            34976, 35398, 34977, 34209, 37751, 39276, 38218, 41138, 37435, 37966, 42345, 35864, 34506, 40105, 38470, 34418, 37141, 39362, 33676, 36674, # sim 07
                            36748, 36059, 35158, 40735, 35483, 42198, 34433, 41390, 39229, 40044, 37740, 40122, 36364, 35113, 38793, 40560, 36857, 37553, 41271, 39981, # sim 08 
                            41439, 38171, 39183, 41890, 37925, 37824, 38002, 35649, 41579, 38806, 37520, 40430, 33822, 39202, 37863, 41253, 33571, 35332, 35748, 39340, # sim 09
                            33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806, # sim 10
                            ]
            target_pop='VPL_TC'
            plotFig=False
        elif int(sys.argv[1])==5:
            print ("Running analysis for a single RT cells")
            cellNumbers = [
                            30550, 30541, 29582, 29174, 29908, 33109, 31521, 32579, 32893, 32954, 32696, 32933, 33334, 31927, 30299, 29934, 30694, 31191, 31989, 32369, # sim 11
                            30242, 30823, 29379, 31241, 31793, 31492, 32974, 30653, 29993, 30022, 29770, 32501, 29195, 29892, 30730, 30655, 32740, 32640, 28671, 28831, # sim 12
                            28660, 29828, 31704, 28988, 29183, 29690, 31254, 30838, 31637, 30922, 30182, 33200, 28663, 31412, 31625, 31778, 29791, 31120, 30543, 29184, # sim 13
                            28612, 30652, 32453, 32047, 29522, 32049, 29342, 31907, 30072, 32729, 29735, 32221, 30986, 33224, 31309, 30551, 31296, 29803, 29007, 30947, # sim 14
                            28805, 30849, 33463, 29657, 30946, 32631, 31840, 30892, 31646, 31738, 31315, 29086, 29040, 28852, 29608, 30025, 31528, 32662, 32781, 31170, # sim 15
                            32479, 33190, 31420, 28785, 30084, 31972, 30225, 30872, 30506, 32036, 33089, 33362, 32299, 32620, 29371, 32292, 32978, 32313, 32267, 30174, # sim 16
                            33014, 30007, 31239, 28733, 32470, 31044, 28694, 29087, 29476, 29687, 30990, 29126, 31800, 28834, 31881, 28925, 30252, 29621, 29094, 29304, # sim 17
                            31400, 29526, 31674, 32147, 31113, 29861, 32413, 29052, 30152, 29731, 29205, 31864, 31393, 33031, 30772, 28731, 30090, 33325, 30891, 29863, # sim 18
                            30403, 31638, 32406, 33043, 30905, 32926, 30014, 30813, 30854, 29679, 29049, 31751, 31816, 29689, 32540, 29846, 28833, 32411, 32730, 29805, # sim 19
                            32846, 29328, 30216, 32641, 29663, 30936, 32371, 29722, 31923, 30609, 32591, 30670, 31012, 31181, 33204, 31924, 32040, 28873, 33230, 31602, # sim 20
                            ]
            target_pop='Rt_RC'
            plotFig=False
        else:
            print ("Running analysis for a single TC cell")
            cellNumbers = [35103]
            plotFig=True

        folderPath = '../stats/validation/simulator_validation/traces_correlation/'
        files = [file for file in os.listdir(folderPath) if ('pearson_correlation.json' in file) and (int(file.split('_')[1]) in cellNumbers)]

        save_figure_folder='../stats/validation/_figures/correlation_statistics/'

        # --- Pearson correlation for the NetPyNE vs NEURON simulation
        store_corr_coef                 = AnalyzeValidation.getPearsonCorrelation(files)
        AnalyzeValidation.plotPearsonCorrelationHistogram(store_corr_coef,save_figure_folder+target_pop+'cummulative_corr_coeff_histogram.png')

        cummulative_corr_coeff,cummulative_corr_coeff_dict = AnalyzeValidation.getCummulativeCoeff(store_corr_coef)
        
        highest_cummulative_corr_coeff  = AnalyzeValidation.indices_of_top_percent_values(cummulative_corr_coeff,100)
        best_cells,best_cells_dict      = AnalyzeValidation.getBestCells(highest_cummulative_corr_coeff,store_corr_coef)

        # --- Get spike output from BBP data
        store_bbp_spiking_dict = AnalyzeValidation.getBBPSpiking(store_corr_coef, target_pop)

        # --- Pearson correlation for the NetPyNE validation vs BBP data
        folder_netpyne_bbp_spiking = '../stats/validation/init_sims/'
        spk_output_files = [file for file in os.listdir(folder_netpyne_bbp_spiking) if ('simplified_single_cell_traces_data__' in file) and (int(file.split('__')[1]) in cellNumbers)]
        store_cell_IO, store_correlation_coeff = AnalyzeValidation.characterizeBBPSpiking(spk_output_files)

        # save_correlation_data = {   'target_pop':target_pop,
        #                             'netpyne_neuron':{  'individual':store_corr_coef,
        #                                                 'sum':cummulative_corr_coeff_dict,
        #                                                 },
        #                             'netpyne_bbp'   :{  'vals':store_correlation_coeff
        #                                                 },
        #                             }
        # save_correlation_data_obj = json.dumps(save_correlation_data)
        # save_json_path = save_figure_folder+target_pop+'_traces_correlation_data.json'
        # with open(save_json_path, "w") as outfile: outfile.write(save_correlation_data_obj)

        store_firing_rate = AnalyzeValidation.get_BBP_NetPyNE_firingRate(spk_output_files)
        save_firing_rate_data = {   'target_pop':target_pop,
                                    'firing_rate':store_firing_rate
                                    }
        save_firing_rate_data_obj = json.dumps(save_firing_rate_data)
        save_firing_rate_json_path = save_figure_folder+target_pop+'_firing_rate_data.json'
        with open(save_firing_rate_json_path, "w") as outfile: outfile.write(save_firing_rate_data_obj)

        # --- Compare with output of the correlation coefficient
        AnalyzeValidation.compareSims(cummulative_corr_coeff_dict,store_bbp_spiking_dict,  filename = save_figure_folder+target_pop+'_pearson_validation_simulators_vs_bbp_spiking.png')
        
        AnalyzeValidation.compareSims(cummulative_corr_coeff_dict,store_correlation_coeff, filename = save_figure_folder+target_pop+'_pearson_validation_simulators_vs_peason_validation_netpyne.png')
        
        norm_cummulative_corr_coeff_dict = AnalyzeValidation.normalize_dict_values(cummulative_corr_coeff_dict, original_min=-6, original_max=6, new_min=-1, new_max=1)
        
        # norm_store_bbp_spiking_dict      = AnalyzeValidation.normalize_dict_values(store_bbp_spiking_dict,      original_min=0, original_max=max(store_bbp_spiking_dict.values()), new_min=0, new_max=1)
        
        norm_store_correlation_coeff     = store_correlation_coeff # no need to normalize

        min_norm_cumm_corr = 0.97
        min_spiking = 3
        store_candidates=[]
        for key in store_bbp_spiking_dict.keys():
            if (store_bbp_spiking_dict[key]>min_spiking) and (norm_cummulative_corr_coeff_dict[key]>min_norm_cumm_corr):
                store_candidates.append(key)
        print('Percentage of candidates compared to all dataset: ', (len(store_candidates)/len(store_bbp_spiking_dict.keys()))*100, ' percent')
        for candidate in store_candidates:
            print(candidate,': \tpearson NEURON_NetPyNE', norm_cummulative_corr_coeff_dict[candidate],': \t spiking',store_bbp_spiking_dict[candidate],': \tpearson simulation',norm_store_correlation_coeff[candidate])

        selected_keys = AnalyzeValidation.select_keys_with_highest_combined_values(norm_cummulative_corr_coeff_dict, store_bbp_spiking_dict, norm_store_correlation_coeff, min_b_value=3)

        for key in selected_keys:
            print(
                '\n norm_cummulative_corr_coeff_dict:   ', norm_cummulative_corr_coeff_dict[key],
                '\n store_bbp_spiking_dict:             ', store_bbp_spiking_dict[key],
                '\n norm_store_correlation_coeff:       ', norm_store_correlation_coeff[key],
                )
    
    ########################################################################################################################################################################################
    if plotFig:
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

        filenames = ['../stats/validation/init_sims/simplified_single_cell_traces_data__'+str(cellNumber)+'__dict.json' for cellNumber in selected_keys]
        from analyze_init import AnalyzeOriginal
        for filename in filenames:
            AnalyzeOriginal.plotTraces(filename=filename)
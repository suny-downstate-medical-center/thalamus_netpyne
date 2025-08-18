#! /Users/joao/opt/anaconda3/bin/ipython
'''
Class to build the stimulation dataset
Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''
########################################################################################################################
import numpy as np
import json
import math
import pandas as pd
########################################################################################################################
from GenerateStimDataset import LoadPSTH
from GenerateStimDataset import PlotFigs
from GenerateStimDataset import SampleData
from GenerateStimDataset import FitData
########################################################################################################################
NetPyNE_rootFolder              = '/Users/joao/Research/Models/BBP/thalamus_netpyne'
stim_source_code_folder         = '/stims/Minnery_2003/fig1'
stim_dataset_folder             = NetPyNE_rootFolder+stim_source_code_folder+'/stim_datasets/'      # Folder with the digitized PSTH
save_figs_folder                = NetPyNE_rootFolder+stim_source_code_folder+'/figs_stim/'          # Folder to save the output figures
save_deflection_model_folder    = NetPyNE_rootFolder+stim_source_code_folder+'/deflection_model/'   # Folder to save the output deflection model and sampled raster

generate_ct_spikes=False # (default: False)
# generate_ct_spikes=True # (default: False)

MLe_angular_tuning_source = 'minnery' # (default)

plot_paper_figure=False
if plot_paper_figure:
    # --- Data save filename
    file_header_dataset = 'mleSpikes_'
    file_header_deflection_dict = 'deflectionModel_'
    n_cells = 180
    angle_range=[0,360]
    MLe_angular_tuning_source = 'minnery'

else:
    if generate_ct_spikes:
        # # --- Data save filename
        # file_header_dataset = 'ctSpikes2_'
        # file_header_deflection_dict = 'ctSpikesModel_'
        # n_cells = 818
        # angle_range=[0,818]
        # MLe_angular_tuning_source = 'minnery'

        # # --- Data save filename
        # file_header_dataset = 'ctSpikes3_'
        # file_header_deflection_dict = 'ctSpikesModel_'
        # n_cells = 180
        # angle_range=[0,360]
        # MLe_angular_tuning_source = 'minnery'

        # # --- Data save filename
        # file_header_dataset = 'ctSpikes4_'
        # file_header_deflection_dict = 'ctSpikesModel4_'
        # n_cells = 360
        # angle_range=[0,360]
        # MLe_angular_tuning_source = 'minnery'
        
        # # --- Data save filename - Current CT feedback
        # file_header_dataset = 'ctSpikes41_'
        # file_header_deflection_dict = 'ctSpikesModel41_'
        # n_cells = 360
        # angle_range=[0,360]
        # MLe_angular_tuning_source = 'hartings_compensated2'


        ######################################################################## 
        
        # 10% / 25% / 50% / 100% ratios
        # --- Data save filename - test
        file_header_dataset = 'ctSpikes60_'
        file_header_deflection_dict = 'ctSpikesModel60_'
        n_cells = 82
        angle_range=[0,82]
        MLe_angular_tuning_source = 'hartings_compensated2'

        # # --- Data save filename - test
        # file_header_dataset = 'ctSpikes61_'
        # file_header_deflection_dict = 'ctSpikesModel61_'
        # n_cells = 204
        # angle_range=[0,204]
        # MLe_angular_tuning_source = 'hartings_compensated2'

        # # --- Data save filename - test
        # file_header_dataset = 'ctSpikes62_'
        # file_header_deflection_dict = 'ctSpikesModel62_'
        # n_cells = 409
        # angle_range=[0,409]
        # MLe_angular_tuning_source = 'hartings_compensated2'

        # # --- Data save filename - test
        # file_header_dataset = 'ctSpikes63_'
        # file_header_deflection_dict = 'ctSpikesModel63_'
        # n_cells = 818
        # angle_range=[0,818]
        # MLe_angular_tuning_source = 'hartings_compensated2'

        ######################################################################## 


        # # --- Data save filename
        # file_header_dataset = 'ctSpikes5_'
        # file_header_deflection_dict = 'ctSpikesModel5_'
        # n_cells = 409
        # angle_range=[0,409]
        # MLe_angular_tuning_source = 'minnery'

        # # --- Data save filename
        # file_header_dataset = 'ctSpikes6_'
        # file_header_deflection_dict = 'ctSpikesModel6_'
        # n_cells = 409
        # angle_range=[0,409]
        # MLe_angular_tuning_source = 'hartings_compensated2'
    else:
        # # --- Data save filename - Minnery MLe angular tuning
        # file_header_dataset = 'mleSpikes_'
        # file_header_deflection_dict = 'deflectionModel_'
        # n_cells = 180
        # angle_range=[0,360]
        # MLe_angular_tuning_source = 'minnery'

        # --- Data save filename - Hartings MLe angular tuning
        file_header_dataset = 'mleSpikes2_'
        file_header_deflection_dict = 'deflectionModel2_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings'

        # --- Data save filename - Hartings MLe angular tuning
        file_header_dataset = 'mleSpikes3_'
        file_header_deflection_dict = 'deflectionModel3_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings_compensated'
        
        # --- Data save filename - Hartings MLe angular tuning
        file_header_dataset = 'mleSpikes4_'
        file_header_deflection_dict = 'deflectionModel4_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings_compensated2'

        # # # --- Data save filename - Hartings MLe angular tuning - 
        # # file_header_dataset = 'mleSpikes5_'
        # # file_header_deflection_dict = 'deflectionModel5_'
        # # n_cells = 180
        # # angle_range=[0,360]
        # # MLe_angular_tuning_source = 'hartings_compensated3'
        # # --- Data save filename - Hartings MLe angular tuning - 
        # file_header_dataset = 'mleSpikes55_'
        # file_header_deflection_dict = 'deflectionModel55_'
        # n_cells = 180
        # angle_range=[0,360]
        # MLe_angular_tuning_source = 'hartings'
        
        file_header_dataset = 'mleSpikes6_'
        file_header_deflection_dict = 'deflectionModel6_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings'
        
        file_header_dataset = 'mleSpikes7_'
        file_header_deflection_dict = 'deflectionModel7_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings'

        file_header_dataset = 'mleSpikes71_'
        file_header_deflection_dict = 'deflectionModel71_'
        n_cells = 180
        angle_range=[0,360]
        MLe_angular_tuning_source = 'hartings_scaled'

        # # --- Data save filename - commands to create a dataset to stimulate the TRN_ring pop using VPM baseline firing
        # file_header_dataset = 'baselineVPMSpikes_'
        # file_header_deflection_dict = 'baselineVPMModel_'
        # n_cells = 2880
        # angle_range=[0,360]

########################################################################################################################
UpdateDeflectionModel       = True

if UpdateDeflectionModel: 
    SampleRaster            = True
    UpdateSampledRasterData = True
else:
    SampleRaster = True
    UpdateSampledRasterData = False

# UpdateSampledRasterFigs  = False
UpdateSampledRasterFigs  = True
saveSimplifiedRasterData = True # Saves an aditional dictionary dataset with only the spike times organized by GID - lighter, to be synchronized in Github
save_format              = 'pkl'
########################################################################################################################
# --- Plot figures
plotFigs = True
########################################################################################################################
# --- scaling probability to adjust for a higher re-sampling rate
rescale_sampling = True            # (True: 0.025 ms) / (False: 1.0 ms)
# rescale_sampling = False            # (True: 0.025 ms) / (False: 1.0 ms)

dt_step_PSTH    = 1.0    # (ms) - extracted from paper methods (measured at 0.1 ms resolution, and binned at 1.0 ms bins)
dt_step_NetPyNE = 0.025  # (ms) - from netpyne simulation cfg file

if rescale_sampling: dt_step_resampling = dt_step_NetPyNE
else:                dt_step_resampling = dt_step_PSTH

########################################################################################################################
# --- Angle to sample rasters from
# Test angle to make paper figure
if plot_paper_figure: 
    target_angles = [148]
    deflection_events = [[2000,2250]]
    last_spike_time=[3500,3520]
else:
        
    # target_angles = [180]
    target_angles = [0,45,90,135,225,270,315,None]
    # target_angles = [45, 180]
    # target_angles = [0,45,90,135,180,225,270,315,None]
    # target_angles = [180,225,270,315,None]
    # target_angles = [0,45,90,135]
    # target_angles = [0,45,90,135,225,270,315,None]
    # target_angles = [180]
    # target_angles = [None]

    # deflection_events = [[500,650],[1000,1250],[1500,1850]]
    # deflection_events = [[2000,2150],[3000,3150],[4000,4250],[5000,5250],[6000,6350],[7000,7350]]
    deflection_events = [[2000,2250],[3000,3250],[4000,4250],[5000,5250],[6000,6250],[7000,7250],[8000,8250],[9000,9250]] # 2024_08_07 - changed all deflection times to 250 ms (like original paper) to allow statistical analysis

    last_spike_time=[10000,10020]

    # Transforms target_angles into the [0-818] equivalent of [0-360]
    if generate_ct_spikes:
        if None in target_angles:   addNone = True
        else:                       addNone = False
        target_angles=[round(((angle_range[1]-angle_range[0])*target_angle)/360) for target_angle in target_angles if target_angle is not None]
        if addNone:target_angles.append(None)

        ct_delay = round(((2/3)*13.1),2) # (Hirai, 2018) - Whiskers->L6A:  13.1 ± 2.5 ms   //  Assuming 2/3 of the time, since inputs come from brainstem, 1 synapse above in the hierarchy

        deflection_events = [[deflection_on+ct_delay,deflection_off+ct_delay] for [deflection_on,deflection_off] in deflection_events]

if generate_ct_spikes:  save_full_dataset = False
else:                   save_full_dataset = True

########################################################################################################################

if UpdateDeflectionModel:
    # Specify the path to your CSV file
    folderPath = stim_dataset_folder
    fileNames = ['PrV.csv','VPM.csv','0stim.csv']
    hist_dict = LoadPSTH.load(folderPath, fileNames)
    if plotFigs: PlotFigs.plotAllTraces(hist_dict, save_path = save_figs_folder)
    ########################################################################################################################
    # # --- Skips the samples for the rise time, until I find a better model to fit the curve (poisson, maybe?)
    # target_nuclei = 'PrV'

    # ON_start  = 62
    # OFF_start = 202

    # skip_ON  = 3 # samples
    # skip_OFF = 2 # samples

    # # --- Selecting and formating the data from PrV
    # PrV_data = LoadPSTH.formatData(hist_dict,target_nuclei='PrV', ON_start=ON_start, OFF_start=OFF_start, skip_ON=skip_ON, skip_OFF=skip_OFF)

    # # --- Calculating the x_data shifted to x=0
    # x_data_ON_zero  = [x-min(PrV_data['ON']['x'])  for x in PrV_data['ON']['x']]
    # x_data_OFF_zero = [x-min(PrV_data['OFF']['x']) for x in PrV_data['OFF']['x']]

    # # --- Updating the dataset with a zero aligned data
    # PrV_data.update({'ON_zero' :{'x':x_data_ON_zero, 'y':PrV_data['ON']['y']}})
    # PrV_data.update({'OFF_zero':{'x':x_data_OFF_zero,'y':PrV_data['OFF']['y']}})

    # #######################################################################################################################################
    # # --- Fitting the baseline
    # baseline_mean = np.mean(PrV_data['baseline']['y'])
    # baseline_std  = np.std( PrV_data['baseline']['y'])

    # fitDegree = 25

    # # --- DATA ON ORIGINAL TIMESTAMPS
    # # --- Fits polynomial to data
    # yPoly_ON       = FitData.fitPolynomial(x_data=PrV_data['ON']['x'],       y_data=PrV_data['ON']['y'], degree=fitDegree)
    # yPoly_OFF      = FitData.fitPolynomial(x_data=PrV_data['OFF']['x'],      y_data=PrV_data['OFF']['y'], degree=fitDegree)
    # # --- Creates polynomial Object to sample from
    # yFit_ON       = FitData.getPolynomialObj(yPoly_ON)
    # yFit_OFF      = FitData.getPolynomialObj(yPoly_OFF)

    # # --- DATA ALIGNED TO ZERO (x0=0)
    # # --- Fits polynomial to data
    # yPoly_ON_zero  = FitData.fitPolynomial(x_data=PrV_data['ON_zero']['x'],  y_data=PrV_data['ON']['y'], degree=fitDegree)
    # yPoly_OFF_zero = FitData.fitPolynomial(x_data=PrV_data['OFF_zero']['x'], y_data=PrV_data['OFF']['y'], degree=fitDegree)
    # # --- Creates polynomial Object to sample from
    # yFit_ON_zero   = FitData.getPolynomialObj(yPoly_ON_zero)
    # yFit_OFF_zero  = FitData.getPolynomialObj(yPoly_OFF_zero)

    # if plotFigs: PlotFigs.plotFits(PrV_data,yFit_ON,yFit_OFF,yFit_ON_zero,yFit_OFF_zero,nSamples=250, save_path = save_figs_folder)

    #######################################################################################################################################
    # === Updated model
    
    # --- Skips the samples for the rise time, until I find a better model to fit the curve (poisson, maybe?)
    # --- Solution: skip samples, starting from peak value for each curve, and use a triple and double exponential decays to fit the ON and OFF events, respectively
    target_nuclei = 'PrV'

    ON_start  = 62
    OFF_start = 202

    skip_ON  = 3 # samples
    skip_OFF = 2 # samples

    # --- Selecting and formating the data from PrV
    PrV_data = LoadPSTH.formatData(hist_dict,target_nuclei='PrV', ON_start=ON_start, OFF_start=OFF_start, skip_ON=skip_ON, skip_OFF=skip_OFF)

    # --- Calculating the x_data shifted to x=0
    x_data_ON_zero  = [x-min(PrV_data['ON']['x'])  for x in PrV_data['ON']['x']]
    x_data_OFF_zero = [x-min(PrV_data['OFF']['x']) for x in PrV_data['OFF']['x']]

    # --- Updating the dataset with a zero aligned data
    PrV_data.update({'ON_zero' :{'x':x_data_ON_zero, 'y':PrV_data['ON']['y']}})
    PrV_data.update({'OFF_zero':{'x':x_data_OFF_zero,'y':PrV_data['OFF']['y']}})

    # --- Fitting the baseline
    baseline_mean = np.mean(PrV_data['baseline']['y'])
    baseline_std  = np.std( PrV_data['baseline']['y'])

    if generate_ct_spikes:
        '''
        import pandas as pd
        # --- Import a dataset of baseline firing generated based on the mean firing of <baseline_mean = np.mean(PrV_data['baseline']['y'])>
        deflection_dataset_path = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
        # --- Load the Pickle file into a Python dictionary
        with open(deflection_dataset_path, 'rb') as f: spikes_dict = pd.read_pickle(f)
        group_spikes = [spkt for cell_id in spikes_dict['spkts'].keys() for spkt in spikes_dict['spkts'][cell_id] if (spkt<10000)]
        group_spikes.sort()
        time_s=10;num_cells = len(spikes_dict['spkts'].keys())
        mean_firing = len(group_spikes)/(time_s*num_cells)
        print('Mean firing: ', mean_firing, ' Hz')
        ct_target_baseline_frequency = 5 # (Hz)
        ct_probability_rescaling = mean_firing/ct_target_baseline_frequency
        print('CT baseline probability rescaling factor: ', ct_probability_rescaling)

        # 
        Mean firing:                                21.57444987775061  Hz
        CT baseline probability rescaling factor:   4.3148899755501215

        '''
        # --- Fitting the baseline
        baseline_mean = np.mean(PrV_data['baseline']['y'])/4.3148899755501215
        baseline_std  = np.std( PrV_data['baseline']['y'])/4.3148899755501215
    else:
        # --- Fitting the baseline
        baseline_mean = np.mean(PrV_data['baseline']['y'])
        baseline_std  = np.std( PrV_data['baseline']['y'])
    
    print('Baseline mean/std: ', baseline_mean, '+/-', baseline_std)

    # --- DATA ALIGNED TO ZERO (x0=0)
    # --- Fits polynomial to data

    # yExpTriple_ON  = FitData.fit_double_exp(x_data=PrV_data['ON']['x'],  y_data=PrV_data['ON']['y'])
    # yExpDouble_OFF = FitData.fit_double_exp(x_data=PrV_data['ON']['x'],  y_data=PrV_data['ON']['y'])
    # # --- Creates polynomial Object to sample from
    # yFit_ON  = FitData.getFittedModel(params=yExpTriple_ON,  model_type='triple_exp')
    # yFit_OFF = FitData.getFittedModel(params=yExpDouble_OFF, model_type='double_exp')
    
    if generate_ct_spikes:
        # --- Reducing probability evenly
        raw_baseline = np.mean(PrV_data['baseline']['y'])
        y_probs_ON  = [(y_prob/4.3148899755501215) for y_prob in PrV_data['ON']['y']]
        y_probs_OFF = [(y_prob/4.3148899755501215) for y_prob in PrV_data['OFF']['y']]
        # # --- Reducing only the baseline
        # raw_baseline = np.mean(PrV_data['baseline']['y'])
        # y_probs_ON  = [(y_prob-raw_baseline+baseline_mean) for y_prob in PrV_data['ON']['y']]
        # y_probs_OFF = [(y_prob-raw_baseline+baseline_mean) for y_prob in PrV_data['OFF']['y']]
    else:
        y_probs_ON  = PrV_data['ON']['y']
        y_probs_OFF = PrV_data['OFF']['y']
    
    print('y deflection mean - ON : ', np.mean(y_probs_ON), ' | min value: ', min(y_probs_ON))
    print('y deflection mean - OFF: ', np.mean(y_probs_OFF), ' | min value: ', min(y_probs_OFF))

    yExpTriple_ON_zero  = FitData.fit_triple_exp(x_data=PrV_data['ON_zero']['x'],  y_data=y_probs_ON)
    yExpDouble_OFF_zero = FitData.fit_double_exp(x_data=PrV_data['OFF_zero']['x'], y_data=y_probs_OFF)
    # --- Creates polynomial Object to sample from
    yFit_ON_zero  = FitData.getFittedModel(params=yExpTriple_ON_zero,  model_type='triple_exp')
    yFit_OFF_zero = FitData.getFittedModel(params=yExpDouble_OFF_zero, model_type='double_exp')
    
    if plotFigs: PlotFigs.plotFits_zero(PrV_data,yFit_ON_zero,yFit_OFF_zero,nSamples=250, save_path = save_figs_folder)
    print(save_figs_folder)
    # import sys;sys.exit()

    #######################################################################################################################################
    # --- Tunning function - Extracted from Fig 3A

    if      MLe_angular_tuning_source=='hartings':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[  1.0,     0.64055063,         0.32471083,         0.20137985,          0.19764487038033318]}
    
    elif      MLe_angular_tuning_source=='hartings_scaled':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[   1.0,
                                    0.64055063 * 0.85,
                                    0.32471083 * 0.85,
                                    0.20137985,
                                    0.19764487038033318]}
    
    
    elif      MLe_angular_tuning_source=='hartings_compensated':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[   1.0                     * 1.0,     
                                    0.64055063              * 0.9131733746489962,         
                                    0.32471083              * 0.7872142060514346,         
                                    0.20137985              * 0.8124547209392146,          
                                    0.19764487038033318     * 1.1504130324564774,
                                    ]} # adds a compensation factor to approximate the measured vs the observed values of angular tunning for each angle
    elif      MLe_angular_tuning_source=='hartings_compensated2':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[   1.0                     * 1.0,     
                                    0.64055063              * 0.9131733746489962,         
                                    0.32471083              * 0.75,         
                                    0.20137985              * 0.8124547209392146,          
                                    0.19764487038033318     * 1.45,
                                    ]} # adds a compensation factor to approximate the measured vs the observed values of angular tunning for each angle
    elif      MLe_angular_tuning_source=='hartings_compensated3':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[   1.0                     * 1.0,     
                                    0.64055063              * 0.9131733746489962,         
                                    0.32471083              * 0.75,         
                                    0.20137985              * 0.95,     # recalibrated to increase response in the opposing direction          
                                    0.19764487038033318     * 2.0,      # recalibrated to increase response in the opposing direction
                                    ]} # adds a compensation factor to approximate the measured vs the observed values of angular tunning for each angle
    elif      MLe_angular_tuning_source=='hartings_compensated4':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[   1.0                     * 1.0,     
                                    0.64055063              * 0.9131733746489962,         
                                    0.32471083              * 1.0,
                                    0.20137985              * 1.0,     # recalibrated to increase response in the opposing direction          
                                    0.19764487038033318     * 1.0,      # recalibrated to increase response in the opposing direction
                                    ]} # adds a compensation factor to approximate the measured vs the observed values of angular tunning for each angle
    elif    MLe_angular_tuning_source=='minnery':
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[  1.0,     0.6982148208637337, 0.5676691729323308, 0.44418423989406364, 0.4282115869017633]}
    else:
        angular_tuning = {  'x':[  0,       1,                  2,                  3,                   4],
                            'y':[  1.0,     0.6982148208637337, 0.5676691729323308, 0.44418423989406364, 0.4282115869017633]}

    print('Angular tuning from - ', MLe_angular_tuning_source)
    print('Angular tuning values: ', angular_tuning)

    coords=[(angular_tuning['x'][i],angular_tuning['y'][i]) for i in range(len(angular_tuning['x']))]

    # --- Adding automatic probability compensation so that angular tuning doesn't result in a reduction of spikes, and replacing manual prob adjustment
    # Calculating the area of the polygon defined by the angular tuning curve
    x_coords = np.array(angular_tuning['x'])
    y_coords = np.array(angular_tuning['y'])
    y_ref    = np.array([1 for i in range(len(angular_tuning['x']))])
    # Area under the curve
    auc     = np.trapz(y_coords, x_coords)
    auc_ref = np.trapz(y_ref,    x_coords)
    auc_ratio = auc_ref/auc
    # Rescaled angular tuning compensation
    angular_tuning.update({'y_scaled': [y_val*auc_ratio for y_val in angular_tuning['y']]})

    # yPoly_angTuning = FitData.fitPolynomial(x_data=angular_tuning['x'], y_data=angular_tuning['y'], degree=3)
    yPoly_angTuning = FitData.fitPolynomial(x_data=angular_tuning['x'], y_data=angular_tuning['y_scaled'], degree=3)

    poly = np.poly1d(yPoly_angTuning)  
    from scipy.integrate import quad
    auc_fit, auc_fit_error = quad(poly, angular_tuning['x'][0], angular_tuning['x'][-1])

    # yFit_angTuning  = FitData.getPolynomialObj(yPoly_angTuning)
    yFit_angTuning  = FitData.getFittedModel(params=yPoly_angTuning, model_type='polynomial')

    if plotFigs: PlotFigs.plotAngularTuning(angular_tuning,nSamples=250,yFit_angTuning=yFit_angTuning, save_path = save_figs_folder, area_params={'auc_ref':auc_ref,'auc':auc,'auc_fit':auc_fit})

    #######################################################################################################################################
    # --- Make a histogram of spikes

    # --- Dictionary to store the Sampled Data
    sampledProb={}

    # --- Sampling probabilities from the baseline
    baseline_duration = math.floor(max(PrV_data['baseline']['x']))
    sampledProb.update({'baseline':SampleData.sampleNormal(mean=baseline_mean,std=baseline_std,start=0,stop=baseline_duration,dt=dt_step_resampling)})

    # # --- Sampling probabilities from the polynomial fit
    # sampledProb.update({'ON' :SampleData.sampleFitFunction(yFit_ON , PrV_data['ON']['x'] , dt_step_resampling)})
    # sampledProb.update({'OFF':SampleData.sampleFitFunction(yFit_OFF, PrV_data['OFF']['x'], dt_step_resampling)})

    # --- Sampling probabilities from the polynomial fit - dataset with time set to ZERO
    sampledProb.update({'ON_zero' :SampleData.sampleFitFunction(yFit_ON_zero , PrV_data['ON_zero']['x'] , dt_step_resampling)})
    sampledProb.update({'OFF_zero':SampleData.sampleFitFunction(yFit_OFF_zero, PrV_data['OFF_zero']['x'], dt_step_resampling)})
    zero_times = [PrV_data['ON']['x'][0],PrV_data['OFF']['x'][0]]
    # --- Plotting figs
    # if plotFigs: PlotFigs.plotSampledDataBar(sampledProb,save_path = save_figs_folder) # Fig 8
    # if plotFigs: PlotFigs.plotFitAndSampled(PrV_data, yFit_ON, yFit_OFF, nSamples=250, sampledData=sampledProb, save_path = save_figs_folder) # Fig 81
    if plotFigs: PlotFigs.plotFitAndSampled_zero(PrV_data, yFit_ON_zero, yFit_OFF_zero, nSamples=250, sampledData=sampledProb, save_path = save_figs_folder,zero_times=zero_times, baseline_data=[baseline_mean,baseline_std]) # Fig 81
    # import sys;sys.exit()

    #######################################################################################################################################
    # --- Generating whisker deflection dictionary
    for target_angle in target_angles:
        if target_angle is not None:    deflection_times = deflection_events+[last_spike_time]
        else:                           deflection_times = [last_spike_time]

        print('Updating whisker deflection Model')
        if target_angle is not None:
            deflection_dict={   
                                # 'deflection_times':             [[500,650],[1000,1250],[1500,1850],[10000,10020]],
                                'deflection_times':             deflection_times,
                                'baseline_prob':                [baseline_mean,baseline_std],
                                'ON_polyFit':                   list(yExpTriple_ON_zero),
                                'OFF_polyFit':                  list(yExpDouble_OFF_zero),
                                # 'ON_polyFit':                   list(yPoly_ON_zero),
                                # 'OFF_polyFit':                  list(yPoly_OFF_zero),
                                'resetTime':                    True,
                                'dt':                           dt_step_resampling,
                                'plotFigs':                     True,
                                'target_angle':                 target_angle,
                                'n_cells':                      n_cells,
                                'angle_range':                  angle_range,
                                # 'prob_compensation':            1.25,
                                'prob_compensation':            1,
                                'yPoly_angTuning':              list(yPoly_angTuning),
                                'maximize_target_angle':        True,
                                'save_stim_data':               True,
                                'analyze_stim_data':            True,
                                'paper_sampling':               dt_step_PSTH,
                                'addMouseHead':                 True,
                                'stim_dataset_folder':          stim_dataset_folder,
                                'save_figs_folder':             save_figs_folder,
                                'save_deflection_model_folder': save_deflection_model_folder,
                                }
            stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in deflection_times)
        else:
            # --- Generates a baseline stimulation dataset
            deflection_dict={   
                                'deflection_times':             deflection_times,
                                'baseline_prob':                [baseline_mean,baseline_std],
                                'ON_polyFit':                   list(yExpTriple_ON_zero),
                                'OFF_polyFit':                  list(yExpDouble_OFF_zero),
                                # 'ON_polyFit':                   list(yPoly_ON_zero),
                                # 'OFF_polyFit':                  list(yPoly_OFF_zero),
                                'resetTime':                    True,
                                'dt':                           dt_step_resampling,
                                'plotFigs':                     True,
                                'target_angle':                 target_angle,
                                'n_cells':                      n_cells,
                                'angle_range':                  angle_range,
                                'prob_compensation':            1,
                                'yPoly_angTuning':              None,
                                'maximize_target_angle':        False,
                                'save_stim_data':               True,
                                'analyze_stim_data':            True,
                                'paper_sampling':               dt_step_PSTH,
                                'addMouseHead':                 True,
                                'stim_dataset_folder':          stim_dataset_folder,
                                'save_figs_folder':             save_figs_folder,
                                'save_deflection_model_folder': save_deflection_model_folder,
                                }
            stims_string = '_stims|None'

        # --- Saves the stim_string into the deflection dictionary for retrieval
        deflection_dict.update({'stims_string':stims_string})

        # --- Saves the sampling interval (dt) into the deflection dictionary for retrieval
        dt_string = '_samplingBins|'+str(int(dt_step_resampling*1000))+'us'

        # --- Save the deflection dictionary
        with open(save_deflection_model_folder+file_header_deflection_dict+str(deflection_dict['target_angle'])+'|deg'+dt_string+stims_string+'.json', 'w') as fp: json.dump(deflection_dict, fp, indent=4)

if SampleRaster:
    #######################################################################################################################################
    # --- Sample whisker deflection from deflection model
    for target_angle in target_angles:
        if target_angle is not None:
            deflection_times = deflection_events+[last_spike_time]
            stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in deflection_times)
        else:                           
            deflection_times = [last_spike_time]
            stims_string = '_stims|None'
        
        # --- Saves the sampling interval (dt) into the deflection dictionary for retrieval
        dt_string = '_samplingBins|'+str(int(dt_step_resampling*1000))+'us'


        deflection_dict = SampleData.LoadJSON(file_path=save_deflection_model_folder+file_header_deflection_dict+str(target_angle)+'|deg'+dt_string+stims_string+'.json')
        # deflection_dict = SampleData.LoadJSON(file_path=save_deflection_model_folder+file_header_deflection_dict+str(target_angle)+'|deg.json')
        print(' \n\n\n##################################################################')
        print(' ---> Target Angle: ', target_angle)
        

        # --- Make a new raster
        if UpdateSampledRasterData:
            print('Generating dataset')
            spikes_dict = SampleData.deflectionEvent(deflection_dict)
            if save_full_dataset:
                if save_format=='json':
                    with open(deflection_dict['save_deflection_model_folder']+'spike_dicts/'+file_header_dataset+'deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.json', 'w') as fp: json.dump(spikes_dict, fp, indent=4)
                else:
                    f_name=deflection_dict['save_deflection_model_folder']+'spike_dicts/'+file_header_dataset+'deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.pkl'
                    with open(f_name, 'wb') as f: pd.to_pickle(spikes_dict, f)

            if UpdateSampledRasterFigs: SampleData.analyzeDeflectionEvent(spikes_dict)
        
        # --- Load saved raster
        else:
            try:
                print('Loading dataset')
                # obs: dt_step_resampling is redefined, so that when loading you can skip updating the deflection model
                file_path = save_deflection_model_folder+'spike_dicts/'+file_header_dataset+'deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.pkl'
                # spikes_dict = SampleData.LoadJSON(file_path)
                spikes_dict = SampleData.LoadPickle(file_path)
                if UpdateSampledRasterFigs: SampleData.analyzeDeflectionEvent(spikes_dict)
                print('Loading dataset - worked')
                print('Loaded spike dict - target angle: ',spikes_dict['target_angle'])
            except:
                print('Loading dataset - failed')

        # Saves an aditional dictionary dataset with only the spike times organized by GID - lighter, to be synchronized in Github
        try:
            if saveSimplifiedRasterData:
                print('Saving simplified dataset')
                simplified_spikes_dict={'spkts':spikes_dict['spkts']}
                # obs: dt_step_resampling is redefined, so that when loading you can skip updating the deflection model
                f_name=deflection_dict['save_deflection_model_folder']+'spike_dicts/'+file_header_dataset+'deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'_simplifiedDataset'+'.pkl'
                with open(f_name, 'wb') as f: pd.to_pickle(simplified_spikes_dict, f)
                print('Saving simplified dataset - worked')
        except:
            print('Saving simplified dataset - failed')

        if plot_paper_figure:
            # spkt_min=1750;spk_max=10020 # for ct spikes
            spkt_min=1750;spk_max=2500
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 100})
            figsize=(50,25)
            plt.figure(figsize=figsize)
            store_scatter_spikes=[]
            for spk_gid in spikes_dict['spkts'].keys():
                scatter_spikes = [spkt-spkt_min for spkt in spikes_dict['spkts'][spk_gid] if (spkt>spkt_min) and (spkt<=spk_max)]
                plt.scatter(scatter_spikes,[spk_gid for i in range(len(scatter_spikes))],c='k')
                store_scatter_spikes+=scatter_spikes # store spikes to make histogram plot
            tick_angles=[i*45 for i in range(5)]
            tick_labels=[str(int(2*i*45)) for i in range(5)]
            plt.ylabel('Brainstem cell angle (deg)')
            plt.yticks(tick_angles,tick_labels)
            # plt.xlabel('Time (ms)')
            # plt.xticks([i*250 for i in range(4)])
            plt.xlim([-10,750])
            plt.xticks([])
            # --- Remove corners
            ax=plt.gca()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            # plt.tight_layout()
            plt.savefig('../paper_figs/deflection_model/00_raster_'+file_header_dataset+'_'+str(deflection_dict['target_angle'])+'.png')

            
            store_scatter_spikes.sort()
            bin_size_ms=1
            plt.figure(figsize=figsize)
            plt.hist(store_scatter_spikes,range=[0,spk_max-spkt_min],bins=int((spk_max-spkt_min)/bin_size_ms),color='k')
            plt.xticks([i*250 for i in range(4)])
            plt.ylim([0,85])
            plt.yticks([0,25,50,75])
            plt.ylabel('Spikes (1/ms)')
            plt.xlim([-10,750])
            plt.xlabel('Time (ms)')
            # --- Remove corners
            ax=plt.gca()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            # plt.tight_layout()
            plt.savefig('../paper_figs/deflection_model/01_PSTH_'+file_header_dataset+'_'+str(deflection_dict['target_angle'])+'.png')



        #######################################################################################################################################

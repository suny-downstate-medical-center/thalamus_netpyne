#!/bin/bash

# Check if the user provided the folder path
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_folder>"
    exit 1
fi

# Assign the first argument to folder_path
folder_path="$1"

# Run the commands with predefined subfolder names and match strings
./copy_matching_files.sh "$folder_path" sim_output/angular_tuning_all angularTuning3_spikes_all
./copy_matching_files.sh "$folder_path" sim_output/angular_tuning_ON angularTuning3_spikes_ON
./copy_matching_files.sh "$folder_path" sim_output/spike_hist spikeHist__data.json
./copy_matching_files.sh "$folder_path" sim_output/raster raster.png
./copy_matching_files.sh "$folder_path" sim_output/rasterCustom _rasterCustom_1.png
./copy_matching_files.sh "$folder_path" sim_output/spikes spikeHist__spikes.json
./copy_matching_files.sh "$folder_path" sim_output/mean_voltage/data _meanVoltage__data.json
./copy_matching_files.sh "$folder_path" sim_output/mean_voltage/figs_VPM mean_traces__VPM__pop.png
./copy_matching_files.sh "$folder_path" sim_output/mean_voltage/figs_TRN traces__TRN__pop.png
./copy_matching_files.sh "$folder_path" sim_output/spikePSD spikePSD.png
./copy_matching_files.sh "$folder_path" sim_output/spikePSD_custom/figs _spikePSDPlot_
./copy_matching_files.sh "$folder_path" sim_output/spikePSD_custom/hist _spikePSDData_
./copy_matching_files.sh "$folder_path" sim_output/spikePSD_custom/data ratePSD_values.json
./copy_matching_files.sh "$folder_path" sim_output/lfpPSD/data _lfpPSD_values.json
./copy_matching_files.sh "$folder_path" sim_output/lfpPSD_windows/data _lfpPSD_values_windows.json
./copy_matching_files.sh "$folder_path" sim_output/lfpTraces/data _LFP_data.json
./copy_matching_files.sh "$folder_path" sim_output/lfpTraces/figs _LFP_timeSeries_
./copy_matching_files.sh "$folder_path" sim_output/spectrogram/figs _LFPSpectrogram.png
./copy_matching_files.sh "$folder_path" sim_output/spectrogram/data_windows _spectrogram_values_windows.pkl

echo "All copy operations completed."



import matplotlib.pyplot as plt
import numpy as np

class PlotSpikingHist():
    def plot_hist(pop_name, netParams, spk_window, pop_color='k'):
        print('Plotting Raster and Histogram for pop ', pop_name)

        cellsList = netParams.popParams[pop_name]['cellsList']

        cells_spikes=[]
        for cell_ind,cell in enumerate(cellsList):
            if spk_window is not None:
                cell_spikes = []
                for spkt in cellsList[cell_ind]['spkTimes']:
                    if (spkt>=spk_window[0]) and (spkt<spk_window[1]):
                        cell_spikes.append(spkt)
            else: cell_spikes = cellsList[cell_ind]['spkTimes']
            cells_spikes.append(cell_spikes)

        cells_spikes_hist = []
        for cell_spikes in cells_spikes: cells_spikes_hist.extend(cell_spikes)
        cells_spikes_hist.sort()

        # plt.figure(figsize=(20,40))
        plt.figure(figsize=(70,20))
        plt.rcParams.update({'font.size': 40})
        # plt.subplot(2,1,1)
        for cell_spikes_ind, cell_spikes in enumerate(cells_spikes):
            plt.scatter(cell_spikes,[cell_spikes_ind for i in range(len(cell_spikes))],c=pop_color)
        # plt.ylim([-5,830])
        # plt.xticks([8000,8250],[0,250])

        # if 'MLe' in pop_name:   plt.yticks([0,90,180])
        # else:                   plt.yticks([0,409,818])

        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        ax1.invert_yaxis()

        plt.savefig('../paper_figs/mle_ct_spiking/raster_'+pop_name+'.png')

        # plt.figure(figsize=(40,20))
        # plt.rcParams.update({'font.size': 40})
        # # plt.subplot(2,1,1)
        # for cell_spikes_ind, cell_spikes in enumerate(cells_spikes):
        #     plt.scatter(cell_spikes,[cell_spikes_ind for i in range(len(cell_spikes))],c=pop_color)
        # # plt.ylim([-5,830])
        # # plt.xticks([8000,8250],[0,250])

        # # if 'MLe' in pop_name:   plt.yticks([0,90,180])
        # # else:                   plt.yticks([0,409,818])

        # ax1=plt.gca()
        # ax1.spines['top'].set_visible(False)
        # ax1.spines['right'].set_visible(False)
        # ax1.spines['bottom'].set_visible(False)
        # ax1.spines['left'].set_visible(False)

        # ax1.invert_yaxis()

        # plt.savefig('../paper_figs/mle_ct_spiking/raster_'+pop_name+'_fullRaster.png')
        
        if spk_window is not None:
            # plt.figure(figsize=(20,10))
            plt.figure(figsize=(70,20))
            plt.rcParams.update({'font.size': 40})
            # window_size=1
            window_size=2.5
            # plt.subplot(2,1,2)
            plt.hist(cells_spikes_hist,bins=np.arange(spk_window[0],spk_window[1]+window_size,window_size), color=pop_color)
            plt.ylim([0,150])
            plt.xticks([8000,8250],[0,250])

            ax1=plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)

            plt.savefig('../paper_figs/mle_ct_spiking/spkHist_'+pop_name+'.png')

        # --- Plots colored raster for MLe_pop
        if 'MLe' in pop_name:
            # plt.figure(figsize=(20,40))
            plt.figure(figsize=(70,20))
            plt.rcParams.update({'font.size': 40})
            cell_plot_colors = ['blue','darkviolet','mediumvioletred','red','gold','darkkhaki','forestgreen','dodgerblue']
            cell_color_dict={}
            for i in np.arange(0,360,2):
            # for i in range(len(cells_spikes)):
                if i>315+22.5: cell_color_dict.update({int(i/2):'blue'})
                else:
                    for bin_i, bin in enumerate([0,45,90,135,180,225,270,315]):
                        if (bin-22.5<i<bin+22.5): cell_color_dict.update({int(i/2):cell_plot_colors[bin_i]})

            # plt.subplot(2,1,1)
            for cell_spikes_ind, cell_spikes in enumerate(cells_spikes):
                plt.scatter(cell_spikes,[cell_spikes_ind for i in range(len(cell_spikes))],c=cell_color_dict[cell_spikes_ind])
            plt.ylim([-5,830])
            plt.xticks([8000,8250],[0,250])

            if 'MLe' in pop_name:   plt.yticks([0,90,180])
            else:                   plt.yticks([0,409,818])

            ax1=plt.gca()
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)

            ax1.invert_yaxis()

            plt.savefig('../paper_figs/mle_ct_spiking/paperFigure_raster_'+pop_name+'.png')
        

# if __name__ == "__main__":


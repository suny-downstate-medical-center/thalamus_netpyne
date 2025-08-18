import matplotlib.pyplot as plt
import numpy as np

class MembranePotentialHeatmap:
    def __init__(self, data):
        """
        Initialize the class with the dataset.
        
        :param data: List of lists with float values representing the membrane potential of neurons over time.
        """
        self.data = np.array(data)
    
    def plot_heatmap(self, cmap='viridis', vmin=None, vmax=None, xlabel='Time', ylabel='Neuron Index', title='Resting Membrane Potential Heatmap', savefig=False, figsize=(10, 8)):
        """
        Plot the heatmap using the provided dataset.
        
        :param cmap: Colormap for the heatmap (default is 'viridis').
        :param vmin: Minimum value for the color scale (default is None, which uses the min value from the data).
        :param vmax: Maximum value for the color scale (default is None, which uses the max value from the data).
        :param xlabel: Label for the x-axis (default is 'Time').
        :param ylabel: Label for the y-axis (default is 'Neuron Index').
        :param title: Title of the heatmap (default is 'Resting Membrane Potential Heatmap').
        """
        plt.figure(figsize=figsize)
        plt.imshow(self.data, aspect='auto', cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)

        # Adding color bar
        plt.colorbar(label='Membrane Potential (mV)')

        # Labeling the axes
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Adding ticks and labels if necessary
        plt.xticks(ticks=np.arange(self.data.shape[1]))
        plt.yticks(ticks=np.arange(self.data.shape[0]))
        plt.xlim([0,self.data.shape[1]])
        # plt.xticks(ticks=np.arange(self.data.shape[1]), labels=[f'Time {i+1}' for i in range(self.data.shape[1])])
        # plt.yticks(ticks=np.arange(self.data.shape[0]), labels=[f'Neuron {i+1}' for i in range(self.data.shape[0])])

        # Display the heatmap
        plt.title(title)
        
        if savefig:
            if savefig==True:
                plt.savefig('heatmap_fig.png')
            elif type(savefig)==str:
                if ('.png' not in savefig) or ('.jpg' not in savefig): savefig+='.png'
                plt.savefig(savefig)

# Example usage
if __name__ == "__main__":
    # Example dataset: List of lists with float values
    data = [
        [ -70.5, -71.0, -70.2, -69.8, -70.0 ],
        [ -65.2, -66.0, -65.8, -66.1, -65.5 ],
        [ -60.3, -61.2, -60.8, -61.0, -60.5 ],
        # Add more data as needed
    ]
    
    heatmap = MembranePotentialHeatmap(data)
    # You can specify the colormap and vmin/vmax values as needed
    heatmap.plot_heatmap(cmap='inferno', vmin=-75, vmax=-55, savefig=False, figsize=(10, 8))

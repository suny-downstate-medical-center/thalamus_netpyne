import random
from neuron import h
import json

class NeuronalSpikingGenerator:
    def __init__(self, num_cells, interval, seed=100000, down_intervals=[],frequency=10,dt=0.025,frequency_std=None):
        self.num_cells = num_cells
        self.interval = interval
        self.seed = seed
        self.frequency=frequency
        self.down_intervals = down_intervals
        self.spiking_dict = {}
        self.random = h.Random(seed)
        self.dt = dt
        self.frequency_std = frequency_std

    def generate_spikes(self):
        for cell in range(self.num_cells):
            print('Generating spikes for cell #: ', cell, '\t- ', self.num_cells-cell, ' cells left')
            spiking_times = []
            t = self.interval[0]
            while t <= self.interval[1]:
                if not self.is_down_interval(t):
                    spike = self.generate_spike(t)
                    if spike:
                        spiking_times.append(spike)
                t += 0.025  # Adjusted for 0.025 ms time step
            spiking_times.sort()  # Sort spike times
            self.spiking_dict[cell] = spiking_times

    def generate_spike(self, t):
        # Use NEURON's random number generator
        spike_chance = self.random.uniform(0, 1)
        if spike_chance < self.calc_spike_probability():
            return t
        else:
            return None

    def calc_spike_probability(self):
        if self.frequency_std is not None:
            # Sample from a normal distribution with mean 'frequency' and std 'frequency_std'
            # frequency_sample = random.gauss(self.frequency, self.frequency_std)
            frequency_sample = self.random.normal(self.frequency, self.frequency_std**2)
            # Ensure the frequency is non-negative
            frequency_sample = max(frequency_sample, 0)
        else:
            frequency_sample = self.frequency
        # Calculate spike probability for the given dt
        return self.dt * frequency_sample / 1000

    # def calc_spike_probability(self):
    #     # Calculate spike probability for 0.025 ms time step
    #     return self.dt * self.frequency / 1000

    def is_down_interval(self, t):
        for interval in self.down_intervals:
            if t >= interval[0] and t <= interval[1]:
                return True
        return False

    def get_spiking_dict(self):
        return self.spiking_dict

# Example usage:
num_cells       = 356
interval        = (0, 10000)    # ms
dt              = 0.025         # ms
seed            = 100000        # Seed for random number generator
# down_intervals  = []  # Example of intervals with no spikes
down_intervals  = [(1000, 1500), (2000, 2500), (3000, 3500)]  # Example of intervals with no spikes

# freq = (1.5, 1.7) # (Dash, 2022)
# freq = (0.000001, None) # Very low frequency to simulate no spikes
# freq = (  4, None) # (Wake CT background from Iavarone, 2023)
# freq = ( 10, None) # (Sleep ML input from Iavarone, 2023)

# num_cells = 356; freq = (0.226, 0.55)       # (Dash, 2022) CT sparse
# # num_cells = 235; freq = (0.000001, None)    # (Dash, 2022) CT silent
# # num_cells = 92;  freq = (6.3, 5.2)          # (Dash, 2022) CT suppressed
# # num_cells = 135; freq = (1.5, 1.7)          # (Dash, 2022) CT activated

# num_cells = 356; freq = (0.226, None)       # (Dash, 2022) CT sparse
# num_cells = 235; freq = (0.000001, None)    # (Dash, 2022) CT silent
# num_cells = 92;  freq = (6.3, None)         # (Dash, 2022) CT suppressed
# num_cells = 135; freq = (1.5, None)         # (Dash, 2022) CT activated

#########################################################################################################
# --- Tentative rescaling of cells based on no silent cells in L6A
# num_cells = 359; freq = (0.226, 0.55)       # (Dash, 2022) CT sparse
# num_cells = 187; freq = (6.3, 5.2)          # (Dash, 2022) CT suppressed
# num_cells = 272; freq = (1.5, 1.7)          # (Dash, 2022) CT activated

# num_cells = 359; freq = (0.226, None)       # (Dash, 2022) CT sparse
# num_cells = 187; freq = (6.3, None)         # (Dash, 2022) CT suppressed
# num_cells = 272; freq = (1.5, None)         # (Dash, 2022) CT activated
#########################################################################################################
# num_cells = 818; freq = (5, None)          # (Dash, 2022) CT activated
# num_cells = 359; freq = (5, None)          # (Dash, 2022) CT sparse
# num_cells = 187; freq = (5, None)          # (Dash, 2022) CT suppressed
num_cells = 272; freq = (5, None)          # (Dash, 2022) CT activated
#########################################################################################################

# (Hirai, 2018) - L6A deflection response latency - (n = 30, 13.1 ± 2.5 ms)

# freq = (0.5, None)
# freq = (  1, None)
# freq = (  2, None)
# freq = (  5, None)
# freq = ( 20, None)

spiking_generator = NeuronalSpikingGenerator(num_cells, interval, seed, down_intervals,freq[0],dt,freq[1])
spiking_generator.generate_spikes()
spiking_dict = spiking_generator.get_spiking_dict()

print(spiking_dict)

if len(down_intervals)>0:
    down_intervals_str_ = ['_'+str(i1)+'|'+str(i2) for (i1,i2) in down_intervals]
    down_intervals_str = ''.join(down_intervals_str_)
else:
    down_intervals_str='_None'

if freq[1] is not None:     freq_str = str(freq[0])+'@'+str(freq[1])
else:                       freq_str = str(freq[0])

filename = '../stims/Artificial/CT/spiking_dict__nCells_'+str(num_cells)+'__sampling_'+str(int(dt*1000))+'_uS__freq_'+freq_str+'_Hz__downTimes_'+down_intervals_str+'__'

print('\t\tSaving artificial spikes in: \t', filename+'.json')
with open(filename+'.json', "w") as json_file: json.dump(spiking_dict, json_file)
    

import matplotlib.pyplot as plt
plt.figure(figsize=(20, 10))
for k in spiking_dict.keys():
    plt.scatter(spiking_dict[k], [int(k) for i in range(len(spiking_dict[k]))], s=0.5)
# plt.xlim([0, 2000])
plt.savefig(filename+'.png', dpi=500)

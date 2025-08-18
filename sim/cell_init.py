"""
init.py

Starting script to run NetPyNE-based thalamus model for thesis.

Usage:
    python cell_init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 8 nrniv -python -mpi cell_init.py

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""


# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
from    matplotlib  import  pyplot  as plt
# print("Matplotlib backend (default): %s" %plt.get_backend())
# modules = []
# for module in sys.modules:
#     if module.startswith('matplotlib'):
#         modules.append(module)
# for module in modules:
#     sys.modules.pop(module)
# import matplotlib
# matplotlib.use("MacOSX")
# from    matplotlib  import  pyplot  as plt
# print("Matplotlib backend (dynamic): %s" %plt.get_backend())

import matplotlib; matplotlib.use('agg')  # to avoid graphics error in servers
from netpyne import sim

############################################################################################################
# --- Running simulation
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cell_cfg.py', netParamsDefault='cell_netParams.py')
createOnly          = False
addNoise            = True

############################################################################################################
# --- Custom package imports
from CurrentStim            import CurrentStim          as CS

############################################################################################################
# --- Creating network using NetPyNE high-level specifications
sim.initialize(
    simConfig = cfg,     
    netParams = netParams)                  # create network object and set cfg and net params
sim.net.createPops()                        # instantiate network populations
sim.net.createCells()                       # instantiate network cells based on defined populations
sim.net.connectCells()                      # create connections between cells based on params
sim.net.addStims()                          # add network stimulation

if addNoise and sim.cfg.addNoiseIClamp:     sim, vecs_dict = CS.addNoiseIClamp(sim)

############################################################################################################
# --- Running the network
if createOnly: print('\n\n\n ---- Running CREATE ONLY mode to inspect for bugs during network creation --- \n\n\n')
else:
    sim.setupRecording()                    # setup variables to record for each cell (spikes, V traces, etc)
    sim.runSim()                            # run parallel Neuron simulation  
    sim.gatherData()                        # gather spiking data and cell info from each node
    sim.analyze()
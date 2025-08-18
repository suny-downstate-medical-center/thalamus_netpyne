# Computational Model of the Whisker Thalamocortical Pathway

## Description

Detailed computational model of the whisker thalamocortical pathway of the mouse developed in NetPyNE ([www.netpyne.org](https://www.netpyne.org/)).

The model is described in the following publication:

- Moreira JVS, Borges FS, Atherton Z, Crandall S, Varela C, Dura-Bernal S (2025) TITLE. JOURNAL.


## Overview

This repository contains all the source code and data required to generate the paper results. Each subfolder has a separate README.md with instructions:

- `/sim`: Python source code of the thalamocortical circuit model. Implemented in NetPyNE/NEURON.
- `/analysis`: Python source code to perform the analyses used to generate the paper results.

## System requirements
The computational model has the following software dependencies:

- [Python](http://python.org) (v 3.9.7) - Programming language
- [NetPyNE](https://www.netpyne.org/) (v 1.0.7) - Framework used for model development
- [NEURON](https://neuron.yale.edu) (v 7.8.2) - Simulation engine
- [MPI](https://pypi.org/project/mpi4py-mpich/) - Package for supporting parallel simulations
- [snap / bluepysnap](https://github.com/BlueBrain/snap) (v 0.13.0) - Package for loading Blue Brain circuit properties
- [SciPy](https://scipy.org) (v 1.10.1) - Package for running statistical tests

Follow this [link](https://neuron.yale.edu/ftp/neuron/2019umn/neuron-quickstart.pdf) for a quickstart with NEURON/MPI installation.

Additionaly, the following packages are required to run the full model locally (not recommended) or in an HPC environment:  

- bluepysnap (v 0.13.0)  
- wheel (v 0.38.4)
- numpy (v 1.23.5)
- matplotlib (v 3.7.1)  
- future (v 0.18.3)  
- pandas (v 1.4.1)  
- scipy (v 1.10.1)  
- matplotlib_scalebar (v 0.8.0)  
- bokeh (v 2.4.2)  
- schema (v 0.7.5)  
- libsonata (v 0.1.14)  
- mpi4py-mpich   
- h5py  

For a detailed description of tin installation process, see the file `README_Installation.md`

The computational model and analysis code has been tested on Mac OS 15.3.2 (Sequoia) and Rocky Linux 8.7 (SUNY Downstate HPC) and 8.9 (Expanse HPC).

## Instructions for use
Instructions on how to install and run the full thalamocortical model, and analyze the model output can be found in the  
`README_installation.md`, `sim/README_sim.md`, `sim_analysis/README_analysis.md` files.


For further information please contact: salvador.dura-bernal@downstate.edu and/or joao.moreira@downstate.edu.
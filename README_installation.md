
# Installation (Computational Model of the Whisker Thalamocortical Pathway)

## Description
Detailed computational model of the whisker thalamocortical pathway of the mouse developed in NetPyNE ([www.netpyne.org](https://www.netpyne.org/)).

The model is described in the following publication:

- Moreira JVS, Borges FS, Atherton Z, Crandall S, Varela C, Dura-Bernal S (2025) TITLE. JOURNAL.

## Overview
The following is a detailed description how to install the necessary software to run the "Detailed computational model of the whisker thalamocortical pathway of the mouse". 

The computational model and analysis code has been tested on Mac OS 15.3.2 - Sequoia (local machine) and Rocky Linux 8.7 (SUNY Downstate HPC) and 8.9 (Expanse HPC).

## Setting up a Conda Environment for NetPyNE

We recommend you install Anaconda toolkit as your package manager.  
Download the miniconda installer:  
- `curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh`  

Run the miniconda installer:  
- `sh ./miniconda.sh -b`  
- `source miniconda3/bin/activate`  
- `conda init`  

Activate miniconda
 - `source ~/.bashrc`


Setting up a new environment using the Anaconda distribution:  
- `conda create -n py397_neuron782 python=3.9.7`  
Move to your python environment:
- `conda activate py397_neuron782`  
Add a source line to your .bashrc so that it starts in `py397_neuron782`  
- `echo -e '\nconda activate py397_neuron782\n' >> ~/.bashrc`  
Creating a virtual environment prevents you from messing up the OS in case something goes wrong during the installation.  
You can always delete the `py397_neuron782` virtual environment and start over.

## Installing MPI, NEURON and NetPyNE

Install MPI:
- `pip install mpi4py-mpich`

Install the NEURON simulator:
- `pip install neuron`  

Install NetPyNE via one of the following options:
- Install the latest released version of NetPyNE via pip (Recommended)  
Linux or MacOS: 
    - `pip install netpyne`  
    or
    - `python -m pip install netpyne`


- Install the development version of NetPyNE via GitHub and pip:  
The NetPyNE package source files, as well as example models, are available via GitHub at: https://github.com/Neurosim-lab/netpyne. The following instructions will install the version in the GitHub “development” branch – it will include some of the latest enhancements and bug fixes, but could also include temporary bugs:
    - `git clone https://github.com/Neurosim-lab/netpyne.git`
    - `cd netpyne`
    - `git checkout development`
    - `pip install -e .`

## Installing remaining packages

All remaining packages can be installed using the following command on the terminal:

- `pip install wheel==0.38.4 numpy==1.23.5 neuron==7.8.2 matplotlib==3.7.1 future==0.18.3 pandas==1.4.1 scipy==1.10.1 matplotlib_scalebar==0.8.0 bokeh==2.4.2 schema==0.7.5 bluepysnap==0.13.0 libsonata==0.1.14 h5py`


## (Extra) Debugging MPI installation on HPCs
The Following steps are required to make MPI work on the SUNY Downstate and Expanse HPCs, and might not be necessary in other cases. If you encounter issues with MPI when running your model, give this a try.

Add the MPI binary path to your `~/.bashrc` file:
- find where your mpich library is:
    - check ~/miniconda3/envs/<ENV>/lib/<PYTHON>/site-packages/
    - `th` (Make sure NOT to proceed) and look for the libmpi.so in the prompt it provides, it will look something like this, then exit the pip interface (do NOT uninstall)  
    `/<PATH_TO_MINICONDA>/envs/dev/lib/python3.9.7/site-packages/mpi4py_mpich.libs/libmpi-1e85b3e0.so.12.1.8`  
    (e.g.: `/HOME_DIR/USERNAME/miniconda3/envs/dev/lib/python3.9.7/site-packages/mpi4py_mpich.libs/libmpi-1e85b3e0.so.12.1.8`)  

- set the LD_LIBRARY_PATH environment to the parent directory <DIRPATH> of your libmpi.so  
`echo -e '\nexport LD_LIBRARY_PATH="<DIRPATH>"\n' >> ~/.bashrc`

At the end of this, the last 5 lines of your ~/.bashrc should look something like this:

    conda activate py397_neuron782
    
    
    
    export LD_LIBRARY_PATH="/HOME_DIR/USERNAME/miniconda3/envs/dev/lib/python3.9.7/site-packages/mpi4py_mpich.libs"


- since nrniv is looking specifically for a file named `libmpi.so`, create a symbolic link at this directory path after re-sourcing your `~/.bashrc`  
    
    `source ~/.bashrc`  
    `cd $LD_LIBRARY_PATH`  
    `ln -s libmpi-1e... libmpi.so`  

    or go to the path of the `libmpi` file (e.g.: `/HOME_DIR/USERNAME/CONDA_VERSION/envs/dev/lib/PYTHON_VERSION/site-packages/mpi4py_mpich.libs/`) and make a copy of the file with `libmpi` in its name, change the name of the copied file to `libmpi.so`.

---

## Verifying installation and testing
- To verify the installation of NetPyNE, you can run the following code on the terminal.  
    - Innitialize the Python interface:  
    run the command `python` or `python3` in the terminal.
    - In the Python terminal:  
    (1) `import netpyne`  
    (2) `netpyne.__file__`  
    expected output: `'/PATH_TO_YOUR_NETPYNE_INSTALLATION/netpyne/netpyne/__init__.py'`  
    (3) `netpyne.__version__`  
    expected output: `'1.0.7'`  

- To verify if the NEURON simulator is working together with NetPyNE, you can run the following code on the terminal.  
    - Innitialize the Python interface:  
    run the command `python` or `python3` in the terminal.
    - In the Python terminal:  
    (1) `import netpyne`  
    (2) `from neuron import h,gui`  
    expected output: a small window should pop up, which is the NEURON gui  
    (3) type `h` in the terminal window  
    expected output: The following should be printed in the terminal `<TopLevelHocInterpreter>`

- To verify the MPI installation, you can run the following command:  
    `mpiexec -n NPPROCS hostname`
    replacing `NPPROCS` by the number of cores available in your machine

- Finnally, to test if the sotware are working to run the model, you can set the varable `createOnly = True` in the script `barreloid_init.py`, and run:
    (1) `python barreloid_init.py`  
    This will run the creation of the network, but stop before the network simulation, which is the step that requires multiple cores and high RAM memory.  
    At this point, you should be able to inspect the network properties, including cell morphology, biophisics, connectivity, stimulation, and recording parameters.

## Running the model
- If you follow these steps, you should be able to clone and run the computational model in your local machine (not recommended), or in an HPC environment.  
- We do not recommend attempting to run the full scale model locally, because it is very computationally expensive. A single run of the model using a integration step `dt = 0.025 ms` for `16 seconds` requires on average 2h35 minutes using 64 cores, and 128Gb of RAM memory.
- For a detailed description of how to setup and run the model, see `/sim/README_sim.md`.
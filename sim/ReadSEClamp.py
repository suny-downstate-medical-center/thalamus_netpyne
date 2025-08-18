from netpyne import sim
import json
import os
import sys
import pickle

class ReadSEClamp():

    def loadSimOutput(filename):
        if      '.pkl'  in filename:    return ReadSEClamp.loadPickle(filename)
        elif    '.json' in filename:    return ReadSEClamp.loadJSON(filename)
        else:                           return ReadSEClamp.loadSimObj(filename)

    def loadJSON(filename):
        print('\t>> Loading JSON file \n')
        with open(filename, 'r') as fileObj: simFile=fileObj.read()
        sim = json.loads(simFile)
        return sim

    def loadPickle(filename):
        print('\t>> Loading Pickle file \n')
        try:
            with open(filename, 'rb') as file:
                sim = pickle.load(file)
            return sim
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
            return None
        except pickle.UnpicklingError as e:
            print(f"Error: Failed to unpickle file '{filename}': {e}")
            return None

    def loadSimObj(filename=None):
        print('\t>> Loading SIM object \n')
        from netpyne import sim
        sim.load(filename)
        return sim
    
    def parseSim(sim):
        try:        
            simData     = sim['simData']
            simCfg   = sim['simConfig']
            simNet      = sim['net']
            print('\t>> Dictionary data parsed')
        except: 
                try:    
                    simData     = sim.simData
                    simCfg      = sim.cfg
                    simNet      = sim.net
                    print('\t>> Object data parsed')
                except: print('\t>> Parsing failed')
        
        return simData, simCfg, simNet
        
    def getSEClampIVResponse(sim):
        simData, simCfg, simNet = ReadSEClamp.parseSim(sim)

        try:    cellList = sim['net']['cells']
        except: cellList = sim.net.cells
        print('\t>> Updating currents from SEClamp protocol')
        store_final_current = {}
        for cell_ind in range(len(cellList)):
            try:    cell_gid = sim['net']['cells'][cell_ind]['gid']
            except: cell_gid = sim.net.cells[cell_ind].gid

            try:
                i_vals = list(simData['I_soma']['cell_'+str(cell_ind)])
                i_ = i_vals[-1]
                v_vals = list(simData['V_soma']['cell_'+str(cell_ind)])
                v_ = v_vals[-1]
            except:
                i_=0
                v_=None
            store_final_current.update({str(cell_gid):{'i':i_,'V':v_}})

        print(store_final_current)

        # Save the dictionary in JSON format
        with open('final_current_data/final_current_data.json', 'w') as json_file: json.dump(store_final_current, json_file)


if __name__ == "__main__":
    
    base_dir        = os.path.expanduser("~")
    folder_path     = base_dir+'/Research/Models/BBP/thalamus_netpyne/data/barreloid_sims/barreloid_network/'
    if len(sys.argv) < 2:
        # sys.exit(1)  # Exit with error
        filename        = folder_path+'barr_net_sim_000166_270_deg_simplifiedGids_USE_0.4029343148532312_data.pkl'
        print("Reading default file ", filename)
    else:
        filename        = folder_path+sys.argv[1]
        print("Reading file ", filename)

    # filename        = folder_path+'barr_net_sim_000166_270_deg_simplifiedGids_USE_0.4029343148532312_data.pkl'
    # filename        = folder_path+'barr_net_sim_000166_270_deg_simplifiedGids_USE_0.4029343148532312_data.json'
    
    save_file_name  = 'pooled_final_current_data.json'

    sim = ReadSEClamp.loadSimOutput(filename=filename)
    # sim = ReadSEClamp.loadJSON(filename=filename)
    # sim = ReadSEClamp.loadPickle(filename=filename)

    ReadSEClamp.getSEClampIVResponse(sim)
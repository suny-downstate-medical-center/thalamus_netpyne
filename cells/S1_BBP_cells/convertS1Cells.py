#####################################################################################################################
#   Code to convert cell morphologyes from NetPyNE-BBP S1 to thalamus model
#####################################################################################################################

import os
import json

from netpyne import specs
netParams = specs.NetParams()

#####################################################################################################################

# --- Recursive function to replace the keys of a dictionary at any nested level
def replace_keys(old_dict, key_dict):
    new_dict = { }
    for key in old_dict.keys():
        new_key = key_dict.get(key, key)
        if isinstance(old_dict[key], dict): new_dict[new_key] = replace_keys(old_dict[key], key_dict)
        else:                               new_dict[new_key] = old_dict[key]
    return new_dict

#####################################################################################################################

def dict_replace_value(d: dict, old: str, new: str) -> dict:
    x = {}
    for k, v in d.items():
        if   isinstance(v, dict):   v = dict_replace_value(v, old, new)
        elif isinstance(v, list):   v = list_replace_value(v, old, new)
        elif isinstance(v, str):    v = v.replace(old, new)
        x[k] = v
    return x

def list_replace_value(l: list, old: str, new: str) -> list:
    x = []
    for e in l:
        if   isinstance(e, list):   e = list_replace_value(e, old, new)
        elif isinstance(e, dict):   e = dict_replace_value(e, old, new)
        elif isinstance(e, str):    e = e.replace(old, new)
        x.append(e)
    return x

#####################################################################################################################

# files = os.listdir()
source_dir='./source_cells/'
files = os.listdir(source_dir)
fileNames = [file for file in files if 'L6_TPC_L4_cADpyr231' in file]

keys=['Im','CaDynamics_E2','SK_E2','Ca_HVA','Ih','NaTs2_t','NaTa_t','K_Pst','K_Tst','Nap_Et2','Ca_LVAst','SKv3_1'] # native NEURON mechanisms should not be changed (e.g. 'pas')

replace_dict={}
for k in keys:replace_dict.update({k:'S1_'+k})
replace_dict.update({'soma':'soma_0'})

print('\t>>>\t Converting cells')
for ind,fileName in enumerate(fileNames):
    print('\t\t',fileName.split('.')[0])
    with open(source_dir+fileName, 'r') as openfile: cell_dict_ = json.load(openfile)

    cell_dict = cell_dict_.copy()
    # cell_dict = replace_keys(cell_dict, {'soma':'soma_0'})
    cell_dict = dict_replace_value(cell_dict,'L6_TPC_L4_cAD','L6A__cell')
    cell_dict = dict_replace_value(cell_dict,'soma','soma_0')
    
    for key in keys: cell_dict = replace_keys(cell_dict,replace_dict)
    # for key in keys: cell_dict = dict_replace_value(cell_dict,key,'S1_'+key)

    try:    cell_index=fileName.split('.')[0].split('_')[-2]
    except: cell_index=1000+ind
    json_object = json.dumps(cell_dict, indent=4)

    # --- Adds 3 zeros to the beginning of the filename for easier indexing
    ind_str=str(ind)
    ind_str2=ind_str.zfill(3)

    # with open('L6A@'+str(ind_str2)+'__cell_cellParams.json', "w") as outfile: outfile.write(json_object)
    
    cellLabel    = 'L6A@'+str(ind_str2)+'__cell'
    
    netParams.cellParams[cellLabel]=cell_dict
    netParams.cellParams[cellLabel]['conds']={'cellType':cell_dict['conds']['cellType']}
    
    netParams.saveCellParams(   label    = cellLabel,
                                fileName = 'L6A@'+str(ind_str2)+'__cell_cellParams.json')
    
    # --- Adding a .pkl format cell
    netParams.saveCellParams(   label    = cellLabel,
                                fileName = 'L6A@'+str(ind_str2)+'__cell_cellParams.pkl')
    
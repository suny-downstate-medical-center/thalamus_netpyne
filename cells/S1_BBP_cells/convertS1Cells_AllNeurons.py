#####################################################################################################################
#   Code to convert cell morphologyes from NetPyNE-BBP S1 to thalamus model
#####################################################################################################################

import os
import json
from netpyne import specs

########################################################################################################################################################################
class ConvertS1Cell():

    # --- Recursive function to replace the keys of a dictionary at any nested level
    def replace_keys(old_dict, key_dict):
        new_dict = { }
        for key in old_dict.keys():
            new_key = key_dict.get(key, key)
            if isinstance(old_dict[key], dict): new_dict[new_key] = ConvertS1Cell.replace_keys(old_dict[key], key_dict)
            else:                               new_dict[new_key] = old_dict[key]
        return new_dict

    ########################################################################################################################################################################

    def dict_replace_value(d: dict, old: str, new: str) -> dict:
        x = {}
        for k, v in d.items():
            if   isinstance(v, dict):   v = ConvertS1Cell.dict_replace_value(v, old, new)
            elif isinstance(v, list):   v = ConvertS1Cell.list_replace_value(v, old, new)
            elif isinstance(v, str):    v = v.replace(old, new)
            x[k] = v
        return x

    def list_replace_value(l: list, old: str, new: str) -> list:
        x = []
        for e in l:
            if   isinstance(e, list):   e = ConvertS1Cell.list_replace_value(e, old, new)
            elif isinstance(e, dict):   e = ConvertS1Cell.dict_replace_value(e, old, new)
            elif isinstance(e, str):    e = e.replace(old, new)
            x.append(e)
        return x

    def convertS1Cell(source_dir, fileNames, S1_cell_dict):
        
        netParams = specs.NetParams()

        replace_dict={}
        for k in S1_cell_dict['mods']:replace_dict.update({k:S1_cell_dict['mod_label']+k})
        replace_dict.update({'soma':'soma_0'})

        for ind,fileName in enumerate(fileNames):
            print('\t\t',fileName.split('.')[0])
            with open(source_dir+fileName, 'r') as openfile: cell_dict_ = json.load(openfile)

            cell_dict = cell_dict_.copy()
        # cell_dict = replace_keys(cell_dict, {'soma':'soma_0'})
            cell_dict = ConvertS1Cell.dict_replace_value(cell_dict,S1_cell_dict['cellName'],S1_cell_dict['newLabel']+'__cell')
            cell_dict = ConvertS1Cell.dict_replace_value(cell_dict,'soma','soma_0')
        
            for key in S1_cell_dict['mods']: cell_dict = ConvertS1Cell.replace_keys(cell_dict,replace_dict)
        # for key in keys: cell_dict = dict_replace_value(cell_dict,key,'S1_'+key)

            try:    cell_index=fileName.split('.')[0].split('_')[-2]
            except: cell_index=1000+ind
            json_object = json.dumps(cell_dict, indent=4)

        # --- Adds 3 zeros to the beginning of the filename for easier indexing
            ind_str=str(ind)
            ind_str2=ind_str.zfill(3)

        # with open('L6A@'+str(ind_str2)+'__cell_cellParams.json', "w") as outfile: outfile.write(json_object)
        
            cellLabel    = S1_cell_dict['newLabel']+'@'+str(ind_str2)+'__cell'
        
            netParams.cellParams[cellLabel]=cell_dict
            netParams.cellParams[cellLabel]['conds']={'cellType':cell_dict['conds']['cellType']}
        
            netParams.saveCellParams(   label    = cellLabel,
                                    fileName = S1_cell_dict['newLabel']+'@'+str(ind_str2)+'__cell_cellParams.json')
        
        # --- Adding a .pkl format cell
            netParams.saveCellParams(   label    = cellLabel,
                                    fileName = S1_cell_dict['newLabel']+'@'+str(ind_str2)+'__cell_cellParams.pkl')

########################################################################################################################################################################

# --- Cells to be converted
            
'''
type    morph                   label           
CT
        L6_TPC_L4_cADpyr231     L6_TPC_L4 cADpyr231     666     666
CC
        L6_TPC_L1_cADpyr231     L6_TPC_L1 cADpyr231     150     150
        L6_UTPC_cADpyr231       L6_UTPC cADpyr231       150     150
        L6_BPC_cADpyr231        L6_BPC cADpyr231        250     250
        L6_IPC_cADpyr231        L6_IPC cADpyr231        300     300
INH
        L6_LBC_bAC217           L6_LBC bAC217           18      36
        L6_LBC_bNAC219          L6_LBC bNAC219          9       36
        L6_LBC_cNAC187          L6_LBC cNAC187          9       36
        L6_MC_bAC217            L6_MC bAC217            9       35
        L6_MC_bNAC219           L6_MC bNAC219           6       35
        L6_MC_cACint209         L6_MC cACint209         20      35
'''


S1_cells    = {
                ####################################################   TC NEURONS   #####################################################################
                'L6_TPC_L4_cADpyr231':{     # --- TC neuron
                    'newLabel': 'L6A',
                    'cellName': 'L6_TPC_L4_cAD',
                    'mods':     ['CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                ####################################################   CC NEURONS   #####################################################################
                'L6_TPC_L1_cADpyr231':{     # --- CC neuron
                    'newLabel': 'L6CC_TPC_L1_cAD',
                    'cellName': 'L6_TPC_L1_cAD',
                    'mods':     ['CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_UTPC_cADpyr231':{     # --- CC neuron
                    'newLabel': 'L6CC_UTPC_cAD',
                    'cellName': 'L6_UTPC_cAD',
                    'mods':     ['CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_BPC_cADpyr231':{     # --- CC neuron
                    'newLabel': 'L6CC_BPC_cAD',
                    'cellName': 'L6_BPC_cAD',
                    'mods':     ['CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_IPC_cADpyr231':{     # --- CC neuron
                    'newLabel': 'L6CC_IPC_cAD',
                    'cellName': 'L6_IPC_cAD',
                    'mods':     ['CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                ####################################################   IN NEURONS   #####################################################################
                'L6_LBC_bAC217':{          # --- Interneuron
                    'newLabel': 'L6IN_LBC_bAC',
                    'cellName': 'L6_LBC_bAC',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_LBC_bNAC219':{          # --- Interneuron
                    'newLabel': 'L6IN_LBC_bNA',
                    'cellName': 'L6_LBC_bNA',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_LBC_cNAC187':{          # --- Interneuron
                    'newLabel': 'L6IN_LBC_cNA',
                    'cellName': 'L6_LBC_cNA',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_MC_bAC217':{          # --- Interneuron
                    'newLabel': 'L6IN_MC_bAC',
                    'cellName': 'L6_MC_bAC',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_MC_bNAC219':{          # --- Interneuron
                    'newLabel': 'L6IN_MC_bNA',
                    'cellName': 'L6_MC_bNA',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                'L6_MC_cACint209':{          # --- Interneuron
                    'newLabel': 'L6IN_MC_cAC',
                    'cellName': 'L6_MC_cAC',
                    'mods':     ['Ca', 'CaDynamics_E2', 'Ca_LVAst', 'Ih', 'Im', 'K_Pst', 'K_Tst', 'NaTa_t', 'NaTs2_t', 'Nap_Et2', 'SK_E2', 'SKv3_1'],
                    'mod_label': 'S1_',
                    },
                ####################################################   OLD NEURONS  #####################################################################
                # 'L6_LBC_bNAC219':{          # --- Old version of a single type of Interneuron
                #     'newLabel': 'L6IN',
                #     'cellName': 'L6_LBC_bNA',
                #     'mods':     ['Ca','Im','CaDynamics_E2','SK_E2','Ca_HVA','Ih','NaTs2_t','NaTa_t','K_Pst','K_Tst','Nap_Et2','Ca_LVAst','SKv3_1'],
                #     'mod_label': 'S1_',
                #     },
                ####################################################      END       #####################################################################
                }

########################################################################################################################################################################
# files = os.listdir()
source_dir='./source_cells/'
files = os.listdir(source_dir)
########################################################################################################################################################################

for S1_cellType in S1_cells.keys():

    S1_cell_dict = S1_cells[S1_cellType]
    fileNames = [file for file in files if S1_cellType in file]

    print('\t>>>\t Converting cells: ', S1_cellType)
    ConvertS1Cell.convertS1Cell(source_dir, fileNames, S1_cell_dict)

########################################################################################################################################################################

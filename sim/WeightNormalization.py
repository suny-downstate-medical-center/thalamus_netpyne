import json
import numpy as np
import pickle

class WeightNormalization():
    def normalizeWeights(sim):
        print('Loading weight normalization dictionaries for each pathway')
        weightNorm={}
        for cellTemplate in sim.net.params.cellParams.keys():
            if ('VPL_TC' in cellTemplate) or ('Rt_RC' in cellTemplate):
                pop,template_id = cellTemplate.split('__')
                if template_id not in weightNorm.keys():weightNorm.update({template_id:{}})
                for synMech in sim.net.params.synMechParams.keys():
                    if 'syn|' not in synMech:continue
                    if (('VPL_TC' in cellTemplate) and synMech.split('|')[2]=='VPM') or (('Rt_RC' in cellTemplate) and synMech.split('|')[2]=='TRN'):
                        if synMech not in weightNorm[template_id].keys():weightNorm[template_id].update({synMech:{}})
                        if '|' in synMech:  synMech_join='_'.join(synMech.split('|'))
                        else:               synMech_join=synMech

                        wNorm_filename = '../cells/calibrated_weightNorm/calibrated_weightNorm__'+cellTemplate+'__'+synMech_join+'.json'
                        with open(wNorm_filename, "r") as json_file: weightNorm_dict=json.load(json_file)
                        weightNorm[template_id][synMech]={wn_keys:weightNorm_dict[wn_keys] for wn_keys in weightNorm_dict.keys()}
                        # weightNorm[template_id][synMech]={wn_keys:np.mean(weightNorm_dict[wn_keys]) for wn_keys in weightNorm_dict.keys()}

        # print('Normalizing conn weights')
        print('Normalizing conn weights - Code changed to normalize <gmax> for each synapse instead of <weight> of each connection, to allow efficient use of vector play on weight')
        for cell_ind,cell in enumerate(sim.net.cells):
            if 'cellModel' in sim.net.cells[cell_ind].tags.keys(): continue # skips vecstim pops because they don't receive inputs and have no morpohology
            cellTemplate = sim.net.cells[cell_ind].tags['label'][0]
            template_id = cellTemplate.split('__')[1]
            for conn_ind,conn in enumerate(sim.net.cells[cell_ind].conns):
                if 'syn|' in sim.net.cells[cell_ind].conns[conn_ind]['synMech']:
                    synMech = sim.net.cells[cell_ind].conns[conn_ind]['synMech']
                    post_sec = sim.net.cells[cell_ind].conns[conn_ind]['sec']
                    wNorm_vals=weightNorm[template_id][synMech][post_sec]
                    num_wNorm_vals = len(wNorm_vals)
                    if      num_wNorm_vals==1: wNorm = wNorm_vals[0]
                    elif    num_wNorm_vals<1:  wNorm = 1 # bug case
                    else:   # find closest value to loc
                        wNorm_nbins = num_wNorm_vals+1
                        wNorm_bins  = [i/num_wNorm_vals for i in range(wNorm_nbins)]
                        post_loc = sim.net.cells[cell_ind].conns[conn_ind]['loc']

                        for i in range(wNorm_nbins-1):
                            wNorm_bin=[wNorm_bins[i],wNorm_bins[i+1]]
                            # print(i, wNorm_bin)
                            if (post_loc>=wNorm_bins[i]) and (post_loc<wNorm_bins[i+1]): 
                                wNorm=wNorm_vals[i]
                            #     print(i, wNorm_bin, wNorm_vals[i] , '<------ selected weight norm for post_loc ', post_loc)
                            # else:
                            #     print(i, wNorm_bin, wNorm_vals[i])
                    
                    # # store_w=sim.net.cells[cell_ind].conns[conn_ind]['weight']
                    # sim.net.cells[cell_ind].conns[conn_ind]['weight']*=wNorm
                    # # store_w2=sim.net.cells[cell_ind].conns[conn_ind]['weight']
                    # # print(cell_ind, conn_ind, conn['label'], store_w,store_w2)
                    
                    # # test - 2024-07-30
                    # sim.net.cells[cell_ind].conns[conn_ind]['hObj'].weight[0]*=wNorm

                    store_gmax=sim.net.cells[cell_ind].conns[conn_ind]['hObj'].syn().gmax
                    sim.net.cells[cell_ind].conns[conn_ind]['hObj'].syn().gmax*=wNorm
                    store_gmax_=sim.net.cells[cell_ind].conns[conn_ind]['hObj'].syn().gmax
                    print(cell_ind,conn_ind,'\tbefore: ', store_gmax,'\tafter: ',store_gmax_, ' \t ratio: ', store_gmax_/store_gmax, '\twNorm: ',wNorm)


        return sim



        # return sim

    #     return weightNorm
    
    # def normalizeWeights(sim):



'''
pathway_cellNumber      = { 'VPL_TC': [36963],
                            'Rt_RC' : [31864]}
pathway_synMechs_dict   = { 'VPL_TC': ['syn|MLe|VPM|exc', 'syn|L6A|VPM|exc', 'syn|TRN|VPM|inh'],
                            'Rt_RC' : ['syn|VPM|TRN|exc', 'syn|L6A|TRN|exc', 'syn|TRN|TRN|inh']}

weightNorm={}
for pathway in pathway_cellNumber.keys():
    if pathway not in weightNorm.keys():weightNorm.update({pathway:{}})
    for cellNumber in pathway_cellNumber[pathway]:
        if cellNumber not in weightNorm[pathway].keys():weightNorm[pathway].update({cellNumber:{}})
        for synMech in pathway_synMechs_dict[pathway]:
            if synMech not in weightNorm[pathway][cellNumber].keys():weightNorm[pathway][cellNumber].update({synMech:{}})
            if '|' in synMech:  synMech_join='_'.join(synMech.split('|'))
            else:               synMech_join=synMech
            wNorm_filename = '../cells/calibrated_weightNorm/calibrated_weightNorm__'+pathway+'__'+str(cellNumber)+'__'+synMech_join+'.json'
            with open(wNorm_filename, "r") as json_file: weightNorm_dict=json.load(json_file)

            weightNorm[pathway][cellNumber][synMech]={wn_keys:np.mean(weightNorm_dict[wn_keys]) for wn_keys in weightNorm_dict.keys()}
'''
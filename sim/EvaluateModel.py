'''
Class to evaluate the barreloid thalamus model on NetPyNE

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''
##################################################################################################################################################################################################################
import numpy as np
class Evaluate():
    def getPopGids(sim):
        # --- Find pop names
        pop_gids={}
        for pop_name_ in sim.net.pops.keys():
            if '@' in pop_name_:    pop=pop_name_.split('__')[0].split('@')[0]
            else:                   pop=pop_name_.split('__')[0]
            # Dictionary of cell gids by pop
            if pop not in pop_gids.keys(): pop_gids.update({pop:[]})
            pop_gids[pop]+=sim.net.pops[pop_name_].cellGids
        return pop_gids
    
    def getPopConvergence(sim,verbose=False):
        # --- Pop parameters
        pop_gids  = Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Convergence to a population
        pop_convergence = {target_pop:{} for target_pop in pop_names}
        for target_pop in pop_names:
            for cell_gid in pop_gids[target_pop]:
                for conn in sim.net.cells[cell_gid].conns:
                    if conn['preGid']=='NetStim': continue # ignore NetStim conns
                    if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                    try: 
                        pre_pop = conn['label'].split('|')[1].split('@')[0]
                    except:
                        if verbose: print('skipping conn ',conn['label'])
                        continue
                    # print(cell_gid, pre_pop)
                    # print(pre_pop,target_pop)
                    if pre_pop not in pop_convergence[target_pop].keys():   pop_convergence[target_pop].update({pre_pop:1})
                    else:                                                   pop_convergence[target_pop][pre_pop]+=1
                # try:
                #     for conn in sim.net.cells[cell_gid].conns:
                #         if conn['preGid']=='NetStim': continue # ignore NetStim conns
                #         if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                #         pre_pop = conn['label'].split('|')[1].split('@')[0]
                #         print(cell_gid, pre_pop)
                #         if pre_pop not in pop_convergence[target_pop].keys():   pop_convergence[target_pop].update({pre_pop:1})
                #         else:                                                   pop_convergence[target_pop][pre_pop]+=1
                # except:
                #     print(conn)
                #     print('cell gid out of range: ', cell_gid)

        return pop_convergence

    def getPopDivergence(sim,verbose=False):
        # --- Pop parameters
        pop_gids=Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Convergence to a population
        pop_divergence  = {source_pop:{} for source_pop in pop_names}
        for source_pop in pop_names:
            for cell_gid in pop_gids[source_pop]:
                for conn in sim.net.cells[cell_gid].conns:
                    if conn['preGid']=='NetStim': continue # ignore NetStim conns
                    if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                    try: 
                        pre_pop = conn['label'].split('|')[1].split('@')[0]
                        post_pop= conn['label'].split('|')[2].split('@')[0]
                    except:
                        if verbose: print('skipping conn ',conn['label'])
                        continue

                    # pre_pop = conn['label'].split('|')[1].split('@')[0]
                    # post_pop= conn['label'].split('|')[2].split('@')[0]
                    if post_pop not in pop_divergence[pre_pop].keys():      pop_divergence[pre_pop].update({post_pop:1})
                    else:                                                   pop_divergence[pre_pop][post_pop]+=1
        return pop_divergence
    
    def getCellConvergence(sim,verbose=False):
        pop_gids        = Evaluate.getPopGids(sim)
        pop_convergence = Evaluate.getPopConvergence(sim,verbose)
        cell_convergence={}
        for target_pop in pop_convergence.keys():
            cell_convergence.update({target_pop:{}})
            for pre_pop in pop_convergence[target_pop].keys():
                cell_convergence[target_pop].update({pre_pop:pop_convergence[target_pop][pre_pop]/len(pop_gids[target_pop])})
        return cell_convergence
    
    def getCellConvergence_v2(sim,verbose=False):
        # --- Pop parameters
        pop_gids  = Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Convergence to a population
        cell_convergence_byCell = {target_pop:{} for target_pop in pop_names}
        for target_pop in pop_names:
            for cell_gid in pop_gids[target_pop]:
                for conn in sim.net.cells[cell_gid].conns:
                    if conn['preGid']=='NetStim': continue # ignore NetStim conns
                    if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                    try: 
                        pre_pop = conn['label'].split('|')[1].split('@')[0]
                    except:
                        if verbose: print('skipping conn ',conn['label'])
                        continue
                    # print(cell_gid, pre_pop)
                    # print(pre_pop,target_pop)
                    if pre_pop not in cell_convergence_byCell[target_pop].keys():  
                        cell_convergence_byCell[target_pop].update({pre_pop:{}})
                    if str(cell_gid) not in cell_convergence_byCell[target_pop][pre_pop].keys():
                        cell_convergence_byCell[target_pop][pre_pop].update({str(cell_gid):1})
                    else:                                                   cell_convergence_byCell[target_pop][pre_pop][str(cell_gid)]+=1
        
        cell_convergence_mean_std = {target_pop:{} for target_pop in cell_convergence_byCell.keys()}
        for target_pop in cell_convergence_byCell.keys():
            for pre_pop in cell_convergence_byCell[target_pop].keys():
                try:    cell_convergence_mean_std[target_pop].update({pre_pop:[np.mean(list(cell_convergence_byCell[target_pop][pre_pop].values())),np.std(list(cell_convergence_byCell[target_pop][pre_pop].values()))]})
                except: pass
        return cell_convergence_mean_std, cell_convergence_byCell

    def getCellDivergence(sim,verbose=False):
        pop_gids        = Evaluate.getPopGids(sim)
        pop_divergence  = Evaluate.getPopDivergence(sim,verbose)
        cell_divergence={}
        for pre_pop in pop_divergence.keys():
            cell_divergence.update({pre_pop:{}})
            for post_pop in pop_divergence[pre_pop].keys():
                cell_divergence[pre_pop].update({post_pop:pop_divergence[pre_pop][post_pop]/len(pop_gids[pre_pop])})
        return cell_divergence
    
    def getConvergenceRatio(sim,verbose=False):
        cell_convergence  = Evaluate.getCellConvergence(sim,verbose)
        convergence_ratio = {target_pop:{} for target_pop in cell_convergence.keys()}
        for target_pop in cell_convergence.keys():
            for source_pop in cell_convergence[target_pop].keys(): convergence_ratio[target_pop].update({source_pop:cell_convergence[target_pop][source_pop]/sum(cell_convergence[target_pop].values())})
        return convergence_ratio
    
    def getDivergenceRatio(sim,verbose=False):
        cell_divergence  = Evaluate.getCellDivergence(sim,verbose)
        divergence_ratio = {source_pop:{} for source_pop in cell_divergence.keys()}
        for source_pop in cell_divergence.keys():
            for target_pop in cell_divergence[source_pop].keys(): divergence_ratio[source_pop].update({target_pop:cell_divergence[source_pop][target_pop]/sum(cell_divergence[source_pop].values())})
        return divergence_ratio



    def runAll(sim,verbose=False):
        # --- Pop parameters
        pop_gids  = Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Pop  Convergence
        pop_convergence  = Evaluate.getPopConvergence( sim,verbose)
        # --- Pop  Divergence
        pop_divergence   = Evaluate.getPopDivergence(  sim,verbose)
        # --- Cell Convergence
        cell_convergence = Evaluate.getCellConvergence(sim,verbose)
        cell_convergence_mean_std, cell_convergence_byCell = Evaluate.getCellConvergence_v2(sim,verbose)
        # --- Cell Divergence
        cell_divergence  = Evaluate.getCellDivergence( sim,verbose)

        return pop_gids,pop_names,pop_convergence,pop_divergence,cell_convergence,cell_divergence, cell_convergence_mean_std, cell_convergence_byCell
    
    def printStats(sim,plotFigs=False,printExtendedStats=True,verbose=False):

        pop_gids,pop_names,pop_convergence,pop_divergence,cell_convergence,cell_divergence, cell_convergence_mean_std, cell_convergence_byCell = Evaluate.runAll(sim,verbose)

        save_connectivity_dict={
            'cell_convergence':  cell_convergence,
            'cell_divergence':   cell_divergence,
            'pop_convergence':   pop_convergence,
            'pop_divergence':    pop_divergence,
            'pop_names':         pop_names,
            'pop_gids':          pop_gids,
            'cell_convergence_mean_std':cell_convergence_mean_std,
            'cell_convergence_byCell':cell_convergence_byCell,
            }

        save_filename = '../conn/conn_validation/connectivity_data__'+sim.cfg.dump_cell_properties.split('/')[-1].split('.json')[0].split('barreloid_network_cell_properties__')[1]+'.json'
        # save_filename = 'figs_conn_validation/connectivity_data__'+sim.cfg.dump_cell_properties.split('/')[-1].split('.json')[0].split('barreloid_network_cell_properties__')[1]+'.json'
        import json
        with open(save_filename, 'w') as f:
            json.dump(save_connectivity_dict, f, indent=4)

        # --- Pathway ratios for cell convergence/divergence
        convergence_ratio = Evaluate.getConvergenceRatio(sim,verbose)
        divergence_ratio  = Evaluate.getDivergenceRatio( sim,verbose)

        print('\n------------------------------------------------------------------------------------------------------------------------------\n')

        print('\t>>\t--- Convergence Data ---')
        print('\t>>\tConvergence Absolute [(number of conns * syns per conn) for each cell]:  ')
        [print('\t\t',key,cell_convergence[key]) for key in cell_convergence.keys()]
        print('\t>>\tRatio - Convergence Absolute [(number of conns * syns per conn) for each cell]:  ')
        [print('\t\t',key,convergence_ratio[key]) for key in convergence_ratio.keys()]

        print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
        print('\t>>\tPathway ratios:  ')
        # --- Calculate convergence of L6A CT pops into VPM
        try:    conv_L6A_VPM=sum([cell_convergence['VPM'][spop] for spop in cell_convergence['VPM'].keys() if ('L6A' in spop) or ('CTvirtual' in spop)])
        except:conv_L6A_VPM=0
        # --- Calculate convergence of L6A CT pops into TRN
        try:    conv_L6A_TRN=sum([cell_convergence['TRN'][spop] for spop in cell_convergence['TRN'].keys() if ('L6A' in spop) or ('CTvirtual' in spop)])
        except: conv_L6A_TRN=0
        # --- Calculate convergence of TRN pops into VPM
        try:    conv_TRN_VPM=sum([cell_convergence['VPM'][spop] for spop in cell_convergence['VPM'].keys() if 'TRN' in spop])
        except: conv_TRN_VPM=0
        
        # --- Conns targeting VPM (BBP data)
        BBP_conns_L6A_to_VPM = sim.cfg.conn_data['VPM']['chem']['L6A']['convergence'][0]*sim.cfg.conn_data['VPM']['chem']['L6A']['synsPerConn'][0] # = 581.4493258129902
        BBP_conns_TRN_to_VPM = sim.cfg.conn_data['VPM']['chem']['TRN']['convergence'][0]*sim.cfg.conn_data['VPM']['chem']['TRN']['synsPerConn'][0] # = 454.1682034044653
        BBP_conns_MLe_to_VPM = sim.cfg.conn_data['VPM']['chem']['MLe']['convergence'][0]*sim.cfg.conn_data['VPM']['chem']['MLe']['synsPerConn'][0] # = 122.97733554449974

        # --- Conns targeting TRN (BBP data)
        BBP_conns_L6A_to_TRN = sim.cfg.conn_data['TRN']['chem']['L6A']['convergence'][0]*sim.cfg.conn_data['TRN']['chem']['L6A']['synsPerConn'][0] # = 255.23635294982
        BBP_conns_VPM_to_TRN = sim.cfg.conn_data['TRN']['chem']['VPM']['convergence'][0]*sim.cfg.conn_data['TRN']['chem']['VPM']['synsPerConn'][0] # = 161.43545983636517
        BBP_conns_TRN_to_TRN = sim.cfg.conn_data['TRN']['chem']['TRN']['convergence'][0]*sim.cfg.conn_data['TRN']['chem']['TRN']['synsPerConn'][0] # = 161.43545983636517

        # --- Ratio of (L6A->VPM)/(MLe->VPM)
        try:    ratio_L6ACT_MLe_to_VPM = conv_L6A_VPM/cell_convergence['VPM']['MLe']
        except: ratio_L6ACT_MLe_to_VPM = 0
        print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', ratio_L6ACT_MLe_to_VPM,           '\t-- target - BBP: ', BBP_conns_L6A_to_VPM/BBP_conns_MLe_to_VPM)

        # --- Ratio of (TRN->VPM)/(MLe->VPM)
        try:    ratio_TRN_MLe_to_VPM = conv_TRN_VPM/cell_convergence['VPM']['MLe']
        except: ratio_TRN_MLe_to_VPM = 0
        print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', ratio_TRN_MLe_to_VPM,             '\t-- target - BBP: ', BBP_conns_TRN_to_VPM/BBP_conns_MLe_to_VPM)

        # --- Ratio of (L6A->VPM)/(L6A->TRN)
        try:    ratio_L6ACT_to_TRN_VPM = conv_L6A_VPM/conv_L6A_TRN
        except: ratio_L6ACT_to_TRN_VPM = 0
        print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', ratio_L6ACT_to_TRN_VPM,           '\t-- target - BBP: ', BBP_conns_L6A_to_VPM/BBP_conns_L6A_to_TRN)
        
        # --- Ratio of (L6A->TRN)/(VPM->TRN)
        try:    ratio_L6ACT_VPM_to_TRN = conv_L6A_TRN/cell_convergence['TRN']['VPM']
        except: ratio_L6ACT_VPM_to_TRN = 0
        try: print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', ratio_L6ACT_VPM_to_TRN,      '\t-- target - BBP: ', BBP_conns_L6A_to_TRN/BBP_conns_VPM_to_TRN)
        except: 0
        
        # --- Ratio of (TRN->TRN)/(VPM->TRN)
        try:    ratio_TRN_VPM_to_TRN = cell_convergence['TRN']['TRN']/cell_convergence['TRN']['VPM']
        except: ratio_TRN_VPM_to_TRN = 0
        try: print('\t>>\t (TRN->TRN)/(VPM->TRN) ratio: ', ratio_TRN_VPM_to_TRN,      '\t-- target - BBP: ', BBP_conns_TRN_to_TRN/BBP_conns_VPM_to_TRN)
        except: 0

        pathway_ratios={
                '(toVPM)_L6A/MLe':  {'BBP':BBP_conns_L6A_to_VPM/BBP_conns_MLe_to_VPM,    'model':ratio_L6ACT_MLe_to_VPM},
                '(toVPM)_TRN/MLe':  {'BBP':BBP_conns_TRN_to_VPM/BBP_conns_MLe_to_VPM,    'model':ratio_TRN_MLe_to_VPM},
                '(fromL6A)_VPM/TRN':{'BBP':BBP_conns_L6A_to_VPM/BBP_conns_L6A_to_TRN,       'model':ratio_L6ACT_to_TRN_VPM},
                '(toTRN)_L6A/VPM':  {'BBP':BBP_conns_L6A_to_TRN/BBP_conns_VPM_to_TRN,      'model':ratio_L6ACT_VPM_to_TRN},
                }

        # --- Plots conn validation figure
        try:    
            name_flag__=sim.cfg.dump_cell_properties
            if '.pkl' in name_flag__:   name_flag_=name_flag__.split('.pkl')[0]
            else:                       name_flag_=name_flag__.split('.json')[0]
            name_flag = name_flag_.split('__')[1]+'__'+name_flag_.split('__')[2]
        except: name_flag=''
        if plotFigs: Evaluate.plotStats(pathway_ratios,name_flag=name_flag)

        print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

        print('\t>>\t--- Divergence Data ---')
        print('\t>>\tDivergence Absolute [(number of conns * syns per conn) for each cell]:  ')
        [print('\t\t',key,cell_divergence[key])  for key in cell_divergence.keys()]
        print('\t>>\tRatio - Divergence Absolute [(number of conns * syns per conn) for each cell]:  ')
        [print('\t\t',key,divergence_ratio[key]) for key in divergence_ratio.keys()]

        print('\n------------------------------------------------------------------------------------------------------------------------------\n')
        
        # --- BBP model stats
        print('\t>>\t--- BBP Convergence Data ---')
        print('\t>>\tBBP Convergence Absolue:  ')
        [print('\t\t',key2,'-to->',key,': ',sim.cfg.conn_data[key]['chem'][key2]['synsPerConn'][0]*sim.cfg.conn_data[key]['chem'][key2]['convergence'][0]) for key in sim.cfg.conn_data.keys() for key2 in sim.cfg.conn_data[key]['chem'].keys()]

        print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
        print('\t>>\tPathway ratios BBP:  ')
        print('\t>>\t (L6A->VPM)/(L6A->MLe) ratio: ', BBP_conns_L6A_to_VPM/BBP_conns_MLe_to_VPM)
        print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', BBP_conns_TRN_to_VPM/BBP_conns_MLe_to_VPM)
        print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', BBP_conns_L6A_to_VPM/BBP_conns_L6A_to_TRN)
        print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', BBP_conns_L6A_to_TRN/BBP_conns_VPM_to_TRN)
        print('\t>>\t (TRN->TRN)/(VPM->TRN) ratio: ', BBP_conns_TRN_to_TRN/BBP_conns_VPM_to_TRN)

        print('\n------------------------------------------------------------------------------------------------------------------------------\n')
    
        if printExtendedStats:

            # --- Ratio of synapses into VB thalamus from (Van Horn and Sherman, 2007)
            driver_RL   = 28
            modul_RS    = 146
            trn_F       = 48
            sum_all = driver_RL+modul_RS+trn_F
            print('\t>>\t--- Convergence Data (Van Horn and Sherman, 2007) ---')
            print('\t>>\tConvergence Proportional:  ')
            print('\t\t MLe-to->VPM: ', driver_RL)
            print('\t\t L6A-to->VPM: ', modul_RS)
            print('\t\t TRN-to->VPM: ', trn_F)

            print('\t>>\tConvergence Ratio:  ')
            print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
            print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
            print('\t\t TRN-to->VPM: ', trn_F/sum_all)

            print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
            print('\t>>\tPathway ratios:  ')
            print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
            print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)
            # print('\t>>\tBBP Convergence ratio:  ')
            # [print('\t\t',key,conn_data[key]) for key in conn_data.keys()]

            print('\n------------------------------------------------------------------------------------------------------------------------------\n')
            # print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

            '''
            (çavdar, 2011) - VB terminals and synapses - Fig 8
                            driver - RL             modulator - RS          trn - F
            terminals       4.91                    86.81818181818181       6.092307692307692
            synapses        12.248803827751196      79.6969696969697        7.569230769230771       (USE THIS VALUE)
            '''
            driver_RL   = 12.248803827751196
            modul_RS    = 79.6969696969697
            trn_F       = 7.569230769230771
            sum_all = driver_RL+modul_RS+trn_F

            print('\t>>\t--- Convergence Data (çavdar, 2007) ---')
            print('\t>>\tConvergence Proportional:  ')
            print('\t\t MLe-to->VPM: ', driver_RL)
            print('\t\t L6A-to->VPM: ', modul_RS)
            print('\t\t TRN-to->VPM: ', trn_F)

            print('\t>>\tConvergence Ratio:  ')
            print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
            print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
            print('\t\t TRN-to->VPM: ', trn_F/sum_all)

            print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
            print('\t>>\tPathway ratios:  ')
            print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
            print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)

            print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    def plotStats(pathway_ratios, name_flag=None):
        # --- Bar Plot with conn values
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 36})
        plt.figure( figsize=(20,20))

        for pathway_ratio_ind,pathway_ratio in enumerate(pathway_ratios.keys()):
            for source in pathway_ratios[pathway_ratio].keys():
                width = 0.25
                if source == 'model':
                    c   = 'green'
                    loc =  pathway_ratio_ind+(width/2)
                elif source == 'BBP':
                    c   = 'purple'
                    loc =  pathway_ratio_ind-(width/2)

                plt.bar(loc,pathway_ratios[pathway_ratio][source],width=width,color=c)

        plt.xticks([0,1,2,3],['L6A/MLe ->VPM','TRN/MLe ->VPM','L6A-> VPM/TRN','L6A/VPM ->TRN'],rotation=10)
        plt.xlabel('Pathway')
        plt.ylabel('Conn Ratio')
        plt.legend(['Thal_NB model + literature','NetPyNE model'])
        if name_flag is not None:
            try:    plt.savefig('../conn/conn_validation/figs_conn_validation/Conn_validation__'+str(name_flag)+'.png',dpi=500)
            except: plt.savefig('../conn/conn_validation/figs_conn_validation/Conn_validation_.png',dpi=500)
        else:
            plt.savefig('../conn/conn_validation/figs_conn_validation/Conn_validation.png',dpi=500)
        # if name_flag is not None:
        #     try:    plt.savefig('figs_conn_validation/Conn_validation__'+str(name_flag)+'.png',dpi=500)
        #     except: plt.savefig('figs_conn_validation/Conn_validation_.png',dpi=500)
        # else:
        #     plt.savefig('figs_conn_validation/Conn_validation.png',dpi=500)
    
    def verifyConns(sim, cellStats=True, popStats=True, plotFigs=False):
        conns_perCell={}
        conns_perPop={}
        for cell_ind in range(len(sim.net.cells)):
            # print('cell ', cell_ind)
            if 'cellModel' in sim.net.cells[cell_ind].tags.keys(): continue # skips VecStims
            secList_proximal    = sim.net.cells[cell_ind].secLists['inputs__proximal']
            secList_intermediate= sim.net.cells[cell_ind].secLists['inputs__intermediate']
            secList_distal      = sim.net.cells[cell_ind].secLists['inputs__distal']
            if len(sim.net.cells[cell_ind].conns)==0: continue
            
            conns_perCell.update(
                            {
                                cell_ind:{
                                            'proximal':     len([conn_ind for conn_ind in range(len(sim.net.cells[cell_ind].conns)) if sim.net.cells[cell_ind].conns[conn_ind]['sec'] in secList_proximal]),
                                            'intermediate': len([conn_ind for conn_ind in range(len(sim.net.cells[cell_ind].conns)) if sim.net.cells[cell_ind].conns[conn_ind]['sec'] in secList_intermediate]),
                                            'distal':       len([conn_ind for conn_ind in range(len(sim.net.cells[cell_ind].conns)) if sim.net.cells[cell_ind].conns[conn_ind]['sec'] in secList_distal]),
                                            }
                                }
                            )
            
            if cellStats: print('Cell # ', cell_ind, '\t - Conns (proximal: ', conns_perCell[cell_ind]['proximal'], ', intermediate: ', conns_perCell[cell_ind]['intermediate'], ', distal: ', conns_perCell[cell_ind]['distal'], ')')

            if 'pop' in sim.net.cells[cell_ind].tags.keys():
                pop = sim.net.cells[cell_ind].tags['pop']
                if pop not in conns_perPop.keys(): conns_perPop.update({pop:{'proximal':[],'intermediate':[],'distal':[]}})
                
                conns_perPop[pop]['proximal'].append(       conns_perCell[cell_ind]['proximal'])
                conns_perPop[pop]['intermediate'].append(   conns_perCell[cell_ind]['intermediate'])
                conns_perPop[pop]['distal'].append(         conns_perCell[cell_ind]['distal'])

        if popStats:
            for pop in conns_perPop.keys():
                print('Pop stats: ', pop, '\t - Avg Conns (proximal: ',     np.mean(conns_perPop[pop]['proximal']),     '+/-',      np.std(conns_perPop[pop]['proximal']), 
                                                        ', intermediate: ', np.mean(conns_perPop[pop]['intermediate']), '+/-',      np.std(conns_perPop[pop]['intermediate']), 
                                                        ', distal: ',       np.mean(conns_perPop[pop]['distal']),       '+/-',      np.std(conns_perPop[pop]['distal']), 
                                                        ')')
        
        if plotFigs: Evaluate.plotConnStats(conns_perPop)

    def plotConnStats(conns_perPop):
        # --- Bar Plot with conn values
        import numpy as np
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 36})

        # Data
        populations         = conns_perPop.keys()
        proximal_means      = [np.mean(conns_perPop[pop]['proximal']) for pop in populations]
        proximal_std        = [np.std(conns_perPop[pop]['proximal']) for pop in populations]
        
        intermediate_means  = [np.mean(conns_perPop[pop]['intermediate']) for pop in populations]
        intermediate_std    = [np.std(conns_perPop[pop]['intermediate']) for pop in populations]
        
        distal_means        = [np.mean(conns_perPop[pop]['distal']) for pop in populations]
        distal_std          = [np.std(conns_perPop[pop]['distal']) for pop in populations]
        
        # Number of groups
        n_groups = len(populations)

        # Set the position of each group of bars
        index = np.arange(n_groups)
        bar_width = 0.25

        # Plotting
        fig, ax = plt.subplots(figsize=(20,20))
        proximal_bar = ax.bar(index, proximal_means, bar_width, label='Proximal', yerr=proximal_std)
        intermediate_bar = ax.bar(index + bar_width, intermediate_means, bar_width, label='Intermediate', yerr=intermediate_std)
        distal_bar = ax.bar(index + 2 * bar_width, distal_means, bar_width, label='Distal', yerr=distal_std)

        # Add labels, title, and legend
        ax.set_xlabel('Populations')
        ax.set_ylabel('Average Connections')
        ax.set_title('Average Connections by Population and Synapse Type')
        ax.set_xticks(index + bar_width)
        ax.set_xticklabels(populations, rotation=45, ha='right')
        ax.legend()

        # Show plot
        plt.tight_layout()
        plt.savefig('../conn/conn_validation/figs_conn_validation/Conn_perPop.png',dpi=500)
        # plt.savefig('figs_conn_validation/Conn_perPop.png',dpi=500)

    def verifyPopConns(sim,prePop,postPop,verbose=False):
        inspect_conns_dict={}
        for cell_ind, cell in enumerate(sim.net.cells):
            conns=0
            for conn in cell.conns:
                try:
                    mech_prePop  = conn['synMech'].split('|')[1]
                    mech_postPop = conn['synMech'].split('|')[2]
                except: continue
                if prePop in conn['synMech'].split('|')[1]:
                    if verbose: print(mech_prePop)
                    conns+=1
            if postPop in cell.tags['pop']:
                if verbose: print(cell, conns)
                inspect_conns_dict.update({cell_ind:conns})
        print('\tAverage number of connections from ',prePop,' to ',postPop,': ', np.mean(list(inspect_conns_dict.values())))
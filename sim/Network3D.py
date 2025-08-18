import math
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from neuron import h
import copy
import json

class Network3D:

    # --- Runs the gap detection within a single method
    def runGapDetection(sim,pops=['TRN__pop','TRN_ring__pop'],max_distance=5,max_appositions_per_section=6,max_gaps_perCell=None,skipSoma=True,seed=100000,plotFig=True,filename='../conn/gap_connecivity/gap_connecivity_template_testingClass', verbose=False):
        ########################################################################################################################################################################
        print('\n\t=Building Gap junction network connectivity')
        print('\t-Mapping spatial coordinates')
        spatial_coordinates                     = Network3D.getSpatialCoordinates(sim,pops=pops,skipSoma=skipSoma)
        print('\t-Measuring distance-based connectivity')
        sections_within_distance                = Network3D.find_points_within_distance(spatial_coordinates, max_distance, verbose)

        ########################################################################################################################################################################
        # --- Stores the bidirectional connectivity in dictionary format, with post_gid as the referencing key
        print('\t-Filtering connection sources')
        filter_sources_bidirectional            = Network3D.getConnDictionary(sections_within_distance,verbose)

        ########################################################################################################################################################################
        # --- Scaling factor to compensate the fact that connections are bidirectional, but are sampled unidirectionally
        rescale_gap_mean                        = 1.6 # approximated from (1.638516179952644)
        gap_mean, gap_std                       = Network3D.getGapData(rescale_gap_mean = rescale_gap_mean)
        
        ########################################################################################################################################################################
        print('\t-Prunning connectivity')
        # --- Prune sources by sampling from all conns using a lognormal distribution
        sources_sampled                         = Network3D.sampleSources(filter_sources_bidirectional,gap_mean, gap_std, seed=seed)
        # --- Recalculating, looking at final bidirectional connectivity
        sources_sampled_bidirectional           = Network3D.getBidirectionalSampledSources(sources_sampled)

        ########################################################################################################################################################################
        # --- Prunning conns based on the defined <max_appositions_per_section> so that sections have a limited number of gap junctions
        prunned_connectivity_dict               = Network3D.pruneConnsUnidirectional(sources_sampled,sections_within_distance,max_appositions_per_section)
        prunned_connectivity_bidirectional_dict = Network3D.pruneConnsBidirectional(sources_sampled_bidirectional,sections_within_distance,max_appositions_per_section)

        ########################################################################################################################################################################
        if max_gaps_perCell is not None:
            print('Prunning Gaps')
            # --- Prunning conns based on the defined <max_gaps_per_cell> so that cells have a limited total number of gap junctions
            pruned_cell_dict = Network3D.pruneCellConnsUnidirectional(prunned_connectivity_dict,max_gaps_perCell,seed=seed)
            # pruned_cell_dict = pruneCellConnsUnidirectional(prunned_connectivity_dict,max_gaps_perCell=40,seed=seed)
            pruned_cell_bidirectional_dict = Network3D.getBidirectionalprunedCellConns(pruned_cell_dict)
        
        ########################################################################################################################################################################
        print('\t-Saving data')
        filename_unidirectional = filename + '_unidirectional'
        filename_bidirectional  = filename + '_bidirectional'
        if '.json' not in filename: 
            filename_unidirectional += '.json'
            filename_bidirectional  += '.json'
        Network3D.savePrunnedConns(prunned_connectivity_dict,               spatial_coordinates,filename=filename_unidirectional, verbose=verbose)
        Network3D.savePrunnedConns(prunned_connectivity_bidirectional_dict, spatial_coordinates,filename=filename_bidirectional,  verbose=verbose)
        print('\t-Plotting connectivity histogram')
        if plotFig: 
            if      plotFig==bool:  saveFileName = '../conn/gap_connecivity/figs/gap_conn_pruned_hist.png'
            elif    plotFig==str:   saveFileName = plotFig
            else:                   saveFileName = '../conn/gap_connecivity/figs/gap_conn_pruned_hist.png'
            Network3D.plotConnsHist(sim,sources_sampled,sources_sampled_bidirectional,
                                    prunned_connectivity_dict,prunned_connectivity_bidirectional_dict,
                                    saveFileName=saveFileName)
        ########################################################################################################################################################################
        if max_gaps_perCell is not None:
            print('Saving prunned Gaps')
            filename_unidirectional_maxGaps = filename + '_unidirectional_maxGaps'
            filename_bidirectional_maxGaps  = filename + '_bidirectional_maxGaps'
            if '.json' not in filename: 
                filename_unidirectional_maxGaps  += '.json'
            Network3D.savePrunnedConns(pruned_cell_dict,                        spatial_coordinates,filename=filename_unidirectional_maxGaps, verbose=verbose)
            Network3D.savePrunnedConns(pruned_cell_bidirectional_dict,          spatial_coordinates,filename=filename_bidirectional_maxGaps,  verbose=verbose)

            if plotFig: 
                if      plotFig==bool:  saveFileName2 = '../conn/gap_connecivity/figs/gap_conn_pruned_hist2.png'
                else:                   saveFileName2 = '../conn/gap_connecivity/figs/gap_conn_pruned_hist2.png'
                Network3D.plotConnsHist(sim,sources_sampled,sources_sampled_bidirectional,
                                        pruned_cell_dict,pruned_cell_bidirectional_dict,
                                        saveFileName=saveFileName2)


        # Network3D.pruneConns(sim, sections_within_distance, spatial_coordinates,max_appositions_per_section = 6,seed=100000,plotFig=False)

    ########################################################################################################################################################################
    # --- Gets all the pt3d for all sections in the cells within the selected pops
    def getSpatialCoordinates(sim,pops=['TRN__pop','TRN_ring__pop'],skipSoma=True):
        cell_gids=[]
        for pop in pops:cell_gids.extend(sim.net.pops[pop].cellGids)
        
        spatial_coordinates={}
        for cell_ind in range(len(sim.net.cells)):
            if cell_ind in cell_gids:
                if 'secs' in sim.net.cells[cell_ind].__dict__.keys():
                    spatial_coordinates.update({cell_ind:{}})
                    for sec in sim.net.cells[cell_ind].__dict__['secs'].keys():
                        # --- Skips gap junctions in the soma
                        if skipSoma and ('soma' in sec): continue
                        # --- Skips gap junctions in the basal dendrites (directly attached to the soma)
                        if skipSoma and ('soma' in sim.net.cells[cell_ind].__dict__['secs'][sec]['topol']['parentSec']): 
                            print('skipping conn in dend ', sec, ' - connected to ', sim.net.cells[cell_ind].__dict__['secs'][sec]['topol']['parentSec'])
                            continue
                        if (type(sim.net.cells[cell_ind].__dict__['secs'][sec]['geom']['pt3d'])==list) and (len(sim.net.cells[cell_ind].__dict__['secs'][sec]['geom']['pt3d'])>0) and (len(sim.net.cells[cell_ind].__dict__['secs'][sec]['geom']['pt3d'][0])==4):
                            spatial_coordinates[cell_ind].update({sec:[[x+sim.net.cells[cell_ind].tags['x'],y+sim.net.cells[cell_ind].tags['y'],z+sim.net.cells[cell_ind].tags['z'],d] for [x,y,z,d] in sim.net.cells[cell_ind].__dict__['secs'][sec]['geom']['pt3d']]})
        return spatial_coordinates

    def getSpatialMapping(spatial_coordinates):
        spatial_mapping=[]
        for cell_ind in spatial_coordinates.keys():
            for sec in spatial_coordinates[cell_ind]:
                for coord in spatial_coordinates[cell_ind][sec]:
                    spatial_mapping.append((coord,(cell_ind,sec)))
                # spatial_mapping.append([(coord,(cell_ind,sec)) for coord in spatial_coordinates[cell_ind][sec]])
        
        return spatial_mapping

    # --- Finds the pairs of intersections based on the distance between 3d points
    def find_points_within_distance(spatial_coordinates, max_distance, verbose=False):
        """
        Find pairs of points within a certain defined distance of each other using a KD-tree.
        
        Args:
        - spatial_coordinates (dict): Dictionary containing spatial coordinates of points.
        - max_distance (float): The maximum allowed distance between points.
        
        Returns:
        - list: List of tuples containing pairs of points within the specified distance.
        """

        all_points = []
        all_points_reference=[]
        for cell_ind in spatial_coordinates.keys():
            for sec in spatial_coordinates[cell_ind].keys():
                all_points.extend(spatial_coordinates[cell_ind][sec])
                refs=[(cell_ind,sec,l) for l in range(len(spatial_coordinates[cell_ind][sec]))]
                all_points_reference.extend(refs)

        tree = KDTree(all_points)
        # pairs_within_distance = []
        sections_within_distance=[]
        for point_ind,point in enumerate(all_points):
            point_reference_cell     = all_points_reference[point_ind][0]
            point_reference_section  = all_points_reference[point_ind][1]
            point_reference_pt3d_ind = all_points_reference[point_ind][1]
            if verbose: print('processing cell/sec: ',point_reference_cell,' - ',point_reference_section, ' \tremaining ',len(all_points)-point_ind)
            # Find points within the specified distance of the current point
            indices = tree.query_ball_point(point, max_distance)
            # Exclude points within the same cell
            indices = [i for i in indices if all_points_reference[i][0]!=point_reference_cell]
            # # Exclude points within the same cell
            # indices = [i for i in indices if all_points[i] not in spatial_coordinates[point[0]].values()]
            # Add pairs of points within distance to the result
            for idx in indices:
                # pairs_within_distance.append((point, all_points[idx]))
                sections_within_distance.append((all_points_reference[point_ind], all_points_reference[idx]))

        return sections_within_distance
        # return pairs_within_distance, sections_within_distance
    
    # --- Stores the bidirectional connectivity in dictionary format, with post_gid as the referencing key
    def getConnDictionary(sections_within_distance,verbose=False):
        filter_sources_bidirectional={}
        for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(sections_within_distance):
            if post_gid not in list(filter_sources_bidirectional.keys()):   
                if verbose: print('Adding ',post_gid)
                filter_sources_bidirectional.update({post_gid:[]})
            if pre_gid not in filter_sources_bidirectional[post_gid]:       filter_sources_bidirectional[post_gid].append(pre_gid)
        count_filtered_sources_bidirectional=[len(filter_sources_bidirectional[post_gid]) for post_gid in filter_sources_bidirectional.keys()]
        neuronal_coupling = [np.mean(count_filtered_sources_bidirectional),np.std(count_filtered_sources_bidirectional)]
        return filter_sources_bidirectional

    def getGapData(rescale_gap_mean=None):
        '''
        # --- Lee 2014
        gap_number = [2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 12, 12, 16, 17, 18, 21, 21]
        np.mean(gap_number) # = 8.625+-5.072905971925756
        # --- Predicted number of gap junctions
        num_of_gaps  = [14.25,              22.75,              31.25,              39.75,              48.25,              56.75,              65.25,              73.75,              82.25]
        neuron_count = [99.66555183946485,  226.75585284280933, 311.03678929765886, 168.56187290969896, 128.42809364548492, 42.809364548494955, 15.384615384615358, 6.688963210702298,  3.3444816053511772]
        '''
        # --- Number of gap junctions (Lee, 2014)
        gap_number = [2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 12, 12, 16, 17, 18, 21, 21]
        gap_mean=np.mean(gap_number)    # = 8.625
        gap_std =np.std(gap_number)     # = +-5.072905971925756

        # --- Scaling factor to compensate the fact that connections are bidirectional, but are sampled unidirectionally
        if rescale_gap_mean is not None: gap_mean = gap_mean/rescale_gap_mean    # worked: bidirectional mean went from ~11.6 to ~8.24, if unidirectional gap_mean is scaled from 8.625 to 5.390625

        return gap_mean, gap_std
    
    def sampleSources(filter_sources_bidirectional,gap_mean,gap_std,seed=100000):
        if seed is not None:    
            rand            = h.Random(seed)
            rand_poisson    = h.Random(seed+1)
            rand_lognormal  = h.Random(seed+2)
        else:                   
            rand            = h.Random(0)
            rand_poisson    = h.Random(1)
            rand_lognormal  = h.Random(2)

        sources_sampled={post_gid:[] for post_gid in filter_sources_bidirectional.keys()}
        for post_gid in filter_sources_bidirectional.keys():
            sample_sources_ = set()
            max_divergence= int(rand_lognormal.lognormal(gap_mean,gap_std**2)) # --- 2 to 20 coupled neurons (Iavarone, 2023 - Fig 3)
            # max_divergence= int(rand_poisson.poisson(gap_mean)) # --- 2 to 20 coupled neurons (Iavarone, 2023 - Fig 3)
            # --- Dont need to resample if there a less indexes that the target number of gap_conns
            if len(filter_sources_bidirectional[post_gid])<=max_divergence:
                sample_sources_=filter_sources_bidirectional[post_gid]
            else:
                while len(sample_sources_) < max_divergence:    
                    sample_index_min=0
                    sample_index_max=len(filter_sources_bidirectional[post_gid])
                    sample_sources_.add(filter_sources_bidirectional[post_gid][int(rand.uniform(sample_index_min,sample_index_max))])
            # print(post_gid, list(sample_sources_))
            sample_sources_=list(sample_sources_)
            sample_sources_.sort()
            sources_sampled[post_gid]=sample_sources_
        return sources_sampled
    
    # --- Recalculating, looking at final bidirectional connectivity
    def getBidirectionalSampledSources(sources_sampled):
        sources_sampled_bidirectional=copy.deepcopy(sources_sampled)
        for post_gid in sources_sampled_bidirectional.keys():
            for pre_gid in sources_sampled_bidirectional[post_gid]:
                if pre_gid not in sources_sampled_bidirectional.keys():sources_sampled_bidirectional.update({pre_gid:[]})
                if post_gid not in sources_sampled_bidirectional[pre_gid]:sources_sampled_bidirectional[pre_gid].append(post_gid)
        return sources_sampled_bidirectional
    
    def pruneConnsUnidirectional(sources_sampled,sections_within_distance,max_appositions_per_section):
        prunned_connectivity_dict={post_gid:[] for post_gid in sources_sampled.keys()}
        track_gids={post_gid:[] for post_gid in sources_sampled.keys()}
        for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(sections_within_distance):
            if (pre_gid in sources_sampled[post_gid]) and track_gids[post_gid].count(pre_gid)<max_appositions_per_section:
                prunned_connectivity_dict[post_gid].append(((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)))
                track_gids[post_gid].append(pre_gid)
        return prunned_connectivity_dict
    
    def pruneConnsBidirectional(sources_sampled_bidirectional,sections_within_distance,max_appositions_per_section):
        prunned_connectivity_bidirectional_dict={post_gid:[] for post_gid in sources_sampled_bidirectional.keys()}
        track_gids_bidirectional={post_gid:[] for post_gid in sources_sampled_bidirectional.keys()}
        for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(sections_within_distance):
            if (pre_gid in sources_sampled_bidirectional[post_gid]) and track_gids_bidirectional[post_gid].count(pre_gid)<max_appositions_per_section:
                prunned_connectivity_bidirectional_dict[post_gid].append(((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)))
                track_gids_bidirectional[post_gid].append(pre_gid)
        return prunned_connectivity_bidirectional_dict
    
    def pruneCellConnsUnidirectional(connectivity_dict,max_gaps_perCell,seed=100000):

        gap_mean, gap_std = Network3D.getGapData(rescale_gap_mean = None)

        if seed is not None:    
            rand            = h.Random(seed)
            rand_poisson    = h.Random(seed+1)
            rand_lognormal  = h.Random(seed+2)
        else:                   
            rand            = h.Random(0)
            rand_poisson    = h.Random(1)
            rand_lognormal  = h.Random(2)

        # --- Stores the index of each conn made by a given GID for sampling
        store_postGid_indexes={}
        for post_gid in connectivity_dict.keys():
            if post_gid not in store_postGid_indexes.keys():    store_postGid_indexes.update({post_gid:[]}) # - for internal processing
            for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(connectivity_dict[post_gid]): 
                store_postGid_indexes[post_gid].append(ind)

        # --- Sampling conns
        store_postGid_indexes_pruned={}
        for post_gid in store_postGid_indexes.keys():
            
            rand_maxGaps = h.Random(post_gid)
            max_gaps_ = int(rand_maxGaps.lognormal(gap_mean,gap_std**2))
            if max_gaps_<= max_gaps_perCell:    random_max_gaps_perCell = max_gaps_
            else:                               random_max_gaps_perCell = gap_mean

            # --- In case the cell has less or the same number of gap conns as the target number
            if len(store_postGid_indexes[post_gid])<=max_gaps_perCell:
                sampled_indexes = store_postGid_indexes[post_gid]
            else:
                sampled_indexes=[]
                while(len(sampled_indexes)<random_max_gaps_perCell):
                    sampled_conn_ind = int(rand.uniform(0,len(store_postGid_indexes[post_gid])))
                    if sampled_conn_ind not in sampled_indexes: 
                        sampled_indexes.append(sampled_conn_ind)
                        print('cell gid: ', post_gid, '\t| sampled_indexes: ' ,len(sampled_indexes), ' | max gaps for this cell: ', random_max_gaps_perCell)
                sampled_indexes.sort()
            store_postGid_indexes_pruned.update({post_gid:sampled_indexes})

        pruned_cell_dict={}
        for post_gid in connectivity_dict.keys():
            if post_gid not in pruned_cell_dict.keys():         pruned_cell_dict.update(     {post_gid:[]}) # - for output
            for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(connectivity_dict[post_gid]):
                if ind in store_postGid_indexes_pruned[post_gid]:
                    pruned_cell_dict[post_gid].append(connectivity_dict[post_gid][ind])

        return pruned_cell_dict

    def getBidirectionalprunedCellConns(pruned_cell_dict):
        pruned_cell_bidirectional_dict=copy.deepcopy(pruned_cell_dict)
        for post_gid in pruned_cell_bidirectional_dict.keys():
            for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(pruned_cell_bidirectional_dict[post_gid]):
                if ((pre_gid,pre_sec,pre_pt3d_ind),(post_gid,post_sec,post_pt3d_ind)) not in pruned_cell_bidirectional_dict[pre_gid]:
                    pruned_cell_bidirectional_dict[pre_gid].append(((pre_gid,pre_sec,pre_pt3d_ind),(post_gid,post_sec,post_pt3d_ind)))
        return pruned_cell_bidirectional_dict
        

    def savePrunnedConns(final_conns_dictionary,spatial_coordinates,filename='../conn/gap_connecivity/gap_connecivity_template.json', verbose=False):
        # --- Calculating the 'loc' value for each point - and removing 0 and 1 locs for 0.001 and 0.999
        final_dict={}
        # for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(prunned_connectivity_bidirectional_dict):
        for post_gid in final_conns_dictionary.keys():
            if post_gid not in final_dict.keys(): final_dict.update({post_gid:[]})

            for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(final_conns_dictionary[post_gid]):
                post_loc = post_pt3d_ind/len(spatial_coordinates[post_gid][post_sec])
                pre_loc  = pre_pt3d_ind/len(spatial_coordinates[pre_gid][pre_sec])

                if   post_loc==1.0: post_loc=0.999
                elif post_loc==0.0: post_loc=0.001
                if   pre_loc==1.0:  pre_loc=0.999
                elif pre_loc==0.0:  pre_loc=0.001

                final_dict[post_gid].append(((post_gid,post_sec,post_pt3d_ind,post_loc),(pre_gid,pre_sec,pre_pt3d_ind,pre_loc)))
        
        final_list=[]
        for post_gid in final_dict.keys(): final_list.extend(final_dict[post_gid])


        # Write the list to a JSON file
        print('\t\tSaving connectivity template in: \t', filename)
        with open(filename, "w") as json_file: json.dump(final_list, json_file)
    
    def plotConnsHist(sim,sources_sampled,sources_sampled_bidirectional,prunned_connectivity_dict,prunned_connectivity_bidirectional_dict,saveFileName='../conn/gap_connecivity/figs/test_prune_hist.png'):

        n_sources_sampled=[len(sources_sampled[post_gid]) for post_gid in sources_sampled.keys()]
        n_sources_sampled_mainPop=[len(sources_sampled[post_gid]) for post_gid in sources_sampled.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_sources_sampled_bidirectional=[len(sources_sampled_bidirectional[post_gid]) for post_gid in sources_sampled_bidirectional.keys()]
        n_sources_sampled_bidirectional_mainPop=[len(sources_sampled_bidirectional[post_gid]) for post_gid in sources_sampled_bidirectional.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_prunned_connectivity_dict=[len(prunned_connectivity_dict[post_gid]) for post_gid in prunned_connectivity_dict.keys()]
        n_prunned_connectivity_dict_mainPop=[len(prunned_connectivity_dict[post_gid]) for post_gid in prunned_connectivity_dict.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_prunned_connectivity_bidirectional_dict=[len(prunned_connectivity_bidirectional_dict[post_gid]) for post_gid in prunned_connectivity_bidirectional_dict.keys()]
        # --- quantifying only cells in the main pop
        n_prunned_connectivity_bidirectional_dict_mainPop=[len(prunned_connectivity_bidirectional_dict[post_gid]) for post_gid in prunned_connectivity_bidirectional_dict.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]


        # --- Predicted number of gap junctions
        predicted_num_gaps  = [14.25,              22.75,              31.25,              39.75,              48.25,              56.75,              65.25,              73.75,              82.25]
        neuron_count        = [99.66555183946485,  226.75585284280933, 311.03678929765886, 168.56187290969896, 128.42809364548492, 42.809364548494955, 15.384615384615358, 6.688963210702298,  3.3444816053511772]

        hist_step=10
        # hist_bins=[0,200,hist_step]
        hist_bins=predicted_num_gaps
        plt.figure(figsize=(20,20))
        plt.rcParams.update({'font.size': 30})
        
        plt.subplot(2,2,1)
        # plt.hist(n_sources_sampled_bidirectional,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='r')
        plt.hist(n_sources_sampled,bins=list(range(0,20,1)),color='b')
        gap_hist_unidirectional = plt.hist(n_sources_sampled_mainPop,bins=list(range(0,20,1)),color='lightblue')
        plt.title('')
        plt.xlabel('Num of neurons coupled\nto target neuron')

        plt.subplot(2,2,3)
        plt.hist(n_prunned_connectivity_dict,bins=hist_bins,color='b',width=6.5)
        plt.hist(n_prunned_connectivity_dict_mainPop,bins=hist_bins,color='lightblue',width=4.5)
        # plt.hist(n_prunned_connectivity_dict,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='b',width=hist_bins[2]*0.8)
        plt.title('')
        plt.xlabel('Predicted number of gap junctions')

        plt.subplot(2,2,2)
        # plt.hist(n_sources_sampled_bidirectional,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='r')
        plt.hist(n_sources_sampled_bidirectional,bins=list(range(0,50,1)),color='b')
        gap_hist_bidirectional = plt.hist(n_sources_sampled_bidirectional_mainPop,bins=list(range(0,50,1)),color='lightblue')
        plt.title('')
        plt.xlabel('Bidirectional - Num of neurons coupled\nto target neuron')
        
        plt.subplot(2,2,4)
        plt.hist(n_prunned_connectivity_bidirectional_dict,bins=hist_bins,color='b',width=6.5)
        plt.hist(n_prunned_connectivity_bidirectional_dict_mainPop,bins=hist_bins,color='lightblue',width=4.5)
        # plt.hist(n_prunned_connectivity_bidirectional_dict,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='b',width=hist_bins[2]*0.8)
        plt.title('')
        plt.xlabel('Bidirectional - Predicted number of gap junctions')

        print('\t\tPlotting Gap Connectivity histogram in: \t', saveFileName)
        plt.savefig(saveFileName,dpi=300)
        
        # PAPER FIG: Gap connectivity validation figure - plotting the number of neurons coupled to each cell
        # NOTE: We used the unilateral connectivity template to create the connections in netpyne, because it already creates the bidirectional conns by default
        #       but in this code we need to calculate the bidirectional connectivity by hand, so that the statistics match the network created in netpyne.
        
        color_lee = 'k'
        color_allTRN = 'cornflowerblue'
        color_mainTRN = 'mediumblue'

        plt.figure(figsize=(15,10))
        plt.rcParams.update({'font.size': 25})
        # plt.subplot(2,1,2)
        gap_hist_bidirectional_allTRN = plt.hist(n_sources_sampled_bidirectional,           bins=list(np.arange(0.5,20,1)), width=0.75,     color=color_allTRN,align='mid',alpha=1.0)
        gap_hist_bidirectional_mainTRN = plt.hist(n_sources_sampled_bidirectional_mainPop,  bins=list(np.arange(0.5,20,1)), width=0.5,    color=color_mainTRN,align='mid',alpha=1.0)

        gap_mean, gap_std = Network3D.getGapData(rescale_gap_mean = None)           # data from Lee, 2014
        model_gap_allTRNmean = np.mean(n_sources_sampled_bidirectional)             # (mean) model data - all TRN cells
        model_gap_allTRNstd  = np.std(n_sources_sampled_bidirectional)              # (std)  model data - all TRN cells
        model_gap_mainTRN_mean = np.mean(n_sources_sampled_bidirectional_mainPop)   # (mean) model data - main TRN cells (center column)
        model_gap_mainTRN_std  = np.std(n_sources_sampled_bidirectional_mainPop)    # (std)  model data - main TRN cells (center column)

        plt.plot(gap_mean,35,color_lee, marker = 'o', markersize=10, label='Lee, 2014')
        plt.plot(model_gap_allTRNmean,33,color_allTRN, marker = '^', markersize=10, label= 'TRN - all cells')
        plt.plot(model_gap_mainTRN_mean,31,color_mainTRN, marker = 'd', markersize=10, label='TRN - central column')

        plt.errorbar(gap_mean,35,xerr=gap_std,color=color_lee,capsize=7)
        plt.errorbar(model_gap_allTRNmean,33,xerr=model_gap_allTRNstd,color=color_allTRN,capsize=7)
        plt.errorbar(model_gap_mainTRN_mean,31,xerr=model_gap_mainTRN_std,color=color_mainTRN,capsize=7)
        plt.legend()
        plt.title('')
        plt.xlim([0,21])
        plt.xticks([0,5,10,15,20])
        plt.xlabel('Number of neurons coupled to target neuron')
        plt.ylim([0,51])
        # plt.ylim([0,31])
        plt.yticks([0,5,10,15,20,25,30])
        plt.ylabel('Neuron count                                      ')
        ax1=plt.gca()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        plt.savefig('../conn/gap_connecivity/figs/paper_fig__gap_validation__bidirectional_connectivity.png',dpi=300)



    # --- Prunes the connectivity based on maximum number of appositions per section, to avoid overconnected sections in crowded regions
    # NOTE: No need to accound for bidirectionality of conns in <sections_within_distance>, because the point estimation already creates symetric A->B and B->A connections
    # NOTE: But, after prunning, bidirectional conns have to be accounted again
    def pruneConns(sim, sections_within_distance, spatial_coordinates,max_appositions_per_section=6,seed=100000,plotFig=False):

        # --- Stores the bidirectional connectivity in dictionary format, with post_gid as the referencing key
        filter_sources_bidirectional = Network3D.getConnDictionary(sections_within_distance)

        # --- Scaling factor to compensate the fact that connections are bidirectional, but are sampled unidirectionally
        rescale_gap_mean = 1.6                  # approximated from (1.638516179952644)
        gap_mean, gap_std = Network3D.getGapData(rescale_gap_mean = rescale_gap_mean)
        
        # --- Prune sources by sampling from all conns using a lognormal distribution
        sources_sampled = Network3D.sampleSources(filter_sources_bidirectional,gap_mean, gap_std, seed=seed)
        # --- Recalculating, looking at final bidirectional connectivity
        sources_sampled_bidirectional = Network3D.getBidirectionalSampledSources(sources_sampled)

        # --- Prunning conns based on the defined <max_appositions_per_section> so that cells have a limited number of gap junctions
        prunned_connectivity_dict               = Network3D.pruneConnsUnidirectional(sources_sampled,sections_within_distance,max_appositions_per_section)
        prunned_connectivity_bidirectional_dict = Network3D.pruneConnsBidirectional(sources_sampled_bidirectional,sections_within_distance,max_appositions_per_section)

        
        ########################################################################################################################################################################
        Network3D.savePrunnedConns(prunned_connectivity_dict,               spatial_coordinates,filename='../conn/gap_connecivity/gap_connecivity_template.json')
        Network3D.savePrunnedConns(prunned_connectivity_bidirectional_dict, spatial_coordinates,filename='../conn/gap_connecivity/gap_connecivity_template_bidirectional.json')
        ########################################################################################################################################################################
        
        # if plotFig:


        n_sources_sampled=[len(sources_sampled[post_gid]) for post_gid in sources_sampled.keys()]
        n_sources_sampled_mainPop=[len(sources_sampled[post_gid]) for post_gid in sources_sampled.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_sources_sampled_bidirectional=[len(sources_sampled_bidirectional[post_gid]) for post_gid in sources_sampled_bidirectional.keys()]
        n_sources_sampled_bidirectional_mainPop=[len(sources_sampled_bidirectional[post_gid]) for post_gid in sources_sampled_bidirectional.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_prunned_connectivity_dict=[len(prunned_connectivity_dict[post_gid]) for post_gid in prunned_connectivity_dict.keys()]
        n_prunned_connectivity_dict_mainPop=[len(prunned_connectivity_dict[post_gid]) for post_gid in prunned_connectivity_dict.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]

        n_prunned_connectivity_bidirectional_dict=[len(prunned_connectivity_bidirectional_dict[post_gid]) for post_gid in prunned_connectivity_bidirectional_dict.keys()]
        # --- quantifying only cells in the main pop
        n_prunned_connectivity_bidirectional_dict_mainPop=[len(prunned_connectivity_bidirectional_dict[post_gid]) for post_gid in prunned_connectivity_bidirectional_dict.keys() if post_gid in sim.net.pops['TRN__pop'].cellGids]


        # --- Predicted number of gap junctions
        predicted_num_gaps  = [14.25,              22.75,              31.25,              39.75,              48.25,              56.75,              65.25,              73.75,              82.25]
        neuron_count        = [99.66555183946485,  226.75585284280933, 311.03678929765886, 168.56187290969896, 128.42809364548492, 42.809364548494955, 15.384615384615358, 6.688963210702298,  3.3444816053511772]

        hist_step=10
        # hist_bins=[0,200,hist_step]
        hist_bins=predicted_num_gaps
        plt.figure(figsize=(20,20))
        plt.rcParams.update({'font.size': 30})
        
        plt.subplot(2,2,1)
        # plt.hist(n_sources_sampled_bidirectional,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='r')
        plt.hist(n_sources_sampled,bins=list(range(0,20,1)),color='b')
        plt.hist(n_sources_sampled_mainPop,bins=list(range(0,20,1)),color='lightblue')
        plt.title('')
        plt.xlabel('Num of neurons coupled\nto target neuron')
        plt.subplot(2,2,3)
        plt.hist(n_prunned_connectivity_dict,bins=hist_bins,color='b',width=6.5)
        plt.hist(n_prunned_connectivity_dict_mainPop,bins=hist_bins,color='lightblue',width=4.5)
        # plt.hist(n_prunned_connectivity_dict,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='b',width=hist_bins[2]*0.8)
        plt.title('')
        plt.xlabel('Predicted number of gap junctions')

        plt.subplot(2,2,2)
        # plt.hist(n_sources_sampled_bidirectional,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='r')
        plt.hist(n_sources_sampled_bidirectional,bins=list(range(0,50,1)),color='b')
        plt.hist(n_sources_sampled_bidirectional_mainPop,bins=list(range(0,50,1)),color='lightblue')
        plt.title('')
        plt.xlabel('Bidirectional - Num of neurons coupled\nto target neuron')
        plt.subplot(2,2,4)
        plt.hist(n_prunned_connectivity_bidirectional_dict,bins=hist_bins,color='b',width=6.5)
        plt.hist(n_prunned_connectivity_bidirectional_dict_mainPop,bins=hist_bins,color='lightblue',width=4.5)
        # plt.hist(n_prunned_connectivity_bidirectional_dict,bins=list(range(hist_bins[0],hist_bins[1],hist_bins[2])),color='b',width=hist_bins[2]*0.8)
        plt.title('')
        plt.xlabel('Bidirectional - Predicted number of gap junctions')

        plt.savefig('../conn/gap_connecivity/figs/test_prune_hist.png',dpi=300)


        ####################################################################################
        

        ####################################################################################
        

            


    def analyzeConnPairs(sections_within_distance):
        conn_pairs_dict={}
        for ind, ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in enumerate(sections_within_distance):
            if post_gid not in conn_pairs_dict.keys(): conn_pairs_dict.update({post_gid:[]})
            conn_pairs_dict[post_gid].append(((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)))
        num_conns=[len(conn_pairs_dict[post_gid]) for post_gid in conn_pairs_dict.keys()]
        return num_conns

    def calculate_absolute_distance(point1, point2):
        """
        Calculate the absolute distance between two points in 3D space.
        
        Args:
        - point1 (tuple): A tuple containing the XYZ coordinates of the first point (x1, y1, z1).
        - point2 (tuple): A tuple containing the XYZ coordinates of the second point (x2, y2, z2).
        
        Returns:
        - float: The absolute distance between the two points.
        """
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        
        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        
        return distance

    def validateDistances(sim,sections_within_distance,spatial_coordinates,max_distance=1):
        validated_conns=0
        for ((post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)) in sections_within_distance:
            
            post_point = spatial_coordinates[post_gid][post_sec][post_pt3d_ind]
            pre_point  = spatial_coordinates[pre_gid][pre_sec][pre_pt3d_ind]
            
            # if calculate_absolute_distance(post_point[0:3], pre_point[0:3])>max_distance:
            if Network3D.calculate_absolute_distance(post_point[0:3], pre_point[0:3])>max_distance:
                print('Failed: ', post_gid,post_sec,post_pt3d_ind,' \t ',pre_gid,pre_sec,pre_pt3d_ind)
            else:validated_conns+=1
            if len(sections_within_distance)==validated_conns:print('All conns are valid')

    def plotConn(sim,sec_pair,spatial_coordinates):
        
        (post_gid,post_sec,post_pt3d_ind),(pre_gid,pre_sec,pre_pt3d_ind)=sec_pair
        # cells = [{'secs':spatial_coordinates[cell]} for cell in [post_gid,pre_gid]]
        cell_gids=[post_gid,pre_gid]
        cell_secs=[post_sec,pre_sec]
        cells = [sim.net.cells[cell].__dict__ for cell in cell_gids]

        post_point = spatial_coordinates[post_gid][post_sec][post_pt3d_ind]
        pre_point  = spatial_coordinates[pre_gid][pre_sec][pre_pt3d_ind]

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(projection='3d')
        colors=['b','grey']
        for ind,cell in enumerate(cells):
            # --- Plot Soma center - SWC files are loaded with cells centered at (0,0,0)
            # ax.scatter(0,0,0, marker='o',color='r',linewidths=0, alpha=1)

            # --- Plot 3d locations of cell morphology sections
            x3d=[];     y3d=[];     z3d=[];     d3d=[]
            for sec in cell['secs'].keys():
                if sec == cell_secs[ind]:   color='r'
                else:                       color = colors[ind]
                # if 'soma' in sec: 
                #     pt3d_diams = [pt3d[3] for ind_pt3d, pt3d in enumerate(cell['secs'][sec]['geom']['pt3d'])]
                #     ax.scatter(0,0,0, s= max(pt3d_diams), marker='o',color='k', alpha=1)
                #     continue
                if 'soma_0' in sec: 
                    pt3d_diams = [pt3d[3] for ind_pt3d, pt3d in enumerate(cell['secs'][sec]['geom']['pt3d'])]
                    pt3d_ = [sim.net.cells[cell_gids[ind]].tags[axis]  for axis_ind, axis in enumerate(['x','y','z'])]
                    ax.scatter(pt3d_[0],pt3d_[1],pt3d_[2], s= max(pt3d_diams), marker='o',color='k', alpha=1)
                    continue
                if 'pt3d' not in cell['secs'][sec]['geom'].keys():continue
                if len(cell['secs'][sec]['geom']['pt3d'])>0 and type(cell['secs'][sec]['geom']['pt3d'])==list:
                    for ind_pt3d, pt3d in enumerate(cell['secs'][sec]['geom']['pt3d']):
                        if ind_pt3d==0:
                            parent_sec = cell['secs'][sec]['topol']['parentSec']
                            previous_pt3d = cell['secs'][parent_sec]['geom']['pt3d'][-1]
                        else:
                            previous_pt3d = cell['secs'][sec]['geom']['pt3d'][ind_pt3d-1]

                        pt3d_            = [pt3d[axis_ind]         + sim.net.cells[cell_gids[ind]].tags[axis]  for axis_ind, axis in enumerate(['x','y','z'])]
                        previous_pt3d_   = [previous_pt3d[axis_ind]+ sim.net.cells[cell_gids[ind]].tags[axis]  for axis_ind, axis in enumerate(['x','y','z'])]
                        
                        ax.plot([previous_pt3d_[0],pt3d_[0]],[previous_pt3d_[1],pt3d_[1]],[previous_pt3d_[2],pt3d_[2]],linewidth=pt3d[3],color=color, alpha=0.7)

        ax.plot([post_point[0],pre_point[0]],[post_point[1],pre_point[1]],[post_point[2],pre_point[2]],color='limegreen',linewidth=1)

        ax.set_xlabel('X Label');ax.set_ylabel('Y Label');ax.set_zlabel('Z Label')
        # ax.view_init(90, 0)
        # ax.set_xlim(-150,150)
        # ax.set_ylim(-200,200)
        ax.set_axis_off()

        plt.savefig('../conn/gap_connecivity/figs/test_cell_spatial_conn.png',dpi=500)

    '''
    # Example usage:
    spatial_coordinates = { ... }  # Your dataset here
    max_distance = 5  # Define the maximum distance allowed between points
    # max_distance = 1.2  # Define the maximum distance allowed between points

    # Find pairs of points within the specified distance
    pairs_within_distance, sections_within_distance = find_points_within_distance(spatial_coordinates, max_distance)

    # Print the pairs of points within the specified distance
    print("Pairs of points within distance:")
    for pair in pairs_within_distance:
        print(pair)
    '''


    # def filter_points_within_distance(points, target_point, max_distance):
    #     """
    #     Filter points within a certain distance of a target point in 3D space using a KD-tree.
        
    #     Args:
    #     - points (list): A list of tuples, each containing the XYZ coordinates of a point.
    #     - target_point (tuple): A tuple containing the XYZ coordinates of the target point (x, y, z).
    #     - max_distance (float): The maximum distance allowed between points.
        
    #     Returns:
    #     - list: A list of tuples containing the XYZ coordinates of points within the specified distance.
    #     """
    #     tree = KDTree(points)
    #     indices = tree.query_ball_point(target_point, max_distance)
    #     filtered_points = [points[i] for i in indices]
    #     return filtered_points

    # # Example usage:
    # points = [(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)] # Example points
    # target_point = (3, 4, 5)
    # max_distance = 5
    # filtered_points = filter_points_within_distance(points, target_point, max_distance)
    # print("Points within the specified distance:", filtered_points)



''''
spatial_coordinates={
                        0:{
                            'soma_0':[
                                [489.71476598986976, 3898.3813194746554, 522.1686024262397, 6.03573751449585],
                                [495.75050326594703, 3898.3813194746554, 522.1686024262397, 6.03573751449585]
                                ],
                            'dend_0':[
                                [502.135255301088,   3882.8267979139864, 531.5313577248542, 2.6500000953674316],
                                [504.4782356859147,  3877.009128045517,  533.0508274628608, 2.2799999713897705],
                                [505.1860255838334,  3875.2282280439913, 533.5109877182929, 2.430000066757202],
                                [505.74948545703285, 3873.689528893906,  533.9193777634589, 2.430000066757202]
                                ],
                            'dend_1':[
                                [505.79562565097206, 3872.553729485947,  533.7544374062506, 2.6500000953674316],
                                [505.82886549243324, 3871.7355284208834, 533.6356277062384, 2.5],
                                ]
                            },
                        1:{
                            'soma_0':[
                                [492.7326345086991,  3898.3813194746554, 522.1686024262397, 6.03573751449585],
                                [495.75050326594703, 3898.3813194746554, 522.1686024262397, 6.03573751449585]
                                ],
                            'dend_4':[
                                [502.135255301088,   3882.8267979139864, 531.5313577248542, 2.6500000953674316],
                                [502.38873573550575, 3882.266868066269,  531.692357499787,  2.6500000953674316],
                                [502.84988542804115, 3881.1064981932223, 531.9921774460761, 2.5],
                                [503.50111576327674, 3879.467828225571,  532.4155673577277, 2.359999895095825],
                                ],
                            'dend_6':[
                                [505.74948545703285, 3873.689528893906,  533.9193777634589, 2.430000066757202],
                                [505.74948545703285, 3873.689528893906,  533.9193777634589, 2.6500000953674316],
                                [505.82886549243324, 3871.7355284208834, 533.6356277062384, 2.5],
                                ]
                            },
                        2:{
                            'soma_0':[
                                [489.71476598986976, 3898.3813194746554, 522.1686024262397, 6.03573751449585],
                                [492.7326345086991,  3898.3813194746554, 522.1686024262397, 6.03573751449585],
                                ],
                            'dend_8':[
                                [502.135255301088,   3882.8267979139864, 531.5313577248542, 2.6500000953674316],
                                [502.38873573550575, 3882.266868066269,  531.692357499787,  2.6500000953674316],
                                [504.4782356859147,  3877.009128045517,  533.0508274628608, 2.2799999713897705],
                                [505.1860255838334,  3875.2282280439913, 533.5109877182929, 2.430000066757202],
                                ],
                            'dend_9':[
                                [505.74948545703285, 3873.689528893906,  533.9193777634589, 2.430000066757202],
                                [505.74948545703285, 3873.689528893906,  533.9193777634589, 2.6500000953674316],
                                ]
                            },
                        }
'''
import math
import numpy as np
from NetPyNE_BBP import Utils as NPUtils

class SynapticPruning:

    def pruneGaps(sim,prune_pops=['TRN__pop'],prune_vals=[8],verbose=False):
        for prune_pop_ind,prune_pop in enumerate(prune_pops):
            prune_val = prune_vals[prune_pop_ind]
            print('\t>> Prunning Gap Junctions in pop ', prune_pop)
            for cell_ind in range(len(sim.net.cells)):
                if sim.net.cells[cell_ind].tags['pop'] == prune_pop:
                    cell_conns = sim.net.cells[cell_ind].conns
                    if verbose: print('\t\t>Prunning cell #: ', cell_ind, ' - ', len(cell_conns), ' total connections')
                    gap_conns=0
                    conns_left=0
                    pruned_conns=[]
                    for cell_conn_ind,cell_conn in enumerate(cell_conns):
                        if '|gap' in cell_conn['label'] and 'esyn' in cell_conn['synMech']: 
                        # if '|elec' in cell_conn['label'] and 'gap_nr' in cell_conn['synMech']: 
                            gap_conns+=1
                            if conns_left<prune_val:
                                conns_left+=1
                                pruned_conns.append(cell_conn)
                            else:
                                if verbose: print('skipping conn #: ', cell_conn_ind,'\t',cell_conn)
                                continue
                        else:
                            pruned_conns.append(cell_conn)
                    sim.net.cells[cell_ind].conns=pruned_conns
                    if verbose: print('\tleft ', len(sim.net.cells[cell_ind].conns), ' connections | of which ', conns_left, ' are gap conns | previous ', gap_conns)
        return sim

    def pruneGaps_byDistance(sim,prune_pops=['TRN__pop'],prune_vals=[8],verbose=False):
        for prune_pop_ind,prune_pop in enumerate(prune_pops):
            prune_val = prune_vals[prune_pop_ind]
            print('\t>> Prunning Gap Junctions in pop ', prune_pop)
            cell_path_dists={}

            cell_path_dists={cell_ind:NPUtils.pathDistance(sim.net.cells[cell_ind].__dict__) for cell_ind in range(len(sim.net.cells)) if 'secs' in sim.net.cells[cell_ind].__dict__.keys()}


            for cell_ind in range(len(sim.net.cells)):
                cell_path_dists.update({cell_ind:{}})
                if sim.net.cells[cell_ind].tags['pop'] == prune_pop:
                    cell_conns = sim.net.cells[cell_ind].conns
                    pre_gids = [cell_conn['preGid'] for cell_conn in cell_conns]

                    cell_path_dists[cell_ind].update({pre_gid:NPUtils.pathDistance(sim.net.cells[pre_gid].__dict__) for pre_gid in pre_gids})

                    NPUtils.pathDistance



                    if verbose: print('\t\t>Prunning cell #: ', cell_ind, ' - ', len(cell_conns), ' connections')
                    conns_left=0
                    pruned_conns=[]
                    for cell_conn_ind,cell_conn in enumerate(cell_conns):
                        if '|elec' in cell_conn['label'] and 'gap_nr' in cell_conn['synMech']: 
                            if conns_left<prune_val:
                                conns_left+=1
                                pruned_conns.append(cell_conn)
                            else:
                                if verbose: print('skipping conn #: ', cell_conn_ind,'\t',cell_conn)
                                continue
                        else:
                            pruned_conns.append(cell_conn)
                    sim.net.cells[cell_ind].conns=pruned_conns
                    if verbose: print('\tleft ', len(sim.net.cells[cell_ind].conns), ' connections')
        return sim
    
class QuantifyESyns:
    def quantifyGap(sim,plotFig=False):
        gap_mech = 'esyn'
        # gap_mech = sim.cfg.TRNe_TRNe__syn_mech
        TRN_pops    = [pop for pop in sim.net.pops.keys() if 'TRN' in pop]
        TRN_cells   = [gid for pop in TRN_pops for gid in sim.net.pops[pop].cellGids]

        # --- Gets the presynaptic cells that connect using gap junctions to a given cell
        pre_gids_dict_unidirectional={gid:[] for gid in TRN_cells}
        for gid in TRN_cells:
            for conn in sim.net.cells[gid].conns:
                if gap_mech == conn['synMech']:
                    pre_gids_dict_unidirectional[gid].append(conn['preGid'])
        
        # --- Matches the postsynaptic cells that are connected to the given cell (once netpyne only accounts for the presynaptic cells// A->B, but no B->, being B the GID cell)
        #  -  Creates a copy of <pre_gids_dict_unidirectional> that will store the bidirectional connections
        pre_gids_dict={gid:pre_gids_dict_unidirectional[gid] for gid in pre_gids_dict_unidirectional.keys()}
        
        for gid in pre_gids_dict.keys():
            for pre_gid in pre_gids_dict[gid]:
                if pre_gid not in pre_gids_dict.keys():  pre_gids_dict.update({pre_gid:[]}) # --- Making sure all <pre_gid> values are also <gid> keys
                
                if gid not in pre_gids_dict[pre_gid]:    pre_gids_dict[pre_gid].append(gid) # --- Adding <gid> to <pre_gid> if not present, to create bidirectional matching
        
        if plotFig: QuantifyESyns.plotGapConns(sim,pre_gids_dict)
    
        print([len(pre_gids_dict[gid]) for gid in pre_gids_dict.keys()])

        return pre_gids_dict
    
    def plotGapConns(sim,pre_gids_dict):
        # --- Scatter plot of pre and post cell positions
        pre_cell_xyz=[]
        for cell_gid in pre_gids_dict.keys():      pre_cell_xyz.append((sim.net.cells[cell_gid].tags['x'],sim.net.cells[cell_gid].tags['y'],sim.net.cells[cell_gid].tags['z']))
        
        # --- Extracting first and last elements of each tuple
        cell_x_values  = [item[0] for item in pre_cell_xyz]
        cell_y_values  = [item[1] for item in pre_cell_xyz]
        cell_z_values  = [item[2] for item in pre_cell_xyz]
        
        store_coords=[]
        store_coords_dict={}
        for gid in pre_gids_dict.keys():            
            store_coords_dict.update({gid:[]})
            cell_xyz = [sim.net.cells[gid].tags['x'],sim.net.cells[gid].tags['y'],sim.net.cells[gid].tags['z']]
            for pre_gid in pre_gids_dict[gid]:
                if gid in pre_gids_dict[pre_gid]:   linecolor = 'r' # checks if present in both lists
                else:                               linecolor = 'b'
                pre_cell_xyz = [sim.net.cells[pre_gid].tags['x'],sim.net.cells[pre_gid].tags['y'],sim.net.cells[pre_gid].tags['z']]
                x_coords = [cell_xyz[0],pre_cell_xyz[0]]
                y_coords = [cell_xyz[1],pre_cell_xyz[1]]
                z_coords = [cell_xyz[2],pre_cell_xyz[2]]
                store_coords.append([x_coords,y_coords,z_coords,linecolor])
                store_coords_dict[gid].append([x_coords,y_coords,z_coords,linecolor])
        
        linewidth = 1
        linealpha = 0.05

        import matplotlib.pyplot as plt
        # --- Creating a scatter plot of cell positions
        plt.figure(figsize=(30, 30))

        plt.subplot(2,2,1)
        plt.title('Scatter plot of Gap junction connections (x-y plane)')
        for x_coords,y_coords,z_coords,linecolor in store_coords: plt.plot(x_coords,y_coords,linewidth=linewidth,alpha=linealpha,color=linecolor)
        plt.scatter(cell_x_values, cell_y_values,   color='k', alpha=1.0, s=20)
        plt.xlim([375,625])
        plt.ylim([2500,2850])
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.gca().invert_yaxis()


        plt.subplot(2,2,2)
        plt.title('Scatter plot of Gap junction connections (x-z plane)')
        for x_coords,y_coords,z_coords,linecolor in store_coords: plt.plot(x_coords,z_coords,linewidth=linewidth,alpha=linealpha,color=linecolor)
        plt.scatter(cell_x_values, cell_z_values,   color='k', alpha=1.0, s=20)
        plt.xlim([375,625])
        plt.ylim([375,625])
        plt.xlabel('X-axis')
        plt.ylabel('Z-axis')
        plt.gca().invert_yaxis()

        plt.subplot(2,2,3)
        plt.title('Histogram of connection distances (in 3d)')
        lengths = Utils.calculate_vector_lengths(store_coords)
        plt.hist(lengths,bins=16)
        plt.xlabel('Count')
        plt.ylabel('Soma-Soma distance')

        # plt.subplot(2,2,4)
        # plt.title('Histogram of connection distances per cell (in 3d)')
        # store_val=[]
        # for gid in store_coords_dict.keys(): 
        #     gid_lengths = Utils.calculate_vector_lengths(store_coords_dict[gid])
        #     hist=plt.hist(gid_lengths,bins=16)
        #     store_val.append(val)

        plt.hist(lengths,bins=16)
        plt.xlabel('Count')
        plt.ylabel('Soma-Soma distance')

        plt.savefig('test_scatter_gap.png', dpi=300)

class Utils:

    def calculate_vector_lengths(store_coords):
        """
        Calculate the length of a 3D vector given two XYZ coordinate points.
        
        Args:
        - point1 (tuple): A tuple containing the XYZ coordinates of the first point (x1, y1, z1).
        - point2 (tuple): A tuple containing the XYZ coordinates of the second point (x2, y2, z2).
        
        Returns:
        - float: The length of the 3D vector between the two points.
        """
        
        return [math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2) for [[x1,x2],[y1,y2],[z1,z2],c] in store_coords]
        
        # for [[x1,x2],[y1,y2],[z1,z2],c] in store_coords:
        #     lengths.append(math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2))

        # x1, y1, z1 = point1
        # x2, y2, z2 = point2
        
        # length = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        
        # return length


    def distance(point1, point2):
        return sum((p1 - p2) ** 2 for p1, p2 in zip(point1, point2)) ** 0.5

    def lengths(store_coords):
        return [sum(Utils.distance(coord, store_coords[i][-1]) for coord in store_coords[i][:-1]) for i in range(len(store_coords))]

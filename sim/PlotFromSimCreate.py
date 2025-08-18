import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

class PlotFromSimCreate():
    def make_segments(x, y):
        '''
        Create list of line segments from x and y coordinates, in the correct format for LineCollection:
        an array of the form   numlines x (points per line) x 2 (x and y) array
        '''

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        return segments

    def colorline(x, y, z=None, cmap_name='viridis', norm_range=[0.0,1.0], linewidth=3, alpha=1.0):
        '''
        Plot a colored line with coordinates x and y
        Optionally specify colors in the array z
        Optionally specify a colormap, a norm function and a line width
        '''
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(cmap_name)
        
        norm=plt.Normalize(norm_range[0],norm_range[1])

        # Default colors equally spaced on [0,1]:
        if z is None:
            z = np.linspace(0.0, 1.0, len(x))
            
        # Special case if a single number:
        if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
            z = np.array([z])
            
        z = np.asarray(z)
        
        segments = PlotFromSimCreate.make_segments(x, y)
        lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
        
        ax = plt.gca()
        ax.add_collection(lc)

        return lc
    
    def colorline_byAlpha(x, y, z=None, cmap_name='viridis', linewidth=3, alpha=1.0):
        '''
        Plot a colored line with coordinates x and y
        Optionally specify colors in the array z
        Optionally specify a colormap, a norm function and a line width
        '''
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(cmap_name)

        rgba = cmap(alpha)
        return rgba
    
    def getCellProperties(sim,show_pops=None):
        x_=[]
        y_=[]
        z_=[]
        c_=[]
        for cell_ind in range(len(sim.net.cells)):
            # --- get cell tags
            cell_tags=sim.net.cells[cell_ind].tags

            # --- skips if the cell is from a pop that should not be plotted in the scatter
            if show_pops is not None:
                # if cell_tags['pop'].split('__pop')[0] not in show_pops: continue 
                if cell_tags['pop'].split('__pop')[0] not in show_pops: 
                    pop_name = cell_tags['pop'].split('__pop')[0]
                    if pop_name.split('@')[0] not in show_pops: continue 
            
            x_.append(cell_tags['x'])
            y_.append(cell_tags['y'])
            z_.append(cell_tags['z'])
            
            if   'L6A_activated'    in cell_tags['pop']: c_.append('firebrick')
            elif 'L6A_silent'       in cell_tags['pop']: c_.append('grey')
            elif 'L6A_sparse'       in cell_tags['pop']: c_.append('red')
            elif 'L6A_suppressed'   in cell_tags['pop']: c_.append('lightcoral')
            # if   'L6A_activated'    in cell_tags['pop']: c_.append('cyan')
            # elif 'L6A_silent'       in cell_tags['pop']: c_.append('limegreen')
            # elif 'L6A_sparse'       in cell_tags['pop']: c_.append('magenta')
            # elif 'L6A_suppressed'   in cell_tags['pop']: c_.append('y')
            elif 'L6CC'             in cell_tags['pop']: c_.append('mediumturquoise')
            elif 'L6IN'             in cell_tags['pop']: c_.append('gold')
            elif 'VPM'              in cell_tags['pop']: c_.append('g')
            elif 'TRN'              in cell_tags['pop']: 
                if cell_tags['pop']=='TRN__pop': c_.append('b')
                else                      : c_.append('lightskyblue')
            elif 'MLe'              in cell_tags['pop']: c_.append('k')
            elif 'CTvirtual'        in cell_tags['pop']: c_.append('red')
            else                                       : c_.append('grey')
        return x_, y_, z_, c_

    def plotScatter(sim):
        import matplotlib.pyplot as plt
        x_, y_, z_, c_ = PlotFromSimCreate.getCellProperties(sim)
        
        # create an axis instance with subplots
        fig, ax = plt.subplots(1,1,figsize=(10,25))
        scatter = ax.scatter(x_, y_, c=c_, alpha=1) #here I need to invert
        ax.invert_yaxis()   # here how you can invert
        plt.savefig('scatter_L6A_diversity.png')

    def plot2dNet(sim,select_sources=None,select_targets=None,plotEvery=1,singleCell=False,singleConn=False,show_pops='all',fig_projection='3d',show_range=None, plotConns=True, cells_alpha=0.25, fig_size=(10,10), dpi=100):
        
        plotted_conns=0
        import matplotlib.pyplot as plt
        if show_pops == 'all': show_pops_=[pop.split('__pop')[0] for pop in sim.net.pops.keys()]
        else: show_pops_=show_pops
        x_, y_, z_, c_ = PlotFromSimCreate.getCellProperties(sim,show_pops=show_pops_)

        if select_sources is not None:  prePops     = 'pre_'+ '_'.join(select_sources)
        else:                           prePops     = 'all_'
        if select_sources is not None:  postPops    = 'post_'+ '_'.join(select_targets)
        else:                           postPops    = 'all_'

        print('Plotting ', fig_projection, ' net from ', prePops, ' to ', postPops)

        ######### plot conn from sim obj
        # fig, ax = plt.subplots(1,1,figsize=(10,25))
        if      fig_projection == '3d':     
            ax = plt.figure(figsize=fig_size).add_subplot(projection='3d')
            ax.invert_zaxis()   # here how you can invert --- obs: this is the PYPLOT axis reference, where Z=height
            ax.scatter(x_, z_, y_, c=c_, alpha=cells_alpha) #here I need to invert

            # --- Invisible y point to fix aspect ratio of plot
            size_y    = (max(y_)-min(y_))
            x1_point = 500+(size_y/2)
            x2_point = 500-(size_y/2)
            z1_point = x1_point
            z2_point = x2_point
            ax.scatter(x1_point, z1_point, max(y_), c='k', alpha=cells_alpha) #here I need to invert
            ax.scatter(x2_point, z2_point, min(y_), c='k', alpha=cells_alpha) #here I need to invert

            view = '3d'
        if      fig_projection == '3d_special':     
            ax = plt.figure(figsize=fig_size).add_subplot(projection='3d')
            ax.invert_zaxis()   # here how you can invert --- obs: this is the PYPLOT axis reference, where Z=height
            ax.scatter(x_, z_, y_, c=c_, alpha=cells_alpha) #here I need to invert

            # --- Invisible y point to fix aspect ratio of plot
            size_y    = (max(y_)-min(y_))
            x1_point = 500+(size_y/2)
            x2_point = 500-(size_y/2)
            z1_point = x1_point
            z2_point = x2_point
            ax.scatter(x1_point, z1_point, max(y_), c='k', alpha=0) #here I need to invert
            ax.scatter(x2_point, z2_point, min(y_), c='k', alpha=0) #here I need to invert

            # # make the panes transparent
            # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.set_axis_off()

            view = '3d'
        elif    fig_projection == '2d_top': 
            ax = plt.figure(figsize=fig_size).add_subplot(projection='3d')
            ax.invert_zaxis()   # here how you can invert --- obs: this is the PYPLOT axis reference, where Z=height
            # ax.view_init(azim=0, elev=90)
            ax.scatter(x_, z_, y_, c=c_, alpha=cells_alpha) #here I need to invert
            ax.view_init(azim=-90, elev=90)
            view = '3d'
        elif    fig_projection == '2d':                       
            ax = plt.figure(figsize=fig_size)
            plt.gca().invert_yaxis()
            plt.scatter(x_, y_, c=c_, alpha=cells_alpha) #here I need to invert
            view = '2d'
        else:                               
            ax = plt.figure(figsize=fig_size)
            plt.gca().invert_yaxis()
            plt.scatter(x_, y_, c=c_, alpha=cells_alpha) #here I need to invert
            view = '2d'
        
        # scatter = ax.scatter(x_, y_, z_, c=c_, alpha=0.25) #here I need to invert

        # scatter = ax.scatter(x_, y_, c=c_, alpha=1) #here I need to invert
        # ax.invert_yaxis()   # here how you can invert --- obs: this is the PYPLOT axis reference, where Z=height
        
        if plotConns:
            store_weights=[]
            for cell in sim.net.cells:
                conns=cell.conns
                for conn in conns:
                    store_weights.append(conn['weight'])
            if len(store_weights)>0: max_weight=max(store_weights)

            basewidth=2

            if select_targets is not None:
                selected_pops=[]
                for pop in sim.net.pops.keys():
                    for pathway in select_targets:
                        if pathway in pop:
                            selected_pops.append(pop)
            else:selected_pops=sim.net.pops.keys() # selects all pops

            cellGids=[]
            for pop in selected_pops:cellGids+=sim.net.pops[pop].cellGids
            # for pop in selected_pops:cellGids.append(sim.net.pops[pop].cellGids)

            gid_depth=[]
            for gid in cellGids:
                post_y   = sim.net.cells[gid].tags['y']
                gid_depth.append((post_y,gid))

            for (depth,gid) in gid_depth:
                if gid%plotEvery!=0: continue
                post_x   = sim.net.cells[gid].tags['x']
                post_y   = sim.net.cells[gid].tags['y']
                post_z   = sim.net.cells[gid].tags['z']
                postPop = sim.net.cells[gid].tags['pop']
                conns    = sim.net.cells[gid].conns
                for conn in conns:
                    if type(conn['preGid'])!=int: continue
                    else:preCell=sim.net.cells[conn['preGid']]

                    prePop=sim.net.cells[conn['preGid']].tags['pop']
                    # --- skips unselected pops
                    if select_sources is not None:
                        for source in select_sources:
                            if source in prePop:

                                pre_x   = preCell.tags['x']
                                pre_y   = preCell.tags['y']
                                pre_z   = preCell.tags['z']
                                # print('pre pop:', prePop,'\t post pop:', postPop)

                                if 'MLe' in prePop:
                                    pre_x_range = [500-55,500+55]
                                    pre_y_range = [4500,5700]
                                    pre_z_range = [500-55,500+55]
                                else:
                                    pre_x_range = sim.net.pops[prePop].tags['xRange']
                                    pre_y_range = sim.net.pops[prePop].tags['yRange']
                                    pre_z_range = sim.net.pops[prePop].tags['zRange']

                                center_point=500
                                if 'L6A' in prePop: 
                                    angle=np.remainder((((np.arctan2(pre_z - center_point, pre_x - center_point))*(180/np.pi))+360),360)#-((%s-%f)/(%f-%f))
                                    alpha_gradient = angle/360
                                    # if 'activated' in prePop:   
                                    #     angle=np.remainder((((np.arctan2(pre_z - center_point, pre_x - center_point))*(180/np.pi))+360),360)#-((%s-%f)/(%f-%f))
                                    #     alpha_gradient = angle/360
                                    # else:                       
                                    #     alpha_gradient = (pre_y-pre_y_range[0])/(pre_y_range[1]-pre_y_range[0]) # --- gradient over the y-axis
                                else:               alpha_gradient = (pre_y-pre_y_range[0])/(pre_y_range[1]-pre_y_range[0])

                                if show_range is not None:
                                    if alpha_gradient<show_range[0] or alpha_gradient>show_range[1]: 
                                        # print(alpha_gradient)
                                        continue # --- skips plotting if the cell is outside of the selected range

                                if   'exc' in conn['synMech']:  cmap_name = 'hsv'
                                elif 'inh' in conn['synMech']:  cmap_name = 'twilight'
                                else:                           cmap_name = 'Greys'
                                # if   'exc' in conn['synMech']:  cmap_name = 'Reds'
                                # elif 'inh' in conn['synMech']:  cmap_name = 'Blues'
                                # else:                           cmap_name = 'Greys'

                                n_dots=100
                                x = np.linspace(pre_x, post_x, n_dots)
                                y = np.linspace(pre_y, post_y, n_dots)
                                z = np.linspace(pre_z, post_z, n_dots)
                                
                                cmap = plt.get_cmap(cmap_name)
                                rgba = cmap(alpha_gradient)
                                if len(store_weights)>0: 
                                    if view == '3d':    plt.plot([pre_x,post_x],[pre_z,post_z],[pre_y,post_y],c=rgba[0:3],linewidth=basewidth*(conn['weight']/max_weight))
                                    elif view == '2d':  plt.plot([pre_x,post_x],[pre_y,post_y],c=rgba[0:3],linewidth=basewidth*(conn['weight']/max_weight))
                                    # plt.plot([pre_x,post_x],[pre_y,post_y],[pre_z,post_z],c=rgba[0:3],linewidth=basewidth*(conn['weight']/max_weight))
                                
                                plotted_conns+=1

                                if singleConn:continue
                if singleCell:continue

        if plotConns: 
            conns_num = PlotFromSimCreate.countAllConns(sim)
            plt.title('Showing \n%f percent \nof all conns'%((plotted_conns/conns_num)*100))
        # center_point = 500
        # ax.set_xlim([center_point-(6000/2),center_point+(6000/2)])
        # ax.set_ylim([center_point-(6000/2),center_point+(6000/2)])
        # ax.set_zlim([0,6000])
        saveFigName = '2dnet_'+prePops+'_'+postPops
        if singleCell:                  saveFigName+='_singleCell'
        if singleConn:                  saveFigName+='_singleConn'
        if plotEvery!=1:                saveFigName+='_every|'+str(plotEvery)
        if type(fig_projection) is str: saveFigName+='_'+fig_projection
        if show_range is not None:      saveFigName+='_'+str(int(show_range[0]*100))+'to'+str(int(show_range[1]*100))+'percent'
        plt.xlabel('x axis')
        # plt.ylabel('y axis')
        # if fig_projection == '3d_special': plt.savefig('figs_network/'+saveFigName+'.png',dpi=1000)
        # else: plt.savefig('figs_network/'+saveFigName+'.png',dpi=500)
        plt.savefig('figs_network/'+saveFigName+'.png',dpi=dpi)

    def countAllConns(sim):
        conns_num=0
        for cell in sim.net.cells: conns_num+=len(cell.conns)
        return conns_num
    
    def plotNetPyNEConnMatrix(sim):
        if not sim.cfg.removeConns:
            sim.net.allCells    = [cell.__getstate__() for cell in sim.net.cells]
            sim.net.allPops     = {label: pop.__getstate__() for label, pop in sim.net.pops.items()}

            conn_features=['strength','weight','convergence','divergence']
            for feature in conn_features:
                all_pops=[pop for pop in sim.net.params.popParams.keys() if 'MLe' not in pop]
                print('Plotting all pops conn mat')
                sim.analysis.plotConn(includePre=all_pops,includePost=all_pops,feature=feature,saveFig='conn_mat_test_allPops__'+feature+'.png',figSize=(20,20))
                if sim.cfg.simulateL6only:
                    print('Plotting CTX conn mat')
                    pops_ctx=[pop for pop in sim.net.params.popParams.keys() if 'L6' in pop]
                    sim.analysis.plotConn(includePre=pops_ctx,includePost=pops_ctx,feature=feature,saveFig='conn_mat_test_ctxPops__'+feature+'.png',figSize=(20,20))
                if sim.cfg.simulateThalonly:
                    print('Plotting Thal conn mat')
                    pops_thal=[pop for pop in sim.net.params.popParams.keys() if ('VPM' in pop) or ('TRN' in pop)]
                    sim.analysis.plotConn(includePre=pops_thal,includePost=pops_thal,feature=feature,saveFig='conn_mat_test_thalPops__'+feature+'.png',figSize=(20,20))
        else: print('Skipping conn mat plot because removeConns is', sim.cfg.removeConns)
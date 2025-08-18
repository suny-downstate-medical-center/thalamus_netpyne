"""
Script: ring_population_script.py

Description:
This script defines a class 'CreateCellPositions' that generates 3D coordinates
for a ring-shaped population of neurons. The ring is divided into multiple regions,
each specified by intervals from the center and relative density. The class also
provides a method for example usage, visualizing the generated cell positions in 3D
and 2D views.

Author: Joao Moreira - joao.moreira@downstate.edu
Date:   January 2024
"""

from neuron import h
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math

import matplotlib
matplotlib.rcParams.update({'font.size': 18})

class CreateCellPositions():

    @staticmethod
    def ring_population(num_cells, regions, relative_density, center_point=500, y_interval=[890, 1090], seed=None, save_fig=None, map_color='viridis',map_color_range=[0,1]):
        """
        Generate 3D coordinates for a ring-shaped population of neurons with variable density.

        Parameters:
        - num_cells (int):          Total number of cells in the population.
        - regions (list of lists):  Intervals describing the distance from the center_point for each region.
        - relative_density (list):  Relative density for each region. Should add up to 1.
        - center_point (int):       Center point for the ring along the x and z axes.
        - y_interval (list):        Interval of y values for the cell positions.
        - save_fig (str):           File name to save the generated figure (optional).
        - seed (int):               Seed value for the random number generator.
        - map_color (str):          Colormap for coloring cells along the y-axis.

        Returns:
        - list of tuples:           Generated 3D coordinates for the cells.
        """

        # Initialize NEURON random number generator
        if seed is not None: rand = h.Random(seed)

        # Validate input lengths
        if len(regions) != len(relative_density):
            raise ValueError("Length of 'regions' and 'relative_density' must be the same.")

        # Normalize relative_density to ensure it adds up to 1
        relative_density = np.array(relative_density) / sum(relative_density)

        # Initialize arrays to store coordinates
        x = np.empty(0)
        y = np.empty(0)
        z = np.empty(0)

        # Creates an array to store the indexes of the extra cells added to prevent rounding error
        extra_cell_indexes=[0]
        
        # Loop through each region to generate positions
        for i in range(len(regions)):

            num_cells_in_region = int(num_cells * relative_density[i]) + 1 # --- adding an extra cell (+1) to each region to avoid rounding issues (should not be used if # of regions is too large)
            # num_cells_in_region = int(num_cells * relative_density[i])
            extra_cell_indexes.append(num_cells_in_region+extra_cell_indexes[-1])
            # Generate positions for the current region using NEURON's Random
            theta_region = [rand.uniform(0, 2 * np.pi) for _ in range(num_cells_in_region)]

            # Use h.Random().uniform to generate random numbers from a uniform distribution
            min_radius_region = regions[i][0]
            max_radius_region = regions[i][1]
            
            # Compensate for the scaling of the area as the radius increases
            radius_region = [rand.uniform(min_radius_region**2, max_radius_region**2)**0.5 for _ in range(num_cells_in_region)]

            x_region = [center_point + radius * h.cos(theta) for radius, theta in zip(radius_region, theta_region)]
            y_region = [rand.uniform(y_interval[0], y_interval[1]) for _ in range(num_cells_in_region)]
            z_region = [center_point + radius * h.sin(theta) for radius, theta in zip(radius_region, theta_region)]

            # Concatenate positions for the current region
            x = np.concatenate([x, x_region])
            y = np.concatenate([y, y_region])
            z = np.concatenate([z, z_region])
        
        del extra_cell_indexes[0]
        # print('cells_in_region ', extra_cell_indexes)

        # Converts to list so that elements can be deleted below
        x = x.tolist()
        y = y.tolist()
        z = z.tolist()

        # Trim or duplicate cells to ensure the total number of cells is exactly equal to num_cells - removes the excess cells added to prevent 
        excess_cells = len(x) - num_cells
        if excess_cells > 0:
            # Removes extra cells -  one from each region, until the length of x matches num_cells
            while excess_cells > 0:
                del x[extra_cell_indexes[-1]-1]
                del y[extra_cell_indexes[-1]-1]
                del z[extra_cell_indexes[-1]-1]
                del extra_cell_indexes[-1]
                excess_cells-=1
        elif excess_cells < 0:
            print('Missing cells - please check the code')
            print('The number of cells is less than the specified number in num_cells: ', num_cells)
            return
            # x = np.concatenate([x, x[:excess_cells]])
            # y = np.concatenate([y, y[:excess_cells]])
            # z = np.concatenate([z, z[:excess_cells]])

        # # Trim or duplicate cells to ensure the total number of cells is exactly equal to num_cells - removes the excess cells added to prevent 
        # excess_cells = len(x) - num_cells
        # if excess_cells > 0:
        #     x = x[:-excess_cells]
        #     y = y[:-excess_cells]
        #     z = z[:-excess_cells]
        # elif excess_cells < 0:
        #     x = np.concatenate([x, x[:excess_cells]])
        #     y = np.concatenate([y, y[:excess_cells]])
        #     z = np.concatenate([z, z[:excess_cells]])

        cell_positions_3D = list(zip(x, y, z))

        if save_fig is not None: CreateCellPositions.plot_ring_population(center_point, regions, y_interval, cell_positions_3D, save_fig, map_color, map_color_range=[0,1], xy_view = 'XZ_distance')

        return cell_positions_3D

    # Extract x, y, and z coordinates for the updated case
    @staticmethod
    def plot_ring_population(center_point, regions, y_interval, cell_positions_3D, save_fig=None, map_color='heat', distance_plane='xz', map_color_range=[0,1], xy_view = 'trueXY'):
        x_updated_params, y_updated_params, z_updated_params = zip(*cell_positions_3D)

        # Calculates the maximum border defined by the regions
        max_radius = max([region[1] for region in regions])

        # Create a figure with 4 subplots
        fig = plt.figure(figsize=(20, 20))

        # Subplot 1: 3D view
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.scatter(x_updated_params, y_updated_params, z_updated_params, c=y_updated_params, cmap=map_color, vmin=map_color_range[0], vmax=map_color_range[1], marker='o',
                    label='Cell Positions')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')
        ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.set_title('Cell Positions in Ring')
        ax1.legend()
        ax1.invert_yaxis()


        # Subplot 2: 2D view from above
        ax2 = fig.add_subplot(222)
        ax2.scatter(x_updated_params, z_updated_params, c=y_updated_params, cmap=map_color, vmin=map_color_range[0], vmax=map_color_range[1], marker='o',
                    label='Cell Positions')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Z')
        ax2.set_title('2D View from Above')
        ax2.set_xlim([center_point-(max_radius*1.3), center_point+(max_radius*1.3)])
        ax2.set_ylim([center_point-(max_radius*1.3), center_point+(max_radius*1.3)])
        # Adding a circle to define the barrel diameter
        circle1 = plt.Circle((center_point, center_point), max_radius, edgecolor='k', facecolor='none', linestyle='dotted') # Create a circle with center (center_point, center_point) and radius max_radius
        circle2 = plt.Circle((center_point, center_point), 140, edgecolor='grey', facecolor='none', linestyle='dotted')        # Create a circle with center (center_point, center_point) and radius max_radius
        ax2.add_patch(circle1) # Add the circle to the axis
        ax2.add_patch(circle2) # Add the circle to the axis
        ax2.invert_yaxis()

        # Calculates the distance from each cell respective to center_point
        distance_from_center = CreateCellPositions.calculate_distance_from_center(cell_positions_3D, center_point, distance_plane)

        if xy_view == 'trueXY':
            # Subplot 3: 2D view from X-Y plane
            ax3 = fig.add_subplot(223)
            ax3.set_xlabel('X')
            ax3.set_ylabel('Y')
            ax3.set_title('2D View from X-Y Plane')
            ax3.set_xlim([center_point-(max_radius*1.3),    center_point+(max_radius*1.3)])
            line11 = plt.Line2D([center_point-(max_radius), center_point-(max_radius)], y_interval, linestyle='dotted', color='k')
            line12 = plt.Line2D([center_point+(max_radius), center_point+(max_radius)], y_interval, linestyle='dotted', color='k')
            line21 = plt.Line2D([center_point-140,          center_point-140],          y_interval, linestyle='dotted', color='grey')
            line22 = plt.Line2D([center_point+140,          center_point+140],          y_interval, linestyle='dotted', color='grey')
            ax3.add_line(line11)
            ax3.add_line(line12)
            ax3.add_line(line21)
            ax3.add_line(line22)
            ax3.scatter(x_updated_params, y_updated_params, c=y_updated_params, cmap=map_color, vmin=map_color_range[0], vmax=map_color_range[1], marker='o',
                        label='Cell Positions')
            ax3.invert_yaxis()
        else:
            # Subplot 3: 2D view from X-Y plane
            ax3 = fig.add_subplot(223)

            # Create an Arrow with a dotted line
            from matplotlib.patches import FancyArrow
            for ind,x in enumerate(distance_from_center):
                arrow = FancyArrow(0, y_updated_params[ind], x, 0, 
                                #    linestyle='dotted', 
                                   color='k', width=0.001,
                                   alpha=0.05)
                # Add the Arrow to the axis
                ax3.add_patch(arrow)
            # ax3.scatter(x_updated_params, y_updated_params, c=y_updated_params, cmap=map_color, marker='o',
            #             label='Cell Positions')
            ax3.set_xlabel('X')
            ax3.set_ylabel('Y')
            ax3.set_title('2D View cell distances to the center point')
            ax3.set_xlim([-10,    (max_radius*1.3)])
            line12 = plt.Line2D([(max_radius), (max_radius)], y_interval, linestyle='dotted', color='k')
            line22 = plt.Line2D([140,          140],          y_interval, linestyle='dotted', color='grey')
            ax3.add_line(line12)
            ax3.add_line(line22)
            ax3.scatter(distance_from_center, y_updated_params, c=y_updated_params, cmap=map_color, vmin=map_color_range[0], vmax=map_color_range[1], marker='o',
                        label='Cell Positions')
            ax3.invert_yaxis()
        
        # Subplot 4: Histogram in X-Z plane
        ax4 = fig.add_subplot(224)
        CreateCellPositions.plot_histogram(ax4, distance_from_center, max_radius, distance_plane)

        plt.tight_layout()
        if save_fig is None:
            save_fig = 'ring_network.png'
        plt.savefig(save_fig, dpi=500)

    @staticmethod
    def calculate_distance_from_center(cell_positions_3D, center_point, distance_plane):
        x_updated_params, _, z_updated_params = zip(*cell_positions_3D)

        if distance_plane == 'xz':
            x_dist = [x - center_point for x in x_updated_params]
            z_dist = [z - center_point for z in z_updated_params]
            distance_from_center = [np.sqrt(x_dist[ind] ** 2 + z_dist[ind] ** 2) for ind, x in enumerate(x_dist)]
        elif distance_plane == 'x':
            distance_from_center = [x - center_point for x in x_updated_params]
        elif distance_plane == 'z':
            distance_from_center = [z - center_point for z in z_updated_params]

        return distance_from_center

    @staticmethod
    def plot_histogram(ax, distance_from_center, max_radius, distance_plane, norm_hist=True):
        bins = 10
        hist, bin_edges = np.histogram(distance_from_center, bins=bins, range=(0, max_radius))
        # hist, bin_edges = np.histogram(distance_from_center, bins=bins, range=(-max_radius, max_radius))
        bin_width = bin_edges[1] - bin_edges[0]

        if norm_hist:
            # Normalize histogram
            norm_hist = []
            for col_ind, col in enumerate(hist):
                inner_radius = bin_edges[col_ind]
                outer_radius = bin_edges[col_ind+1]
                result_area = CreateCellPositions.ring_area(inner_radius, outer_radius)
                norm_hist.append(col/(result_area*1e-6))
            hist_data = norm_hist
            ax.set_ylabel('Cells/mm2')
        else:
            # Convert counts to percentages
            hist_percentage = (hist / len(distance_from_center)) * 100
            hist_data = hist_percentage
            ax.set_ylabel('Percentage of Total Cells')

        # Plot the histogram with percentages
        ax.bar(bin_edges[:-1], hist_data, width=bin_width, color='grey', edgecolor='black', alpha=0.7)
        if distance_plane == 'xz':  
            ax.set_xlim([-bin_width, max_radius*1.3])
            ax.set_title('Histogram of Cells in X-Z Plane')
        else:                       
            ax.set_xlim([-max_radius*1.3, max_radius*1.3])
            ax.set_title('Histogram of Cells in '+str(distance_plane)+' Plane')


    @staticmethod
    def ring_area(inner_radius, outer_radius):
        """
        Calculate the area of a 2D ring defined by two radii.

        Parameters:
        - inner_radius (float): Inner radius of the ring.
        - outer_radius (float): Outer radius of the ring.

        Returns:
        - float: Area of the ring.
        """
        area = np.pi * (outer_radius**2 - inner_radius**2)
        return area


if __name__ == '__main__':

    print('\t>> Creating cell positions')
    # Example usage with updated parameters
    num_cells = 500
    center_point = 500
    seed = 100000  # Specify your desired seed value

    ring_pop_dicts={
                # 'CT_pop':{
                #     'num_cells':            818,
                #     # 'num_cells':            num_cells,
                #     'regions':              [[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]],
                #     'relative_density':     [0.10, 0.15, 0.25, 0.25, 0.25],
                #     'center_point':         center_point,
                #     'y_interval':           [890, 1090],
                #     'seed':                 seed,
                #     'colormap':             'autumn_r',
                #     'colormap_vmin':       0.5,
                #     'colormap_vmax':       0.5,
                # },
                # 'CC_pop':{
                #     'num_cells':            634,
                #     # 'num_cells':            num_cells,
                #     'regions':              [[0.0, 25.0], [25.0, 50.0], [50.0, 75.0], [75.0, 125.0], [125.0, 175.0]],
                #     'relative_density':     [0.01, 0.025, 0.065, 0.3, 0.6],
                #     'center_point':         center_point,
                #     'y_interval':           [890, 1090],
                #     'seed':                 seed+1,
                #     'colormap':             'winter_r',
                #     'colormap_vmin':       0.5,
                #     'colormap_vmax':       0.5,
                # },
                # 'IN_pop':{
                #     'num_cells':            131,
                #     # 'num_cells':            num_cells,
                #     'regions':              [[0.0, 175.0]],
                #     'relative_density':     [1.0],
                #     'center_point':         center_point,
                #     'y_interval':           [890, 1090],
                #     'seed':                 seed+1,
                #     'colormap':             'copper',
                #     'colormap_vmin':       0.5,
                #     'colormap_vmax':       0.5,
                # },
                # 'UniformCT_pop':{
                #     'num_cells':            num_cells,
                #     'regions':              [[0, 350]],
                #     'relative_density':     [1],
                #     'center_point':         center_point,
                #     'y_interval':           [890, 1090],
                #     'seed':                 seed,
                #     'colormap':             'gist_heat',
                # },
                # 'UniformCC_pop':{
                #     'num_cells':            num_cells,
                #     'regions':              [[0, 350]],
                #     'relative_density':     [1],
                #     'center_point':         center_point,
                #     'y_interval':           [890, 1090],
                #     'seed':                 seed+1,
                #     'colormap':             'gist_heat',
                # },
                'UniformTRN_pop':{
                    'num_cells':            num_cells,
                    'regions':              [[100, 200]],
                    'relative_density':     [1],
                    'center_point':         center_point,
                    'y_interval':           [2550,2800],
                    'seed':                 seed+1,
                    'colormap':             'autumn_r',
                },
                'UniformTRN_pop2':{
                    'num_cells':            num_cells,
                    'regions':              [[50, 100]],
                    'relative_density':     [1],
                    'center_point':         center_point,
                    'y_interval':           [2550,2800],
                    'seed':                 seed+1,
                    'colormap':             'autumn_r',
                },
    }

    for pop_name_ind, pop_name in enumerate(ring_pop_dicts.keys()):

        print('\t\t>     ', pop_name)
        
        num_cells           = ring_pop_dicts[pop_name]['num_cells']
        regions             = ring_pop_dicts[pop_name]['regions']
        relative_density    = ring_pop_dicts[pop_name]['relative_density']
        center_point        = ring_pop_dicts[pop_name]['center_point']
        y_interval          = ring_pop_dicts[pop_name]['y_interval']
        seed                = ring_pop_dicts[pop_name]['seed']
        colormap            = ring_pop_dicts[pop_name]['colormap']
        try:    vmin        = ring_pop_dicts[pop_name]['colormap_vmin']
        except: vmin        = min(min(regions))
        try:    vmax        = ring_pop_dicts[pop_name]['colormap_vmax']
        except: vmax        = max(max(regions))

        ring_pop_dicts[pop_name].update({'cell_positions_3D':CreateCellPositions.ring_population(num_cells, regions, relative_density, center_point, y_interval, seed, save_fig='CreateCellPositions_'+str(pop_name_ind)+'_'+pop_name+'_plot_ring_network_seed_'+str(seed)+'.png',map_color=colormap,map_color_range=[vmin,vmax])})

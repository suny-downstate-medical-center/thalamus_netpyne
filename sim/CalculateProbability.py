import numpy as np
class CalculateProbability:
    def calculate_probability_of_connection(diameter, height, avg_connections, max_distance, numCells=551, density=None, seed = 100000):
        import numpy as np
        from scipy.spatial.distance import cdist
        from neuron import h

        # Initialize NEURON random number generator
        if seed is not None: rand = h.Random(seed)

        # Convert diameter and height to mm
        diameter_mm = diameter / 1000  # Convert micrometers to millimeters
        height_mm = height / 1000      # Convert micrometers to millimeters
        max_distance_mm = max_distance/1000

        if numCells is not None:
            total_cells = numCells
        else:        
            # Calculate the total number of cells
            volume = np.pi * (diameter_mm / 2)**2 * height_mm
            total_cells = density * volume

        # Generate random 3D coordinates for cells in the volume
        cell_coordinates_=[]
        for i in range(total_cells): cell_coordinates_.append([diameter_mm*rand.uniform(0,1), diameter_mm*rand.uniform(0,1), height_mm*rand.uniform(0,1)])
        cell_coordinates = np.array(cell_coordinates_)

        # cell_coordinates = np.random.rand(int(total_cells), 3) * np.array([diameter_mm, diameter_mm, height_mm])

        # Initialize connections count for each cell
        connections_count = np.zeros(int(total_cells))

        # Calculate distances between all pairs of cells
        distances = cdist(cell_coordinates, cell_coordinates)

        # Check distances and establish connections within the specified distance
        for i in range(int(total_cells)):
            for j in range(i + 1, int(total_cells)):
                if distances[i, j] <= max_distance_mm:
                # if distances[i, j] <= max_distance:
                    connections_count[i] += 1
                    connections_count[j] += 1

        # Calculate the average number of connections within the specified distance
        avg_connections_within_distance = np.mean(connections_count)

        # Calculate the probability of connection to achieve the desired average
        probability_of_connection = avg_connections / avg_connections_within_distance

        return probability_of_connection
    
    def verifyConns(sim,synsPerConn=1):
        TRN_cells = sum([sim.net.pops[pop].tags['numCells'] for pop in sim.net.pops.keys() if 'TRN' in pop])
        if TRN_cells==0:return
        count_conns=0
        store_conns_per_cell={}
        for cell in sim.net.cells:
            conns_per_cell=0
            for conn in cell.conns:
                if '|gap' in conn['label']:
                    count_conns+=1
                    conns_per_cell+=1
            
            if conns_per_cell>0:store_conns_per_cell.update({cell.gid:conns_per_cell})
            

        mean_conns=count_conns/TRN_cells
        mean_conns_perCell = mean_conns/synsPerConn
        print('mean_conns: ', mean_conns, '\t-\tmean_conns_perCell: ', mean_conns_perCell)
        try:
            print('max_conns: ', max(list(store_conns_per_cell.values())), '\t-\tmin_conns: ', min(list(store_conns_per_cell.values())))
            print('median_conns: ', np.median(list(store_conns_per_cell.values())))
        except:
            print('empty conns')
            pass

        # return store_conns_per_cell
    
if __name__ == '__main__':
    print('\t>> Creating connection probability using NEURON rng')
    # Example usage with updated values
    density         = 56160.40  # cells/mm^3
    diameter        = 200      # um
    height          = 250        # um
    avg_connections = 8.625
    max_distance    = 80   # μm

    probability = CalculateProbability.calculate_probability_of_connection(diameter, height, avg_connections, max_distance, numCells=551, density=None, seed = 100000)
    print(f"The probability of connection is approximately {probability:.10f}")
    
    # for i in range(100):
    #     probability = CalculateProbability.calculate_probability_of_connection(diameter, height, avg_connections, max_distance, numCells=100, density=None, seed = 100000)
    #     print(f"The probability of connection is approximately {probability:.10f}")

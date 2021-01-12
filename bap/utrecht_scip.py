import networkx as nx
import random

import numpy as np


def MaternProcess(lambdaParent = 10, lambdaDaughter = 5, radiusCluster = 5, plot=False):
    """
    Simulate a Matern point process on a rectangle.
    Author: H. Paul Keeler, 2018.
    https://github.com/hpaulkeeler/posts/blob/master/MaternClusterRectangle/MaternClusterRectangle.py
    https://hpaulkeeler.com/simulating-a-matern-cluster-point-process/
    
    PARAMETERS:
    ----------
    lambdaParent   # density of parent Poisson point process
    lambdaDaughter   # mean number of points in each cluster
    radiusCluster   # radius of cluster disk (for daughter points)
    """
    # Simulation window parameters
    xMin = 0;
    xMax = 100;
    yMin = 0;
    yMax = 100;

    # Extended simulation windows parameters
    rExt = radiusCluster;  # extension parameter -- use cluster radius
    xMinExt = xMin - rExt;
    xMaxExt = xMax + rExt;
    yMinExt = yMin - rExt;
    yMaxExt = yMax + rExt;
    # rectangle dimensions
    xDeltaExt = xMaxExt - xMinExt;
    yDeltaExt = yMaxExt - yMinExt;


    # Simulate Poisson point process for the parents
    numbPointsParent = np.random.poisson(lambdaParent);  # Poisson number of points
    
    # x and y coordinates of Poisson points for the parent
    xxParent = xMinExt + xDeltaExt * np.random.uniform(0, 1, numbPointsParent);
    yyParent = yMinExt + yDeltaExt * np.random.uniform(0, 1, numbPointsParent);

    # Simulate Poisson point process for the daughters (ie final poiint process)
    numbPointsDaughter = np.random.poisson(lambdaDaughter, numbPointsParent);
    numbPoints = sum(numbPointsDaughter);  # total number of points
    print(f'===================> {numbPoints} zones created')

    # Generate the (relative) locations in polar coordinates by
    # simulating independent variables.
    theta = 2 * np.pi * np.random.uniform(0, 1, numbPoints);  # angular coordinates
    rho = radiusCluster * np.sqrt(np.random.uniform(0, 1, numbPoints));  # radial coordinates

    # Convert from polar to Cartesian coordinates
    xx0 = rho * np.cos(theta);
    yy0 = rho * np.sin(theta);

    # replicate parent points (ie centres of disks/clusters)
    xx = np.repeat(xxParent, numbPointsDaughter);
    yy = np.repeat(yyParent, numbPointsDaughter);

    # translate points (ie parents points are the centres of cluster disks)
    xx = xx + xx0;
    yy = yy + yy0;

    # thin points if outside the simulation window
    booleInside = ((xx >= xMin) & (xx <= xMax) & (yy >= yMin) & (yy <= yMax));
    # retain points inside simulation window
    xx = xx[booleInside];  
    yy = yy[booleInside];

    # Plotting
    if plot:
        plt.scatter(xx, yy, c = 'blue')#, edgecolor='b', facecolor='none', alpha=0.5);
        plt.xlabel('x');
        plt.ylabel('y');
        plt.axis('equal');
    
    return xx, yy

class Utrecht:
    def __init__(self, randomize, seed):
        
        # The graph is directed since it may not take the same time to go from
        # A to B than from B to A.
        self.G = nx.DiGraph()
        
        # The data identifies vertices with postal coldes.
        # We map these postal codes to vertex ids 0, 1, ..., |V|-1.
        self.pc = []
        
        # A vertex is a list of [vertex id, population density, x coordinate,
        # y coordinate].
        # The sum of inhabitant densities of all the vertices is equal to 1.
        self.V = []
        
        # The edges are weighted by the number of seconds it takes to go from
        # one vertex to another.
        self.E = []
        
        # Location of the ambulance bases.
        self.B = []
        
        # Location of the hospitals.
        self.H = []
        
        # Number of ambulances (as stated by the authors of "Improving fairness
        # in ambulance planning by time sharing").
        self.num_ambulances = 19

        self.randomize = randomize
        random.seed(seed)
        np.random.seed(seed)
        
        self.load()
        self.build_graph()        

    def load(self):
        """Loads the data from the data/ directory and converts the postal codes
           to numbers from 0 to |V|-1.
        """
        if not self.randomize:
            # Load vertices.
            with open('../data/demand_nodes.txt', 'r') as f:
                for i in f.readlines()[1:]:
                    vertex = [float(j) for j in i.strip('\n').split('\t')]
                    vertex[0] = int(vertex[0])
                    self.V.append(vertex)
                # Construct mapping.
                self.pc = [j[0] for j in self.V]
                # Replace postal code by mapped index.
                for j in self.V:
                    j[0] = self.pc.index(j[0])
            f.close()

            # Load edges.
            with open('../data/rijtijden_ravu_rivm2009.txt', 'r') as f:
                def randomize_dist(di):
                    """adds or removes 0-20% of the initial distance"""
                    diff = 0.2
                    modif = 1+random.uniform(-diff, diff)
                    return int(int(di)*modif)

                raw_edges = [randomize_dist(i) if self.randomize else int(i) for i in f.readline().strip('\n').split(',')[:-1]]
                raw_edges = [raw_edges[i:i+len(self.V)] for i in
                             range(0, len(raw_edges), len(self.V))]
            f.close()
            for i in range(len(self.V)):
                for j in range(len(self.V)):
                    edge = (i, j, raw_edges[i][j])
                    self.E.append(edge)

            # Load bases.
            bases = [3436, 3447, 3561, 3582, 3608, 3645, 3707, 3811, 3823,
                     3911, 3941, 3958, 3743, 3991, 3417, 3769, 3648, 3931]
            for i in bases:
                self.B.append(self.pc.index(i))

            # Load hospitals.
            hospitals = [3584, 3435, 3543, 3582, 3813]
            for i in hospitals:
                self.H.append(self.pc.index(i))
        else:
            xx,yy = MaternProcess() # Generate clustered process
            nV = len(xx)
            locs = [i for i in range(nV)]
            den = np.random.exponential(5, (nV))
            den = den/sum(den)
            self.V = [[i, den[i], xx[i], yy[i]] for i in range(nV)]
            def dist (x1, y1, x2, y2):
                """
                Distance metric
                """
                # Manhattan distance
                # return abs(x1-x2)+abs(y1-y2)
                # Asymmetric Manhattan distance
                return abs(x1-x2)+abs(y1-y2) + 0.05*(x1-x2)
                # Euclidean distance
                # return np.sqrt((x1-x2)*(x1_x2) + (y1-y2)*(y1-y2))
                # Asymmetric Euclidean distance
                # return np.sqrt((x1-x2)*(x1_x2) + (y1-y2)*(y1-y2)) + 0.05*(x1-x2)

            self.E = [(i, j, int(10*dist(xx[i], yy[i], xx[j], yy[j]))) for i in range(nV) for j in range(nV)]

            # Bases
            nBases = 18 # Or use a smart way to decide this
            self.B = [int((i+0.5)*nV/nBases) for i in range(nBases)]
            self.pc = locs

            
    def build_graph(self):
        """Build graph from the data loaded."""
        self.G.add_nodes_from([i[0] for i in self.V])
        self.G.add_weighted_edges_from(self.E)

        
    def reduce_graph(self, seconds):
        """Removes all edges whose traveling time is
           greater than 'seconds' seconds.
        """
        for i in self.E:
            if i[2] > seconds:
                self.G.remove_edge(i[0], i[1])

                
    def get_population_densities(self):
        densities = [x[1] for x in self.V]
        #random.Random(0).shuffle(densities)
        return densities


    def get_sufficient_coverage(self):
        """Returns a list of the number of ambulances needed in the zones
           to ensure sufficient coverage.
        """
        pop_density = self.get_population_densities()
        sorted_pd = sorted(pop_density)
        n = len(self.V)
        q1 = sorted_pd[n//4]
        q2 = sorted_pd[n//2]
        q3 = sorted_pd[3*n//4]

        def pop(x):
            if x <= q1:
                return 1
            elif x <= q2:
                return 2
            elif x <= q3:
                return 3
            else:
                return 4

        sufficient = [x for x in map(pop, pop_density)]
        return sufficient


    def get_binary_adjacency(self):
        """Returns a binary adjacency matrix of G."""
        adj = nx.to_numpy_matrix(self.G)
        n = len(self.V)
        for i in range(n):
            for j in range(n):
                if (i == j):
                    adj[i, j] = 1
                elif (adj[i, j] != 0):
                    adj[i, j] = 1
        adj = adj.astype(int)
        return adj


    def allocation_coverage(self, allocation):
        """Get the sufficient coverage vector associated with an allocation."""
        sufficient = self.get_sufficient_coverage()
        adj = self.get_binary_adjacency()
        n = len(self.V)
        coverages = [0 for _ in range(n)]
        for i in range(len(allocation)):
            base = self.B[i]
            num_amb = allocation[i]
            coverages = [coverages[j]+num_amb*adj[base, j] for j in range(n)]
        # Normalize ambulance coverage to sufficient coverage
        coverages = [1 if coverages[i] >= sufficient[i] else 0 for i in range(n)]
        return coverages


    def allocations_coverage(self, allocations):
        """Get the cumulative sufficient coverage vector associated with a list of allocations."""
        coverages = []
        for allocation in allocations:
            coverages.append(self.allocation_coverage(allocation))
        return [sum(x) for x in zip(*coverages)]

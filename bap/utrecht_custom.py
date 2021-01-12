import networkx as nx
import random
import math

import numpy as np

class Utrecht:
    def __init__(self, gpickle_file, bases_file):
        
        # The graph is directed since it may not take the same time to go from
        # A to B than from B to A.
        self.G = nx.DiGraph()
        
        # A vertex is a tuple of (x,y).
        self.V = []
        
        # The edges are weighted by the number of seconds it takes to go from
        # one vertex to another.
        self.E = []
        
        # Location of the ambulance bases.
        self.B = []
        
        self.load(gpickle_file, bases_file)

    def load(self, gpickle_file, bases_file):
        """Loads the data from the gpickle file and from the bases.txt file.
        """
        # Load graph (gpickle file)
        self.G = nx.read_gpickle(gpickle_file)
        self.V = list(self.G.nodes)
        
        # Load bases
        self.B = []
        with open(bases_file, 'r') as f:
            for line in f.readlines():
                self.B.append(eval(line.strip('\n')))
        self.B = [self.V.index(i) for i in self.B] # bases are the index of the target vertex
        
        # if not self.randomize:
        #     # Load vertices.
        #     with open('../data/demand_nodes.txt', 'r') as f:
        #         for i in f.readlines()[1:]:
        #             vertex = [float(j) for j in i.strip('\n').split('\t')]
        #             vertex[0] = int(vertex[0])
        #             self.V.append(vertex)
        #         # Construct mapping.
        #         self.pc = [j[0] for j in self.V]
        #         # Replace postal code by mapped index.
        #         for j in self.V:
        #             j[0] = self.pc.index(j[0])
        #     f.close()

        #     # Load edges.
        #     with open('../data/rijtijden_ravu_rivm2009.txt', 'r') as f:
        #         def randomize_dist(di):
        #             """adds or removes 0-20% of the initial distance"""
        #             diff = 0.2
        #             modif = 1+random.uniform(-diff, diff)
        #             return int(int(di)*modif)

        #         raw_edges = [randomize_dist(i) if self.randomize else int(i) for i in f.readline().strip('\n').split(',')[:-1]]
        #         raw_edges = [raw_edges[i:i+len(self.V)] for i in
        #                      range(0, len(raw_edges), len(self.V))]
        #     f.close()
        #     for i in range(len(self.V)):
        #         for j in range(len(self.V)):
        #             edge = (i, j, raw_edges[i][j])
        #             self.E.append(edge)

        #     # Load bases.
        #     bases = [3436, 3447, 3561, 3582, 3608, 3645, 3707, 3811, 3823,
        #              3911, 3941, 3958, 3743, 3991, 3417, 3769, 3648, 3931]
        #     for i in bases:
        #         self.B.append(self.pc.index(i))



        # self.G.add_nodes_from([i[0] for i in self.V])
        # self.G.add_weighted_edges_from(self.E)


    def get_population_densities(self):
        centroid = (sum(x[0] for x in self.V)/len(self.V), sum(x[1] for x in self.V)/len(self.V))
        def dist_2_pts(x1, y1, x2, y2):
            return math.sqrt((x1-x2)**2 + (y1-y2)**2)
        #densities = [1000-dist_2_pts(x[0], x[1], centroid[0], centroid[1])*random.Random(0).uniform(0.0001, 1) for x in self.V] # 1000000-dist because we want the largest number to indicate a lower distance (a higher pop density)
        densities = [dist_2_pts(x[0], x[1], centroid[0], centroid[1])*random.Random(int(x[0])).uniform(1, 1.2) for x in self.V] # <--------------- random.uniform!
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
        #random.Random(0).shuffle(sufficient)
        return sufficient
        

    # def get_sufficient_coverage(self):
    #     """Returns a list of the number of ambulances needed in the zones
    #        to ensure sufficient coverage.
    #     """
    #     # Spread the demand randomly-ish
    #     return [((i+4)%4)+1 for i in range(len(self.V))]


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

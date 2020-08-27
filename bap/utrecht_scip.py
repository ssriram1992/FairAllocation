import networkx as nx


class Utrecht:
    def __init__(self):
        
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
        
        self.load()
        self.build_graph()

    def load(self):
        """Loads the data from the data/ directory and converts the postal codes
           to numbers from 0 to |V|-1.
        """
        
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
            raw_edges = [int(i) for i in f.readline().strip('\n').split(',')[:-1]]
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
        return [x[1] for x in self.V]


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

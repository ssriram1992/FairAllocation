import networkx as nx

class Utrecht:
    def __init__(self):
        # The graph is directed since it may not take the same time to go from A to B than from B to A
        self.G = nx.DiGraph()
        # The data identifies vertices with postal coldes. We map these postal codes to vertex ids 0, 1, ...
        self.pc = []
        # A vertex is a list of [vertex id, population density, x coordinate, y coordinate] (the sum of inhabitant densities of all the vertices is equal to 1)
        self.V = []
        # The edges are weighted by the number of seconds it takes to go from one vertex to another
        self.E = []
        # Location of the ambulance bases
        self.B = []
        # Location of the hospitals
        self.H = []
        # Number of ambulances (as stated by the authors of "Improving fairness in ambulance planning by time sharing")
        self.num_ambulances = 19
        self.load()
        self.build_graph()

    def load(self):
        """Loads the data from the data/ directory and converts the postal codes to numbers from 0 to |V|-1."""
        # Load vertices
        with open('data/demand_nodes.txt', 'r') as f:
            for i in f.readlines()[1:]:
                vertex = [float(j) for j in i.strip('\n').split('\t')]
                vertex[0] = int(vertex[0])
                self.V.append(vertex)
            # Construct mapping
            self.pc = [j[0] for j in self.V]
            # Replace postal code by mapped index
            for j in self.V:
                j[0] = self.pc.index(j[0])
        f.close()
        # Load edges
        with open('data/rijtijden_ravu_rivm2009.txt', 'r') as f:
            raw_edges = [int(i) for i in f.readline().strip('\n').split(',')[:-1]]
            raw_edges = [raw_edges[i:i+len(self.V)] for i in range(0, len(raw_edges), len(self.V))]
        f.close()
        for i in range(len(self.V)):
            for j in range(len(self.V)):
                edge = (i, j, raw_edges[i][j])
                self.E.append(edge)
        # Load bases
        bases = [3436, 3447, 3561, 3582, 3608, 3645, 3707, 3811, 3823, 3911, 3941, 3958, 3743, 3991, 3417, 3769, 3648, 3931]
        for i in bases:
            self.B.append(self.pc.index(i))
        # Load hospitals
        hospitals = [3584, 3435, 3543, 3582, 3813]
        for i in hospitals:
            self.H.append(self.pc.index(i))

    def build_graph(self):
        """Build graph from the data loaded."""
        self.G.add_nodes_from([i[0] for i in self.V])
        self.G.add_weighted_edges_from(self.E)

    def reduce_graph(self, seconds):
        """Removes all edges whose traveling time is greater than 'seconds' seconds."""
        for i in self.E:
            if i[2] > seconds:
                self.G.remove_edge(i[0], i[1])

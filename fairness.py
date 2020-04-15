import gurobipy as gp
import networkx as nx
import numpy as np

def makeGraph(numVert, sparsity=0.5):
    """
    Makes an adjacency matrix of a graph with given sparsity
    
    Parameters
    ----------
    numVert: int
        Number of vertices in the target graph
    sparsity: double
        Probability that an edge exists
    
    Returns
    -------
    numpy.ndarray
        Adjacency matrix
    """
    Adj = np.zeros((numVert, numVert),dtype=int)
    for i in range(numVert):
        for j in range(i):
            val = (np.random.rand() < sparsity)*1
            Adj[i,j] = val
            Adj[j,i] = val
    return Adj
def drawGraph(Adj):
    """
    Draws a graph using networkx
    
    Given the adjacency matrix of a graph, draws it using networkx. 
    Warning: Do not try for huge graphs. 
    
    Parameters
    ----------
    Adj: ndArray
        Adjacency matrix 
    
    Returns
    -------
    networkx.Graph
        An object containing the graph defined by Adj
    """
    G = nx.Graph()
    nV = Adj.shape[0]
    G.add_nodes_from([v for v in range(nV)])
    for i in range(nV):
        for j in range(i):
            if Adj[i,j]:
                G.add_edge(i,j)
    nx.draw_shell(G, with_labels=True, font_weight='bold', font_color='white')
    return G
def getGurobiBasicModel_Diff(Adj, f):
    """
    Creates a gurobi basic model with max-min fairness metric
    
    Given the adjacency matrix, creates the gurobi model that decides the 
    number of allocations to be done to each of the vertices so that the 
    total fairness measured according to the max-min fairness metric is 
    maximized. 
    
    Parameters
    ----------
    Adj: ndArray
        Adjacency matrix 
    f: ndArray
        Benefit of v when an ambulance is in u is in f[u,v]
    
    Returns
    -------
    gurobipy.Model
        The required gurobi model
        
    """
    nV = Adj.shape[0]
    M = gp.Model()
    x = M.addVars(range(nV),lb=0,vtype=gp.GRB.INTEGER, name = "x")
    tau = M.addVars(range(nV),lb=1,vtype=gp.GRB.CONTINUOUS, name = "tau")
    tauMax = M.addVar(obj=1)
    tauMin = M.addVar(obj=-1)
    for v in range(nV):
        M.addConstr(x[v] + gp.quicksum([f[u,v]*Adj[v,u]*x[u] for u in range(nV)]) == tau[v])
    M.addGenConstrMax(tauMax, [tau[v] for v in range(nV)], name="Max")
    M.addGenConstrMin(tauMin, [tau[v] for v in range(nV)], name="Min")
    # (0,0,0, ..., 0) should not be a solution - no good cut
    M.addConstr(gp.quicksum(x[v] for v in range(nV)) >= 1)
    M.update()
    return M
def getGurobiBasicModel_L1(Adj,f):
    """
    Creates a gurobi basic model with L1 fairness metric
    
    Given the adjacency matrix, creates the gurobi model that decides the 
    number of allocations to be done to each of the vertices so that the 
    total fairness measured according to the L1 fairness metric is 
    maximized. 
    
    Parameters
    ----------
    Adj: ndArray
        Adjacency matrix 
    f: ndArray
        Benefit of v when an ambulance is in u is in f[u,v]
    
    Returns
    -------
    gurobipy.Model
        The required gurobi model
        
    """
    nV = Adj.shape[0]
    M = gp.Model()
    M = gp.Model()
    x = M.addVars(range(nV),lb=0,vtype=gp.GRB.INTEGER, name = "x")
    tau = M.addVars(range(nV),lb=1,vtype=gp.GRB.CONTINUOUS, name = "tau")
    tauAvg = M.addVar(name='tauAvg',vtype=gp.GRB.CONTINUOUS)
    tauDiffAbs = M.addVars(range(nV),lb=0,vtype=gp.GRB.CONTINUOUS, name = "tauDiffAbs")
    tauDiff = M.addVars(range(nV),lb=-gp.GRB.INFINITY,vtype=gp.GRB.CONTINUOUS, name = "taudiff")
    M.addConstrs(tau[v]-tauAvg == tauDiff[v] for v in range(nV))
    for v in range(nV):
        M.addConstr(x[v] + gp.quicksum([f[u,v]*Adj[v,u]*x[u] for u in range(nV)]) == tau[v])
        M.addConstr(tauDiffAbs[v]==gp.abs_(tauDiff[v]))
    M.addConstr(nV*tauAvg == gp.quicksum(tau[v] for v in range(nV)))
    M.setObjective(gp.quicksum(tauDiffAbs[v] for v in range(nV)))
    M.update()
    return M
def getMinimalGurobi(M, Adj):
    """
    Finds the _least_ number of allocations that achieves the fairness level
    promised in the model M
    
    Parameters
    ----------
    M:   gurobipy.Model
        Solved gurobipy.Model that identified the maximum possible fairness 
        and containing all the equations defining the said fairness metric.
    Adj: ndArray
        Adjacency matrix 
    
    Returns
    -------
    gurobipy.Model
        The required gurobi model
        
    """
    nV = Adj.shape[0]
    M2 = M.copy()
    obj = M.getAttr("objval")
    M2.addConstr(M2.getObjective() == obj)
    M2.update()
    M2.setObjective(gp.quicksum(M2.getVarByName("x["+str(v)+"]") for v in range(nV)))
    M2.update()
    return M2
def greedDiff(nV, Adj, f, M2, k=1):
    """
    Greedily maximizes fairness (according to max-min metrix) 
    over multiple rounds as informed by getMinimalGurobi()
    
    Parameters
    ----------
    nV: int
        number of vertices
    Adj: ndArray
        Adjacency matrix 
    f:  ndArray
        Benefit of v when an ambulance is in u is in f[u,v]
    M2: gurobipy.Model
        The model returned by getMinimalGurobi() after being solved. 
        This should have details on the number of times ambulance has 
        to be placed in each node v
    k:  int
        Number of ambulance vehicles
    
    Returns
    -------
    Allocs: list
        A list where i-th element has the list of locations where 
        ambulance has to be placed at the i-th round
    Fairnesses: list
        A list where i-th element denotes the unfairness at the end 
        of i-th round.
    """
    # The maximum number of times ambulance can be placed at v
    xLim = [round(M2.getVarByName("x["+str(v)+"]").X) for v in range(nV)]
    # Welfare obtained by v so far.
    tauSoFar = [0 for i in range(nV)]
    # Number of times ambulance placed at v so far 
    # xSoFar should be < xLim componentwise
    xSoFar = [0 for i in range(nV)]
    # Total number of roubds
    nRounds = int(np.ceil(round(M2.objval)/k))
    
    Allocs = []
    Fairnesses = []
    for rnd in range(nRounds):
        # gurobiModel
#         print("Round:",rnd, "of", nRounds)
        Mgreed = gp.Model()
        Mgreed.setParam("LogToConsole",0)
        strict = not (rnd==(nRounds-1))
        strict = False

        xt = Mgreed.addVars(range(nV), vtype=gp.GRB.BINARY, name="xt")
        tau = Mgreed.addVars(range(nV),lb=0,vtype=gp.GRB.CONTINUOUS, name = "tau")

        tauMax = Mgreed.addVar(obj=1)
        tauMin = Mgreed.addVar(obj=-1)

        for v in range(nV):
            Mgreed.addConstr(xt[v] + gp.quicksum([f[u,v]*Adj[v,u]*xt[u] for u in range(nV)]) +tauSoFar[v] == tau[v])
            if xSoFar[v] >= xLim[v]:
                xt[v].ub = 0
        Mgreed.addGenConstrMax(tauMax, [tau[v] for v in range(nV)], name="Max")
        Mgreed.addGenConstrMin(tauMin, [tau[v] for v in range(nV)], name="Min")
        # (0,0,0, ..., 0) should not be a solution - no good cut
        Mgreed.addConstr(gp.quicksum(xt[v] for v in range(nV)) >= 1)

        # Maximum k ambulance
        if strict:
            Mgreed.addConstr(gp.quicksum(xt[v] for v in range(nV)) == k)
        else:
            Mgreed.addConstr(gp.quicksum(xt[v] for v in range(nV)) <= k)
        

        Mgreed.update()
        Mgreed.optimize()
        alloc = [v for v in range(nV) if xt[v].X>0.5]
        fairness = Mgreed.objval
        Allocs.append(alloc)
        Fairnesses.append(fairness)
#         print("Ambulance at: ", alloc)
        xSoFar = [xSoFar[v] + xt[v].X for v in range(nV)]
        tauSoFar =[tau[i].X for i in range(nV)]
    return Allocs, Fairnesses
def extractAllocAgg(M, nV):
    """
    Extracts the x-solution from a gurobipy.Model
    
    Given a solved gurobipy.Model and the number of vertices,
    extracts the value of every variable of the form x[v]
    
    """
    x = np.zeros(nV, dtype=int)
    for i in range(nV):
        x[i] = round(M.getVarByName("x["+str(i)+"]").X)
    return x

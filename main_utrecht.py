import gurobipy as gp
from math import ceil
import networkx as nx
import numpy as np
from utrecht import *

def basic_model():
    # Model
    M = gp.Model()
    x = M.addVars(num_rounds, len(bases), vtype=gp.GRB.BINARY, name='x')

    # Each round, all ambulances must be assigned
    for i in range(num_rounds):
        M.addConstr(gp.quicksum(x[i, j] for j in range(len(bases))) == num_ambulances)

    # A configuration must be "efficient", i.e., >=95% of the zones covered
    tau = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau')
    for i in range(num_rounds):
        for j in range(n):
            if j in bases:
                M.addConstr((tau[i, j] == 1) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 1))
            else:
                M.addConstr((tau[i, j] == 1) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 1))
    for i in range(num_rounds):
        M.addConstr(gp.quicksum(tau[i, j] for j in range(n)) >= ceil(0.95*n))

    # Transition cost: A maximum of max_transition vehicles can be moved when switching between configurations
    remain = M.addVars(num_rounds-1, range(len(bases)), vtype=gp.GRB.BINARY)
    for i in range(num_rounds-1):
        for j in range(len(bases)):
            M.addConstr((remain[i, j] == 1) >> ((x[i, j] + x[i+1, j]) >= 2))
    for i in range(num_rounds-1):
        M.addConstr(gp.quicksum(remain[i, j] for j in range(len(bases))) >= ceil((1-max_transition)*num_ambulances))
        
    return M, tau

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes".
utrecht = Utrecht()
utrecht.reduce_graph(900)

# Bases are the only vertices that can be allocated ambulances
bases = utrecht.B
n = len(utrecht.V)

# Replace the seconds by a boolean representing if the zones are close enough
# to cover each other or not.
adj = nx.to_numpy_matrix(utrecht.G)
for i in range(n):
    for j in range(n):
        if (adj[i, j] != 0):
            adj[i, j] = 1
adj = adj.astype(int)

# Model parameters
num_ambulances = 7
num_rounds = 5
max_transition = 0.25

# Model phase 1
M1, tau = basic_model()
    
# Objective
tau_max = M1.addVar(obj=1)
tau_min = M1.addVar(obj=-1)
tau_col = M1.addVars(range(n), vtype=gp.GRB.INTEGER, name='tau_col')
for i in range(n):
    M1.addConstr(tau_col[i] == gp.quicksum([tau[j, i] for j in range(num_rounds)]))
M1.addGenConstrMax(tau_max, [tau_col[i] for i in range(n)], name='tau_max')
M1.addGenConstrMin(tau_min, [tau_col[i] for i in range(n)], name='tau_min')

# Solve
M1.update()
M1.optimize()
phi = M1.getObjective().getValue()

# Model phase 2
M2, tau = basic_model()
    
# Objective
tau_max = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
tau_min = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
tau_col = M2.addVars(num_rounds, n, vtype=gp.GRB.INTEGER, name='tau_col')
for i in range(num_rounds):
    for j in range(n):
        M2.addConstr(tau_col[i, j] == gp.quicksum([tau[k, j] for k in range(i+1)]))
for i in range(num_rounds):
    M2.addGenConstrMax(tau_max[i], [tau_col[i, j] for j in range(n)], name='tau_max')
    M2.addGenConstrMin(tau_min[i], [tau_col[i, j] for j in range(n)], name='tau_min')
M2.setObjective(gp.quicksum([tau_max[i]-tau_min[i] for i in range(num_rounds)]), gp.GRB.MINIMIZE)

# Ensure that the final path leads to the optimal solution
M2.addConstr(tau_max[num_rounds-1]-tau_min[num_rounds-1] == phi)

# Solve
M2.update()
M2.optimize()

print('Status:', M2.Status)
sol = M2.getVars()
for i in sol:
    if ('x' in i.Varname) and (i.X > 0):
        print(i.Varname, '---', i.X)

import gurobipy as gp
from math import ceil
import networkx as nx
import numpy as np
from utrecht import *

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
num_ambulances = 8
num_rounds = 28
max_transition = 0.25

# Model
M = gp.Model()
x = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='x')

# Ambulances can only be assigned to bases
for i in range(num_rounds):
    for j in range(n):
        if (j not in bases):
            M.addConstr(x[i, j] == 0)

# Each round, all ambulances must be assigned
for i in range(num_rounds):
    M.addConstr(gp.quicksum(x[i, j] for j in range(n)) == num_ambulances)

# A configuration must be "efficient", i.e., >=95% of the zones covered
tau = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau')
for i in range(num_rounds):
    for j in range(n):
        M.addConstr((tau[i, j] == 1) >> (x[i, j]+gp.quicksum(x[i, k] * adj[k, j] for k in range(n)) >= 1))
for i in range(num_rounds):
    M.addConstr(gp.quicksum(tau[i, j] for j in range(n)) >= ceil(0.95*n))

# Transition cost: A maximum of max_transition vehicles can be moved when switching between configurations
remain = M.addVars(num_rounds-1, range(n), vtype=gp.GRB.BINARY)
for i in range(num_rounds-1):
    for j in range(n):
        M.addConstr((remain[i, j] == 1) >> ((x[i, j] + x[i+1, j]) >= 2))
for i in range(num_rounds-1):
    M.addConstr(gp.quicksum(remain[i, j] for j in range(n)) >= ceil((1-max_transition)*num_ambulances))
    
# Objective
tau_max = M.addVar(obj=1)
tau_min = M.addVar(obj=-1)
tau_col = M.addVars(range(n), vtype=gp.GRB.INTEGER, name='tau_col')
for i in range(n):
    M.addConstr(tau_col[i] == gp.quicksum([tau[j, i] for j in range(num_rounds)]))
M.addGenConstrMax(tau_max, [tau_col[i] for i in range(n)], name='tau_max')
M.addGenConstrMin(tau_min, [tau_col[i] for i in range(n)], name='tau_min')

# Solve
M.update()
M.optimize()

# print('Status:', M.Status)
# sol = M.getVars()
# for i in sol:
#     if ('x' in i.Varname) and (i.X > 0):
#         print(i.Varname, '---', i.X)

# TODO: Take into account zone density when considering if a zone is covered or not

import gurobipy as gp
from math import ceil
import networkx as nx
import numpy as np
import re
from utrecht import *


def basic_model():
    M = gp.Model()
    x = M.addVars(num_rounds, len(bases), vtype=gp.GRB.BINARY, name='x')

    # Each round, all ambulances must be assigned.
    for i in range(num_rounds):
        M.addConstr(gp.quicksum(x[i, j] for j in range(len(bases))) == num_ambulances)

    # A configuration must be "efficient", i.e., >=95% of the zones covered.
    tau = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau')
    tau2 = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau2')
    for i in range(num_rounds):
        for j in range(n):
            if j in bases:
                M.addConstr((tau[i, j] == 1) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 1))
                M.addConstr((tau[i, j] == 0) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) == 0))
                M.addConstr((tau2[i, j] == 1) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 2))
                M.addConstr((tau2[i, j] == 0) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= 1))
            else:
                M.addConstr((tau[i, j] == 1) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 1))
                M.addConstr((tau[i, j] == 0) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) == 1))
                M.addConstr((tau2[i, j] == 1) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= 2))
                M.addConstr((tau2[i, j] == 0) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= 1))
    for i in range(num_rounds):
        M.addConstr(gp.quicksum(tau[i, j] for j in range(n)) >= ceil(0.95*n))

    # Transition cost: A maximum of max_transition vehicles can be moved when switching between configurations.
    remain = M.addVars(num_rounds-1, len(bases), vtype=gp.GRB.BINARY)
    for i in range(num_rounds-1):
        for j in range(len(bases)):
            M.addConstr((remain[i, j] == 1) >> ((x[i, j] + x[i+1, j]) >= 2))
            M.addConstr((remain[i, j] == 0) >> ((x[i, j] + x[i+1, j]) <= 1))
    for i in range(num_rounds-1):
        M.addConstr(gp.quicksum(remain[i, j] for j in range(len(bases))) >= ceil((1-max_transition)*num_ambulances))
        
    return M, tau, tau2

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes (900 seconds)".
utrecht = Utrecht()
utrecht.reduce_graph(900)

# Bases are the only vertices that can be allocated ambulances.
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
num_rounds = 5
max_transition = 0.25

# Model phase 1
M1, tau, tau2 = basic_model()
    
# Objective
tau2_max = M1.addVar(obj=1)
tau2_min = M1.addVar(obj=-1)
tau2_col = M1.addVars(n, vtype=gp.GRB.INTEGER, name='tau2_col')
for i in range(n):
    M1.addConstr(tau2_col[i] == gp.quicksum([tau2[j, i] for j in range(num_rounds)]))
M1.addGenConstrMax(tau2_max, [tau2_col[i] for i in range(n)], name='tau2_max')
M1.addGenConstrMin(tau2_min, [tau2_col[i] for i in range(n)], name='tau2_min')

# Solve
M1.update()
M1.optimize()
phi = M1.getObjective().getValue()

print('Status (2 is optimal):', M1.Status)
sol = M1.getVars()
out_x = [[] for _ in range(num_rounds)]
out_tau = [[0 for _ in range(n)] for _ in range(num_rounds)]
out_tau2 = [[0 for _ in range(n)] for _ in range(num_rounds)]
for i in sol:
    search_x = re.search('^x\[([0-9]*),([0-9]*)\]$', i.varname)
    search_tau = re.search('^tau\[([0-9]*),([0-9]*)\]$', i.varname)
    search_tau2 = re.search('^tau2\[([0-9]*),([0-9]*)\]$', i.varname)
    if search_x and i.x > 0:
        out_x[int(search_x.group(1))].append(search_x.group(2))
    if search_tau:
        out_tau[int(search_tau.group(1))][int(search_tau.group(2))] = i.x
    if search_tau2:
        out_tau2[int(search_tau2.group(1))][int(search_tau2.group(2))] = i.x

print('---x-----------------------------')
for i in out_x:
    print([int(j) for j in i])
# print('---tau-----------------------------')
# for i in out_tau:
#     print([int(j) for j in i])
# print('---tau2-----------------------------')
# for i in out_tau2:
#     print([int(j) for j in i])
print('------------------------------------')
print('Single coverages:', [int(sum(i)) for i in out_tau])
print('Double coverages:', [int(sum(i)) for i in out_tau2])

# # Model phase 2
# M2, tau, tau2 = basic_model()
    
# # Objective
# # tau_max = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
# # tau_min = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
# # tau_col = M2.addVars(num_rounds, n, vtype=gp.GRB.INTEGER, name='tau_col')
# # for i in range(num_rounds):
# #     for j in range(n):
# #         M2.addConstr(tau_col[i, j] == gp.quicksum([tau[k, j] for k in range(i+1)]))
# # for i in range(num_rounds):
# #     M2.addGenConstrMax(tau_max[i], [tau_col[i, j] for j in range(n)], name='tau_max')
# #     M2.addGenConstrMin(tau_min[i], [tau_col[i, j] for j in range(n)], name='tau_min')
# # M2.setObjective(gp.quicksum([tau_max[i]-tau_min[i] for i in range(num_rounds)]), gp.GRB.MINIMIZE)
# tau_max = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
# tau_min = M2.addVars(num_rounds, vtype=gp.GRB.INTEGER)
# tau_col = M2.addVars(num_rounds, n, vtype=gp.GRB.INTEGER, name='tau_col')
# for i in range(num_rounds):
#     for j in range(n):
#         M2.addConstr(tau_col[i, j] == gp.quicksum([tau2[k, j] for k in range(i+1)]))
# for i in range(num_rounds):
#     M2.addGenConstrMax(tau_max[i], [tau_col[i, j] for j in range(n)], name='tau_max')
#     M2.addGenConstrMin(tau_min[i], [tau_col[i, j] for j in range(n)], name='tau_min')
# M2.setObjective(gp.quicksum([tau_max[i]-tau_min[i] for i in range(num_rounds)]), gp.GRB.MINIMIZE)

# # Ensure that the final path leads to the optimal solution
# M2.addConstr(tau_max[num_rounds-1]-tau_min[num_rounds-1] == phi)

# # Solve
# M2.update()
# M2.optimize()

# print('Status:', M2.Status)
# sol = M2.getVars()
# for i in sol:
#     if ('2' in i.Varname):# and (i.X > 0):
#         print(i.Varname, '---', i.X)
#     # if ('x' in i.Varname) and (i.X > 0):
#     #     print(i.Varname, '---', i.X)

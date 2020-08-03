import gurobipy as gp
from math import ceil, floor
import networkx as nx
import numpy as np
import re
from utrecht import *


def basic_model():
    M = gp.Model()

    # x[i, j]: Number of ambulances allocated to base j on round i
    x = M.addVars(num_rounds, len(bases), vtype=gp.GRB.INTEGER, lb=0, name='x')

    # Each round, all ambulances must be assigned.
    for i in range(num_rounds):
        M.addConstr(gp.quicksum(x[i, j] for j in range(len(bases))) == num_ambulances)

    # A configuration must be "efficient", i.e., >=95% of the zones covered.
    tau = M.addVars(num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau')
    for i in range(num_rounds):
        for j in range(n):
            if j in bases:
                M.addConstr((tau[i, j] == 1) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= sufficient[j]))
                M.addConstr((tau[i, j] == 0) >> (x[i, bases.index(j)]+gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= sufficient[j]-1))
            else:
                M.addConstr((tau[i, j] == 1) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= sufficient[j]))
                M.addConstr((tau[i, j] == 0) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= sufficient[j]-1))

    for i in range(num_rounds):
        M.addConstr(gp.quicksum(tau[i, j] for j in range(n)) >= ceil(min_sufficient*n))

    # Transition cost: A maximum of max_transition vehicles can be moved when switching between configurations.
    # delta[i, j]: Difference in the number of ambulances in zone j between rounds i and i+1
    delta = M.addVars(num_rounds-1, len(bases), lb=-gp.GRB.INFINITY, vtype=gp.GRB.INTEGER)
    # delta_abs[i, j]: Same as delta, but absolute value of the difference
    delta_abs = M.addVars(num_rounds-1, len(bases), lb=0, vtype=gp.GRB.INTEGER)
    for i in range(num_rounds-1):
        for j in range(len(bases)):
            M.addConstr(delta[i, j] == (x[i, j] - x[i+1, j]))
            M.addGenConstrAbs(delta_abs[i, j], delta[i, j])
    for i in range(num_rounds-1):
        # /2 because if A sends an ambulance to B, only one ambulance has moved, but A and B both show a delta of 1 (1+1=2, so we need to divide by 2)
        M.addConstr(gp.quicksum(delta_abs[i, j] for j in range(len(bases)))/2 <= floor(max_transition*num_ambulances))
        
    return M, tau

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

# Sufficient coverage required for various opulation densities.
pop_density = utrecht.get_population_densities()
sorted_pd = sorted(pop_density)
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

# Model parameters
num_ambulances = 20
num_rounds = 1
max_transition = 0.5 # 0 = no transition allowed, 1 = unlimited transitions
min_sufficient = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered

# Objective (coverage = sufficient coverage)
#   - efficiency: Maximize the sum of coverages.
#   - difference: Minimize the difference between the zone the most often covered, and the zone the least often covered.
#   - minmax: Minimize the coverage of the zone the most covered.
#   - maxmin: Maximize the coverage of the zone the least covered.
#   - min_uncovered: Minimize the number of zones which are never covered.
objective = 'maxmin'

# Model
mainM, tau = basic_model()

# tau_col[i]: The number of times zone i was sufficiently covered over interval T.
tau_col = mainM.addVars(n, lb=0, ub=num_rounds, vtype=gp.GRB.INTEGER, name='tau_col')
for i in range(n):
    mainM.addConstr(tau_col[i] == gp.quicksum([tau[j, i] for j in range(num_rounds)]))

# Maximum and minimum values of tau_col.
tau_max = mainM.addVar(name='tau_max')
tau_min = mainM.addVar(name='tau_min')
mainM.addGenConstrMax(tau_max, tau_col)
mainM.addGenConstrMin(tau_min, tau_col)

if objective == 'efficiency':
    efficiency = mainM.addVar(obj=-1, name='efficiency')
    mainM.addConstr(efficiency == gp.quicksum(tau_col))
elif objective == 'difference':
    tau_max.setAttr('obj', 1)
    tau_min.setAttr('obj', -1)
elif objective == 'minmax':
    tau_max.setAttr('obj', -1)
elif objective == 'maxmin':
    tau_min.setAttr('obj', -1)
elif objective == 'min_uncovered':
    uncovered = mainM.addVars(n, vtype=gp.GRB.BINARY, name='uncovered')
    for i in range(n):
        mainM.addConstr((uncovered[i] == 1) >> (tau_col[i] == 0))
        mainM.addConstr((uncovered[i] == 0) >> (tau_col[i] >= 1))
    sum_uncovered = mainM.addVar(obj=1, name='sum_uncovered')
    mainM.addConstr(sum_uncovered == gp.quicksum(uncovered))
    
# Solve
mainM.update()
mainM.optimize()
phi = mainM.getObjective().getValue()

print('----------------------------------------')
print('Status (2 is optimal):', mainM.Status)
print('Objective:', objective)
sol = mainM.getVars()
out_x = [[] for _ in range(num_rounds)]
out_tau = [[0 for _ in range(n)] for _ in range(num_rounds)]
out_tau_max = None
out_tau_min = None
for i in sol:
    search_x = re.search('^x\[([0-9]*),([0-9]*)\]$', i.varname)
    search_tau = re.search('^tau\[([0-9]*),([0-9]*)\]$', i.varname)
    search_tau_max = re.search('^tau_max$', i.varname)
    search_tau_min = re.search('^tau_min$', i.varname)
    if search_x and i.x > 0:
        out_x[int(search_x.group(1))].append((search_x.group(2), i.x))
    if search_tau:
        out_tau[int(search_tau.group(1))][int(search_tau.group(2))] = i.x
    if search_tau_max:
        out_tau_max = i.x
    if search_tau_min:
        out_tau_min = i.x

print('-----x----------------------------------')
for i in out_x:
    for j in i:
        print(j[0]+'('+str(int(j[1]))+') ', end='')
    print()
        #print([int(j) for j in i])
print('-----tau--------------------------------')
for i in out_tau:
    print([int(j) for j in i])
print('----------------------------------------')
print('Sufficient coverages:', [int(sum(i)) for i in out_tau])
print('tau_max:', out_tau_max)
print('tau_min:', out_tau_min)

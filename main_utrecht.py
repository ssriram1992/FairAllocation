import gurobipy as gp
from math import ceil, floor
import networkx as nx
import numpy as np
import re
from utrecht import *


def basic_model(fix_x0, warm_xs, warm_tau_adjust, warm_num_rounds):
    """
    fix_x0: Boolean to decide if we fix warm_xs[0] (for greedy).
    warm_xs: List of lists of the xs to warm start.
    warm_tau_adjust: List of tau adjustments (for greedy).
    warm_num_rounds: Number of rounds to consider.
    """
    M = gp.Model()

    # x[i, j]: Number of ambulances allocated to base j on round i.
    x = M.addVars(warm_num_rounds, len(bases), vtype=gp.GRB.INTEGER, name='x')
    
    # Warm start.
    if warm_xs is not None:
        for i in range(len(warm_xs)):
            for j in range(len(bases)):
                if fix_x0 and i == 0:
                    M.addConstr(x[i, j] == warm_xs[i][j])
                else:
                    x[i, j].start = warm_xs[i][j]
    if warm_tau_adjust is None:
        warm_tau_adjust = [0 for _ in range(n)]
    
    # Each round, all ambulances must be assigned.
    for i in range(warm_num_rounds):
        M.addConstr(gp.quicksum(x[i, j] for j in range(len(bases))) == num_ambulances)

    # tau[i, j]: Binary indication of if zone j is sufficiently covered on round i (this takes into account population density).
    tau = M.addVars(warm_num_rounds, range(n), vtype=gp.GRB.BINARY, name='tau')
    for i in range(warm_num_rounds):
        for j in range(n):
            if j in bases:
                M.addConstr((tau[i, j] == 1) >> (x[i, bases.index(j)] + gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= sufficient[j]))
                M.addConstr((tau[i, j] == 0) >> (x[i, bases.index(j)] + gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= sufficient[j]-1))
            else:
                M.addConstr((tau[i, j] == 1) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) >= sufficient[j]))
                M.addConstr((tau[i, j] == 0) >> (gp.quicksum(x[i, k] * adj[bases[k], j] for k in range(len(bases))) <= sufficient[j]-1))

    # A configuration must be "efficient", i.e., must satisfy a minimum coverage threshold.
    for i in range(warm_num_rounds):
        M.addConstr(gp.quicksum(tau[i, j] for j in range(n)) >= ceil(min_coverage*n))

    # delta[i, j]: Difference in the number of ambulances in zone j between rounds i and i+1.
    delta = M.addVars(warm_num_rounds-1, len(bases), lb=-gp.GRB.INFINITY, vtype=gp.GRB.INTEGER)
    # delta_abs[i, j]: Same as delta, but absolute value of the difference.
    delta_abs = M.addVars(warm_num_rounds-1, len(bases), lb=0, vtype=gp.GRB.INTEGER)
    for i in range(warm_num_rounds-1):
        for j in range(len(bases)):
            M.addConstr(delta[i, j] == (x[i, j] - x[i+1, j]))
            M.addGenConstrAbs(delta_abs[i, j], delta[i, j])
    # Transition cost: A maximum of max_transition vehicles can be moved when switching between configurations.
    for i in range(warm_num_rounds-1):
        # /2 because if A sends an ambulance to B, only one ambulance has moved, but A and B both show a delta of 1 (1+1=2, so we need to divide by 2).
        M.addConstr(gp.quicksum(delta_abs[i, j] for j in range(len(bases)))/2 <= floor(max_transition*num_ambulances))

    # tau_col[i]: The number of times zone i was sufficiently covered over interval T.
    tau_col = M.addVars(n, lb=0, ub=warm_num_rounds, vtype=gp.GRB.INTEGER, name='tau_col')
    for i in range(n):
        M.addConstr(tau_col[i] + warm_tau_adjust[i] == gp.quicksum([tau[j, i] for j in range(warm_num_rounds)]))

    # Maximum and minimum values of tau_col.
    tau_max = M.addVar(name='tau_max')
    tau_min = M.addVar(name='tau_min')
    M.addGenConstrMax(tau_max, tau_col)
    M.addGenConstrMin(tau_min, tau_col)

    if objective == 'efficiency':
        efficiency = M.addVar(obj=-1, name='efficiency')
        M.addConstr(efficiency == gp.quicksum(tau_col))
    elif objective == 'difference':
        tau_max.setAttr('obj', 1)
        tau_min.setAttr('obj', -1)
    elif objective == 'minmax':
        tau_max.setAttr('obj', -1)
    elif objective == 'maxmin':
        tau_min.setAttr('obj', -1)
    elif objective == 'min_uncovered':
        uncovered = M.addVars(n, vtype=gp.GRB.BINARY, name='uncovered')
        for i in range(n):
            M.addConstr((uncovered[i] == 1) >> (tau_col[i] == 0))
            M.addConstr((uncovered[i] == 0) >> (tau_col[i] >= 1))
        sum_uncovered = M.addVar(obj=1, name='sum_uncovered')
        M.addConstr(sum_uncovered == gp.quicksum(uncovered))
        
    return M

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

# Sufficient coverage required for various population densities.
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
num_rounds = 3
max_transition = 1 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered

# Objective (coverage = sufficient coverage)
#   - efficiency: Maximize the sum of coverages.
#   - difference: Minimize the difference between the zone the most often covered, and the zone the least often covered.
#   - minmax: Minimize the coverage of the zone the most covered.
#   - maxmin: Maximize the coverage of the zone the least covered.
#   - min_uncovered: Minimize the number of zones which are never covered.
objective = 'maxmin'

def tau_adjuster(x):
    """Returns an array of tau adjustments for the x array in input."""
    coverings = [0 for _ in range(n)]
    for base_idx in range(len(bases)): # Zone/x index
        # print(len(x))
        # print(len(bases))
        # print(bases)
        # print(ambulance)
        base_zone = bases[base_idx] # Zone number
        base_ambulances = x[base_idx] # Number of ambulances in that base
        coverings[base_zone] = base_ambulances
        for zone in range(n):
            if (zone != base_idx):
                coverings[zone] += adj[base_zone, zone]*base_ambulances

    adjustments = [0 for _ in range(n)]
    for i in range(n):
        if coverings[i] >= sufficient[i]:
            adjustments[i] += 1

    return adjustments

def tau_adjusters(xs):
    """Returns an array of tau adjustments for the xs array in input."""
    adjustments = [0 for _ in range(n)]
    for i in xs:
        tau_for_x = tau_adjuster(i)
        for j in range(n):
            adjustments[j] += tau_for_x[j]
    return adjustments

def get_x(model):
    """Get the values of the xs of the various rounds (list of lists)."""
    sol = model.getVars()
    out_x = [[0 for _ in range(len(bases))] for _ in range(greedy_rounds)]
    for i in sol:
        search_x = re.search('^x\[([0-9]*),([0-9]*)\]$', i.varname)
        if search_x:
            out_x[int(search_x.group(1))][int(search_x.group(2))] = int(i.x)
    return out_x
    
# Greedy solution
greedy_x = []
greedy_rounds = 4
greedy_tau_adjustments = []

# First iteration has no fixed x
M = basic_model(False, None, None, greedy_rounds)
M.update()
M.optimize()
greedy_x.append(get_x(M)[0])
# Middle iterations keep the first x
for i in range(num_rounds-greedy_rounds-1):
    print("----", i)
    # print(greedy_x[-1])
    # print(tau_adjusters(greedy_x[:-1]))
    M = basic_model(True, [greedy_x[-1]], tau_adjusters(greedy_x[:-1]), greedy_rounds)
    M.update()
    M.optimize()
    greedy_x.append(get_x(M)[1])
    print(greedy_x)
print("======================OUT OF LOOP")
# Last iteration keeps the remaining xs
M = basic_model(True, [greedy_x[-1]], tau_adjusters(greedy_x[:-1]), greedy_rounds)
M.update()
M.optimize()
for i in range(greedy_rounds-1):
    greedy_x.append(get_x(M)[i+1])

print("===================START REAL MODEL")
    
# Model
mainM = basic_model(False, greedy_x, None, num_rounds)

# lets say the interval is 3:
# - solve for 3 with no warm start, add x[0] to the greedy list
# - for all the rest: solve for 3 with x[0] warm start, (include tau adjust (compute all the taus with the tau_adjuster function, then compile them into a single array)), and add x[1] to the greedy list
# IMPORTANT: do not do warm start for the greedy, but rather fix x[0] to some values

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
# print('-----tau--------------------------------')
# for i in out_tau:
#     print([int(j) for j in i])
print('----------------------------------------')
print('Sufficient coverages:', [int(sum(i)) for i in out_tau])
print('tau_max:', out_tau_max)
print('tau_min:', out_tau_min)

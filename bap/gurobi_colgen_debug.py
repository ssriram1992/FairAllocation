import networkx
import numpy
from utrecht_scip import *
import math
import gurobipy as gp

EPS = 1e-6

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes (900 seconds)".
utrecht = Utrecht()
utrecht.reduce_graph(900)

# Bases are the only vertices that can be allocated ambulances.
bases = utrecht.B
n = len(utrecht.V)

# Replace the seconds by a boolean representing if the zones are close enough
# to cover each other or not.
adj = utrecht.get_binary_adjacency()

# Sufficient coverage required for various population densities.
sufficient = utrecht.get_sufficient_coverage()

# Model parameters
num_ambulances = 20
num_rounds = 3
max_transition = 1 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = 4 # Maximum practical number of ambulances in a zone (doesn't seem to make much of a difference). This is 4 because a base with 4 ambulances can sufficiently cover any zone no matter what its population density


# Find one feasible configuration for the initial column.
def initial_feasible_configuration():
    """Returns one feasible allocation, along with its coverage."""
    m = gp.Model('M')
    x = [m.addVar(lb=0, ub=max_practical_ambulances, vtype='I', name='x'+str(i)) for i in range(len(bases))]

    # Constraint: Don't use more than the available number of ambulances.
    m.addConstr(gp.quicksum(x) <= num_ambulances)

    # Constraint: The configuration must sufficiently cover 95% of the zones.
    c = [m.addVar(vtype='B', name='c'+str(i)) for i in range(n)]
    for i in range(n):
        m.addGenConstrIndicator(c[i], True, gp.quicksum(x[j] * adj[bases[j], i] for j in range(len(bases))) >= sufficient[i])
    m.addConstr(gp.quicksum(c) >= math.ceil(min_coverage*n))
    m.update()
    m.optimize()

    allocation = [m.getVarByName('x'+str(i)).X for i in range(len(x))]
    coverage = utrecht.allocation_coverage(allocation)

    return allocation, coverage

# Initial feasible column.
allocation, coverage = initial_feasible_configuration()

# Keeping track of allocations and coverages.
allocations = [allocation]
coverages = [coverage]

# DEBUG: Optimal solution for 20 ambulances, 95% coverage and 3 rounds (solution = 2)
# DEBUG: Uncomment the following lines to start with the optimal solution
# allocations = [[0, 1, 0, 2, 1, 2, 0, 0, 4, 0, 2, 4, 0, 2, 2, 0, 0, 0],
#                [2, 3, 0, 1, 0, 4, 0, 0, 1, 1, 3, 0, 0, 1, 1, 3, 0, 0],
#                [2, 1, 1, 1, 0, 3, 0, 0, 2, 0, 0, 4, 2, 0, 3, 0, 1, 0]]
# coverages = [utrecht.allocation_coverage(x) for x in allocations]

while True:

    ##### MASTER PROBLEM
    mp = gp.Model('MP')
    mp.Params.OutputFlag = 0

    # Objective.
    phi = mp.addVar(name='phi', vtype='C')
    mp.setObjective(phi, gp.GRB.MAXIMIZE)

    x_vars = [mp.addVar(vtype='C') for _ in allocations]

    # The n constraints constraining the value of phi.
    phi_cons = []
    for i in range(n):
        phi_cons.append(mp.addConstr(gp.quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) - phi >= 0, name='phi'+str(i)))

    # Time horizon constraint: We must use T configurations.
    time_cons = mp.addConstr(gp.quicksum(x_vars) == num_rounds, name='time')
    
    mp.update()
    mp.optimize()
#     print('==========>>> MP objective:', mp.getObjective().getValue())
    duals = [mp.getConstrByName(c).Pi for c in ['phi'+str(i) for i in range(n)]]
#     print('==========>>> Duals:', [i for i in duals if i != 0])
    mu = mp.getConstrByName('time').Pi
#     print('==========>>> mu:', mu)


    ##### PRICING PROBLEM
    
    pp = gp.Model('PP')
    pp.Params.OutputFlag = 0
    c_vars = [pp.addVar(vtype='B', name='c'+str(i)) for i in range(n)]
    u_vars = [pp.addVar(vtype='I', lb=0, ub=max_practical_ambulances, name='u'+str(i)) for i in range(len(bases))]

    # Big M constraints (link c and u).
    bigM = num_ambulances*2
    for i in range(n):
        pp.addConstr(gp.quicksum(adj[i, bases[j]]*u_vars[j] for j in range(len(bases))) - sufficient[i] + bigM*(1-c_vars[i]) >= 0)
        pp.addConstr(gp.quicksum(adj[i, bases[j]]*u_vars[j] for j in range(len(bases))) - sufficient[i] + 1 -bigM*c_vars[i] <= 0)

    # Constraint: Don't use more than the available number of ambulances.
    pp.addConstr(gp.quicksum(u_vars) <= num_ambulances)
    
    # Constraint: The configuration must sufficiently cover 95% of the zones.
    pp.addConstr(gp.quicksum(c_vars) >= math.ceil(min_coverage*n))

    # Objective
    pp.setObjective(gp.quicksum(c_vars[i]*duals[i] for i in range(n)), gp.GRB.MINIMIZE)
    pp.update()
    pp.optimize()
#     print('==========>>> PP objective:', pp.getObjective().getValue())
    violation = - (pp.getObjective().getValue() + mu)
    print('=============>>> Dual constraint violation: ',violation)
    if violation < EPS:
        print("Violation of Dual constraint (11b) is very small. ")
    

    # Get the newly generated allocation and coverage
    cov = [int(pp.getVarByName('c'+str(i)).X) for i in range(n)]
    alloc = [int(pp.getVarByName('u'+str(i)).X) for i in range(len(bases))]
    if alloc in allocations:
        print('==========>>> Generation of an already existing column. Exiting.')
        break
#     print('==========>>> Adding new column')
    coverages.append(cov)
    allocations.append(alloc)
    print()

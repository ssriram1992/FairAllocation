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
num_rounds = 30
max_transition = 0.4 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = 4 # Maximum practical number of ambulances in a zone (doesn't seem to make much of a difference). This is 4 because a base with 4 ambulances can sufficiently cover any zone no matter what its population density
cuts = []

"""
TODO:
- The CP model finds infeasible solutions both because the graph is disconnected and because there is no hamiltonian path in a connected graph. We should be able to find a hamiltonian path in a connected graph
"""


def initial_feasible_configuration():
    """Generates one feasible allocation, along with its coverage."""
    m = gp.Model('M')
    x = [m.addVar(lb=0, ub=max_practical_ambulances, vtype='I', name='x'+str(i)) for i in range(len(bases))]

    # Constraint: Don't use more than the available number of ambulances.
    m.addConstr(gp.quicksum(x) == num_ambulances)

    # Constraint: The configuration must sufficiently cover a portion of the zones.
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



def pricing_problem(duals, mu):
    # INPUT: duals, etc
    pp = gp.Model('PP')
    pp.Params.OutputFlag = 0
    c_vars = [pp.addVar(vtype='B', name='c'+str(i)) for i in range(n)]
    u_vars = [pp.addVar(vtype='I', lb=0, ub=max_practical_ambulances, name='u'+str(i)) for i in range(len(bases))]

    # Big M constraints (link c and u). TODO: put indicator constraints instead
    bigM = num_ambulances*2
    for i in range(n):
        pp.addConstr(gp.quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) - sufficient[i] + bigM*(1-c_vars[i]) >= 0)
        pp.addConstr(gp.quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) - sufficient[i] + 1 -bigM*c_vars[i] <= 0)

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
    # print('=============>>> Dual constraint violation: ',violation)
    # if violation < EPS:
    #     print("Violation of Dual constraint (11b) is very small. ")
    

    # Get the newly generated allocation and coverage
    cov = [int(pp.getVarByName('c'+str(i)).X) for i in range(n)]
    alloc = [int(pp.getVarByName('u'+str(i)).X) for i in range(len(bases))]

    return cov, alloc



def master_problem(cuts, integrality=False):
    """Creates and solves the master problem.

    Returns obj, duals, mu
    """
    ##### MASTER PROBLEM
    mp = gp.Model('MP')
    mp.Params.OutputFlag = 0
    
    var_type='I' if integrality else 'C'
    
    # Objective.
    phi = mp.addVar(name='phi', vtype=var_type)
    mp.setObjective(phi, gp.GRB.MAXIMIZE)
    x_vars = mp.addVars(len(allocations), vtype=var_type, name='x')

    # The n constraints constraining the value of phi.
    phi_cons = []
    for i in range(n):
        phi_cons.append(mp.addConstr(gp.quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) - phi >= 0, name='phi'+str(i)))

    # Add cuts
    print('num cuts:', len(cuts))
    for cut in cuts:
        print(f'adding cut {cut}')
        mp.addConstr(gp.quicksum([x_vars[i] for i in cut]) <= num_rounds-1)
        
    # Time horizon constraint: We must use T configurations.
    time_cons = mp.addConstr(gp.quicksum(x_vars) == num_rounds, name='time')    
    mp.update()
    mp.optimize()
    
    obj = mp.getObjective().getValue()
    if integrality:
        duals = -1
        mu = -1
    if not integrality:
        duals = [mp.getConstrByName(c).Pi for c in ['phi'+str(i) for i in range(n)]]
        mu = mp.getConstrByName('time').Pi
    x_res = [mp.getVarByName('x['+str(i)+']').X for i in range(len(allocations))]

    return obj, x_res, duals, mu





optimal = False
relaxed_obj = -1
while True:

    if optimal:
        obj, x_res, duals, mu = master_problem(cuts, True)
    else:
        obj, x_res, duals, mu = master_problem(cuts, False)
    
    if not optimal:
        relaxed_obj = obj
    print('==========>>> MP objective:', obj)
    if optimal:
        x_res = [(i, int(x_res[i])) for i in range(len(x_res))]
        V = [i[0] for i in x_res if i[1] > 0]
        V0 = [i for i in range(len(V))]
        E = []
        for i in range(len(V0)-1):
            for j in range(i+1, len(V0)):
                a1 = allocations[V[V0[i]]]
                a2 = allocations[V[V0[j]]]
                max_amb = math.floor(num_ambulances*max_transition) # max number of ambulances that can transition between two configurations
                a_diff = [abs(a1[x]-a2[x]) for x in range(len(a1))]
                if sum(a_diff)/2 < max_amb:
                    E.append((V0[i], V0[j]))
                    E.append((V0[j], V0[i]))
        c = [coverages[i] for i in range(len(coverage)) if i in V]
        c = [list(x) for x in zip(*c)]

        import cp
        cp_obj, cp_sol = cp.cp_solve(V0, E, c, num_rounds, obj)
        print('CP OBJ:', cp_obj)
        if cp_obj < 0:
            print('cp infeasible for solution:', [i for i in x_res if i[1] > 0])
            optimal = False
            cuts.append(V)
            continue
        else:
            print('================================================================')
            print(f'CP SOL is {cp_obj} and allocation is {cp_sol}')
            #print('Variables used:', x_res)
            print('OBJ=', obj)
            # if obj == math.floor(relaxed_obj):
            #     print('Root node is integer optimal')
            # else:
            #     print('Root node IS NOT integer optimal')
            break



    ##### PRICING PROBLEM
    
    cov, alloc = pricing_problem(duals, mu)
    if alloc in allocations:
        print('==========>>> Generation of an already existing column. Exiting.')
        optimal = True
#     print('==========>>> Adding new column')
    coverages.append(cov)
    allocations.append(alloc)
    print()





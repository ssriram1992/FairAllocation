# TODO:
# - Add the 1-config graph inequalities proposed by Sriram.
# - Add the cuts as lazy cuts instead of starting from scratch each time (actually, this is not possible? See https://support.gurobi.com/hc/en-us/articles/360013197972-How-do-I-implement-lazy-constraints-in-Gurobi-).
# - Replace OR-Tools by MiniZinc: It allows a CP-based pricing problem, and GCC for the other CP model.
# - Find a better way to determine what is the sufficient coverage associated with a population density?


import gurobipy as gp
import math
import networkx
import numpy
from ortools.sat.python import cp_model
from utrecht_scip import *
import sys


### PARAMETERS #################################################################

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

# Other parameters.
num_ambulances = 20
num_rounds = 30
max_transition = 0.45 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = max(sufficient)
mp_integer = True # True: The once the LP is optimal, the MP is solved with integer constraints on these columns, and this solution is given to the CP model. False: The LP solution is given to the CP model.
verbose = True

################################################################################


def initial_feasible_configuration():
    """Generates one feasible allocation, along with its coverage.

    Returns:
      A feasible allocation and its coverage.
    """
    if verbose:
        print('Generating the initial column')
        
    m = gp.Model()
    x = m.addVars(len(bases), lb=0, ub=max_practical_ambulances, vtype=gp.GRB.INTEGER, name='x')

    # Constraint: Don't use more than the available number of ambulances.
    m.addConstr(gp.quicksum(x) <= num_ambulances)

    # Constraint: The configuration must sufficiently cover a portion of the zones.
    c = m.addVars(n, vtype=gp.GRB.BINARY)
    for i in range(n):
        m.addGenConstrIndicator(c[i], True, gp.quicksum(x[j]*adj[bases[j], i] for j in range(len(bases))) >= sufficient[i])
        m.addGenConstrIndicator(c[i], False, gp.quicksum(x[j]*adj[bases[j], i] for j in range(len(bases))) <= sufficient[i] - 1)
    m.addConstr(gp.quicksum(c) >= math.ceil(min_coverage*n))
    
    m.update()
    m.optimize()
    allocation = [m.getVarByName(f'x[{i}]').X for i in range(len(x))]
    coverage = utrecht.allocation_coverage(allocation)

    return allocation, coverage


# def pricing_problem_cp(duals, mu):
#     model = cp_model.CpModel()
#     u_vars = [model.NewIntVar(0, max_practical_ambulances, 'u'+str(i)) for i in range(len(bases))]
#     c_vars = [model.NewBoolVar('') for i in range(n)]

#     for i in range(n):
#         model.Add(sum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) >= sufficient[i]).OnlyEnforceIf(c_vars[i])
#         model.Add(sum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) <= sufficient[i]-1).OnlyEnforceIf(c_vars[i].Not())

#     # Constraint: Don't use more than the available number of ambulances.
#     model.Add(sum(u_vars) <= num_ambulances)
    
#     # Constraint: The configuration must sufficiently cover 95% of the zones.
#     model.Add(sum(c_vars) >= math.ceil(min_coverage*n))

#     # Don't generate a column which already exists.
#     model.AddForbiddenAssignments(u_vars, [tuple(int(j) for j in i) for i in allocations])
    
#     # Objective
#     duals = [int(i*100) for i in duals]
#     model.Add(sum(c_vars[i]*duals[i] for i in range(n)) < -1)

#     solver = cp_model.CpSolver()
#     status = solver.Solve(model)
#     s = {cp_model.OPTIMAL: 'optimal',
#          cp_model.FEASIBLE: 'feasible',
#          cp_model.INFEASIBLE: 'infeasible',
#          cp_model.MODEL_INVALID: 'model invalid',
#          cp_model.UNKNOWN: 'unknown'}
#     print(s[status])

#     # Get the newly generated allocation and coverage
#     cov = [int(solver.Value(c_vars[i])) for i in range(n)]
#     alloc = [int(solver.Value(u_vars[i])) for i in range(len(bases))]
#     print(alloc)
#     return alloc, cov, []


def pricing_problem(duals, mu):
    """Generates a new column.
    
    Args:
      duals: List of the duals of the master problem.
      mu: Dual value of mu.

    Returns:
      A new allocation and its coverage.
    """
    pp = gp.Model()
    # pp.setParam(gp.GRB.Param.IntFeasTol, 1e-09)
    # pp.setParam(gp.GRB.Param.FeasibilityTol, 1e-09)
    # pp.setParam(gp.GRB.Param.OptimalityTol, 1e-09)
    pp.Params.OutputFlag = 0
    c_vars = pp.addVars(n, vtype=gp.GRB.BINARY, name='c')
    u_vars = pp.addVars(len(bases), vtype=gp.GRB.INTEGER, lb=0, ub=max_practical_ambulances, name='u')

    # Constraints: Link c and u.
    for i in range(n):
        pp.addGenConstrIndicator(c_vars[i], True, gp.quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) >= sufficient[i])
        pp.addGenConstrIndicator(c_vars[i], False, gp.quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) <= sufficient[i]-1)

    # Constraint: Don't use more than the available number of ambulances.
    pp.addConstr(gp.quicksum(u_vars) <= num_ambulances)
    
    # Constraint: The configuration must sufficiently cover the zones.
    pp.addConstr(gp.quicksum(c_vars) >= math.ceil(min_coverage*n))

    # Objective.
    pp.setObjective(gp.quicksum(c_vars[i]*duals[i] for i in range(n)), gp.GRB.MINIMIZE)

    pp.update()
    pp.optimize()    
    allocation = [int(pp.getVarByName(f'u[{i}]').X) for i in range(len(bases))]
    coverage = [int(pp.getVarByName(f'c[{i}]').X) for i in range(n)]
    c_temp = [pp.getVarByName(f'c[{i}]').X for i in range(n)]
    a_temp = [pp.getVarByName(f'u[{i}]').X for i in range(len(bases))]

    violation = -(pp.getObjective().getValue()+mu)
    if verbose:
        print(f'Generated a column (violation {violation})')
    
    return allocation, coverage, c_temp, a_temp


def master_problem(cuts, integrality=False):
    """Creates and solves the master problem.

    Args:
      cuts: The cuts generated by the CP model for the master problem.
      integrality: Enforce integrality constraints.

    Returns:
      Objective value, solution, duals, mu
    """
    mp = gp.Model()
    mp.Params.OutputFlag = 0

    var_type=gp.GRB.INTEGER if integrality else gp.GRB.CONTINUOUS
    
    # Objective.
    phi = mp.addVar(name='phi', vtype=var_type)
    mp.setObjective(phi, gp.GRB.MAXIMIZE)
    x_vars = mp.addVars(len(allocations), vtype=var_type, name='x')

    # The n constraints constraining the value of phi.
    mp.addConstrs((gp.quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) - phi >= 0 for i in range(n)), name='phi')

    # Add the CP-generated cuts.
    for cut in cuts:
        mp.addConstr(gp.quicksum([x_vars[i] for i in cut]) <= num_rounds-1)
        
    # Time horizon constraint: We must use T configurations.
    mp.addConstr(gp.quicksum(x_vars) == num_rounds, name='time')

    mp.update()
    mp.optimize()

    if integrality:
        duals = None
        mu = None
    if not integrality:
        duals = [mp.getConstrByName(c).Pi for c in [f'phi[{i}]' for i in range(n)]]
        mu = mp.getConstrByName('time').Pi

    obj = mp.getObjective().getValue()
    x_res = [mp.getVarByName('x['+str(i)+']').X for i in range(len(allocations))]
    
    return obj, x_res, duals, mu


def cp_solve(V, E, col_cov, cuts=[]):
    """Solves a partial problem with a CP model.
    
    Args:
      V: List of vertices (columns).
      E: List of edges (if a transition between two columns is allowed).
      col_cov: Matrix of the zone coverages of the columns (c[i][j] == 1 if
               zone i is covered by column j).

    Returns:
      - Objective value of the best Hamiltonian path, -1 if there is no
        Hamiltonian path, -2 if the graph is not connected.
      - A feasible solution for this objective value.
    """
    num_cols = len(V)
    num_zones = len(col_cov)
    
    # First, check if the graph is disconnected (in which case no
    # Hamiltonian path exists).
    G = networkx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    # If the graph is not connected, no Hamiltonian path can exist.
    if not networkx.is_connected(G):
        return -2, []

    # Variables.
    model = cp_model.CpModel()
    x = [model.NewIntVar(0, num_cols-1, 'x'+str(i)) for i in range(num_rounds)]

    # Alternative for GCC, since the constraint is not available in OR-Tools.
    x_occs = []
    for i in range(num_cols):
        occs = []
        for j in range(num_rounds):
            boolvar = model.NewBoolVar('')
            model.Add(x[j] == i).OnlyEnforceIf(boolvar)
            model.Add(x[j] != i).OnlyEnforceIf(boolvar.Not())
            occs.append(boolvar)
        x_occs.append(sum(occs))
        # if mp_integer:
        #     model.AddLinearConstraint(x_occs[i], 1, num_rounds-num_cols+1)

    # Add the CP cuts.
    for cut in cuts:
        model.Add(sum(x_occs[i] for i in range(num_cols) if i in cut) <= num_rounds-1)

    # Objective.
    phi = model.NewIntVar(int(lb), math.floor(ub), 'phi')
    coverages = [model.NewIntVar(0, num_rounds, 'c'+str(i))
                 for i in range(num_zones)]
    for i in range(num_zones):
        model.Add(cp_model.LinearExpr.ScalProd(x_occs, col_cov[i]) == coverages[i])
    model.AddMinEquality(phi, coverages)
    model.Maximize(phi)

    # Regular constraint (Hamiltonian path).
    # For the initial state, we use a dummy node which is connected to
    # all other nodes.
    dummy = max(V)+1
    start = dummy
    end = V
    arcs = [(dummy, i, i) for i in V]
    for e in E:
        arcs.append((e[0], e[1], e[1]))
    model.AddAutomaton(x, start, end, arcs)

    # Solve the model.
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    assert status == cp_model.OPTIMAL or status == cp_model.INFEASIBLE

    if status == cp_model.OPTIMAL:
        solution = [solver.Value(x[i]) for i in range(num_rounds)]
        return solver.ObjectiveValue(), solution
    elif status == cp_model.INFEASIBLE:
        return -1, []


# Initial feasible column.
allocation, coverage = initial_feasible_configuration()

# Keep track of allocations and coverages.
allocations = [allocation]
coverages = [coverage]

mp_cuts = []
lb = 0
ub = math.inf
lp_optimal = False
lp_obj = -1
while True:
    if lp_optimal and mp_integer:
        obj, x_res, duals, mu = master_problem(mp_cuts, True)
    else:
        obj, x_res, duals, mu = master_problem(mp_cuts, False)
        lp_obj = obj
    if lp_optimal:
        ub = math.floor(obj)

    if verbose:
        print(f'LP obj: {round(lp_obj, 2)},\tLB: {lb},\tUB: {ub},\t{len(x_res)} columns,\t{len(mp_cuts)} cuts')

    if lp_optimal:
        x_res = [(i, x_res[i]) for i in range(len(x_res))]
        V = [i[0] for i in x_res if i[1] > 0]
        if verbose:
            print(f'Solving the CP model with {len(V)} variables.')

        # Create the CP cuts.
        V0 = [i for i in range(len(V))] # States of the automaton must start at 0.
        cp_cuts = [] # CP cuts must be adjusted to start form 0, like V0.
        for cut in mp_cuts:
            intersection = list(set(cut) & set(V))
            if intersection == []:
                continue
            cp_cut = []
            for node in intersection:
                cp_cut.append(V0[V.index(node)])
            cp_cuts.append(cp_cut)

        # Create the graph for the CP regular constraint.
        E = []
        for i in range(len(V0)-1):
            for j in range(i+1, len(V0)):
                a1 = allocations[V[V0[i]]]
                a2 = allocations[V[V0[j]]]
                # Maximum number of ambulances that can transition between two configurations.
                max_amb = math.floor(num_ambulances*max_transition)
                a_diff = [abs(a1[x]-a2[x]) for x in range(len(a1))]
                if sum(a_diff)/2 < max_amb:
                    E.append((V0[i], V0[j]))
                    E.append((V0[j], V0[i]))

        # Get the coverages of the subset of columns which take nonzero values.
        c = [coverages[i] for i in range(len(coverage)) if i in V]
        c = [list(x) for x in zip(*c)]

        # Solve the CP model.
        cp_obj, cp_sol = cp_solve(V0, E, c, cp_cuts)
        lb = max(lb, cp_obj)
        if verbose:
            if cp_obj == -2:
                print(f'CP solution infeasible: Graph is not connected.')
            elif cp_obj == -1:
                print(f'CP solution infeasible: Graph is connected, but contains no Hamiltonian path.')
            else:
                print(f'CP solution is feasible with value {cp_obj}.')

        if cp_obj < math.floor(ub):
            lp_optimal = False
            mp_cuts.append(V)
            continue
        else:
            print('Optimality is proven.')
            check_allocations = []
            for i in x_res:
                for j in range(int(i[1])):
                    check_allocations.append(allocations[i[0]])
            print(f'there are {len(check_allocations)} allocations')
            check_value = min(utrecht.allocations_coverage(check_allocations))
            print(f'Checking value of optimal solution: {check_value}')
            assert check_value == ub
            break
    
    allocation, coverage, c_temp, a_temp = pricing_problem(duals, mu)
    # assert coverage == utrecht.allocation_coverage(allocation)
    if not (coverage == utrecht.allocation_coverage(allocation)):
        print('Real coverage:', utrecht.allocation_coverage(allocation))
        print('PP coverage:', coverage)
        print('PP coverage (raw):', c_temp)
        print('PP allocation:', allocation)
        print('PP allocation (raw):', a_temp)
        print('Pricing problem outputs wrong coverage!')
        exit(1)
    if allocation in allocations:
        print('Generated column already exists.')
        lp_optimal = True
    coverages.append(coverage)
    allocations.append(allocation)


        # check_allocations = []
        # for i in range(len(x_res)):
        #     for j in range(int(x_res[i])):
        #         check_allocations.append(allocations[i])
        # print(f'there are {len(check_allocations)} allocations')
        # check_value = min(utrecht.allocations_coverage(check_allocations))
        # print(f'Checking value of optimal solution: {check_value}')
        # exit()




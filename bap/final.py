# TODO:
# - Add the 1-config graph inequalities proposed by Sriram.
# - Add the cuts as lazy cuts instead of starting from scratch each time (actually, this is not possible? See https://support.gurobi.com/hc/en-us/articles/360013197972-How-do-I-implement-lazy-constraints-in-Gurobi-).
# - Replace OR-Tools by MiniZinc: It allows a CP-based pricing problem (although this may not help much), and GCC for the other CP model.
# - Find a better way to determine what is the sufficient coverage associated with a population density?


import gurobipy as gp
import math
import networkx
import numpy
from ortools.sat.python import cp_model
from utrecht_scip import *
import sys
import time
# for i in {20..1}; do mt=`bc -l <<< $i/20`; echo $mt; python final.py $mt > 0$mt.txt; done


### PARAMETERS #################################################################

EPS = 1e-6

randomize_graph = False
seed = 0

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes (900 seconds)".
utrecht = Utrecht(randomize_graph, seed)
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
cl = True # Accept arguments from the command line
if cl:
    max_transition = float(sys.argv[1])
else:
    max_transition = 0.5 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = max(sufficient)
mp_integer = False # True: The once the LP is optimal, the MP is solved with integer constraints on these columns, and this solution is given to the CP model. False: The LP solution is given to the CP model.
verbose = True
time_limit = 300

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
    allocation = [round(m.getVarByName(f'x[{i}]').X) for i in range(len(x))]
    coverage = utrecht.allocation_coverage(allocation)

    return allocation, coverage


def pricing_problem(duals, mu):
    """Generates a new column.
    
    Args:
      duals: List of the duals of the master problem.
      mu: Dual value of mu.

    Returns:
      A new allocation and its coverage.
    """
    pp = gp.Model()
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
    allocation = [round(pp.getVarByName(f'u[{i}]').X) for i in range(len(bases))]
    coverage = [round(pp.getVarByName(f'c[{i}]').X) for i in range(n)]
    assert coverage == utrecht.allocation_coverage(allocation)
    
    violation = -(pp.getObjective().getValue()+mu)
    if verbose:
        print(f'Generated a column (violation {violation})')
    
    return allocation, coverage


def master_problem(cuts, lb=0, ub=num_rounds, integrality=False):
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
    phi_high = mp.addVar(name='phi_high', vtype=gp.GRB.CONTINUOUS)
    phi_low = mp.addVar(name='phi_low', vtype=gp.GRB.CONTINUOUS)
    phi = mp.addVar(name='phi', vtype=gp.GRB.CONTINUOUS)
    mp.addConstr(phi == phi_high-phi_low)
    mp.setObjective(phi, gp.GRB.MINIMIZE)
    x_vars = mp.addVars(len(allocations), vtype=var_type, name='x')

    # The n constraints constraining the value of phi.
    mp.addConstrs((gp.quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) <= phi_high for i in range(n)), name='phi_high_cstr')
    mp.addConstrs((gp.quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) >= phi_low for i in range(n)), name='phi_low_cstr')

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
        duals_low = [mp.getConstrByName(c).Pi for c in [f'phi_low_cstr[{i}]' for i in range(n)]]
        duals_high = [mp.getConstrByName(c).Pi for c in [f'phi_high_cstr[{i}]' for i in range(n)]]
        mu = mp.getConstrByName('time').Pi

    obj = mp.getObjective().getValue()
    x_res = [mp.getVarByName('x['+str(i)+']').X for i in range(len(allocations))]
    
    return obj, x_res, [sum([-x[0], -x[1]]) for x in zip(duals_low, duals_high)], mu


def cp_solve(V, E, col_cov, cuts=[]):
    """Solves a partial problem with a CP model.
    
    Args:
      V: List of vertices (columns).
      E: List of edges (if a transition between two columns is allowed).
      col_cov: Matrix of the zone coverages of the columns (c[i][j] == 1 if
               zone i is covered by column j).

    Returns:
      - Objective value of the best Hamiltonian path, -1 if there is no
        Hamiltonian path within the LB/UB limits, -2 if the graph is not
        connected (this latter case has been removed).
      - A feasible solution for this objective value.
    """
    num_cols = len(V)
    num_zones = len(col_cov)
    
    # First, check if the graph is disconnected (in which case no
    # Hamiltonian path exists).
    G = networkx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    # # If the graph is not connected, no Hamiltonian path can exist.
    # if not networkx.is_connected(G):
    #     return -2, []

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
    phi = model.NewIntVar(int(lb), math.floor(ub)-1, 'phi')
    coverages = [model.NewIntVar(0, num_rounds, 'c'+str(i))
                 for i in range(num_zones)]
    for i in range(num_zones):
        model.Add(cp_model.LinearExpr.ScalProd(x_occs, col_cov[i]) == coverages[i])
    phi_low = model.NewIntVar(0, num_rounds, 'phi_low')
    phi_high = model.NewIntVar(0, num_rounds, 'phi_high')
    model.AddMinEquality(phi_low, coverages)
    model.AddMaxEquality(phi_high, coverages)
    model.Add(phi == phi_high-phi_low)
    model.Minimize(phi)

    # Regular constraint (Hamiltonian path).
    # For the initial state, we use a dummy node which is connected to
    # all other nodes.
    dummy = max(V)+1
    start = dummy
    end = V
    arcs = [(dummy, i, i) for i in V]
    for e in E:
        arcs.append((e[0], e[1], e[1]))
    # If there is only one vertex then a Hamiltonian path exists.
    if len(V) > 1:
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


start_time = time.time()
    
# Initial feasible column.
allocation, coverage = initial_feasible_configuration()

# Keep track of allocations and coverages.
allocations = [allocation]
coverages = [coverage]

mp_cuts = []
lb = 0
ub = num_rounds+1
lp_optimal = False
lp_obj = -1
while True:
    if time.time() - start_time > time_limit:
        break
    if lp_optimal and mp_integer:
        obj, x_res, duals, mu = master_problem(mp_cuts, lb, ub, True)
    else:
        obj, x_res, duals, mu = master_problem(mp_cuts, lb, ub, False)
        lp_obj = obj
    if lp_optimal:
        lb = math.ceil(obj)

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
        c = [coverages[i] for i in range(len(coverages)) if i in V]
        c = [list(x) for x in zip(*c)]

        # Solve the CP model.
        cp_obj, cp_sol = cp_solve(V0, E, c, cp_cuts)
        if verbose:
            if cp_obj == -2:
                print(f'CP solution infeasible: Graph is not connected.')
            elif cp_obj == -1:
                print(f'CP solution infeasible: Graph is connected, but contains no Hamiltonian path within the LB/UB limits.')
            else:
                print(f'CP solution is feasible with value {cp_obj}.')
                lb = min(lb, cp_obj)

        if cp_obj > math.ceil(lb):
            lp_optimal = False
            mp_cuts.append(V)
            ub = min(cp_obj, ub)
            continue
        else:
            # If the CP model outputs an infeasible solution, add the cut and continue
            if cp_obj == -1:
                lp_optimal = False
                mp_cuts.append(V)
                continue
            check_allocations = []
            for i in cp_sol:
                check_allocations.append(allocations[V[V0[i]]])
            check_value = max(utrecht.allocations_coverage(check_allocations))-min(utrecht.allocations_coverage(check_allocations))
            assert check_value == lb, f'{check_value} vs {lb}'
            print('Optimality is proven.')
            break
    
    allocation, coverage = pricing_problem(duals, mu)
    if allocation in allocations:
        print('Generated column already exists.')
        lp_optimal = True
    coverages.append(coverage)
    allocations.append(allocation)

print(f'{len(allocations)},{len(mp_cuts)},{lb},{ub}')
print(time.time()-start_time)

# This CP model returns an integer indicating the value of the best Hamiltonian
# path for the pattern passed as input, or -1 if no such path exists.


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from ortools.sat.python import cp_model
import math
import networkx


# # sys.argv to get the file where the args are stored
# #num_cols = 4
# #num_zones = 5
# num_rounds = 30
# ip_ub = 999 #(the ub of the ip model)
# # IN THIS EXAMPLE: the graph is x0-x1-x2
# #                                  |
# #                                  x3
# V = [i for i in range(4)]
# E = [(0, 1), (1, 0),
#      (1, 2), (2, 1),
#      (3, 1), (1, 3),
#      (0, 0), (1, 1), (2, 2), (3, 3)]
# import random
# import sys
# random.seed(13)#int(sys.argv[1]))
# cov = [[random.randint(0, 1) for _ in range(4)] for _ in range(5)]

# print(cov)

def cp_solve(V, E, cov, num_rounds, ip_ub):
    """Solves the problem with a CP model.
    
    Args:
      V: List of vertices (columns).
      E: List of edges (if a transition between two columns is allowed).
      cov: Matrix of the zone coverages of the columns (c[i][j] == 1 if zone i
           is covered by column j).
      num_rounds: Number of rounds (time horizon).
      ip_ub: Upper bound of the IP node that the CP model is solving.

    Returns:
      Value of the best Hamiltonian path, -1 if no path exists.
    """
    num_cols = len(V)
    num_zones = len(cov)
    
    # First, check if the graph is disconnected (in which case no
    # Hamiltonian path exists).
    G = networkx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    # If the graph is not connected, no Hamiltonian path can exist.
    if not networkx.is_connected(G):
        return -1

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
        model.AddLinearConstraint(x_occs[i], 1, num_rounds-num_cols+1)

    # Objective (it cannot be greater than the UB of the IP model).
    phi = model.NewIntVar(0, math.floor(ip_ub), 'phi')
    coverages = [model.NewIntVar(0, math.floor(ip_ub), 'c'+str(i))
                 for i in range(num_zones)]
    for i in range(num_zones):
        model.Add(cp_model.LinearExpr.ScalProd(x_occs, cov[i]) == coverages[i])
    model.AddMinEquality(phi, coverages)
    model.Maximize(phi)

    # Regular constraint (Hamiltonian path).
    # For the initial state, we use a dummy node which is connected to
    # all other nodes.
    dummy = num_cols
    start = dummy
    end = V
    arcs = [(dummy, i, i) for i in V]
    for e in E:
        arcs.append((e[0], e[1], e[1]))
    model.AddAutomaton(x, start, end, arcs)

    # Solve the model.
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    s = {cp_model.OPTIMAL: 'optimal',
         cp_model.FEASIBLE: 'feasible',
         cp_model.INFEASIBLE: 'infeasible',
         cp_model.MODEL_INVALID: 'model invalid',
         cp_model.UNKNOWN: 'unknown'}

    assert status == cp_model.OPTIMAL or status == cp_model.INFEASIBLE

    # for i in range(num_rounds):
    #     print(f'x{i} = {solver.Value(x[i])}')
    # for i in range(num_cols):
    #     print(f'occ{i} = {solver.Value(x_occs[i])}')
    
    if status == cp_model.OPTIMAL:
        return solver.ObjectiveValue()
    elif status == cp_model.INFEASIBLE:
        return -1


print(cp_solve(V, E, cov, num_rounds, ip_ub))

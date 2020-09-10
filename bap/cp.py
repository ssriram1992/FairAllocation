from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from ortools.sat.python import cp_model
import networkx

"""
Suppose we have 100 columns and 30 rounds.

Naive approach:
- Graph of 100*30=3000 nodes since in theory we could visit the same configuration 30 times, and we need a graph that can accomodate a Hamiltonian circuit.

TSP/succ approach:
- 30 variables with domains of 100 values each (columns).
- successor dynamic to make sure a path exists.

IP->diversity->CP
- Find an optimal solution with a certain pattern (eg: columns X, Y, and Z are used, we don't care how many of each column is used)
- Using CP, check for an optimal solution with a hamiltonian path on this small problem
- Add a cut for this pattern and add it to the MP
- Goto 1
- Let's say that at some point

IP/CP method:

Hi everyone,

I thought of an alternative way to solve our problem:

1. Explore the branch-and-price tree normally. When an integer solution is found: Take the "pattern" of that solution and give it to the CP model of step 2. The pattern is the subset of columns which are used at least once by the integer solution (here, we don't care exactly how many times each one is used).
2. The CP model finds the Hamiltonian path with the best objective value for this pattern (each column in the pattern must be used at least once by the CP model). If this objective value is equal to the upper bound, the problem is solved. If it is lower, update the lower bound. If no Hamiltonian path exists, the pattern is infeasible.
3. In all those cases, add a cut for this pattern to the IP. This cut can be propagated to the whole branch-and-price tree since the pattern has been fully solved by the CP model. Go to 1.

I've looked at a sample of the solution pool for integer optimal solutions, and these solutions generally use about 5 to 10 columns, regardless of the number of rounds. This means that the CP model would deal with a smaller and more manageable number of variables.

Let me know what you think.

Philippe

"""


def SolveWithTimeLimitSampleSat():
    """Minimal CP-SAT example to showcase calling the solver."""
    # Creates the model.
    model = cp_model.CpModel()
    # Creates the variables.
    num_cols = 3
    num_zones = 5
    num_rounds = 5
    x = [model.NewIntVar(0, num_cols - 1, 'x'+str(i)) for i in range(num_rounds)]

    # Alternative for GCC, since the constraint is not available in OR-Tools.
    for i in range(num_cols):
        occs = []
        for j in range(num_rounds):
            boolvar = model.NewBoolVar('')
            model.Add(x[j] == i).OnlyEnforceIf(boolvar)
            model.Add(x[j] != i).OnlyEnforceIf(boolvar.Not())
            occs.append(boolvar)
        model.AddLinearConstraint(sum(occs), 1, num_rounds-num_cols+1)

    
        

    model.Maximize(sum(x))
    
    # Adds an all-different constraint.
    #model.Add(x[0] != x[1])
    #model.AddAllDifferent(x)

    # IN THIS EXAMPLE: the graph is x0-x1-x2



    # V = [0, 1, 2]
    # E = [(0, 1), (1, 0),
    #      (1, 2), (2, 1)]

    # # If the graph is not connected, there is definitely no Hamiltonian path.
    # G_test_conn = networkx.Graph()
    # G_test_conn.add_nodes_from(V)
    # G_test_conn.add_edges_from(E)
    # if not networkx.is_connected(G_test_conn):
    #     print('Graph not connected: No Hamiltonian path!')
    #     exit(0)

    # # We know that the graph is connected. We add a dummy node which is connected
    # # to every node, since the constraint is looking for a Hamiltonian circuit.
    # V.append(len(V))
    # for i in range(len(V)-1):
    #     dummy = len(V)-1
    #     E.append((i, dummy))
    #     E.append((dummy, i))

    # # The Hamiltonian circuit constraint.
    # lits = [model.NewBoolVar(f'b_{E[i][0]},{E[i][1]}') for i in range(len(E))]
    # arcs = []
    # for i in range(len(E)):
    #     arcs.append((E[i][0], E[i][1], lits[i]))
    # model.AddCircuit(arcs)

    # Creates a solver and solves the model.
    solver = cp_model.CpSolver()

    # Sets a time limit of 10 seconds.
    #solver.parameters.max_time_in_seconds = 10.0

    status = solver.Solve(model)
    # https://developers.google.com/optimization/cp/cp_solver#cp-sat-return-values
    s = {cp_model.OPTIMAL: 'optimal',
         cp_model.FEASIBLE: 'feasible',
         cp_model.INFEASIBLE: 'infeasible',
         cp_model.MODEL_INVALID: 'model invalid',
         cp_model.UNKNOWN: 'unknown'}
    print("status:", s[status])

    if status == cp_model.OPTIMAL:# or status ==:
        print('x0 = %i' % solver.Value(x[0]))
        print('x1 = %i' % solver.Value(x[1]))
        print('x2 = %i' % solver.Value(x[2]))
        print('x3 = %i' % solver.Value(x[3]))
        print('x4 = %i' % solver.Value(x[4]))
        # for i in range(len(E)):
        #     print(arcs[i], solver.Value(lits[i]))



SolveWithTimeLimitSampleSat()

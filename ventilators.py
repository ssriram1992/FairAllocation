import gurobipy as gp
import ventutils


########
# DATA #
########

# Q[i][t]: Demand at location i at time t.
Q = [list(x) for x in zip(*ventutils.ode_data())] # ventutils.ihme_data()

# n: Number of locations.
n = len(Q)

# T: Time horizon.
T = len(Q[0])

# D[i][j]: Delay in transfer from location i to location j.
# TODO: Make this relative to the distance between locations (1000 miles = 1 day?)
D = [[3 for _ in range(n)] for _ in range(n)]


##############
# PARAMETERS #
##############

# V[i]: Initial supply of ventilators at location i.
V = ventutils.initial_supply()

# R: Recovery time.
R = 10


#############
# VARIABLES #
#############

M = gp.Model()

# x[i, t]: Number of ventilators at location i at time t.
x = M.addVars(n, T, lb=0, vtype=gp.GRB.INTEGER, name='x')

# N[t]: Number of ventilators injected into the system at time t.
N = M.addVars(T, lb=0, vtype=gp.GRB.INTEGER, name='N')

# z[i, t]: How many of the new ventilators N[t] are sent to location i at time t.
z = M.addVars(n, T, lb=0, vtype=gp.GRB.INTEGER, name='z')

# f[i, j, t]: How many ventilators location i sends to location j at time t.
D_max = max([max(i) for i in D])
f = M.addVars(n, n, range(-D_max, T),
              lb=0, vtype=gp.GRB.INTEGER, name='f')
# Negative indices are required for checking t-D[i][j] when t is small.
# For these indices f is zero.
for t in range(-D_max, 0):
    for i in range(n):
        for j in range(n):
            M.addConstr(f[i, j, t] == 0)

# tau[i, t]: Shortage of ventilators at location i at time t.
tau = M.addVars(n, T, lb=0, vtype=gp.GRB.INTEGER, name='tau')

# phi: Objective.
phi = M.addVar(obj=1, lb=0, vtype=gp.GRB.INTEGER, name='phi')


###############
# CONSTRAINTS #
###############

# The number of ventilators in the system at time t should not be greater than
# the combination of the initial supply of ventilators and what has been
# injected into the system since the beginning.
for t in range(T):
    M.addConstr(gp.quicksum(x[i, t] for i in range(n)) <=
                sum(V) + gp.quicksum(N[tp] for tp in range(t+1)))

# The number of new ventilators assigned at time t should not be greater than
# the number of new ventilators injected into the system at time t.
for t in range(T):
    M.addConstr(gp.quicksum(x[i, t] for i in range(n)) <= N[t])

# The number of ventilators in the various locations should be consistent over
# time, including consistency with the new ventilators injected into the system
# and with the transfer of ventilators between locations.
for i in range(n):
    for t in range(1, T):
        M.addConstr(x[i, t] ==
                    x[i, t-1] +
                    z[i, t] +
                    gp.quicksum(f[j, i, t-D[j][i]] for j in range(n)) -
                    gp.quicksum(f[i, j, t] for j in range(n)))

# A ventilator cannot be moved before its recovery time.
for i in range(n):
    for t in range(1, T):
        M.addConstr(x[i, t] >=
                    gp.quicksum(z[i, tp]+gp.quicksum(f[j, i, tp-D[j][i]] for j in range(n)) for tp in range(max(1, t-R), t)))

# Keep track of the shortages of ventilators.
for i in range(n):
    for t in range(T):
        M.addConstr(tau[i, t] == gp.max_(0, Q[i][t]-x[i, t]))
        # M.addConstr(tau[i, t] == Q[i][t]-x[i, t])

        
#############
# OBJECTIVE #
#############

M.addConstr(phi == gp.quicksum(gp.max_(tau[i, t] for i in range(n)) for t in range(T)))


###########
# SOLVING #
###########

M.update()
M.optimize()
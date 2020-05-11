import gurobipy as gp
import matplotlib.pyplot as plt
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
# TODO: Make this relative to the distance between locations (1000 miles = 1 day?).
# Bertsimas chooses a lead time of 3 days.
D = [[1 for _ in range(n)] for _ in range(n)]


##############
# PARAMETERS #
##############

# V[i]: Initial supply of ventilators at location i.
V = ventutils.initial_supply()
# TODO: Playing with values, can later be changed.
V = [int(vv/10) for vv in V]

# R: Recovery time.
R = 10

# N_max: Number of ventilators available to be injected into the system each day.
# Bertsimas estimates 450 ventilators per day.
N_max = 65000/T


#############
# VARIABLES #
#############

M = gp.Model()

# x[i, t]: Number of ventilators at location i at time t.
x = M.addVars(n, T, lb=0, vtype=gp.GRB.INTEGER, name='x')

# N[t]: Number of ventilators injected into the system at time t.
# The lower bound forces at least most of the ventilators to be used.
N = M.addVars(T, lb=0.9*N_max, ub=N_max, vtype=gp.GRB.INTEGER, name='N')

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
            f[i, j, t].lb = 0
            f[i, j, t].ub = 0
# Prevent transfers from a location to itself.
for t in range(0, T):
    for i in range(n):
        f[i, i, t].lb = 0
        f[i, i, t].ub = 0

# tau[i, t]: Shortage of ventilators at location i at time t.
tau = M.addVars(n, T, lb=-gp.GRB.INFINITY, vtype=gp.GRB.INTEGER, name='tau')
tau_pos = M.addVars(n, T, lb=0, vtype=gp.GRB.INTEGER, name='tau_pos')

# phi: Objective.
phi = M.addVar(obj=1, lb=0, vtype=gp.GRB.INTEGER, name='phi')
phi_ = M.addVars(n, lb=0, vtype=gp.GRB.INTEGER, name='phi_')


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
    M.addConstr(gp.quicksum(z[i, t] for i in range(n)) == N[t])

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
        M.addConstr(tau[i, t] == Q[i][t]-x[i, t])
        M.addGenConstrMax(tau_pos[i, t], [tau[i, t]], 0)
        

#############
# OBJECTIVE #
#############

for i in range(n):
    M.addConstr(phi_[i] == gp.quicksum(tau_pos[i, t] for t in range(T)))
    
M.addGenConstrMax(phi, phi_)


###########
# SOLVING #
###########

M.update()
M.optimize()


############
# PLOTTING #
############

def plotState(loc):
    xval = [t for t in range(T)]
    yval = [x[loc, t].X for t in range(T)]
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1.2, 1.2])
    ax.plot(xval, yval, 'bo--')
    ax.plot(xval, [Q[loc][t] for t in range(T)], 'rs--')
    ax.legend(labels = ('Quantity supplied', 'Quantity demanded'), loc = 'upper right')
    ax.set_title("State "+str(loc))
    ax.set_xlabel('Day')
    ax.set_ylabel('Number of ventilators')
    plt.show()

def plot_philippe(loc):
    """For some reason, plotState() doesn't render well for me."""
    xval = [t for t in range(T)]
    yval = [x[loc, t].X for t in range(T)]
    fig, ax = plt.subplots()
    ax.plot(xval, yval, 'bo--')
    ax.plot(xval, [Q[loc][t] for t in range(T)], 'rs--')
    ax.set(xlabel='Day', ylabel='Number of ventilators', title='State'+str(loc))
    ax.grid()
    ax.legend(labels=('Quantity supplied', 'Quantity demanded'), loc = 'upper right')
    plt.show()
    
for loc in range(50):
    plot_philippe(loc)
    # plotState(loc)

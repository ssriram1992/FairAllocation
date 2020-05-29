from gurobipy import *
import gurobipy as gp
import numpy as np
import ventutils


def printVars(z, string = "", tol = 1e-7):
    for zz in z:
        if abs(z[zz].X) >= tol :
            print(string,zz,"----", z[zz].x)
            
def printParams(z, string = "", tol = 1e-7):
    for zz in z:
        if abs(z[zz]) >= tol :
            print(string,zz,"----", z[zz])

def greedy_iter(n,
                time_start,
                time_end,
                V,
                Q,
                D,
                N_max,
                eff,
                R,
                x_past,
                N_past,
                z_past,
                f_past,
                phi_past):
    """Performs an iteration of greedy algorithm between."""
    
    # Modified T for the RTH interval.
    T_ = [t for t in range(time_start, time_end)]

    # Greedy model
    M_greedy = gp.Model()

    
    #############
    # Variables #
    #############

    # x[i, t]: Number of ventilators at location i at time t.
    x = M_greedy.addVars(n, T_, lb=0, vtype=gp.GRB.INTEGER, name='x')

    # N[t]: Number of ventilators injected into the system at time t.
    # The lower bound forces at least most of the ventilators to be used.
    N = M_greedy.addVars(T_, lb=eff*N_max[time_start:time_end], ub=N_max[time_start:time_end], vtype=gp.GRB.INTEGER, name='N')
    
    # z[i, t]: How many of the new ventilators N[t] are sent to location i at time t.
    z = M_greedy.addVars (n, T_, lb=0, vtype=gp.GRB.INTEGER, name='z')
    
    # f[i, j, t]: How many ventilators location i sends to location j at time t.
    f = M_greedy.addVars(n, n, range(time_start, time_end), obj=0.001, lb=0, vtype=gp.GRB.INTEGER, name='f')

    # tau[i, t]: Shortage of ventilators at location i at time t.
    tau = M_greedy.addVars(n, T_, lb=-gp.GRB.INFINITY, vtype=gp.GRB.INTEGER, name='tau')
    tau_pos = M_greedy.addVars(n, T_, lb=0, vtype=gp.GRB.INTEGER, name='tau_pos')

    
    ###############
    # Constraints #
    ###############

    for t in T_:
        # States should not transfer ventilators to themselves.
        for i in range(n):
            f[i, i, t].ub = 0
            
        # Two states should not transfer ventilators back and forth.
        for i in range(n):
            for j in range(i):
                M_greedy.addSOS(gp.GRB.SOS_TYPE1, [f[i, j, t], f[j, i, t]])
                
        # The total number of ventilators everywhere should not be greater than the number injected.
        M_greedy.addConstr(gp.quicksum(x[i, t] for i in range(n)) <=
                           sum(V) +
                           sum(N_past[tp] for tp in range(time_start)) +
                           gp.quicksum(N[tp] for tp in range(time_start, t+1)))

        # New ventilators should necessarily be allocated. 
        M_greedy.addConstr(gp.quicksum(z[i, t] for i in range(n)) == N[t])
        
        # Ventilator shorages (modified Q for the RTH interval).
        Q_ = {t:[Q[i][t] for i in range(n)] for t in range(time_start, time_end)}
        for i in range(n):
            M_greedy.addConstr(tau[i, t] == Q_[t][i]-x[i, t])
            M_greedy.addGenConstrMax(tau_pos[i, t], [tau[i, t]], 0)

    # inbound[i, t]: The number of incoming ventilators for location i at time t.
    inbound = {}
    for i in range(n):
        for t in T_:
            inbound[(i,t)] = 0
            
            # Get the information from f or f_past, depending.
            for j in range(n):
                if t-D[j][i] >= time_start:
                    inbound[(i,t)] += f[j, i, t-D[j][i]]
                elif t-D[j][i] >= 0:
                    inbound[(i,t)] += f_past[j, i, t-D[j][i]]

            # x_prev: The x value of a state at the previous time point
            x_prev = 0
            if t == time_start and t > 0:
                x_prev += x_past[i, t-1]
            elif t > time_start:
                x_prev += x[i, t-1]

            M_greedy.addConstr(x[i, t] ==
                               x_prev +
                               z[i, t] +
                               inbound[(i,t)] -
                               gp.quicksum(f[i, j, t] for j in range(n)))
            
    # A ventilator cannot be moved before its recovery time.
    # Philippe: I haven't checked if this block is correct.
    for i in range(n):
        for t in T_:
            eqn = 0
            for tp in range(max(1, t-R), t):
                if tp in T_:
                    eqn += (z[i, tp] + inbound[(i, tp)])
                else:
                    eqn += z_past[i, tp] + sum([f_past[j, i, tp-D[j][i]] for j in range(n) if tp-D[j][i] >= 0])
            M_greedy.addConstr(x[i, t] >= eqn)

    # Objective
    # 17:30 normalization considering the phi, instead of considering the shortage, we consider the percentage of unfulfilled demand
    # 20:55
    
#     # Minimizing final disparity
#     phi_max = M_greedy.addVar(lb=0, vtype=gp.GRB.INTEGER, name='phi_max')
#     phi_min = M_greedy.addVar(lb=0,vtype=gp.GRB.INTEGER, name="phi_min")
#     phi_ = M_greedy.addVars(n, lb=0, vtype=gp.GRB.INTEGER, name='phi_')
#     for i in range(n):
#         M_greedy.addConstr(phi_[i] == phi_past[i] + gp.quicksum(tau_pos[i, t] for t in T_))
#     M_greedy.addGenConstrMin(phi_min, phi_)
#     M_greedy.addGenConstrMax(phi_max, phi_)
#     M_greedy.setObjective(phi_max-phi_min)
    # Minimizing day wise disparity

    # phi_[i, t]: Sum of the shortages in location i from the beginning up to time t.
    phi_ = M_greedy.addVars(n, T_, lb=0, vtype=gp.GRB.INTEGER, name='phi_')
    for tt in T_:
        for i in range(n):
            M_greedy.addConstr(phi_[i, tt] == phi_past[i] + gp.quicksum(tau_pos[i, t] for t in T_ if t <= tt))

    # phi_max/min[t]: max/min phi of all locations on time t.
    phi_max = M_greedy.addVars(T_, lb=0, obj=1, vtype=gp.GRB.INTEGER, name='phi_max')
    phi_min = M_greedy.addVars(T_, lb=0, obj=-1, vtype=gp.GRB.INTEGER, name="phi_min")    
    for tt in T_:
        M_greedy.addGenConstrMin(phi_min[tt], [phi_[i,tt] for i in range(n)])
        M_greedy.addGenConstrMax(phi_max[tt], [phi_[i,tt] for i in range(n)])

    # Solve.
    M_greedy.update()
    M_greedy.params.Threads = 2
    M_greedy.optimize()
    
    return M_greedy, x, N, z, f, phi_


########
# DATA #
########

# Q[i][t]: Demand at location i at time t.
Q = np.array([list(x) for x in zip(*ventutils.ode_data())]) # ventutils.ihme_data()

# n: Number of locations.
n = len(Q)

# T: Time horizon.
T = len(Q[0])

# D[i][j]: Delay in transfer from location i to location j.
# TODO: Make this relative to the distance between locations (1000 miles = 1 day?)
D = [[1 for _ in range(n)] for _ in range(n)]


##############
# PARAMETERS #
##############

# V[i]: Initial supply of ventilators at location i.
V = ventutils.initial_supply()
V = [int(v/10) for v in V]

# R: Recovery time.
R = 10

# N_max: Number of ventilators available to be injected into the system each day.
N_max = np.sum(Q, axis=0)/20 # 65000/T

# Efficiency: Fraction of the available ventilators which are required to be injected into the system each day.
eff = 0.9

# Length of a rolling time horizon (RTH) interval.
period = 5


#####################
# GREEDY ITERATIONS #
#####################

# Values of the variables in past RTH intervals.
x_past = {}
N_past = {}
z_past = {}
f_past = {}
phi_past = {i:0 for i in range(n)}

# # Initial RTH interval. (might be needed for taupast, etc)
# time_start = 0
# time_end = period

# phi_past = {(i,t):0 for (i,t) in itertools.product(range(n),range(0,time_start))}
# taupast = {(i,t):0 for (i,t) in itertools.product(range(n),range(0,time_start))}
# tau_pos_past = {(i,t):0 for (i,t) in itertools.product(range(n),range(0,time_start))}

for time in range(0, T, period):
    # Current RTH interval.
    time_start = time
    time_end = min(time+period, T)

    # Solve one greedy iteration.
    M_greedy, x, N, z, f, phi = greedy_iter(n,
                                            time_start,
                                            time_end,
                                            V,
                                            Q,
                                            D,
                                            N_max,
                                            eff,
                                            R,
                                            x_past,
                                            N_past,
                                            z_past,
                                            f_past,
                                            phi_past)

    # Save the results of the current RTH interval's solution.
    for t in range(time_start, time_end):
        for i in range(n):
            x_past[(i, t)] = x[i, t].X
            z_past[(i, t)] = z[i, t].X
            for j in range(n):
                f_past[(i, j, t)] = f[i, j, t].X
        N_past[t] = N[t].X
    for i in range(n):
        phi_past[i] = phi[i, time_end-1].X


print(phi_past)

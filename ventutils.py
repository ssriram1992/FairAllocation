# Note: We are using the IHME and ODE projections of April 9th 2020, so that
# they cover the same range of days.


import datetime
import itertools
import gurobipy as gp


def state_indices():
    """ Lists the states in alphabetical order.

    Returns:
      A list of the states in alphabetical order.
    """
    with open('data/bertsimas/state_list.csv', 'r') as f:
        data = [x.strip('"\n') for x in f.readlines()[1:]]
    f.close()
    data.sort()
    return data


def ihme_clean_data():
    """ Basic clean-up of IHME data.

    Returns:
      Cleaned-up data.
    """
    with open('data/bertsimas/projected_demand_ihme.csv', 'r') as f:
        data = [[x[0].strip('"\n'), x[1], int(x[4].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    return data


def ode_clean_data():
    """ Basic clean-up of ODE data.

    Returns:
      Cleaned-up data.
    """
    with open('data/bertsimas/projected_demand_ode.csv', 'r') as f:
        data = [[x[0].strip('"\n'), x[1], int(x[7].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    return data


def parse_data(data):
    """ Parses cleaned IHME or ODE data.

    Arguments:
      data: Clean IHME or ODE data.
    
    Returns:
      Parsed data.
    """
    dates = list(set([x[1] for x in data]))
    dates.sort(key=lambda date: datetime.datetime.strptime(date, '%Y-%m-%d'))
    states = state_indices()
    projection = [[None for _ in range(len(states))] for _ in range(len(dates))]
    for entry in data:
        projection[dates.index(entry[1])][states.index(entry[0])] = entry[2]
    return projection


def ihme_data():
    """ Get the IHME ventilator projection.

    Returns:
      For each day, a list of ventilator demand for the states,
      indexed in alphabetical order.
    """
    return parse_data(ihme_clean_data())


def ode_data():
    """ Get the ODE ventilator projection.

    Returns:
      For each day, a list of ventilator demand for the states,
      indexed in alphabetical order.
    """
    return parse_data(ode_clean_data())


def initial_supply():
    """ Get the initial supply of ventilators.

    Returns:
      A list of the initial supply of ventilators for each state,
      indexed in alphabetical order.
    """
    with open('data/bertsimas/initial_supply.csv', 'r') as f:
        data = [[x[0].strip('"\n'), int(x[5].strip('\n'))] for x in
                [y.split(',') for y in f.readlines()[1:]]]
    f.close()
    data.sort(key=lambda entry: entry[0])
    return [x[1] for x in data]







def greedyIter(n, tStart, tEnd, V, Q, D, N_max, R, xpast = None, Npast= None, zpast= None, fpast= None, phi_past= None):
    """
    Performs an iteration of greedy algorithm between times tStart and tEnd. All input data
    and past data have to be given starting from second iteration.
    """
    TT = [tt for tt in range(tStart, tEnd)]
    Qcurr = {tt:[Q[i][tt] for i in range(50)] for tt in range(tStart, tEnd)}
    # Past data initialization
    if xpast is None:
        xpast = {(i,t):0 for (i,t) in itertools.product(range(n),range(0,tStart))}
    if Npast is None:
        Npast = {t:0 for t in range(0,tStart)}
    if zpast is None:
        zpast = {(i,t):0 for (i,t) in itertools.product(range(n),range(0,tStart))}
    if fpast is None:
        fpast = {(i,j,t):0 for (i,j,t) in itertools.product(range(n),range(n),range(0,tStart))}
    if phi_past is None:
        phi_past = {i:0 for i in range(n)}
    ###
    Mgreed = gp.Model()
    # Variables
    x = Mgreed.addVars(n, TT, lb = 0, vtype = gp.GRB.INTEGER, name = "x" )
    N = Mgreed.addVars(TT, lb = 0.9*N_max[tStart:tEnd], ub =  N_max[tStart:tEnd], vtype = gp.GRB.INTEGER, name = "N")
    z = Mgreed.addVars (n, TT, lb = 0,  vtype=gp.GRB.INTEGER, name='z')
    f = Mgreed.addVars(n, n, range(tStart, tEnd), obj = 0.001, lb=0, vtype=gp.GRB.INTEGER, name='f')
    #
    tau = Mgreed.addVars(n, TT, lb=-gp.GRB.INFINITY, vtype=gp.GRB.INTEGER, name='tau')
    tauPos = Mgreed.addVars(n, TT, lb=0, vtype=gp.GRB.INTEGER, name='tauPos')
    ###
    # Constraints
    for t in TT:
        # No transfer to self
        for i in range(n):
            f[i, i, t].lb = 0
            f[i, i, t].ub = 0
        # No transfer back and forth
        for i in range(n):
            for j in range(i):
                Mgreed.addSOS(gp.GRB.SOS_TYPE1, [ f[i,j,t], f[j,i,t] ])
        # Total number of ventilators everywhere should be lesser than number injected
        Mgreed.addConstr(gp.quicksum(x[i, t] for i in range(n)) <= sum(V) + 
                    sum(Npast[tp] for tp in range(tStart)) +
                    gp.quicksum(N[tp] for tp in range(tStart,t+1)))
        #
        # New ventilators should necessarily be allocated. 
        Mgreed.addConstr(gp.quicksum(z[i, t] for i in range(n)) == N[t])
        # Ventilator shorages
        for i in range(n):
            Mgreed.addConstr(tau[i, t] == Qcurr[t][i]-x[i, t])
            Mgreed.addGenConstrMax(tauPos[i,t],[tau[i,t]],0)
    inCome = {}
    for i in range(n):
        for t in TT:
            inCome[(i,t)] = 0
            xeqn = 0
            # Handling transfer in the past
            for j in range(n):
                if t-D[j][i] >= tStart:
                    inCome[(i,t)] += f[j, i, t-D[j][i]]
                elif t-D[j][i] >= 0:
                    inCome[(i,t)] += fpast[j, i, t-D[j][i]]

            # How to handle x[i, t-1]
            if t == tStart and  t > 0:
                xeqn += xpast[i, t-1]
            elif t > tStart:
                xeqn += x[i, t-1]

            Mgreed.addConstr(x[i, t] == xeqn +
                        z[i, t] +
                        inCome[(i,t)] -
                        gp.quicksum(f[i, j, t] for j in range(n)))
    # Cant remove before transfer delay
    for i in range(n):
        for t in TT:
            eqn = 0
            for tp in range(max(1,t-R), t):
                if tp in TT:
                    eqn +=  (z[i,tp] + inCome[(i, tp)])
                else:
                    eqn += zpast[i,tp]  + sum([fpast[j, i, tp-D[j][i]] for j in range(n) if tp - D[j][i] >= 0 ])
                        
            Mgreed.addConstr(x[i, t] >= eqn )
    ###
    # Objective
#     # Minimizing final disparity
#     phiMax = Mgreed.addVar(lb=0, vtype=gp.GRB.INTEGER, name='phiMax')
#     phiMin = Mgreed.addVar(lb=0,vtype=gp.GRB.INTEGER, name="phiMin")
#     phi_ = Mgreed.addVars(n, lb=0, vtype=gp.GRB.INTEGER, name='phi_')
#     for i in range(n):
#         Mgreed.addConstr(phi_[i] == phi_past[i] + gp.quicksum(tauPos[i, t] for t in TT))
#     Mgreed.addGenConstrMin(phiMin, phi_)
#     Mgreed.addGenConstrMax(phiMax, phi_)
#     Mgreed.setObjective(phiMax-phiMin)
    # Minimizing day wise disparity
    phi_ = Mgreed.addVars(n, TT, lb=0, vtype=gp.GRB.INTEGER, name='phi_')
    phiMax = Mgreed.addVars(TT,lb=0, obj=1, vtype=gp.GRB.INTEGER, name='phiMax')
    phiMin = Mgreed.addVars(TT,lb=0, obj=-1, vtype=gp.GRB.INTEGER, name="phiMin")
    for tt in TT:
        for i in range(n):
            Mgreed.addConstr(phi_[i, tt] == phi_past[i] + gp.quicksum(tauPos[i, t] for t in TT if t <= tt))
    
    for tt in TT:
        Mgreed.addGenConstrMin(phiMin[tt], [phi_[i,tt] for i in range(n)])
        Mgreed.addGenConstrMax(phiMax[tt], [phi_[i,tt] for i in range(n)])
    Mgreed.update()
    ###
    Mgreed.params.Threads = 2
    Mgreed.optimize()
    ###
    return Mgreed, x,N,z, f, phi_
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
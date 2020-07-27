import gurobipy as gp

M = gp.Model()
x = M.addVars(range(2), lb=0, vtype=gp.GRB.INTEGER, name='x')
tau = M.addVars(range(2),lb=0,vtype=gp.GRB.CONTINUOUS, name='tau')
M.addConstr(gp.quicksum(x) >= 1)
M.addConstr(tau[0] == x[0]+x[1]/2)
M.addConstr(tau[1] == x[1]+x[0]/2)

M.setObjective(tau[0]-tau[1])
M.update()
M.optimize()

sol = M.getVars()
for var in sol:
    print(var.varname, var.x)

from gurobipy import *
import numpy as np

nScen = 15
nD = 10

nodes = ["n"+str(i+1) for i in range(nD)]
scenario = ["xi"+str(i+1) for i in range(nScen)]

demands = dict()
cf = dict()
cv = dict()
transp = dict()
for nn in nodes:
        cf[nn] = np.random.randint(11, 20)*3
        cv[nn] = np.random.randint(1,5)*2
        for xi in scenario:
                demands[nn, xi] = np.random.randint(1, 21)
        for nn2 in nodes:
                transp[nn, nn2] = np.random.randint(0,3)

prob = dict()
for xi in scenario:
        prob[xi] = 1.0/nScen

M = Model()
x = M.addVars(nodes, name = "x", ub = 100, obj = cv)
u = M.addVars(nodes, name = "u", vtype =GRB.BINARY, obj = cf)

M.addConstrs((x[nn] <= 100*u[nn] for nn in nodes))

q = M.addVars(nodes, nodes, scenario, name = "q")

M.addConstrs((quicksum(q[nn, nn2, xi] for nn2 in nodes) <= x[nn] for nn in
nodes for xi in scenario))

M.addConstrs((quicksum(q[nn, nn2, xi] for nn in nodes) >= demands[nn2, xi]
for nn2 in nodes for xi in scenario))

M.update()
oo = M.getObjective()

expr = quicksum(prob[xi]*q[nn, nn2, xi]*transp[nn, nn2] for nn in nodes for
nn2 in nodes for xi in scenario)

M.setObjective(expr+oo)
M.write('sriram.mps')
M.optimize()

import gurobipy as gp
import networkx as nx
import numpy as np
from fairness import *

nV = 10
Adj = np.zeros((nV,nV), dtype=int)

# Adjacency matrix for the example in manuscript
for i in range(nV-1):
    Adj[i,i+1] = 1
    Adj[i+1,i] = 1
Adj[0,nV-1] = 1
Adj[nV-1,0] = 1
Adj[2,7]=1
Adj[7,2]=1
Adj[2,9]=1
Adj[9,2]=1

# Random adjacency matrix 
# Adj = np.array([[0,0,0,1],[0,0,0,1],[0,0,0,1],[1,1,1,0]],dtype=int)


f=0.25*Adj
drawGraph(Adj)

# Max-min metric of fairness
M = getGurobiBasicModel_Diff(Adj,f)
M.optimize()
M2 = getMinimalGurobi(M, Adj)
M2.optimize()


# L1 metric of fairness
N = getGurobiBasicModel_L1(Adj,f)
N.optimize()
N2 = getMinimalGurobi(N, Adj)
N2.optimize()


# solutions
print("With Diff Metric")
var = M2.getVars()
for x in var:
    if("x" in x.Varname):
        print(x.VarName, "-----", round(x.X,2))
print("Total rounds: ",round(M2.objval))
TotalWelfare = sum([M2.getVarByName("tau["+str(v)+"]").X for v in range(nV)])
print("Total Welfare in fair solution: ", TotalWelfare)

print("With L1 Metric")
var = N2.getVars()
for x in var:
    if("x" in x.Varname):
        print(x.VarName, "-----", round(x.X,2))
print("Total rounds: ",round(N2.objval))
TotalWelfare = sum([N2.getVarByName("tau["+str(v)+"]").X for v in range(nV)])
print("Total Welfare in fair solution: ", TotalWelfare)

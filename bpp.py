import gurobipy as gp


w = [497, 497, 495, 485, 480, 478, 474, 473, 472, 470, 466, 450, 446, 445, 445, 444, 439, 434, 430, 420, 419, 414, 412, 410, 407, 405, 400, 397, 395, 376, 372, 370, 366, 366, 366, 366, 366, 363, 363, 362, 361, 357, 357, 356, 356, 355, 352, 351, 350, 350, 350, 347, 336, 333, 329, 325, 320, 315, 314, 313, 307, 303, 302, 301, 299, 298, 298, 298, 295, 294, 292, 290, 288, 287, 283, 282, 282, 276, 275, 275, 274, 273, 273, 272, 272, 271, 271, 269, 269, 268, 267, 267, 266, 263, 263, 262, 262, 261, 260, 259, 259, 259, 258, 256, 255, 254, 254, 254, 253, 253, 253, 253, 252, 252, 252, 252, 251, 251, 250, 250]
c = 1000
m = len(w) # Number of items
n = 50 # Max number of bins

M = gp.Model()
y = M.addVars(n, vtype=gp.GRB.BINARY, obj=1) # y[j] = 1 if bin j is used, 0 otherwise
x = M.addVars(m, n, vtype=gp.GRB.BINARY) # x[i][j] = 1 if item i is packed into bin j, 0 otherwise

# Each item is packed once
for i in range(m):
    M.addConstr(gp.quicksum(x[i, j] for j in range(n)) == 1)

# If a bin contains an item, that bin is used
for i in range(m):
    for j in range(n):
        M.addConstr(x[i, j] <= y[j])

# Bin capacity constraints
for j in range(n):
    M.addConstr(gp.quicksum(x[i, j]*w[i] for i in range(m)) <= c*y[j])

M.write('bpp.lp')

# M.update()
# M.optimize()

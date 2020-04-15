import copy
import networkx as nx
import random


#   B
#   |
# A-D-C
nodes = [_ for _ in range(4)]
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from([(0, 3),
                  (1, 3),
                  (2, 3)])

# Display the graph
# import matplotlib.pyplot as plt
# nx.draw_networkx(G)
# plt.show()

X = {key:[] for key in nodes}
p = 1/3
r = 1
allocations = [3, 1, 2]


def pre_allocate(node):
    """Checks the local profit distribution of a potential resource allocation.

    Args:
      node: The node to be allocated a resource.

    Returns:
      The profit distribution of this allocation.
    """
    alloc = [0 for _ in nodes]
    alloc[node] += r
    for i in G.neighbors(node):
        alloc[i] += p
    return alloc


def allocate(node, some_X):
    """Allocates the resource to a node and updates some_X.

    Args:
      node: The node to which the resource is allocated.
    """
    alloc = pre_allocate(node)
    for i in range(len(nodes)):
        some_X[i].append(alloc[i])    


def fairness(some_X):
    """Computes the fairness of some resource allocation distribution.

    Args:
      some_X: A global resource distribution.

    Returns:
      The fairness of this resource distribution.
    """
    total = 0
    for _, v in some_X.items():
        total += sum(v)
    mean = total/len(nodes)
    balance = 0
    for _, v in some_X.items():
        balance += abs(mean - sum(v))
    return balance
    # val_min = 999999999
    # val_max = 0
    # for _, v in some_X.items():
    #     val_min = min(val_min, sum(v))
    #     val_max = max(val_max, sum(v))
    # return val_max-val_min


def greedy():
    """Picks the best node to be greedy allocated.

    If multiple allocations are equivalent, chooses randomly.

    Returns:
      A greedy local allocation.
    """
    best_fairness = 999999
    best_node = []
    for i in range(len(nodes)):
        some_X = copy.deepcopy(X)
        allocate(i, some_X)
        bal = fairness(some_X)
        if bal <= best_fairness+0.00001 and bal >= best_fairness-0.00001:
            best_fairness = bal
            best_node.append(i)
        elif bal < best_fairness:
            best_fairness = bal
            best_node = [i]
    return random.choice(best_node)



template = '{0:<3} {1:<3} {2:<30} {3:<30} {4:<10}'
print(template.format('t', 'a', 'x(t)', 'X_0,t', 'f(X_0,t)'))
time = 0
while True:
    best_node = greedy()
    line = [time]
    line.append(best_node)
    line.append(','.join([str(round(j, 2)) for j in pre_allocate(best_node)]))
    allocate(best_node, X)
    line.append(','.join([str(round(sum(v), 2)) for _, v in X.items()]))
    line.append(fairness(X))
    print(template.format(*line))
    time += 1
    if time == 50:
        break

# allocations = [3, 1, 2, 0, 1, 2, 0, 1, 2, 0]
# template = '{0:<3} {1:<3} {2:<30} {3:<30} {4:<10}'
# print(template.format('t', 'a', 'x(t)', 'X_0,t', 'f(X_0,t)'))
# for i in range(len(allocations)):
#     line = [i]
#     line.append(allocations[i])
#     line.append(','.join([str(round(j, 2)) for j in pre_allocate(allocations[i])]))
#     allocate(allocations[i], X)
#     line.append(','.join([str(round(sum(v), 2)) for _, v in X.items()]))
#     line.append(fairness(X))
#     print(template.format(*line))




"""
Hi,

If we consider that T goes to infinity, here is what happens.

__For the greedy approach__

As the first choice, the greedy allocation inevitably starts with D, for a fairness value of 1.

As the second choice, all allocations are equivalent as they all result in a fairness of 2. Now, only two paths can be taken from here:
Path 1: From the second choice onwards, a continuous cycle of choices of A, B, C (in any order), which will see the fairness value oscillate between 1, 2, 1.67, 1, 2, 1.67, ...
Path 2: A second choice of D, after which will follow the same continuous cycle as the previous option, but now the fairness value oscillates between 2, 2, 2.67, 2, 2, 2.67, ... Considering that we're not allowed to look into the future to see where these paths lead, we're bound to pick this option at some point (the option of picking Path 2 is available any time we're on Path 1 and our fairness value is 1. But once we're on Path 2, we can't go back to Path 1.).

__For exact optimization__

The continuous cycle of choices is A, B, C (in any order), and the fairness value oscillates between 1.33, 1.33, 0, 1.33, 1.33, 0, ...

__Limits__

For the objective of max f(X_{0,T}), we thus have liminf=2 and limsup=2.67 for the greedy approach, and liminf=0 and limsup=1.33 for the exact approach.

About what you said about F(T) = \sum _{t=1}^T f(X_{0,t}), I'm not sure I understand what you mean. As F(T) represents the area under the curve (which we want to minimize), it will keep increasing for both methods as T approaches infinity since there is always some slight unfairness at each time t.

Philippe

"""

"""
I'm not sure I understand what you mean by F(T) = \sum _{t=1}^T f(X_{0,t}). My understanding is that we are looking at the fairness of cumulative allocations, i.e., f(X_{0, T}), but not at the cumulation of fairness values themselves over time.
"""

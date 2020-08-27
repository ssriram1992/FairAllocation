import networkx
import numpy
from pyscipopt import Model, Pricer, SCIP_RESULT, SCIP_PARAMSETTING, quicksum
from utrecht_scip import *
import math

EPS = 1e-6

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes (900 seconds)".
utrecht = Utrecht()
utrecht.reduce_graph(900)

# Bases are the only vertices that can be allocated ambulances.
bases = utrecht.B
n = len(utrecht.V)

# Replace the seconds by a boolean representing if the zones are close enough
# to cover each other or not.
adj = utrecht.get_binary_adjacency()

# Sufficient coverage required for various population densities.
sufficient = utrecht.get_sufficient_coverage()

# Model parameters
num_ambulances = 20
num_rounds = 3
max_transition = 1 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = 4 # Maximum practical number of ambulances in a zone (doesn't seem to make much of a difference). This is 4 because a base with 4 ambulances can sufficiently cover any zone no matter what its population density


class AmbulancePricer(Pricer):
    def __init__(self):
        self.data = {}
    
    def pricerredcost(self):
        # Retrieving the dual solutions
        dual_solutions = []
        for i, c in enumerate(self.data['phi_cons']):
            dual_solutions.append(self.model.getDualsolLinear(c))
        print('==========>>> Dual list:', [i for i in dual_solutions if i > EPS or i < -EPS])
        time = self.data['time_cons']
        mu = self.model.getDualsolLinear(time)
        print('==========>>> mu:', mu)
        
        # Pricing problem.
        pp = Model("pp")

        # Turning off presolve and disable verbosity.
        pp.setPresolve(SCIP_PARAMSETTING.OFF)
        pp.hideOutput()

        # Variables.
        c_vars = [pp.addVar(vtype='B', obj=-dual_solutions[i]) for i in range(n)] # DEBUG: obj
        u_vars = [pp.addVar(vtype='I', lb=0, ub=max_practical_ambulances) for _ in range(len(bases))] # ub is 4 since that's the best coverage that the base can offer in any case

        # # Forbid the generation of existing columns.
        # binary_bases = [[pp.addVar(vtype='B') for _ in range(num_ambulances+1)] for _ in range(len(bases))]
        # for i in range(len(bases)):
        #     for j in range(num_ambulances+1):
        #         pp.addConsIndicator(u_vars[i] <= j, binary_bases[i][j])
        #         pp.addConsIndicator(u_vars[i] >= j, binary_bases[i][j])
        # for i in self.data['allocations']:
        #     pp.addCons(quicksum(binary_bases[j][i[j]] for j in range(len(i))) <= sum(i) -1)
        
        # Big M constraints (link c and u).
        bigM = num_ambulances*2
        for i in range(n):
            pp.addCons(quicksum(adj[i, bases[j]]*u_vars[j] for j in range(len(bases))) - sufficient[i] + bigM*(1-c_vars[i]) >= 0)
            pp.addCons(quicksum(adj[i, bases[j]]*u_vars[j] for j in range(len(bases))) - sufficient[i] + 1 -bigM*c_vars[i] <= 0)

        # Constraint: Don't use more than the available number of ambulances.
        pp.addCons(quicksum(u_vars) <= num_ambulances)

        # Constraint: The configuration must sufficiently cover 95% of the zones.
        pp.addCons(quicksum(c_vars) >= math.ceil(min_coverage*n))
        
        pp.optimize()

        objval = pp.getObjVal() # DEBUG: obj
        print('==========>>> PP objective value:', objval)

        # Adding the column to the master problem
        if mu-objval < -EPS: # DEBUG: obj
            print('==========>>> Generating new column')
            #currentNumVar = len(self.data['var'])

            # Creating new var; must set pricedVar to True
            # A new variable is created for the new pattern
            newVar = self.model.addVar("NewPattern_", vtype = "I", pricedVar = True) # pricedVar=True says that this variable was generated by the PP

            # Adding the new variable to the constraints of the master problem
            # The new pattern is created, and the demand constraints are updated
            newPattern = []
            for i, c in enumerate(self.data['phi_cons']):
                # print(len(self.data['phi_cons']))
                coeff = round(pp.getVal(c_vars[i]))
                self.model.addConsCoeff(c, newVar, coeff)
                newPattern.append(coeff)
            
            # Also update the time horizon constraint
            self.model.addConsCoeff(self.data['time_cons'], newVar, 1)

            newAllocation = []
            for i in range(len(bases)):
                coeff = round(pp.getVal(u_vars[i]))
                newAllocation.append(coeff)
            print('==========>>>', newAllocation)
                
            # Storing the new variable in the pricer data.
            self.data['coverages'].append(newPattern)
            self.data['allocations'].append(newAllocation)
            self.data['x_vars'].append(newVar)

            print('==========>>> MP objective value:', self.model.getObjVal())
            
            # print([hash(tuple(i)) for i in self.data['allocations']])

        return {'result':SCIP_RESULT.SUCCESS}

    # The initialisation function for the variable pricer to retrieve the transformed constraints of the problem
    # This function seems to be called only once
    def pricerinit(self):
        for i, c in enumerate(self.data['phi_cons']):
            self.data['phi_cons'][i] = self.model.getTransformedCons(c) # TODO: what is a transformed constraint?
        self.data['time_cons'] = self.model.getTransformedCons(self.data['time_cons'])


def master_problem():
    # Find one feasible configuration for the initial column.
    def initial_feasible_configuration():
        """Returns a feasible configuration."""
        m = Model('M')
        m.setPresolve(SCIP_PARAMSETTING.OFF)
        m.hideOutput()
        x = [m.addVar(lb=0, ub=max_practical_ambulances, vtype='I', name='x'+str(i)) for i in range(len(bases))] # again, ub=4 because that's the max coverage the base can offer
        
        # Constraint: Don't use more than the available number of ambulances.
        m.addCons(quicksum(x) <= num_ambulances)

        # Constraint: The configuration must sufficiently cover 95% of the zones.
        c = [m.addVar(vtype='B', name='c'+str(i)) for i in range(n)]
        for i in range(n):
            m.addConsIndicator(quicksum(x[j] * adj[bases[j], i] for j in range(len(bases))) >= sufficient[i], c[i])
        m.addCons(quicksum(c) >= math.ceil(min_coverage*n))
        m.optimize()

        allocation = [int(m.getVal(x[i])) for i in range(len(x))]        
        coverage = utrecht.allocation_coverage(allocation)
                
        return allocation, coverage
    
    # Master problem.
    mp = Model('MP')
    mp.setPresolve(SCIP_PARAMSETTING.OFF)

    # Initial feasible column.
    allocation, coverage = initial_feasible_configuration()

    # Keeping track of allocations and coverages.
    allocations = [allocation]
    coverages = [coverage]
    
    # DEBUG: Optimal solution for 20 ambulances, 95% coverage and 3 rounds (solution = 2)
    # allocations = [[0, 1, 0, 2, 1, 2, 0, 0, 4, 0, 2, 4, 0, 2, 2, 0, 0, 0],
    #                [2, 3, 0, 1, 0, 4, 0, 0, 1, 1, 3, 0, 0, 1, 1, 3, 0, 0],
    #                [2, 1, 1, 1, 0, 3, 0, 0, 2, 0, 0, 4, 2, 0, 3, 0, 1, 0]]
    # coverages = [utrecht.allocation_coverage(x) for x in allocations]
    
    # Objective.
    phi = mp.addVar(name='phi', vtype='I')
    mp.setObjective(phi, sense='maximize')

    x_vars = [mp.addVar(vtype='I') for _ in allocations]
    
    # The n constraints constraining the value of phi.
    phi_cons = []
    for i in range(n):
        phi_cons.append(mp.addCons(quicksum(coverages[j][i]*x_vars[j] for j in range(len(allocations))) - phi >= 0, separate=False, modifiable=True))

    # Time horizon constraint: We must use T configurations.
    time_cons = mp.addCons(quicksum(x_vars) == num_rounds, separate=False, modifiable=True)
    
    # Pricer.
    ap = AmbulancePricer()
    mp.includePricer(ap, 'AP', '')
    ap.data['x_vars'] = x_vars
    ap.data['phi_cons'] = phi_cons
    ap.data['time_cons'] = time_cons
    ap.data['allocations'] = allocations
    ap.data['coverages'] = coverages
    
    # Branching.
    # TODO

    mp.optimize()

    print('=====================================')
    print(f'Objective: {mp.getVal(phi)}')
    x = []
    for i in range(len(x_vars)):
        x_val = mp.getVal(x_vars[i])
        if x_val - EPS > 0:
            x.append((i, x_val))
            
    print(f'x: {x}')

master_problem()



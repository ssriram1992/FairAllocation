import networkx
import numpy
from pyscipopt import *#Model, Pricer, SCIP_RESULT, SCIP_PARAMSETTING, quicksum, Branchrule
from utrecht_scip import *
import math
from ortools.sat.python import cp_model

EPS = 1e-6

randomize = False
seed = 0

# For efficiency we use the performance target for the Netherlands, i.e.,
# that "95% of all calls should be reached within 15 minutes (900 seconds)".
utrecht = Utrecht(randomize, seed)
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
num_rounds = 30
max_transition = 1 # 0 = no transition allowed, 1 = unlimited transitions
min_coverage = 0.95 # 0 = no one needs to be covered, 1 = everyone has to be covered
max_practical_ambulances = 4 # Maximum practical number of ambulances in a zone (doesn't seem to make much of a difference). This is 4 because a base with 4 ambulances can sufficiently cover any zone no matter what its population density

try:
    from types import SimpleNamespace
except:
    class SimpleNamespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

        def __repr__(self):
            keys = sorted(self.__dict__)
            items = ("{}={!r}".format(k, self.__dict__[k]) for k in keys)
            return "{}({})".format(type(self).__name__, ", ".join(items))

        def __eq__(self, other):
            return self.__dict__ == other.__dict__


def has_hamiltonian_path(V, E, col_cov, num_times):
    """Checks if a Hamiltonian path exists.
    
    Args:
      V: List of vertices (columns).
      E: List of edges (if a transition between two columns is allowed).
      col_cov: Matrix of the zone coverages of the columns (c[i][j] == 1 if
               zone i is covered by column j).
      num_times: Number of times the variables are found in the solution.

    Returns:
      Boolean indicating if a Hamiltonian path exists
    """
    num_cols = len(V)
    num_zones = len(col_cov)
    
    # First, check if the graph is disconnected (in which case no
    # Hamiltonian path exists).
    G = networkx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)
    # If the graph is not connected, no Hamiltonian path can exist.
    if not networkx.is_connected(G):
        return False

    # Variables.
    model = cp_model.CpModel()
    x = [model.NewIntVar(0, num_cols-1, 'x'+str(i)) for i in range(num_rounds)]

    # Alternative for GCC, since the constraint is not available in OR-Tools.
    x_occs = []
    #print('num_times', num_times)
    for i in range(num_cols):
        occs = []
        for j in range(num_rounds):
            boolvar = model.NewBoolVar('')
            model.Add(x[j] == i).OnlyEnforceIf(boolvar)
            model.Add(x[j] != i).OnlyEnforceIf(boolvar.Not())
            occs.append(boolvar)
        x_occs.append(sum(occs))
        # Every var must be used the correct number of times
        model.Add(x_occs[i] == int(num_times[i]))

    # Regular constraint (Hamiltonian path).
    # For the initial state, we use a dummy node which is connected to
    # all other nodes.
    dummy = max(V)+1
    start = dummy
    end = V
    arcs = [(dummy, i, i) for i in V]
    for e in E:
        arcs.append((e[0], e[1], e[1]))
    # If there is only one vertex then a Hamiltonian path exists.
    if len(V) > 1:
        model.AddAutomaton(x, start, end, arcs)

    # Solve the model.
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    assert status == cp_model.OPTIMAL or status == cp_model.INFEASIBLE

    if status == cp_model.OPTIMAL:
        return True
    elif status == cp_model.INFEASIBLE:
        print(V)
        print(E)
        print(num_times)
        exit()
        return False

        

class AmbulanceConshdlr(Conshdlr):
    def __init__(self, model, vars, alloc, cov):
        self.model = model
        self.vars = vars
        self.alloc = alloc
        self.cov = cov
    
    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        allocations = self.alloc
        coverages = self.cov
        '''calls enforcing method of constraint handler for LP solution for all constraints added'''
        # We assume that any LP solution is feasible, since infeasibility is determined only once all variables are integer
        print("===> Entering consenfolp()")
        s = self.model
        vars = s.getVars(transformed=True) # TRANSFORMED=TRUE??
        # Check that all the vars are fixed; if one of them is not return FEASIBLE
        at_least_one_fixed = False # there is one integer var
        at_least_one_unfixed = False # there is one continuous var, whose bounds are not 0 to infinity
        at_least_one_infty = False # there is one var whose bounds are 0 to infinity
        for v in vars:
            if v.getUbLocal()-v.getLbLocal() > 0.5 and v.getUbLocal() < 100000:
                at_least_one_unfixed = True
            elif v.getUbLocal()-v.getLbLocal() <= 0.5:
                at_least_one_fixed = True
            else:
                at_least_one_infty = True

        if at_least_one_unfixed == True or at_least_one_fixed == False: # Solution has at least one unfixed var or all infty vars
            print(f'LBS {[v.getLbLocal() for v in vars]}')
            print(f'UBS {[v.getUbLocal() for v in vars]}')
            #input('Enter to continue')
            return {"result": SCIP_RESULT.FEASIBLE}

        # Now we know all vars are fixed, let's see if there is a hamiltonian path
        values = [v.getLbLocal() for v in vars]
        V = [i for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        values = [values[i] for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        V0 = [i for i in range(len(V))]
        E = []
        for i in range(len(V0)-1):
            for j in range(i+1, len(V0)):
                a1 = allocations[V[V0[i]]]
                a2 = allocations[V[V0[j]]]
                # Maximum number of ambulances that can transition between two configurations.
                max_amb = math.floor(num_ambulances*max_transition)
                a_diff = [abs(a1[x]-a2[x]) for x in range(len(a1))]
                if sum(a_diff)/2 < max_amb:
                    E.append((V0[i], V0[j]))
                    E.append((V0[j], V0[i]))
        # Get the coverages of the subset of columns which take nonzero values.
        c = [coverages[i] for i in range(len(coverages)) if i in V]
        c = [list(x) for x in zip(*c)]
        has_path = has_hamiltonian_path(V0, E, c, values)
        print(f'=======================*****************************>>>>>>>>>>>> {has_path}')
        if has_path:
            return {"result": SCIP_RESULT.FEASIBLE}
        else:
            return {"result": SCIP_RESULT.CUTOFF}

    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        allocations = self.alloc
        coverages = self.cov
        print("===> Entering consenfops()")
        s = self.model
        vars = s.getVars(transformed=True) # TRANSFORMED=TRUE??
        # Check that all the vars are fixed; if one of them is not return FEASIBLE
        # for v in vars:
        #     if v.getUbLocal()-v.getLbLocal() > 0.5 and v.getUbLocal() < 100000:
        #         return {"result": SCIP_RESULT.FEASIBLE}
        # Check that all the vars are fixed; if one of them is not return FEASIBLE
        at_least_one_fixed = False # there is one integer var
        at_least_one_unfixed = False # there is one continuous var, whose bounds are not 0 to infinity
        at_least_one_infty = False # there is one var whose bounds are 0 to infinity
        for v in vars:
            if v.getUbLocal()-v.getLbLocal() > 0.5 and v.getUbLocal() < 100000:
                at_least_one_unfixed = True
            elif v.getUbLocal()-v.getLbLocal() <= 0.5:
                at_least_one_fixed = True
            else:
                at_least_one_infty = True

        if at_least_one_unfixed == True or at_least_one_fixed == False: # Solution has at least one unfixed var or all infty vars
            print(f'LBS {[v.getLbLocal() for v in vars]}')
            print(f'UBS {[v.getUbLocal() for v in vars]}')
            #input('Enter to continue')
            return {"result": SCIP_RESULT.FEASIBLE}

        # Now we know all vars are fixed, let's see if there is a hamiltonian path
        values = [v.getLbLocal() for v in vars]
        V = [i for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        values = [values[i] for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        V0 = [i for i in range(len(V))]
        E = []
        for i in range(len(V0)-1):
            for j in range(i+1, len(V0)):
                a1 = allocations[V[V0[i]]]
                a2 = allocations[V[V0[j]]]
                # Maximum number of ambulances that can transition between two configurations.
                max_amb = math.floor(num_ambulances*max_transition)
                a_diff = [abs(a1[x]-a2[x]) for x in range(len(a1))]
                if sum(a_diff)/2 < max_amb:
                    E.append((V0[i], V0[j]))
                    E.append((V0[j], V0[i]))
        # Get the coverages of the subset of columns which take nonzero values.
        c = [coverages[i] for i in range(len(coverages)) if i in V]
        c = [list(x) for x in zip(*c)]
        has_path = has_hamiltonian_path(V0, E, c, values)
        print(f'=======================*****************************>>>>>>>>>>>> {has_path}')
        if has_path:
            return {"result": SCIP_RESULT.FEASIBLE}
        else:
            return {"result": SCIP_RESULT.CUTOFF}

    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason, completely):
        allocations = self.alloc
        coverages = self.cov
        print("===> Entering conscheck()")
        assert len(constraints) == 1
        s = self.model
        vars = s.getVars(transformed=True) # TRANSFORMED=TRUE??
        # Check that all the vars are fixed; if one of them is not return FEASIBLE
        # for v in vars:
        #     if v.getUbLocal()-v.getLbLocal() > 0.5 and v.getUbLocal() < 100000:
        #         return {"result": SCIP_RESULT.FEASIBLE}
        # Check that all the vars are fixed; if one of them is not return FEASIBLE
        at_least_one_fixed = False # there is one integer var
        at_least_one_unfixed = False # there is one continuous var, whose bounds are not 0 to infinity
        at_least_one_infty = False # there is one var whose bounds are 0 to infinity
        for v in vars:
            if v.getUbLocal()-v.getLbLocal() > 0.5 and v.getUbLocal() < 100000:
                at_least_one_unfixed = True
            elif v.getUbLocal()-v.getLbLocal() <= 0.5:
                at_least_one_fixed = True
            else:
                at_least_one_infty = True

        if at_least_one_unfixed == True or at_least_one_fixed == False: # Solution has at least one unfixed var or all infty vars
            print(f'LBS {[v.getLbLocal() for v in vars]}')
            print(f'UBS {[v.getUbLocal() for v in vars]}')
            #input('Enter to continue')
            return {"result": SCIP_RESULT.FEASIBLE}

        # # Now we know all vars are fixed, let's see if there is a hamiltonian path
        values = [v.getLbLocal() for v in vars]
        V = [i for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        values = [values[i] for i in range(len(values)) if values[i] > 0.1 and values[i] < 100000]
        V0 = [i for i in range(len(V))]
        E = []
        for i in range(len(V0)-1):
            for j in range(i+1, len(V0)):
                a1 = allocations[V[V0[i]]]
                a2 = allocations[V[V0[j]]]
                # Maximum number of ambulances that can transition between two configurations.
                max_amb = math.floor(num_ambulances*max_transition)
                a_diff = [abs(a1[x]-a2[x]) for x in range(len(a1))]
                if sum(a_diff)/2 < max_amb:
                    E.append((V0[i], V0[j]))
                    E.append((V0[j], V0[i]))
        # Get the coverages of the subset of columns which take nonzero values.
        c = [coverages[i] for i in range(len(coverages)) if i in V]
        c = [list(x) for x in zip(*c)]
        has_path = has_hamiltonian_path(V0, E, c, values)
        print(f'=======================*****************************>>>>>>>>>>>> {has_path}')
        if has_path:
            return {"result": SCIP_RESULT.FEASIBLE}
        else:
            return {"result": SCIP_RESULT.INFEASIBLE}

        

    def conspresol(self, constraints, nrounds, presoltiming,
                   nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                   nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, result_dict):
        return result_dict




class AmbulanceBranching(Branchrule):
    def __init__(self, model, variables):
        self.model = model
        self.variables = variables # We can add all the parameeters we need, such as the variables

    def branchexeclp(self, allowaddcons):
        print('===> Entering branchexeclp()')
        # print(f'Choices: {self.model.getLPBranchCands()}')
        # choice = self.model.getLPBranchCands()[0][0]
        # Random branching
        choice = random.choice(self.model.getLPBranchCands()[0])
        down, eq, up = self.model.branchVar(choice)
        # print('Branching on', choice)
        return {"result": SCIP_RESULT.BRANCHED}

    def branchexecps(self, alloaddcons):
        print('===> Entering branchexecps()')
        print(f'Choices: {[(i, i.getLPSol()) for i in self.model.getPseudoBranchCands()[0]]}')
        #choice = self.model.getPseudoBranchCands()[0][0]
        choice = random.choice(self.model.getPseudoBranchCands()[0])
        down, eq, up = self.model.branchVar(choice)
        print('Branching on', choice)
        return {"result": SCIP_RESULT.BRANCHED}
        


class AmbulancePricer(Pricer):
    def __init__(self):
        self.data = {}
    
    def pricerredcost(self):
        x = []
        for i in range(len(self.data['x_vars'])):
            x_val = self.data['mp'].getVal(self.data['x_vars'][i])
            if x_val - EPS > 0:
                x.append((i, x_val))
            
        print(f'MP sol: {x}')
        # Retrieving the dual solutions
        dual_solutions = []
        for i, c in enumerate(self.data['phi_cons']):
            dual_solutions.append(self.model.getDualsolLinear(c))
        # print('==========>>> Dual list:', [i for i in dual_solutions if i > EPS or i < -EPS])
        time = self.data['time_cons']
        mu = self.model.getDualsolLinear(time)
        # print('==========>>> mu:', mu)
        
        # Pricing problem.
        pp = Model("pp")

        # Turning off presolve and disable verbosity.
        pp.setPresolve(SCIP_PARAMSETTING.OFF)
        pp.hideOutput()

        # Variables.
        c_vars = [pp.addVar(vtype='B', obj=-dual_solutions[i]) for i in range(n)] # DEBUG: obj
        u_vars = [pp.addVar(vtype='I', lb=0, ub=max_practical_ambulances) for _ in range(len(bases))] # ub is 4 since that's the best coverage that the base can offer in any case
        
        # Big M constraints (link c and u).
        bigM = num_ambulances*2
        for i in range(n):
            pp.addCons(quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) - sufficient[i] + bigM*(1-c_vars[i]) >= 0)
            pp.addCons(quicksum(adj[bases[j], i]*u_vars[j] for j in range(len(bases))) - sufficient[i] + 1 -bigM*c_vars[i] <= 0)

        # Constraint: Don't use more than the available number of ambulances.
        pp.addCons(quicksum(u_vars) <= num_ambulances)

        # Constraint: The configuration must sufficiently cover 95% of the zones.
        pp.addCons(quicksum(c_vars) >= math.ceil(min_coverage*n))
        
        pp.optimize()

        objval = pp.getObjVal() # DEBUG: obj
        print('==========>>> PP objective value, mu:', objval, mu)

        # Adding the column to the master problem
        if abs(mu)-abs(objval) < -EPS: # DEBUG: obj
            print('==========>>> Generating new column')
            #currentNumVar = len(self.data['var'])

            # Creating new var; must set pricedVar to True
            # A new variable is created for the new pattern
            newVar = self.model.addVar(f"NewPattern_{len(self.data['allocations'])}", vtype = "I", pricedVar = True) # pricedVar=True says that this variable was generated by the PP

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
            # print('==========>>> allocation:', newAllocation)

            if newAllocation in self.data['allocations']:
                print('==========>>> Allocation already exists!')
            
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
    
    # Objective.
    # TODO: put back to integer
    phi = mp.addVar(name='phi', vtype='I') # M = implicit integer
    mp.setObjective(phi, sense='maximize')

    # TODO: put back as integer
    x_vars = [mp.addVar(vtype='I', name=f'NewPattern_{i}') for i in range(len(allocations))]
    
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
    ap.data['mp'] = mp
    ap.data['phi_cons'] = phi_cons
    ap.data['time_cons'] = time_cons
    ap.data['allocations'] = allocations
    ap.data['coverages'] = coverages
    
    # Branching.
    branchrule = AmbulanceBranching(mp, x_vars)
    mp.includeBranchrule(branchrule, '', '', priority=100000, maxdepth=65534, maxbounddist=1.0)

    conshdlr = AmbulanceConshdlr(mp, x_vars, allocations, coverages)
    mp.includeConshdlr(conshdlr, "", "", enfopriority = -10, chckpriority = -10)
    cons = mp.createCons(conshdlr, "", modifiable=True, local=True) # modifiable since new vars will be introduced
    #cons.data = SimpleNamespace()
    # cons.data['vars'] = x_vars
    # cons.data['model'] = mp
    # cons.data['allocations'] = allocations
    # cons.data['coverages'] = coverages
    mp.addPyCons(cons)
    
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



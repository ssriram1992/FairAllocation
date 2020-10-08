from pyscipopt import *#Model, Pricer, SCIP_RESULT, SCIP_PARAMSETTING, quicksum, Branchrule

"""
Getting everything to work:
1. Create a virtual environment in a folder
2. Download pyscipopt from github and extract it inside the virtual environment folder
3. In src/scip.pxd, add the line:

    int SCIPgetNPseudoBranchCands(SCIP* scip)

   before the line:

    SCIP_RETCODE SCIPgetPseudoBranchCands(SCIP* scip, SCIP_VAR*** pseudocands, int* npseudocands, int* npriopseudocands)
4. In src/scip.pyx, after the definition of getLPBranchCands(), add the function:

    def getPseudoBranchCands(self):
        cdef int ncands
        cdef int npseudocands
        cdef int npriopseudocands
        cdef SCIP_VAR** pseudocands
        PY_SCIP_CALL(SCIPgetPseudoBranchCands(self._scip, &pseudocands, &npseudocands, &npriopseudocands))
        return ([Variable.create(pseudocands[i]) for i in range(npseudocands)], npseudocands, npriopseudocands)

5. cd to the pyscipopt folder and:

pip install .

(this will do the cpython/cython/whatever routine to update the files)

6. Everything should work

"""

"""
OBJECTIVE:

This is a simple bin packing problem with no special constraints. For this simple example/test, we want to use Branchrule and Conshdlr to ensure that in the end, the value of each integer variable must end in zero (i.e. must be a multiple of 10)

ISSUES:
- In consprop(), why is the status always "unknown"? (Right now, I disabled consprop so this is irrelevant---see argument propfreq in the MP)
- Why is conscheck() not called when the status is "optimal"?
- We might have to change options on lines:

    s.includeConshdlr()

    and

    cons = s.createCons(conshdlr....)
"""

# Be prompted to choose the branching candidates during the solve process
manual = False

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


class MustEndInZeroConshdlr(Conshdlr):
##########################
# I put exit() in the functions which we have not implemented ourselves, so that we know to implement them if needed
###########################
    def constrans(self, sourceconstraint):
        '''sets method of constraint handler to transform constraint data into data belonging to the transformed problem '''
        print("****************** CONSTRANS")
        # Left as default for now
        return {}

    def consinitlp(self, constraints):
        '''calls LP initialization method of constraint handler to separate all initial active constraints '''
        print("****************** CONSINITLP")
        # Left as default for now
        return {}

    def conssepalp(self, constraints, nusefulconss):
        '''calls separator method of constraint handler to separate LP solution '''
        print("****************** CONSSEPALP")
        exit()
        return {}

    def conssepasol(self, constraints, nusefulconss, solution):
        '''calls separator method of constraint handler to separate given primal solution '''
        print("****************** CONSSEPASOL")
        exit()
        return {}

    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        '''calls enforcing method of constraint handler for LP solution for all constraints added'''
        print("****************** CONSENFOLP")
        # print('solution is:', [constraints[0].data.model.getVal(i) for i in constraints[0].data.vars])
        # We assume that any LP solution is feasible, since infeasibility is determined only once all variables are integer
        return {"result": SCIP_RESULT.FEASIBLE}

    def consenforelax(self, solution, constraints, nusefulconss, solinfeasible):
        '''calls enforcing method of constraint handler for a relaxation solution for all constraints added'''
        print("****************** CONSENFORELAX")
        exit()
        return {}

    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        '''calls enforcing method of constraint handler for pseudo solution for all constraints added'''
        print("****************** CONSENFOPS")
        exit()
        return {}

    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason, completely):
        '''calls feasibility check method of constraint handler '''
        # if the solution is not integer then it is valid by default, since we only check integer solutions
        # check for hamiltonian path
        print("****************** CONSCHECK")
        assert len(constraints) == 1
        status = constraints[0].data.model.getStatus()
        print(f'Solution status: {status}')
        assert status == 'unknown' # Why is the status always unknown??
        # if status == 'optimal':
        #     print('solution is:', [constraints[0].data.model.getVal(i) for i in constraints[0].data.vars])
        #     exit()
        #     #return {"result": SCIP_RESULT.INFEASIBLE}
        # else:
        #     # if the problem is not fully solved, this is infeasible
        #     return {"result": SCIP_RESULT.INFEASIBLE}
        # If we always return INFEASIBLE, then why does SCIP return an optimal solution for this problem?
        return {"result": SCIP_RESULT.INFEASIBLE}

        
    def consprop(self, constraints, nusefulconss, nmarkedconss, proptiming):
        '''calls propagation method of constraint handler '''
        print("****************** CONSPROP")
        exit() # I disabled propfreq for now
        status = constraints[0].data.model.getStatus()
        print(status)
        if status == 'optimal':
            print('solution is:', [constraints[0].data.model.getVal(i) for i in constraints[0].data.vars])
            exit()
        return {}#"result": SCIP_RESULT.DIDNOTRUN}

    def conspresol(self, constraints, nrounds, presoltiming,
                   nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                   nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, result_dict):
        '''calls presolving method of constraint handler '''
        print("****************** CONSPRESOL")
        # Left as default for now
        return result_dict

    def consresprop(self):
        '''sets propagation conflict resolving method of constraint handler '''
        print("****************** CONSRESPROP")
        exit()
        return {}

    def conslock(self, constraint, locktype, nlockspos, nlocksneg):
        '''variable rounding lock method of constraint handler'''
        print("****************** CONSLOCK")
        # Left as default for now
        return {}

    def consgetnvars(self, constraint):
        '''sets constraint variable number getter method of constraint handler '''
        print("****************** CONSGETNVARS")
        # Left as default for now
        return {}


# This is the branching object, and the user must define some of the functions below
class CutBranching(Branchrule):
    def __init__(self, model, variables):
        self.model = model
        self.variables = variables # We can add all the parameeters we need, such as the variables

    # Branching rule: Executes branching rule for fractional LP solution
    # This function must be defined by the user
    def branchexeclp(self, allowaddcons):
        print('********** branchexeclp')
        print(f'Choices: {self.model.getLPBranchCands()}')
        if manual:
            choice = input('Choose branching variable: ')
            choice = int(choice)
            choice2 = self.model.getLPBranchCands()[0][choice]
            down, eq, up = self.model.branchVar(choice2)
        else:
            choice2 = self.model.getLPBranchCands()[0][0]
            down, eq, up = self.model.branchVar(choice2)
        print('==========> Branching on', choice2)
        return {"result": SCIP_RESULT.BRANCHED}

    # Optional: Executes branching rule for external branching candidates
    def branchexecext(self, alloaddcons):
        print('********** branchexecext')
        exit()

    # Optional: Executes branching rule for not completely fixed pseudo solution
    def branchexecps(self, alloaddcons):
        print('********** branchexecps')
        print(f'Choices: {[(i, i.getLPSol()) for i in self.model.getPseudoBranchCands()[0]]}')
        if manual:
            choice = input('Choose branching variable: ')
            choice = int(choice)
            choice2 = self.model.getPseudoBranchCands()[0][choice]
            down, eq, up = self.model.branchVar(choice2)
        else:
            choice2 = self.model.getPseudoBranchCands()[0][0]
            down, eq, up = self.model.branchVar(choice2)
        print('==========> Branching on', choice2)
        return {"result": SCIP_RESULT.BRANCHED}
    

# This is the pricer object, and the user must define some of the functions below
class CutPricer(Pricer):
    def __init__(self):
        self.data = {}

    # Initialisation function for the variable pricer to retrieve the transformed constraints of the problem
    # Note: What is a transformed constraint?
    def pricerinit(self):
        for i, c in enumerate(self.data['cons']):
            self.data['cons'][i] = self.model.getTransformedCons(c)
        
    # The reduced cost function for the variable pricer
    def pricerredcost(self):
        # Retrieving the dual solutions
        dualSolutions = []
        for i, c in enumerate(self.data['cons']):
            dualSolutions.append(self.model.getDualsolLinear(c))

        # Building a MIP to solve the subproblem
        subMIP = Model("CuttingStock-Sub")
        subMIP.setPresolve(SCIP_PARAMSETTING.OFF)
        subMIP.hideOutput()

        # Variables for subMIP
        cutWidthVars = []
        varNames = []
        varBaseName = "CutWidth"
        for i in range(len(dualSolutions)):
            varNames.append(varBaseName + "_" + str(i))
            cutWidthVars.append(subMIP.addVar(varNames[i], vtype = "I", obj = -1.0 * dualSolutions[i]))

        # Adding the knapsack constraint (pricing problem constraints)
        knapsackCons = subMIP.addCons(
            quicksum(w*v for (w,v) in zip(self.data['widths'], cutWidthVars)) <= self.data['rollLength'])

        # Solving the subMIP to generate the most negative reduced cost pattern
        subMIP.optimize()

        # Adding a new column to the master problem
        objval = 1 + subMIP.getObjVal()
        if objval < -1e-08:
            currentNumVar = len(self.data['var'])

            # A new variable is created for the new pattern: the pricedVar flag must be set to True (this tells the master problem that this variable was generated by the pricing problem)
            newVar = self.model.addVar("NewPattern_" + str(currentNumVar), vtype="I", obj=1.0, pricedVar=True)

            # Add the new variable to the constraints of the master problem
            # The new pattern is created, and the demand constraints of the master problem are updated here, to take into account this new variable
            newPattern = []
            for i, c in enumerate(self.data['cons']):
                coeff = round(subMIP.getVal(cutWidthVars[i]))
                self.model.addConsCoeff(c, newVar, coeff)
                newPattern.append(coeff)

            # Storing the new variable in the pricer data.
            self.data['patterns'].append(newPattern)
            self.data['var'].append(newVar)
            print(f'==========> Generated a column: NewPattern_{currentNumVar}')

        return {'result':SCIP_RESULT.SUCCESS}


def cuttingstock():
    # Master problem
    s = Model("CuttingStock")
    s.setPresolve(0)

    # Pricing problem
    pricer = CutPricer()
    s.includePricer(pricer, "CuttingStockPricer", "Pricer to identify new cutting stock patterns")

    # Problem data
    widths = [14, 31, 36, 45]
    demand = [211, 395, 610, 97]
    rollLength = 100

    # MP variables
    cutPatternVars = [] # List of int vars: Number of times variable i (pattern i) is used
    varNames = [] # List of strings: Names of the variables
    varBaseName = "Pattern"
    patterns = [] # List of patterns (a pattern is a list of integers: pattern[i] is how many times weight i is in that pattern)

    # Initial columns: Variables
    for i in range(len(widths)):
        varNames.append(varBaseName + "_" + str(i))
        cutPatternVars.append(s.addVar(varNames[i], obj = 1.0, vtype = "I"))

    # Demand constraints
    demandCons = [] # List of constraints: The demand must be satisfied (these constraints are updated inside of the pricing problem, when a new pattern/variable is generated)
    for i in range(len(widths)):
        numWidthsPerRoll = float(int(rollLength/widths[i]))
        
        # The modifiable flag must be set to True (this tells the msater problem that the pricing problem will automatically modify these constraints to take into account newly generated columns)
        # The separate=False flag is because: "In most cases you should deactivate separators since cutting planes that are added to your master problem may destroy your pricing problem", see https://www.scipopt.org/doc/html/FAQ.php
        demandCons.append(s.addCons(numWidthsPerRoll*cutPatternVars[i] >= demand[i], separate = False, modifiable = True))
        # Initial columns: Patterns
        newPattern = [0]*len(widths)
        newPattern[i] = numWidthsPerRoll
        patterns.append(newPattern)

    # Setting the pricer_data for use in the init and redcost functions
    pricer.data['var'] = cutPatternVars # List of int vars: Number of tiles variable/pattern i is used
    pricer.data['cons'] = demandCons # List of constraints (demand)
    pricer.data['widths'] = widths # List of integers
    pricer.data['demand'] = demand # List of integers
    pricer.data['rollLength'] = rollLength # Capacity: integer
    pricer.data['patterns'] = patterns # List of patterns (pattern: list of integer)

    # Add the branching rule to the MP
    branchrule = CutBranching(s, cutPatternVars)
    s.includeBranchrule(branchrule, '', '', priority=10000000, maxdepth=999, maxbounddist=1.0)

    conshdlr = MustEndInZeroConshdlr()
    s.includeConshdlr(conshdlr, "", "", propfreq = -1, enfopriority = -10, chckpriority = -10) # propfreq = -1 to disable it, it was at 1 before
    cons = s.createCons(conshdlr, "", modifiable=True, local=True) # modifiable since new vars will be introduced, local=True because of: https://www.scipopt.org/doc/html/group__PublicConstraintMethods.php#ga38a9d97e56deea3042bb6348a4e90d26
    cons.data = SimpleNamespace()
    cons.data.vars = cutPatternVars
    cons.data.model = s
    s.addPyCons(cons)
    
    # solve problem
    s.optimize() # This here solves the MP, using the PP and the branching rule

    # print original data
    printWidths = '\t'.join(str(e) for e in widths)
    print('\nInput Data')
    print('==========')
    print('Roll Length:', rollLength)
    print('Widths:\t', printWidths)
    print('Demand:\t', '\t'.join(str(e) for e in demand))

    # print solution
    widthOutput = [0]*len(widths)
    print('\nResult')
    print('======')
    print('\t\tSol Value', '\tWidths\t', printWidths)
    for i in range(len(pricer.data['var'])):
        rollUsage = 0
        solValue = s.getVal(pricer.data['var'][i])
        if solValue > 0:
            outline = 'Pattern_' + str(i) + ':\t' + str(solValue) + '\t\tCuts:\t '
            for j in range(len(widths)):
                rollUsage += pricer.data['patterns'][i][j]*widths[j]
                widthOutput[j] += pricer.data['patterns'][i][j]*solValue
                outline += str(pricer.data['patterns'][i][j]) + '\t'
            outline += 'Usage:' + str(rollUsage)
            print(outline)

    print('\t\t\tTotal Output:\t', '\t'.join(str(e) for e in widthOutput))

    
if __name__ == '__main__':
    cuttingstock()

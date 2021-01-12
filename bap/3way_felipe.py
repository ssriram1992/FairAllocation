from pyscipopt import *
import random

"""
===================
GETTING THIS TO RUN
===================

0. We have to extend the Python interface a little bit first.
1. Create a virtual environment in a folder.
2. Download pyscipopt from github and extract it inside the virtual environment folder.
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

   (this will do the cython routine to update the files)

6. Everything should work.


=======================
OBJECTIVE OF THIS MODEL
=======================

This is a simple bin packing problem with column generation and branch and price (we merged the two examples cutstock.py and test_branch_probing_lp.py).

We want to do two things:
1. When branching reaches an integer node, we want to check if that node is feasible with regard to an external feasibility checker.
2. If the external feasibility checker says that the node is infeasible, we want to do 3-way branching on an integer variable of that node (e.g. if some variable x has value 3, we branch on x<=2, x==3, x>=4).

The external feasibility checker is given an integer solution, and checks the value of each variable in that solution. If any variable has a value that does not end in zero (i.e., is not a multiple of 10), then the solution is deemed infeasible. [We want to get this simple example working in order to use the same principle on the actual problem we're working on, where the external feasibility checker checks to see if a Hamiltonian path exists in a graph that is constructed from the integer solution of a node.]


=================
WHAT WE HAVE DONE
=================

We use Branchrule and Conshdlr. Here is our understanding of the mechanics of branching:

- If a node has some fractional variables, then branchexeclp() is called.
- If a node does not have any fractional variables, then branchexecps() is called.
- When we reach a fractional solution, consenfolp() is called. We will say that all fractional solutions are feasible since we are interested in checking the feasibility of integer solutions only.
- When we branch on a pseudo-solution, which is integral (right?), consenfops() is called. This is the place where we want to make sure that the value of each variable is a multiple of 10.
- We never seem to be able to branch on pseudo-solutions.

"""


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


modulo = 2
var_ub = 200
seed = 0
relative_gap = 0.01 # stop at 1% gap
random.seed(seed)

class MustEndInZeroConshdlr(Conshdlr):
    # def consenfolp(self, constraints, nusefulconss, solinfeasible):
    #     print('===> Entering consenfolp()')
    #     s = self.model
    #     vars = s.getVars()
    #     for v in vars:
    #         if s.getSolVal(None, v) % 10 != 0:
    #             return {"result": SCIP_RESULT.INFEASIBLE}
    #     return {"result": SCIP_RESULT.FEASIBLE}

    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        '''calls enforcing method of constraint handler for LP solution for all constraints added'''
        # We assume that any LP solution is feasible, since infeasibility is determined only once all variables are integer
        s = self.model
        nodeDiscard = True
        var =  s.getVars(transformed=True) # TRANSFORMED=TRUE??
        for v in var:
            if v.getUbLocal() - v.getLbLocal() > 0.5:
                nodeDiscard = False
        # if VERBOSE >= DEBUG:
        #     print("Discardable node found")
        for v in var:
            if s.getSolVal(None, v) % modulo !=0:
                if nodeDiscard:
                    return {"result": SCIP_RESULT.CUTOFF}
                else:
                    return {"result": SCIP_RESULT.INFEASIBLE}
        return {"result": SCIP_RESULT.FEASIBLE}

    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        print("===> Entering consenfops()")
        s = self.model
        vars = s.getVars(transformed=True) # TRANSFORMED=TRUE??
        for v in vars:
            if s.getSolVal(None, v) % modulo != 0:
                return {"result": SCIP_RESULT.INFEASIBLE}
        return {"result": SCIP_RESULT.FEASIBLE}

    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason, completely):
        print("===> Entering conscheck()")
        assert len(constraints) == 1
        s = self.model
        vars = s.getVars(transformed=True) # TRANSFORMED=TRUE??
        for v in vars:
            if s.getSolVal(solution, v) % modulo != 0:
                return {"result": SCIP_RESULT.INFEASIBLE}
        return {"result": SCIP_RESULT.FEASIBLE}

    def conspresol(self, constraints, nrounds, presoltiming,
                   nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                   nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, result_dict):
        return result_dict


class CutBranching(Branchrule):
    def __init__(self, model, variables):
        self.model = model
        self.variables = variables

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
    

class CutPricer(Pricer):
    def __init__(self):
        self.data = {}

    def pricerinit(self):
        for i, c in enumerate(self.data['cons']):
            self.data['cons'][i] = self.model.getTransformedCons(c)
        
    def pricerredcost(self):
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
            newVar = self.model.addVar("NewPattern_" + str(currentNumVar), vtype="I", obj=1.0, pricedVar=True, ub=var_ub)
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
    s.setRealParam('limits/gap', relative_gap)
    s.setPresolve(0)

    # Pricing problem
    pricer = CutPricer()
    s.includePricer(pricer, "CuttingStockPricer", "Pricer to identify new cutting stock patterns")

    # Problem data
    widths = [14, 31, 36, 45]
    demand = [211, 395, 610, 97]
    rollLength = 100

    # MP variables
    cutPatternVars = []
    varNames = []
    varBaseName = "Pattern"
    patterns = []

    # Initial columns: Variables
    for i in range(len(widths)):
        varNames.append(varBaseName + "_" + str(i))
        cutPatternVars.append(s.addVar(varNames[i], obj = 1.0, vtype = "I", ub=310))

    # Demand constraints
    demandCons = []
    for i in range(len(widths)):
        numWidthsPerRoll = float(int(rollLength/widths[i]))
        demandCons.append(s.addCons(numWidthsPerRoll*cutPatternVars[i] >= demand[i], separate = False, modifiable = True))
        newPattern = [0]*len(widths)
        newPattern[i] = numWidthsPerRoll
        patterns.append(newPattern)

    # Setting the pricer_data for use in the init and redcost functions
    pricer.data['var'] = cutPatternVars
    pricer.data['cons'] = demandCons
    pricer.data['widths'] = widths
    pricer.data['demand'] = demand
    pricer.data['rollLength'] = rollLength
    pricer.data['patterns'] = patterns

    # Add the branching rule to the MP
    branchrule = CutBranching(s, cutPatternVars)
    s.includeBranchrule(branchrule, '', '', priority=100000, maxdepth=65534, maxbounddist=1.0)

    conshdlr = MustEndInZeroConshdlr()
    s.includeConshdlr(conshdlr, "", "", enfopriority = -10, chckpriority = -10)
    cons = s.createCons(conshdlr, "", modifiable=True, local=True)
    cons.data = SimpleNamespace()
    cons.data.vars = cutPatternVars
    cons.data.model = s
    s.addPyCons(cons)
    
    # solve problem
    s.optimize()

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

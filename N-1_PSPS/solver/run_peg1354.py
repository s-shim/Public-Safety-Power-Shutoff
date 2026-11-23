"""Run the (FLOW) formulation on one PEG-1354 instance for the EPSR revision.

The model is built exactly as in 3_L_N-1_PSPS_PFlow_Solver_8_default.py; only the
instance selection, Gurobi parameters and output paths come from the command line.

2x2 design (MIPFocus=3 and GomoryPasses=100 fixed in all four cells):
  --families all   : every Gurobi cut family follows Cuts=<level>
  --families nine  : Cuts=0, and the nine selected families are set to <level>
  --level 1        : conservative (moderate) cut generation
  --level 2        : aggressive cut generation (as in the original selected runs)
"""
import argparse, os, resource, json, time
from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket
import datetime

ap = argparse.ArgumentParser()
ap.add_argument('--inst', type=int, required=True)
ap.add_argument('--families', choices=['all', 'nine', 'default'], required=True,
                help="'default' = pure Gurobi defaults (no cut parameters), still MIPFocus as given")
ap.add_argument('--level', type=int, choices=[1, 2], default=2)
ap.add_argument('--mipfocus', type=int, default=3)
ap.add_argument('--gomorypasses', type=int, default=100)
ap.add_argument('--threads', type=int, default=32)
ap.add_argument('--timelimit', type=float, default=172000)      # a little under 48 h
ap.add_argument('--nodefilestart', type=float, default=400)     # GB
ap.add_argument('--nodefiledir', default=os.environ.get('SCRATCH', '.'))
ap.add_argument('--seed', type=int, default=0, help='Gurobi Seed parameter')
ap.add_argument('--data', default='../../matpower/case1354pegase')
ap.add_argument('--nbus', type=int, default=1354, help='bus count in the CSV file names (for small test cases)')
ap.add_argument('--out', default='results')
args = ap.parse_args()

NINE = ['MIRCuts', 'ImpliedCuts', 'FlowCoverCuts', 'ZeroHalfCuts', 'ProjImpliedCuts',
        'RelaxLiftCuts', 'CoverCuts', 'RLTCuts']   # + Gomory via GomoryPasses
if args.families == 'default':
    tag = 'default_mf%d' % args.mipfocus
else:
    tag = '%s_L%d_mf%d' % (args.families, args.level, args.mipfocus)
runname = 'n%d_inst%d_%s' % (args.nbus, args.inst, tag)
os.makedirs(args.out, exist_ok=True)
t_start = time.time()

machineName = socket.gethostname()
laziness = 0
relax = False
spanning = False
print(machineName)
print(datetime.datetime.now())
print(vars(args))

numNodes = args.nbus
sizeArray = []
instArray = []
ratioDemandRequiredArray = []
totalDemandArray = []
totalRiskArray = []
resultRiskArray = []
ratioRiskArray = []
timeArray = []
nodecountArray = []
securityArray = []
machineArray = []
spanningArray = []
inst = args.inst
lines = pd.read_csv('%s/lines_%s_inst%s.csv'%(args.data,numNodes,inst))
nodes = pd.read_csv('%s/nodes_%s.csv'%(args.data,numNodes))
        
dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = md.profile(lines,nodes)

ratioDemandRequired = 0.9
demandRequired = totalDemand * ratioDemandRequired

        
model = Model('N-1 PSPS')
model.setParam('LogFile', '%s/grblog_%s.txt' % (args.out, runname))
model.setParam('Threads', args.threads)
model.setParam('TimeLimit', args.timelimit)
model.setParam('NodefileStart', args.nodefilestart)
model.setParam('NodefileDir', args.nodefiledir)
model.setParam('Seed', args.seed)
model.setParam('MIPFocus', args.mipfocus)
if args.families == 'all':
    model.setParam('Cuts', args.level)
    model.setParam('GomoryPasses', args.gomorypasses)
elif args.families == 'nine':
    model.setParam('Cuts', 0)
    for p in NINE:
        model.setParam(p, args.level)
    model.setParam('GomoryPasses', args.gomorypasses)


x_vars = []
x_names = []
for i in nodeArray:
    x_vars += [(i)]
    x_names += ['X[%s]'%i]
X = model.addVars(x_vars, vtype = GRB.BINARY, name = x_names)
    
z_vars = []
z_names = []
for l in lineArray:
    z_vars += [(l)]
    z_names += ['Z[%s]'%l]
Z = model.addVars(z_vars, vtype = GRB.BINARY, name = z_names)

d_vars = []
d_names = []
for i in nodeArray:
    d_vars += [(i)]
    d_names += ['D[%s]'%i]
D = model.addVars(d_vars, vtype = GRB.CONTINUOUS, name = d_names)

g_vars = []
g_names = []
for i in nodeArray:
    g_vars += [(i)]
    g_names += ['Gen[%s]'%i]
Gen = model.addVars(g_vars, vtype = GRB.CONTINUOUS, name = g_names)

pf_vars = []
pf_names = []
for i in nodeArray:
    if i != 1:
        for l in lineArray:
            pf_vars += [(i,l,sourceNode[l],targetNode[l])]
            pf_names += ['PF[%s,%s,%s,%s]'%(i,l,sourceNode[l],targetNode[l])]
            pf_vars += [(i,l,targetNode[l],sourceNode[l])]
            pf_names += ['PF[%s,%s,%s,%s]'%(i,l,targetNode[l],sourceNode[l])]
PF = model.addVars(pf_vars, vtype = GRB.CONTINUOUS, name = pf_names)

pfc_vars = []
pfc_names = []
for c in theContingent:
    for l in lineArray:
        pfc_vars += [(c,l,sourceNode[l],targetNode[l])]
        pfc_names += ['PFC[%s,%s,%s,%s]'%(c,l,sourceNode[l],targetNode[l])]
        pfc_vars += [(c,l,targetNode[l],sourceNode[l])]
        pfc_names += ['PFC[%s,%s,%s,%s]'%(c,l,targetNode[l],sourceNode[l])]
PFC = model.addVars(pfc_vars, vtype = GRB.CONTINUOUS, name = pfc_names)
    

LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
model.addConstr(LinExpr(LHS)>=demandRequired, name='demandRequired')


LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
    LHS += [(-1,Gen[i])]    
model.addConstr(LinExpr(LHS)==0, name='demand supply balance')


for i in nodeArray:
    LHS = []
    LHS += [(1,D[i])]
    LHS += [(-dBound[i],X[i])]
    model.addConstr(LinExpr(LHS)<=0, name='demand bound (%s)'%i)


for i in nodeArray:
    LHS_low = []
    LHS_up = []

    LHS_low += [(1,Gen[i])]
    LHS_up += [(1,Gen[i])]

    LHS_low += [(-minGen[i],X[i])]
    LHS_up += [(-maxGen[i],X[i])]

    #model.addConstr(LinExpr(LHS_low)>=0, name='generation lower bound (%s)'%i)
    model.addConstr(LinExpr(LHS_up)<=0, name='generation upper bound (%s)'%i)


LHS = [(1,X[1])]
model.addConstr(LinExpr(LHS)==1, name='reference node')


for l in lineArray:
    LHS_source = [(1,Z[l]),(-1,X[sourceNode[l]])]
    LHS_target = [(1,Z[l]),(-1,X[targetNode[l]])]

    model.addConstr(LinExpr(LHS_source)<=0, name='source bound (%s)'%l)
    model.addConstr(LinExpr(LHS_target)<=0, name='target bound (%s)'%l)


# =============================================================================
#     if switchableFunction[l] == 0:
#         LHS_edge = [(-1,Z[l]),(1,X[sourceNode[l]]),(1,X[targetNode[l]])]
#         model.addConstr(LinExpr(LHS_edge)<=1, name='edge inequality (%s)'%l)
# =============================================================================

    if switchableFunction[l] == 0:
        LHS_edge1 = [(-1,Z[l]),(1,X[sourceNode[l]])]
        LHS_edge2 = [(-1,Z[l]),(1,X[targetNode[l]])]
        model.addConstr(LinExpr(LHS_edge1)==0, name='edge inequality1 (%s)'%l)
        model.addConstr(LinExpr(LHS_edge2)==0, name='edge inequality2 (%s)'%l)


for i in nodeArray:
    if i != 1:
        for l in lineArray:
            LHSF1 = [(-1,Z[l]),(1,PF[i,l,sourceNode[l],targetNode[l]]),(1,PF[i,l,targetNode[l],sourceNode[l]])]
            model.addConstr(LinExpr(LHSF1)<=0, name='flow edge bound 1 (%s,%s)'%(i,l))

for i in nodeArray:
    if i != 1:
        LHSLP1 = []
        LHSLM1 = [(-1,X[i])]
        LHSLPi = [(-1,X[i])]
        LHSLMi = []
        for l in lineArray:

            if sourceNode[l] == 1:
                LHSLP1 += [(1,PF[i,l,targetNode[l],sourceNode[l]])]
        
            if targetNode[l] == 1:
                LHSLP1 += [(1,PF[i,l,sourceNode[l],targetNode[l]])]
        
            if sourceNode[l] == 1:
                LHSLM1 += [(1,PF[i,l,sourceNode[l],targetNode[l]])]
        
            if targetNode[l] == 1:
                LHSLM1 += [(1,PF[i,l,targetNode[l],sourceNode[l]])]

            if sourceNode[l] == i:
                LHSLPi += [(1,PF[i,l,targetNode[l],sourceNode[l]])]
        
            if targetNode[l] == i:
                LHSLPi += [(1,PF[i,l,sourceNode[l],targetNode[l]])]

            if sourceNode[l] == i:
                LHSLMi += [(1,PF[i,l,sourceNode[l],targetNode[l]])]
        
            if targetNode[l] == i:
                LHSLMi += [(1,PF[i,l,targetNode[l],sourceNode[l]])]

        model.addConstr(LinExpr(LHSLP1)==0, name='zero inflow to 1 (%s)'%i)#.Lazy = laziness
        model.addConstr(LinExpr(LHSLM1)==0, name='sumX outflow from 1 (%s)'%i)#.Lazy = laziness
        model.addConstr(LinExpr(LHSLPi)==0, name='sumX inflow to i (%s)'%i)#.Lazy = laziness
        model.addConstr(LinExpr(LHSLMi)==0, name='zero outflow from i (%s)'%i)#.Lazy = laziness

        for j in nodeArray:
            if j != i and j != 1:
              
                LHSKCH0 = []
                LHSKCH1 = [(-1,X[i])]
                for l in lineArray:
                    if sourceNode[l] == j:
                        LHSKCH0 += [(1,PF[i,l,targetNode[l],sourceNode[l]])]
                        LHSKCH0 += [(-1,PF[i,l,sourceNode[l],targetNode[l]])]
                        LHSKCH1 += [(1,PF[i,l,targetNode[l],sourceNode[l]])]
                
                    if targetNode[l] == j:
                        LHSKCH0 += [(1,PF[i,l,sourceNode[l],targetNode[l]])]
                        LHSKCH0 += [(-1,PF[i,l,targetNode[l],sourceNode[l]])]
                        LHSKCH1 += [(1,PF[i,l,sourceNode[l],targetNode[l]])]
        
                model.addConstr(LinExpr(LHSKCH0)==0, name='Kirchhoff at (%s,%s)'%(j,i)).Lazy = laziness            
                model.addConstr(LinExpr(LHSKCH1)<=0, name='Kirchhoff Limit at (%s,%s)'%(j,i)).Lazy = laziness            
        


for c in theContingent:
    LHS = [(1,PFC[c,c,sourceNode[c],targetNode[c]]),(1,PFC[c,c,targetNode[c],sourceNode[c]])]
    model.addConstr(LinExpr(LHS)==0, name='Contingency at (%s)'%c).Lazy = laziness            

    for l in lineArray:
        LHSF1 = [(-1,Z[l]),(1,PFC[c,l,sourceNode[l],targetNode[l]]),(1,PFC[c,l,targetNode[l],sourceNode[l]])]
        model.addConstr(LinExpr(LHSF1)<=0, name='flow edge bound 1 (%s,%s)'%(c,l)).Lazy = laziness

    LHSLPs = []
    LHSLMs = []
    LHSLMs += [(-1,Z[c])]
    LHSLPt = []
    LHSLPt += [(-1,Z[c])]
    LHSLMt = []

    for l in lineArray:
    
        if sourceNode[l] == sourceNode[c]:
            LHSLPs += [(1,PFC[c,l,targetNode[l],sourceNode[l]])]
    
        if targetNode[l] == sourceNode[c]:
            LHSLPs += [(1,PFC[c,l,sourceNode[l],targetNode[l]])]
    
        if sourceNode[l] == sourceNode[c]:
            LHSLMs += [(1,PFC[c,l,sourceNode[l],targetNode[l]])]
    
        if targetNode[l] == sourceNode[c]:
            LHSLMs += [(1,PFC[c,l,targetNode[l],sourceNode[l]])]

        if sourceNode[l] == targetNode[c]:
            LHSLPt += [(1,PFC[c,l,targetNode[l],sourceNode[l]])]
    
        if targetNode[l] == targetNode[c]:
            LHSLPt += [(1,PFC[c,l,sourceNode[l],targetNode[l]])]
    
        if sourceNode[l] == targetNode[c]:
            LHSLMt += [(1,PFC[c,l,sourceNode[l],targetNode[l]])]
    
        if targetNode[l] == targetNode[c]:
            LHSLMt += [(1,PFC[c,l,targetNode[l],sourceNode[l]])]

    model.addConstr(LinExpr(LHSLPs)==0, name='zero inflow into sourceNode[%s]'%c)#.Lazy = laziness
    model.addConstr(LinExpr(LHSLMs)==0, name='Z[c] outflow from sourceNode[%s]'%c)#.Lazy = laziness
    model.addConstr(LinExpr(LHSLMt)==0, name='zero outflow from targetNode[%s]'%c)#.Lazy = laziness
    model.addConstr(LinExpr(LHSLPt)==0, name='Z[c] inflow from targetNode[%s]'%c)#.Lazy = laziness
    
    for i in nodeArray:
        if i != sourceNode[c] and i != targetNode[c]:
            
            LHSKCH0 = []
            LHSKCH1 = [(-1,X[i])]
            for l in lineArray:
                if sourceNode[l] == i:
                    LHSKCH0 += [(+1,PFC[c,l,targetNode[l],sourceNode[l]])]
                    LHSKCH0 += [(-1,PFC[c,l,sourceNode[l],targetNode[l]])]
                    LHSKCH1 += [(+1,PFC[c,l,targetNode[l],sourceNode[l]])]
            
                if targetNode[l] == i:
                    LHSKCH0 += [(+1,PFC[c,l,sourceNode[l],targetNode[l]])]
                    LHSKCH0 += [(-1,PFC[c,l,targetNode[l],sourceNode[l]])]
                    LHSKCH1 += [(+1,PFC[c,l,sourceNode[l],targetNode[l]])]
    
            model.addConstr(LinExpr(LHSKCH0)==0, name='Kirchhoff at (%s,%s)'%(c,i)).Lazy = laziness            
            model.addConstr(LinExpr(LHSKCH1)<=0, name='Kirchhoff Limit at (%s,%s)'%(c,i)).Lazy = laziness            
            

    
if spanning == True:
    for i in nodeArray:
        LHS = [(1,X[i])]
        model.addConstr(LinExpr(LHS)==1, name='Fix X (%s)'%i)



    
objTerms = []
for l in lineArray:
    objTerms += [(riskFunction[l],Z[l])]


                                
model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
model.update()
if relax == True:
    model = model.relax()
model.optimize()
if model.SolCount == 0:
    raise SystemExit('no incumbent found')



# read the optimal solution
variableName = ['obj']
variableValue = [model.getObjective().getValue()]
resultRisk = 0.0
resultDemand = 0.0
zValue = {}
for v in model.getVars():
    variableName += [v.varname]
    variableValue += [v.x]
    
    if v.varname[0] == 'D':
        resultDemand += v.x

    if v.varname[0] == 'Z':
        l = int(v.varname[2:-1])
        zValue[l] = 0
        if v.x + 0.0001 > 1:
            resultRisk += riskFunction[l]
            zValue[l] = 1
            # print(v.varname,l,riskFunction[l],v.x)



optSolution = pd.DataFrame(list(zip(variableName, variableValue)),columns =['varName', 'varVal'])
optSolution.to_csv('%s/opt_%s.csv' % (args.out, runname), index = False)#Check

solution = copy.deepcopy(optSolution)
closed = []
for n in solution['varName']:
    if n[0] == 'Z':
        i = int(n[2:-1])
        [val_i] = solution.loc[solution['varName']==n,'varVal']
        if val_i > 1 - 0.0001:
            closed += [i]

G, key, lineFunction = md.multigraph(nodes, lines, lineArray)

subG = nx.MultiGraph()
subContingent = []
for l in closed:
    subG.add_edge(sourceNode[l],targetNode[l],key[l])
    if contingentFunction[l] == 1:
        subContingent += [l]

if relax == False:
    security = nx.is_connected(subG)       
    for c in subContingent:
        copySubG = copy.deepcopy(subG)
        copySubG.remove_edge(sourceNode[c],targetNode[c],key[c])
        if nx.has_path(copySubG,sourceNode[c],targetNode[c]) == False:
            security = False
    
    print(security)
    
    
    sizeArray += [len(nodeArray)]
    instArray += [inst]
    ratioDemandRequiredArray += [ratioDemandRequired]
    totalDemandArray += [totalDemand]
    totalRiskArray += [totalRisk]
    resultRiskArray += [resultRisk]
    ratioRiskArray += [resultRisk/totalRisk]
    timeArray += [model.Runtime]
    nodecountArray += [model.NodeCount]
    securityArray += [security]
    machineArray += [machineName]
    spanningArray += [spanning]
    
    
    resultTable = pd.DataFrame(list(zip(sizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, securityArray, machineArray, spanningArray)),columns =['Nodes', 'inst', 'totalDemand','ratioDemandRequired', 'totalRiskArray', 'minRisk', 'ratioRisk', 'Time', 'B&B Nodes', 'Security', 'Machine', 'Spanning'])
    resultTable.to_csv('%s/result_%s.csv' % (args.out, runname), index = False)#Check




# ---- run summary (peak memory of this process, incl. Gurobi threads) ----
peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
summary = dict(run=runname, inst=args.inst, families=args.families, level=args.level,
               mipfocus=args.mipfocus, gomorypasses=args.gomorypasses, threads=args.threads,
               status=model.Status, obj=model.ObjVal, bound=model.ObjBound, gap=model.MIPGap,
               runtime=model.Runtime, nodes=model.NodeCount, peak_mem_gb=round(peak_gb, 1),
               wall_s=round(time.time() - t_start, 1), host=machineName,
               gurobi='%d.%d.%d' % gurobi.version())
print(summary)
json.dump(summary, open('%s/summary_%s.json' % (args.out, runname), 'w'), indent=1)

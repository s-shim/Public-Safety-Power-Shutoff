from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket

machineName = socket.gethostname()
laziness = 0
relax = True
spanning = False

#numNodes = 30
numNodes = 118
#numNodes = 300
#numNodes = 1354


folderName = 'matpower/case%s'%(numNodes)
#folderName = 'matpower/case%spegase'%(numNodes)


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


inst = 1


lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))
        
dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = md.profile(lines,nodes)

ratioDemandRequired = 0.9
demandRequired = totalDemand * ratioDemandRequired

        
model = Model('N-1 PSPS')
model.setParam('LogFile', 'result%s/grblog_%s_Relax%s_N-1_PSPS_Span%s_Flow_%s_inst%s.txt'%(len(nodeArray),machineName,relax,spanning,len(nodeArray),inst))    

x_vars = []
x_names = []
for i in nodeArray:
    x_vars += [(i)]
    x_names += ['X[%s]'%i]
#X = model.addVars(x_vars, vtype = GRB.BINARY, name = x_names)
X = model.addVars(x_vars, vtype = GRB.CONTINUOUS, name = x_names)
  
z_vars = []
z_names = []
for l in lineArray:
    z_vars += [(l)]
    z_names += ['Z[%s]'%l]
#Z = model.addVars(z_vars, vtype = GRB.BINARY, name = z_names)
Z = model.addVars(z_vars, vtype = GRB.CONTINUOUS, name = z_names)   

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

f_vars = []
f_names = []
for l in lineArray:
    f_vars += [(l,sourceNode[l],targetNode[l])]
    f_names += ['F[%s,%s,%s]'%(l,sourceNode[l],targetNode[l])]
    f_vars += [(l,targetNode[l],sourceNode[l])]
    f_names += ['F[%s,%s,%s]'%(l,targetNode[l],sourceNode[l])]
F = model.addVars(f_vars, vtype = GRB.CONTINUOUS, name = f_names)

fc_vars = []
fc_names = []
for c in theContingent:
    for l in lineArray:
        fc_vars += [(c,l,sourceNode[l],targetNode[l])]
        fc_names += ['FC[%s,%s,%s,%s]'%(c,l,sourceNode[l],targetNode[l])]
        fc_vars += [(c,l,targetNode[l],sourceNode[l])]
        fc_names += ['FC[%s,%s,%s,%s]'%(c,l,targetNode[l],sourceNode[l])]
FC = model.addVars(fc_vars, vtype = GRB.CONTINUOUS, name = fc_names)
    

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


LHSLP1 = []
LHSLM1 = []
for l in lineArray:
    LHSF = [(-len(nodeArray),Z[l]),(1,F[l,sourceNode[l],targetNode[l]]),(1,F[l,targetNode[l],sourceNode[l]])]
    model.addConstr(LinExpr(LHSF)<=0, name='flow edge bound (%s)'%l)

    if sourceNode[l] == 1:
        LHSLP1 += [(1,F[l,targetNode[l],sourceNode[l]])]

    if targetNode[l] == 1:
        LHSLP1 += [(1,F[l,sourceNode[l],targetNode[l]])]

    if sourceNode[l] == 1:
        LHSLM1 += [(1,F[l,sourceNode[l],targetNode[l]])]

    if targetNode[l] == 1:
        LHSLM1 += [(1,F[l,targetNode[l],sourceNode[l]])]

for i in nodeArray:
    if i != 1:
        LHSLM1 += [(-1,X[i])]
        
        LHSKCH = [(-1,X[i])]
        for l in lineArray:
            if sourceNode[l] == i:
                LHSKCH += [(1,F[l,targetNode[l],sourceNode[l]])]
                LHSKCH += [(-1,F[l,sourceNode[l],targetNode[l]])]
        
            if targetNode[l] == i:
                LHSKCH += [(1,F[l,sourceNode[l],targetNode[l]])]
                LHSKCH += [(-1,F[l,targetNode[l],sourceNode[l]])]

        model.addConstr(LinExpr(LHSKCH)==0, name='Kirchhoff at (%s)'%i)            
        
model.addConstr(LinExpr(LHSLP1)==0, name='zero inflow to 1')
model.addConstr(LinExpr(LHSLM1)==0, name='sumX outflow to 1')


for c in theContingent:
    LHS = [(1,FC[c,c,sourceNode[c],targetNode[c]]),(1,FC[c,c,targetNode[c],sourceNode[c]])]
    model.addConstr(LinExpr(LHS)==0, name='Contingency at (%s)'%c).Lazy = laziness            

    LHSLP1 = []
    LHSLM1 = []
    for l in lineArray:
        LHSF = [(-len(nodeArray),Z[l]),(1,FC[c,l,sourceNode[l],targetNode[l]]),(1,FC[c,l,targetNode[l],sourceNode[l]])]
        model.addConstr(LinExpr(LHSF)<=0, name='flow edge bound (%s,%s)'%(c,l)).Lazy = laziness
    
        if sourceNode[l] == 1:
            LHSLP1 += [(1,FC[c,l,targetNode[l],sourceNode[l]])]
    
        if targetNode[l] == 1:
            LHSLP1 += [(1,FC[c,l,sourceNode[l],targetNode[l]])]
    
        if sourceNode[l] == 1:
            LHSLM1 += [(1,FC[c,l,sourceNode[l],targetNode[l]])]
    
        if targetNode[l] == 1:
            LHSLM1 += [(1,FC[c,l,targetNode[l],sourceNode[l]])]
    
    for i in nodeArray:
        if i != 1:
            LHSLM1 += [(-1,X[i])]
            
            LHSKCH = [(-1,X[i])]
            for l in lineArray:
                if sourceNode[l] == i:
                    LHSKCH += [(1,FC[c,l,targetNode[l],sourceNode[l]])]
                    LHSKCH += [(-1,FC[c,l,sourceNode[l],targetNode[l]])]
            
                if targetNode[l] == i:
                    LHSKCH += [(1,FC[c,l,sourceNode[l],targetNode[l]])]
                    LHSKCH += [(-1,FC[c,l,targetNode[l],sourceNode[l]])]
    
            model.addConstr(LinExpr(LHSKCH)==0, name='Kirchhoff at (%s,%s)'%(c,i)).Lazy = laziness            
            
    model.addConstr(LinExpr(LHSLP1)==0, name='zero inflow to 1 (%s)'%c).Lazy = laziness
    model.addConstr(LinExpr(LHSLM1)==0, name='sumX outflow to 1 (%s)'%c).Lazy = laziness

    
if spanning == True:
    for i in nodeArray:
        LHS = [(1,X[i])]
        model.addConstr(LinExpr(LHS)==1, name='Fix X (%s)'%i)



    
objTerms = []
for l in lineArray:
    objTerms += [(riskFunction[l],Z[l])]


                                
model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
model.update()
model.optimize()



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
optSolution.to_csv(r'result%s/opt_N-1_PSPS_Span%s_Flow_Relax%s_%s_inst%s.csv'%(len(nodeArray),spanning,relax,len(nodeArray),inst), index = False)#Check

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
    resultTable.to_csv(r'result%s/result_%s_N-1_PSPS_Span%s_Flow_%s.csv'%(len(nodeArray),machineName,spanning,len(nodeArray)), index = False)#Check



from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket
import time


def relax(model,lineArray,nodeArray,riskFunction,sourceNode,targetNode,minGen,maxGen,dBound):
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
#         if switchableFunction[l] == 0:
#             LHS_edge = [(-1,Z[l]),(1,X[sourceNode[l]]),(1,X[targetNode[l]])]
#             model.addConstr(LinExpr(LHS_edge)<=1, name='edge inequality (%s)'%l)
# =============================================================================

        if switchableFunction[l] == 0:
            LHS_edge1 = [(-1,Z[l]),(1,X[sourceNode[l]])]
            LHS_edge2 = [(-1,Z[l]),(1,X[targetNode[l]])]
            model.addConstr(LinExpr(LHS_edge1)==0, name='edge inequality1 (%s)'%l)
            model.addConstr(LinExpr(LHS_edge2)==0, name='edge inequality2 (%s)'%l)

        
    if spanning == True:
        for i in nodeArray:
            LHS = [(1,X[i])]
            model.addConstr(LinExpr(LHS)==1, name='Fix X (%s)'%i)

            
    objTerms = []
    for l in lineArray:
        objTerms += [(riskFunction[l],Z[l])]
                                
    model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
    
    return model


def BC(model,xVal,zVal,lineArray,nodeArray,sourceNode,targetNode,simpleG):
    cutsAdded = False
    
    for l in lineArray:
        simpleG[sourceNode[l]][targetNode[l]]['capacity'] = 0
    
    for l in lineArray:
        simpleG[sourceNode[l]][targetNode[l]]['capacity'] += zVal[l]
        
    for i in nodeArray:
        if i != 1:
            cut_value, partition = nx.minimum_cut(simpleG, 1, i, capacity = 'capacity')
            reachable, non_reachable = partition
            if cut_value < xVal[i]:
                LHS = [(-1,X[i])]
                for l in lineArray:
                    if sourceNode[l] in reachable and targetNode[l] in non_reachable:
                        LHS += [(1,Z[l])]
                    if sourceNode[l] in non_reachable and targetNode[l] in reachable:
                        LHS += [(1,Z[l])]
                model.addConstr(LinExpr(LHS)>=0)
                cutsAdded = True
                                
    for c in theContingent:
        simpleGC = copy.deepcopy(simpleG)
        simpleGC[sourceNode[c]][targetNode[c]]['capacity'] = simpleGC[sourceNode[c]][targetNode[c]]['capacity'] - zVal[c]
        cut_valueC, partitionC = nx.minimum_cut(simpleGC, sourceNode[c], targetNode[c], capacity = 'capacity')
        reachableC, non_reachableC = partitionC
        if cut_valueC < xVal[sourceNode[c]] + xVal[targetNode[c]] - 1:
            LHSC = [(-1,X[sourceNode[c]]),(-1,X[targetNode[c]])]
            for l in lineArray:
                if l != c:
                    if sourceNode[l] in reachableC and targetNode[l] in non_reachableC:
                        LHSC += [(1,Z[l])]
                    if sourceNode[l] in non_reachableC and targetNode[l] in reachableC:
                        LHSC += [(1,Z[l])]
            model.addConstr(LinExpr(LHSC)>=-1)
            cutsAdded = True

    return model, cutsAdded



machineName = socket.gethostname()
spanning = False


#numNodes = 30
#numNodes = 118
numNodes = 300
#numNodes = 1354


folderName = 'matpower/case%s'%(numNodes)
#folderName = 'matpower/case%spegase'%(numNodes)

inst = 1


lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))
        
dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = md.profile(lines,nodes)
G, key, lineFunction = md.multigraph(nodes, lines, lineArray)
simpleG = nx.Graph(G)

ratioDemandRequired = 0.9
demandRequired = totalDemand * ratioDemandRequired

        
model = Model('N-1 PSPS Relax')
#model.setParam('LogFile', 'result%s/grblog_%s_N-1_PSPS_Span%s_BCRelax_%s_inst%s.txt'%(len(nodeArray),machineName,spanning,len(nodeArray),inst))    
model.setParam('OutputFlag', 0)

x_vars = []
x_names = []
for i in nodeArray:
    x_vars += [(i)]
    x_names += ['X[%s]'%i]
X = model.addVars(x_vars, vtype = GRB.CONTINUOUS, name = x_names)

for i in nodeArray:
    LHS = [(1,X[i])]
    model.addConstr(LinExpr(LHS)<=1, name='X[%s]<=1'%i)
    
z_vars = []
z_names = []
for l in lineArray:
    z_vars += [(l)]
    z_names += ['Z[%s]'%l]
Z = model.addVars(z_vars, vtype = GRB.CONTINUOUS, name = z_names)

for l in lineArray:
    LHS = [(1,Z[l])]
    model.addConstr(LinExpr(LHS)<=1, name='Z[%s]<=1'%l)

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

LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
model.addConstr(LinExpr(LHS)>=demandRequired, name='demandRequired')

LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
    LHS += [(-1,Gen[i])]    
model.addConstr(LinExpr(LHS)==0, name='demand supply balance')

tic = time.time()
model = relax(model,lineArray,nodeArray,riskFunction,sourceNode,targetNode,minGen,maxGen,dBound)
toc = time.time()
CPUTime = toc - tic

bestObj = -1 # model.getObjective().getValue()
  
run = 0
cutsAdded = True
improve = True
while cutsAdded == True and improve == True:
    model.update()
    model.optimize()    

    toc = time.time()
    CPUTime = toc - tic

    obj = model.getObjective()
    objVal = obj.getValue()
    
    if bestObj < objVal:
        print()
        print('### Run =',run)
        print('###CPUTime =',CPUTime)
        print('### objVal =',objVal)

        xVal = {}
        zVal = {}
        for v in model.getVars():
            if v.varname[0] == 'Z':
                l = int(v.varname[2:-1])
                zVal[l] = v.x
            if v.varname[0] == 'X':
                i = int(v.varname[2:-1])
                xVal[i] = v.x

        bestObj = objVal
        improve = True
        run += 1
        model, cutsAdded = BC(model,xVal,zVal,lineArray,nodeArray,sourceNode,targetNode,simpleG)
                                
    else:
        improve = False
    





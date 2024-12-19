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


laziness = 0
machineName = socket.gethostname()
spanning = False


#numNodes = 30
numNodes = 118
#numNodes = 300
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



tic = time.time()
model = relax(model,lineArray,nodeArray,riskFunction,sourceNode,targetNode,minGen,maxGen,dBound)
toc = time.time()
CPUTime = toc - tic

bestObj = -1 # model.getObjective().getValue()
minRiskIncumbent = totalRisk
minRiskRun = -1
minRiskTime = toc - tic
  
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

        def zFunction(l):
            return zVal[l]

        sortedLineArray = copy.deepcopy(lineArray)
        sortedLineArray.sort(key = zFunction, reverse = True)
        
        closed = []
        for l in sortedLineArray:
            if zVal[l] >= 0.5:
                closed += [l]
        minClosed = copy.deepcopy(closed)

        subG = nx.MultiGraph()
        subContingent = []
        riskIncumbent = 0.0
        for l in closed:
            subG.add_edge(sourceNode[l],targetNode[l],key[l])
            riskIncumbent += riskFunction[l]
            if contingentFunction[l] == 1:
                subContingent += [l]

        sumDemand = 0.0
        sumGen = 0.0
        for i in subG.nodes():
            sumDemand += dBound[i]
            sumGen += maxGen[i]
        incumbentDemand = min(sumDemand,sumGen)
        
        #print('### demand feasibility: 0 <=',incumbentDemand - demandRequired)
        
        if incumbentDemand - demandRequired >= 0:
            feasibility = True
        else: 
            feasibility = False
            
        while feasibility == False:
            l = sortedLineArray[len(closed)]
            if sourceNode[l] not in subG.nodes():
                i = sourceNode[l]
                sumDemand += dBound[i]
                sumGen += maxGen[i]
            if targetNode[l] not in subG.nodes():
                i = targetNode[l]
                sumDemand += dBound[i]
                sumGen += maxGen[i]
            incumbentDemand = min(sumDemand,sumGen)

            closed += [l]
            subG.add_edge(sourceNode[l],targetNode[l],key[l])
            riskIncumbent += riskFunction[l]
            if contingentFunction[l] == 1:
                subContingent += [l]

            if incumbentDemand - demandRequired >= 0:
                feasibility = True
            else: 
                feasibility = False

        security = False
        while security == False:
            security = nx.is_connected(subG)       
            if security == True:
                for c in subContingent:
                    copySubG = copy.deepcopy(subG)
                    copySubG.remove_edge(sourceNode[c],targetNode[c],key[c])
                    if nx.has_path(copySubG,sourceNode[c],targetNode[c]) == False:
                        security = False
                        break
            if security == False:
                l = sortedLineArray[len(closed)]
                if sourceNode[l] not in subG.nodes():
                    i = sourceNode[l]
                    sumDemand += dBound[i]
                    sumGen += maxGen[i]
                if targetNode[l] not in subG.nodes():
                    i = targetNode[l]
                    sumDemand += dBound[i]
                    sumGen += maxGen[i]
                incumbentDemand = min(sumDemand,sumGen)
    
                closed += [l]
                subG.add_edge(sourceNode[l],targetNode[l],key[l])
                riskIncumbent += riskFunction[l]
                if contingentFunction[l] == 1:
                    subContingent += [l]
                
        print('### riskIncumbent =',riskIncumbent,'; totalRisk =',totalRisk)      
        
        if minRiskIncumbent > riskIncumbent:
            minRiskIncumbent = riskIncumbent
            minRiskRun = run
            toc = time.time()
            minRiskTime = toc - tic
            minClosed = copy.deepcopy(closed)

        print('### minRiskIncumbent =',minRiskIncumbent,'at minRiskRun =',minRiskRun,'and minRiskTime =',minRiskTime)              
        
        bestObj = objVal
        improve = True
        run += 1
        model, cutsAdded = BC(model,xVal,zVal,lineArray,nodeArray,sourceNode,targetNode,simpleG)
                                
    else:
        improve = False
    





from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket
import datetime

machineName = socket.gethostname()
spanning = False
print(machineName)
print(datetime.datetime.now())


#numNodes = 30
#numNodes = 118
#numNodes = 300
numNodes = 1354


#folderName = 'matpower/case%s'%(numNodes)
folderName = 'matpower/case%spegase'%(numNodes)


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
G, key, lineFunction = md.multigraph(nodes, lines, lineArray)
simpleG = nx.Graph(G)
print('total risk =',totalRisk)

ratioDemandRequired = 0.9
demandRequired = totalDemand * ratioDemandRequired

        
model = Model('N-1 PSPS')
model.setParam('LogFile', 'result%s/grblog_%s_N-1_PSPS_Span%s_BC_%s_inst%s.txt'%(len(nodeArray),machineName,spanning,len(nodeArray),inst))    

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


LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
model.addConstr(LinExpr(LHS)>=demandRequired, name='demandRequired')


LHS = []
for i in nodeArray:
    LHS += [(1,D[i])]
    LHS += [(-1,Gen[i])]    
model.addConstr(LinExpr(LHS)<=0, name='demand supply balance')


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
        print('non-switchable =',l,sourceNode[l],targetNode[l])
        LHS_edge1 = [(-1,Z[l]),(1,X[sourceNode[l]])]
        LHS_edge2 = [(-1,Z[l]),(1,X[targetNode[l]])]
        model.addConstr(LinExpr(LHS_edge1)==0, name='edge inequality1 (%s)'%l)
        model.addConstr(LinExpr(LHS_edge2)==0, name='edge inequality2 (%s)'%l)


model._minRiskIncumbent = totalRisk + 1

def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        xVal = {}
        zVal = {}
        for i in nodeArray:
            xVal[i] = model.cbGetSolution(X[i])

        for l in lineArray:
            simpleG[sourceNode[l]][targetNode[l]]['capacity'] = 0

        for l in lineArray:
            zVal[l] = model.cbGetSolution(Z[l])
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
                    model.cbLazy(LinExpr(LHS)>=0)
                                        
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
                model.cbLazy(LinExpr(LHSC)>=-1)
                
    if where == GRB.callback.MIPNODE:
        status = model.cbGet(GRB.callback.MIPNODE_STATUS)                
        if status == GRB.OPTIMAL:
            # make a list of edges selected in the solution
            xVal = {}
            zVal = {}
            for i in nodeArray:
                xVal[i] = model.cbGetNodeRel(X[i])
            
            for l in lineArray:
                zVal[l] = model.cbGetNodeRel(Z[l])
    
            def zFunction(l):
                return zVal[l]
    
            sortedLineArray = copy.deepcopy(lineArray)
            sortedLineArray.sort(key = zFunction, reverse = True)
            
            closed = []
# =============================================================================
#             for l in sortedLineArray:
#                 if zVal[l] >= 0.5:
#                     closed += [l]
# =============================================================================
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
                                
            if model._minRiskIncumbent > riskIncumbent:
                model._minRiskIncumbent = riskIncumbent
                minClosed = copy.deepcopy(closed)
                zINT = {}
                xINT = {}
                for l in lineArray:
                    zINT[l] = 0
                for l in closed:
                    zINT[l] += 1
                for i in nodeArray:
                    xINT[i] = 0
                for i in subG.nodes():
                    xINT[i] += 1
                for l in lineArray:
                    model.cbSetSolution(Z[l], zINT[l])
                for i in nodeArray:
                    model.cbSetSolution(X[i], xINT[i])                    
                        
                print('### minRiskIncumbent =',model._minRiskIncumbent,'; totalRisk =',totalRisk)              



if spanning == True:
    for i in nodeArray:
        LHS = [(1,X[i])]
        model.addConstr(LinExpr(LHS)==1, name='Fix X (%s)'%i)



    
objTerms = []
for l in lineArray:
    objTerms += [(riskFunction[l],Z[l])]

                                
model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
model.update()
model.Params.lazyConstraints = 1
model.optimize(mycallback)


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
optSolution.to_csv(r'result%s/opt_N-1_PSPS_Span%s_BC_%s_inst%s.csv'%(len(nodeArray),spanning,len(nodeArray),inst), index = False)#Check


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
resultTable.to_csv(r'result%s/result_%s_N-1_PSPS_Span%s_BC_%s.csv'%(len(nodeArray),machineName,spanning,len(nodeArray)), index = False)#Check


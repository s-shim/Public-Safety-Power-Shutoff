#from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket
import datetime
import time
import math
import multiprocessing as mp


def Rounding(zVal,inputPackage):
    augLineArray,augRiskFunction,sourceNode,targetNode,nodeArray,dBound,maxGen,demandRequired = inputPackage

    sortedAugLineArray = sorted(augLineArray, key = lambda l: zVal[l], reverse = True)
    
    # make a list of edges selected in the solution
    is_subNode = {}
    path_to_1 = {}
    
    addDegree = {}
    for i in nodeArray:
        addDegree[i] = 0
        is_subNode[i] = 0
        path_to_1[i] = 0
    path_to_1[1] = 1
    is_subNode[1] = 1
    
    is_subLine = {}
    is_solLine = {}
    for l in augLineArray:
        is_subLine[l] = 0
        is_solLine[l] = 0
    
    subG = nx.MultiGraph()
    subG.add_node(1)
    subKey = {}
    non_contingent = []
    connectedTo1 = {}
    disconnectedTo1 = []
    totalDBound = 0
    totalGenBound = 0
    theSubNodes = []
    theSubLines = []
    for l in sortedAugLineArray:
        if is_subLine[l] == 0:
            theSubLines += [l]
            is_subLine[l] = 1
            i = sourceNode[l]
            j = targetNode[l]
            was_subNode_i = is_subNode[i]
            was_subNode_j = is_subNode[j]
            is_subNode[i] = 1
            is_subNode[j] = 1
            subKey[l] = subG.add_edge(i,j)
            subG[i][j][subKey[l]]['origin'] = l
            
            if was_subNode_i == 0:
                theSubNodes += [i]
                disconnectedTo1 += [i]
                totalDBound += dBound[i]
                totalGenBound += maxGen[i]
                for l_i in augNeighbors[i]:
                    if l_i != l and sourceNode[l_i] == i:
                        j_i = targetNode[l_i]
                        if is_subNode[j_i] == 1:
                            theSubLines += [l_i]
                            is_subLine[l_i] = 1
                            subKey[l_i] = subG.add_edge(i,j_i)
                            subG[i][j_i][subKey[l_i]]['origin'] = l_i
                                
                    if l_i != l and targetNode[l_i] == i:
                        j_i = sourceNode[l_i]
                        if is_subNode[j_i] == 1:
                            theSubLines += [l_i]
                            is_subLine[l_i] = 1
                            subKey[l_i] = subG.add_edge(i,j_i)
                            subG[i][j_i][subKey[l_i]]['origin'] = l_i
                                
            if was_subNode_j == 0:
                theSubNodes += [j]
                disconnectedTo1 += [j]
                totalDBound += dBound[j]
                totalGenBound += maxGen[j]
                for l_j in augNeighbors[j]:
                    if l_j != l and sourceNode[l_j] == j:
                        i_j = targetNode[l_j]
                        if is_subNode[i_j] == 1:
                            theSubLines += [l_j]
                            is_subLine[l_j] = 1
                            subKey[l_j] = subG.add_edge(j,i_j)
                            subG[j][i_j][subKey[l_j]]['origin'] = l_j
                                
                    if l_j != l and targetNode[l_j] == j:
                        i_j = sourceNode[l_j]
                        if is_subNode[i_j] == 1:
                            theSubLines += [l_j]
                            is_subLine[l_j] = 1
                            subKey[l_j] = subG.add_edge(j,i_j)
                            subG[j][i_j][subKey[l_j]]['origin'] = l_j
                
            connectedTo1 = []
            for k in disconnectedTo1:
                if nx.has_path(subG,1,k) == True:
                    connectedTo1 += [k]
            for k in connectedTo1:
                disconnectedTo1.remove(k)
        
            non_contingent += [l]
            contingent = []
            for l_nc in non_contingent:
                tempsubG = copy.deepcopy(subG)
                tempsubG.remove_edge(sourceNode[l_nc],targetNode[l_nc],subKey[l_nc])
                if nx.has_path(tempsubG,sourceNode[l_nc],targetNode[l_nc]) == True:
                    contingent += [l_nc]
            for l_c in contingent:
                non_contingent.remove(l_c)
                
            if len(non_contingent) == 0 and len(disconnectedTo1) == 0:
                if demandRequired <= totalDBound and demandRequired <= totalGenBound:
                    break
    
    simpleSubG = nx.Graph(subG)
    for (i,j) in simpleSubG.edges:
        simpleSubG[i][j]['weight'] = -1
        
    for l in theSubLines:
        if simpleSubG[sourceNode[l]][targetNode[l]]['weight'] < zVal[l]:
            simpleSubG[sourceNode[l]][targetNode[l]]['weight'] = zVal[l]
            simpleSubG[sourceNode[l]][targetNode[l]]['origin'] = l
    
    T = nx.maximum_spanning_tree(simpleSubG,weight='weight')
    solG = nx.MultiGraph(T)
    
    notCovered = []
    covered = {}
    solRisk = 0
    for (i,j) in T.edges:
        l = simpleSubG[i][j]['origin']
        notCovered += [l]
        covered[l] = 0
        solRisk += augRiskFunction[l]
        is_solLine[l] = 1
        theSubLines.remove(l)
    
    sortedSubLines = sorted(theSubLines, key = lambda l: zVal[l], reverse = True)
    remainder = []
    lastLines = []
    fullyCover = False
    for l in sortedSubLines:
        i = sourceNode[l]
        j = targetNode[l]
        if (T.degree(i) + addDegree[i]) % 2 == 1 and (T.degree(j) + addDegree[j]) % 2 == 1:
            sp = nx.shortest_path(T,source=i,target=j)
            coverAnyNew = False
            for k in range(1,len(sp)):
                l_k = simpleSubG[sp[k-1]][sp[k]]['origin']
                if covered[l_k] == 0:
                    notCovered.remove(l_k)
                    covered[l_k] = 1
                    coverAnyNew = True
    
            if coverAnyNew == True:
                is_solLine[l] = 1
                solRisk += augRiskFunction[l]
                addDegree[i] += 1
                addDegree[j] += 1
            else:
                lastLines += [l]
        else:
            remainder += [l]
            
        if len(notCovered) == 0:
            fullyCover = True
            break
    
    if fullyCover == False:
        for l in remainder:
            u = sourceNode[l]
            v = targetNode[l]
    
            sp = nx.shortest_path(T,source=u,target=v)
            coverAnyNew = False
            for i in range(1,len(sp)):
                l_i = simpleSubG[sp[i-1]][sp[i]]['origin']
                if covered[l_i] == 0:
                    notCovered.remove(l_i)
                    covered[l_i] = 1
                    coverAnyNew = True
            if coverAnyNew == True:
                is_solLine[l] = 1
                solRisk += augRiskFunction[l]
                addDegree[u] += 1
                addDegree[v] += 1
    
            if len(notCovered) == 0:
                fullyCover = True
                break
            
    if fullyCover == False:
        for l in lastLines:
            u = sourceNode[l]
            v = targetNode[l]
    
            sp = nx.shortest_path(T,source=u,target=v)
            coverAnyNew = False
            for i in range(1,len(sp)):
                l_i = simpleSubG[sp[i-1]][sp[i]]['origin']
                if covered[l_i] == 0:
                    notCovered.remove(l_i)
                    covered[l_i] = 1
                    coverAnyNew = True
            if coverAnyNew == True:
                is_solLine[l] = 1
                solRisk += augRiskFunction[l]
                addDegree[u] += 1
                addDegree[v] += 1
    
            if len(notCovered) == 0:
                fullyCover = True
                break

    return solRisk, is_solLine, is_subNode, fullyCover



def ARR2(arg):
    iteration, timeLimit, bestSolutionPackage, inputPackage, augInputPackage, inst = arg
    return ARR_INC(iteration, timeLimit, bestSolutionPackage, inputPackage, augInputPackage, inst)

    
def ARR_INC(iteration, timeLimit, bestSolutionPackage, inputPackage, augInputPackage, inst):
    augLineArray,augRiskFunction, sourceNode,targetNode,nodeArray,dBound,maxGen,demandRequired = inputPackage
    augNeighbors, augG, augKey, augLineFunction = augInputPackage

    bestTime, bestTrial, bestSolRisk, best_is_solLine, best_is_solNode, bestMachine, bestIter = bestSolutionPackage
    trial = bestTrial
    tic = time.time()
    toc = time.time()
    timePast = bestTime
    
    initialSolRisk = bestSolRisk
        
    sizeArray = [len(nodeArray)]
    lineSizeArray = [len(lineArray)] ###
    augLineSizeArray = [len(augLineArray)] ###
    instArray = [inst]
    ratioDemandRequiredArray = [ratioDemandRequired] 
    totalDemandArray = [totalDemand]
    totalRiskArray = [totalRisk]
    resultRiskArray = [bestSolRisk]
    ratioRiskArray = [bestSolRisk/totalRisk]
    timeArray = [bestTime]
    nodecountArray = [bestTrial]
    stepArray = ['Incumbent']
    machineArray = [bestMachine]
    spanningArray = ['not spanning']
    iterArray = [bestIter]
    TLArray = [timeLimit]
    
    listHeuristic = list(zip(sizeArray,lineSizeArray,augLineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, stepArray, machineArray, spanningArray,iterArray,TLArray))
    nameColumn = ['Nodes','Lines','augLines','inst','totalDemand','ratioDemandRequired','totalRiskArray','minRisk','ratioRisk','Time','Trial', 'Step','Machine','Spanning','iter','timeLimit']
    heuristicTable = pd.DataFrame(listHeuristic,columns =nameColumn)
    heuristicTable.to_csv(r'parallel20250708/%s/inst%s/process/parallel_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False,machineName,iteration), index = False)#Check
    
# =============================================================================
#     varNameArray = ['risk']
#     varValArray = [bestSolRisk]        
#     varNameArray += ['time']
#     varValArray += [bestTime]        
#     varNameArray += ['trial']
#     varValArray += [bestTrial]        
#     for i in nodeArray:
#         varNameArray += ['X[%s]'%i]
#         varValArray += [best_is_solNode[i]]        
#     for l in augLineArray:
#         varNameArray += ['Z[%s]'%l]
#         varValArray += [best_is_solLine[l]]        
# 
#     intSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
#     intSolution.to_csv(r'relax_int/int_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,False,machineName,iteration), index = False)#Check
# =============================================================================

    print()
    print('bestTrial=',bestTrial)
    print('bestSolRisk =',bestSolRisk)
    
    half_solLine = {}
    seed_solLine = {}
    for l in augLineArray:
        half_solLine[l] = 0.5
        seed_solLine[l] = 0.5
        
    nLocal = 0
    reset = False
    while toc - tic < timeLimit:
        trial += 1
        rxVal = {}
        rzVal = {}
        reset = False
    
        fullyCover = False
        while fullyCover == False and toc - tic < timeLimit:
            for l in augLineArray:
                # rzVal[l] = seed_solLine[l] * random.random()
                rzVal[l] = seed_solLine[l] + (1 - seed_solLine[l]) * random.random()
                # rzVal[l] = random.random()
            
            solRisk, is_solLine, is_solNode, fullyCover = Rounding(rzVal,inputPackage)
            toc = time.time()
        
        if toc - tic < timeLimit:
            same = True
            RMSD = 0
            for l in augLineArray:
                RMSD = RMSD + (seed_solLine[l] - 0.5) ** 2
                if is_solLine[l] != best_is_solLine:
                    same = False
            RMSD = math.sqrt(RMSD / len(augLineArray))
                
            if same == True:
                nLocal += 1
                if random.random() < min(1,nLocal/20) * RMSD:
                    reset = True
                    seed_solLine = copy.deepcopy(half_solLine)
                    nLocal = 0
                    
                    info_sol = pd.read_csv('parallel20250708/%s/inst%s/info/info_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False))                    
                    [info_risk] = info_sol.loc[info_sol['varName']=='risk','varVal']
                    if bestSolRisk > info_risk:
                        inc_sol = pd.read_csv('parallel20250708/%s/inst%s/int/int_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False))                    
                        [inc_risk] = inc_sol.loc[inc_sol['varName']=='risk','varVal']

                        bestSolRisk = inc_risk
                        [bestTime] = inc_sol.loc[inc_sol['varName']=='time','varVal']
                        [bestTrial] = inc_sol.loc[inc_sol['varName']=='trial','varVal']

                        best_is_solLine = {}
                        best_is_solNode = {}
                        for varName in inc_sol['varName']:
                            if varName[0] == 'Z':
                                l = int(varName[2:-1])
                                [is_solLine_l] = inc_sol.loc[inc_sol['varName']==varName,'varVal']
                                is_solLine_l = int(is_solLine_l + 0.0001)
                                best_is_solLine[l] = is_solLine_l
                            if varName[0] == 'X':
                                i = int(varName[2:-1])
                                is_solNode_i = inc_sol.loc[inc_sol['varName']==varName,'varVal']
                                is_solNode_i = int(is_solNode_i + 0.0001)
                                best_is_solNode[i] = is_solNode_i
                            

            if same == False:
                nLocal = 0
        
                if bestSolRisk > solRisk:
                    bestTime = toc - tic + timePast        
                    bestTrial = trial
                    bestSolRisk = solRisk
                    best_is_solLine = copy.deepcopy(is_solLine)
                    best_is_solNode = copy.deepcopy(is_solNode)
                    print()
                    print(datetime.datetime.now())
                    print('bestTrial =',bestTrial)
                    print('bestTime =',bestTime)
                    print('bestSolRisk =',bestSolRisk)
        
                    sizeArray += [len(nodeArray)]
                    lineSizeArray += [len(lineArray)] 
                    augLineSizeArray += [len(augLineArray)] ###
                    instArray += [inst]
                    ratioDemandRequiredArray += [ratioDemandRequired] 
                    totalDemandArray += [totalDemand]
                    totalRiskArray += [totalRisk]
                    resultRiskArray += [bestSolRisk]
                    ratioRiskArray += [bestSolRisk/totalRisk]
                    timeArray += [bestTime]
                    nodecountArray += [bestTrial]
                    stepArray += ['Intermediate']
                    machineArray += [machineName]
                    spanningArray += ['not spanning']
                    iterArray += [iteration]
                    TLArray += [timeLimit]
        
                    listHeuristic = list(zip(sizeArray,lineSizeArray,augLineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, stepArray, machineArray, spanningArray,iterArray,TLArray))
                    nameColumn = ['Nodes','Lines','augLines','inst','totalDemand','ratioDemandRequired','totalRiskArray','minRisk','ratioRisk','Time','Trial', 'Step','Machine','Spanning','iter','timeLimit']
                    heuristicTable = pd.DataFrame(listHeuristic,columns = nameColumn)
                    heuristicTable.to_csv(r'parallel20250708/%s/inst%s/process/parallel_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False,machineName,iteration), index = False)#Check
                    
                    info_sol = pd.read_csv('parallel20250708/%s/inst%s/info/info_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False))                    
                    [info_risk] = info_sol.loc[info_sol['varName']=='risk','varVal']
                    if bestSolRisk < info_risk:                    
                        varNameArray = ['risk']
                        varValArray = [bestSolRisk]        
                        varNameArray += ['time']
                        varValArray += [bestTime]        
                        varNameArray += ['trial']
                        varValArray += [int(bestTrial + 0.0001)]        
                        varNameArray += [machineName]
                        varValArray += [-2]        
                        varNameArray += ['iter']
                        varValArray += [int(iteration + 0.0001)]        
                        infoSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
                        infoSolution.to_csv(r'parallel20250708/%s/inst%s/info/info_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False), index = False)#Check
                        for i in nodeArray:
                            varNameArray += ['X[%s]'%i]
                            varValArray += [int(best_is_solNode[i] + 0.0001)]        
                        for l in augLineArray:
                            varNameArray += ['Z[%s]'%l]
                            varValArray += [int(best_is_solLine[l] + 0.0001)]        
                    
                        intSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
                        intSolution.to_csv(r'parallel20250708/%s/inst%s/int/int_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False), index = False)#Check
                        
            if reset == False:
                alpha = 1 / (1 + math.exp(4 * RMSD))
                for l in augLineArray:
                    seed_solLine[l] = (1 - alpha) * seed_solLine[l] + alpha * best_is_solLine[l]
    
        toc = time.time()
    
    sizeArray += [len(nodeArray)]
    lineSizeArray += [len(lineArray)] 
    augLineSizeArray += [len(augLineArray)] ###
    instArray += [inst]
    ratioDemandRequiredArray += [ratioDemandRequired] 
    totalDemandArray += [totalDemand]
    totalRiskArray += [totalRisk]
    resultRiskArray += [bestSolRisk]
    ratioRiskArray += [bestSolRisk/totalRisk]
    timeArray += [bestTime]
    nodecountArray += [bestTrial]
    stepArray += ['Final']
    machineArray += [machineName]
    spanningArray += ['not spanning']
    iterArray += [iteration]
    TLArray += [timeLimit]
    
    listHeuristic = list(zip(sizeArray,lineSizeArray,augLineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, stepArray, machineArray, spanningArray,iterArray,TLArray))
    nameColumn = ['Nodes','Lines','augLines','inst','totalDemand','ratioDemandRequired','totalRiskArray','minRisk','ratioRisk','Time','Trial', 'Step','Machine','Spanning','iter','timeLimit']
    heuristicTable = pd.DataFrame(listHeuristic,columns = nameColumn)
    heuristicTable.to_csv(r'parallel20250708/%s/inst%s/process/parallel_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False,machineName,iteration), index = False)#Check
            
    return iteration, bestTime, bestTrial, bestSolRisk, best_is_solLine, best_is_solNode, initialSolRisk



### Code Starts Here
machineName = socket.gethostname()
spanning = False
(numNodes,folderName,timeLimit) = (1354,'matpower/case%spegase'%(1354), 3600 * 24 * 6)#[(1354,'matpower/case%spegase'%(1354), 3600 * 24 * 6),(118,'matpower/case%s'%(118),60),(300,'matpower/case%s'%(300),900)]:

print(machineName)
print(datetime.datetime.now())

nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))


sizeArray = []
lineSizeArray = [] 
initialSolArray = []
instArray = []
ratioDemandRequiredArray = [] 
totalDemandArray = []
totalRiskArray = []
resultRiskArray = []
ratioRiskArray = []
timeArray = []
nodecountArray = []
stepArray = [] 
machineArray = []
spanningArray = []
iterArray = []
TLArray = []


for inst in [2]:
        
    lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
            
    dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = md.profile(lines,nodes)
    G, key, lineFunction = md.multigraph(nodes, lines, lineArray)
    simpleG = nx.Graph(G)
    
    ratioDemandRequired = 0.9
    demandRequired = totalDemand * ratioDemandRequired
    
    augG = nx.MultiGraph(G)
    augKey = copy.deepcopy(key)
    augRiskFunction = copy.deepcopy(riskFunction)
    augLineArray = copy.deepcopy(lineArray)
    augLineFunction = copy.deepcopy(lineFunction)
    for l in lineArray:
        if contingentFunction[l] == 0:
            augKey[-l] = augG.add_edge(sourceNode[l],targetNode[l])
            augRiskFunction[-l] = 0
            augLineArray += [-l]
            sourceNode[-l] = sourceNode[l]
            targetNode[-l] = targetNode[l]
            switchableFunction[-l] = 1
            augLineFunction[sourceNode[-l],targetNode[-l],augKey[-l]] = -l
    print(len(G.edges))
    print(len(augG.edges))
    print(datetime.datetime.now())
    print('nodes=',len(G.nodes))
    print('edges=',len(G.edges))
    
    
    augNeighbors = {}
    for i in nodeArray:
        augNeighbors[i] = []    
    for l in augLineArray:
        augNeighbors[sourceNode[l]] += [l]
        augNeighbors[targetNode[l]] += [l]


    inputPackage = augLineArray, augRiskFunction, sourceNode,targetNode,nodeArray,dBound,maxGen,demandRequired
    augInputPackage = augNeighbors, augG, augKey, augLineFunction  



# =============================================================================
#     # construct initial incumbent if no incumbent    
#     varNameArray = ['risk']
#     varValArray = [totalRisk]        
#     varNameArray += ['time']
#     varValArray += [0]        
#     varNameArray += ['trial']
#     varValArray += [0]        
#     varNameArray += [machineName]
#     varValArray += [-2]        
#     varNameArray += ['iter']
#     varValArray += [-1]        
#     for i in nodeArray:
#         varNameArray += ['X[%s]'%i]
#         varValArray += [1]
#     for l in augLineArray:
#         varNameArray += ['Z[%s]'%l]
#         varValArray += [1]
#     intSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
#     intSolution.to_csv(r'parallel20250708/%s/inst%s/int/int_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False), index = False)#Check
# =============================================================================

    
    inc_sol = pd.read_csv('parallel20250708/%s/inst%s/int/int_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False))                    

    best_is_solLine = {}
    best_is_solNode = {}
    
    for i in nodeArray:
        [x_i] = inc_sol.loc[inc_sol['varName']=='X[%s]'%i,'varVal']
        best_is_solNode[i] = int(x_i)
    
    inc_value = 0
    for l in augLineArray:
        [z_l] = inc_sol.loc[inc_sol['varName']=='Z[%s]'%l,'varVal']
        z_l = int(z_l)
        best_is_solLine[l] = z_l
        inc_value += augRiskFunction[l] * z_l

    [bestTime] = inc_sol.loc[inc_sol['varName']=='time','varVal']
    [bestTrial] = inc_sol.loc[inc_sol['varName']=='trial','varVal']
    [bestSolRisk] = inc_sol.loc[inc_sol['varName']=='risk','varVal']

    print(bestSolRisk,'=',inc_value)

    # construct initial sharing info    
    varNameArray = ['risk']
    varValArray = [inc_value]        
    varNameArray += ['time']
    varValArray += [bestTime]        
    varNameArray += ['trial']
    varValArray += [bestTrial]    
    [bestMachine] = inc_sol.loc[inc_sol['varVal']==-2,'varName']    
    varNameArray += [bestMachine]
    varValArray += [-2]        
    [bestIter] = inc_sol.loc[inc_sol['varName']=='iter','varVal'] 
    bestIter = int(bestIter + 0.0001)
    varNameArray += ['iter']
    varValArray += [bestIter]        
    infoSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
    infoSolution.to_csv(r'parallel20250708/%s/inst%s/info/info_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,len(nodeArray),inst,False), index = False)#Check

        
    bestSolutionPackage = bestTime, bestTrial, bestSolRisk, best_is_solLine, best_is_solNode, bestMachine, bestIter
    timePast = bestTime
    
    
    if __name__ == '__main__':
        numCores = mp.cpu_count()
        p = mp.Pool(numCores)
    
        multiArgs = []  
        for iteration in range(numCores):
            multiArgs += [(iteration, timeLimit, bestSolutionPackage, inputPackage, augInputPackage, inst)]  
    
        results = p.map(ARR2, multiArgs)
    
        grandRisk = totalRisk + 1    
        grandTime = timeLimit * 10 + timePast
        
        for bestPackage in results:
            iteration, bestTime, bestTrial, bestSolRisk, best_is_solLine, best_is_solNode, initialSolRisk = bestPackage
    
            if grandRisk == bestSolRisk and grandTime > bestTime:
                grandRisk = bestSolRisk
                grand_is_solLine = copy.deepcopy(best_is_solLine)
                grand_is_solNode = copy.deepcopy(best_is_solNode)
                grandTrial = bestTrial            
                grandTime = bestTime
                grandIteration = iteration
            
            if grandRisk > bestSolRisk:
                grandRisk = bestSolRisk
                grand_is_solLine = copy.deepcopy(best_is_solLine)
                grand_is_solNode = copy.deepcopy(best_is_solNode)
                grandTrial = bestTrial            
                grandTime = bestTime
                grandIteration = iteration
    
        print()
        print(grandRisk,grandTrial,grandTime,grandIteration,timeLimit)
        print(datetime.datetime.now())
        print()
    
        sizeArray += [len(nodeArray)]
        lineSizeArray += [len(lineArray)] 
        initialSolArray += [initialSolRisk] ###
        instArray += [inst]
        ratioDemandRequiredArray += [ratioDemandRequired] 
        totalDemandArray += [totalDemand]
        totalRiskArray += [totalRisk]
        resultRiskArray += [grandRisk]
        ratioRiskArray += [grandRisk/totalRisk]
        timeArray += [grandTime]
        nodecountArray += [grandTrial]
        stepArray += ['ARR'] 
        machineArray += [machineName]
        spanningArray += ['not spanning']
        iterArray += [grandIteration]
        TLArray += [timeLimit]
        
        listGrand = list(zip(sizeArray,lineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, initialSolArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, stepArray, machineArray, spanningArray,iterArray,TLArray))
        nameColumn = ['Nodes','Lines','inst','totalDemand','ratioDemandRequired','totalRiskArray','initialRisk','finalRisk','ratioRisk','Time','Trial', 'Method','Machine','Spanning','iter','timeLimit']
        grandTable = pd.DataFrame(listGrand,columns = nameColumn)
        grandTable.to_csv(r'parallel20250708/%s/inst%s/summary/summary_size%s_span%s_machine%s.csv'%(len(nodeArray),inst,len(nodeArray),False,machineName), index = False)#Check
    
        varNameArray = ['risk']
        varValArray = [grandRisk]        
        varNameArray += ['time']
        varValArray += [grandTime]        
        varNameArray += ['trial']
        varValArray += [grandTrial]        
        for i in nodeArray:
            varNameArray += ['X[%s]'%i]
            varValArray += [grand_is_solNode[i]]        
        for l in augLineArray:
            varNameArray += ['Z[%s]'%l]
            varValArray += [grand_is_solLine[l]]        
    
        grandSolution = pd.DataFrame(list(zip(varNameArray,varValArray)),columns = ['varName','varVal'])
        grandTable.to_csv(r'parallel20250708/%s/inst%s/summary/summary_size%s_span%s_machine%s.csv'%(len(nodeArray),inst,len(nodeArray),False,machineName), index = False)#Check
        
        
        
        
        
        
        
        
        
        
        
        
        
            
        
        
        
    
    

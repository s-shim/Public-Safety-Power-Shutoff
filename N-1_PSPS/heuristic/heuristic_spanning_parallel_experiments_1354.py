#from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
#import myDictionary as md
import socket
import datetime
import time
import math
import multiprocessing as mp


def ARR2(arg):
    lineArray,nodeArray,inst,G,timeLimit,profiling,machineName,iteration = arg
    return ARR(lineArray,nodeArray,inst,G,timeLimit,profiling,machineName,iteration)  
    
    
def ARR(lineArray,nodeArray,inst,G,timeLimit,profiling,machineName,iteration):
    dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = profiling
    demandRequired = totalDemand * ratioDemandRequired

    halfZ = {}
    relZ = {}
    edgeWeight = {}
    for l in lineArray:
        halfZ[l] = 0.5
        relZ[l] = 0.5
        edgeWeight[l] = random.random() * relZ[l]
        
    
    # initial solution    
    tic = time.time()
    solEdges, solFunction, solRisk = rounding(edgeWeight,lineArray,G,riskFunction)
    
    trial = 0
    bestTrial = trial 
    bestSolEdges = copy.deepcopy(solEdges)
    bestSolFunction = copy.deepcopy(solFunction)
    bestSolRisk = solRisk   
    print(trial,bestSolRisk,0)
        
    toc = time.time()
    bestTime = toc - tic
    # timeLimit = bestTime * 20
    nLocal = 0
    
    sizeArray = [len(nodeArray)]
    lineSizeArray = [len(lineArray)] ###
    instArray = [inst]
    ratioDemandRequiredArray = [ratioDemandRequired] 
    totalDemandArray = [totalDemand]
    totalRiskArray = [totalRisk]
    resultRiskArray = [bestSolRisk]
    ratioRiskArray = [bestSolRisk/totalRisk]
    timeArray = [bestTime]
    nodecountArray = [bestTrial]
    #securityArray = [security]
    machineArray = [machineName]
    spanningArray = ['spanning']
    iterArray = [iteration] ###
    TLArray = [timeLimit] ###
    stepArray = ['Initial']
    
    while toc - tic < timeLimit:
        reset = False
        same = False
        trial += 1
        edgeWeight = {}
        RMSD = 0
        for l in lineArray:
            edgeWeight[l] = random.random() * relZ[l]
            RMSD += (relZ[l] - 0.5) ** 2
        RMSD = RMSD / len(lineArray)
        RMSD = math.sqrt(RMSD)
    
        solEdges, solFunction, solRisk = rounding(edgeWeight,lineArray,G,riskFunction)
    
        if bestSolRisk > solRisk:
            bestSolEdges = copy.deepcopy(solEdges)
            bestSolFunction = copy.deepcopy(solFunction)
            bestSolRisk = solRisk   
            bestTrial = trial
            nLocal = 0
            toc = time.time()
            bestTime = toc - tic
            print(trial,bestSolRisk,toc - tic)
            # timeLimit = max(timeLimit,bestTime * 2)
    
            sizeArray += [len(nodeArray)]
            lineSizeArray += [len(lineArray)] ###
            instArray += [inst]
            ratioDemandRequiredArray += [ratioDemandRequired] 
            totalDemandArray += [totalDemand]
            totalRiskArray += [totalRisk]
            resultRiskArray += [bestSolRisk]
            ratioRiskArray += [bestSolRisk/totalRisk]
            timeArray += [bestTime]
            nodecountArray += [bestTrial]
            #securityArray += [security]
            machineArray += [machineName]
            spanningArray += ['spanning']
            iterArray += [iteration] ###
            TLArray += [timeLimit] ###
            stepArray += ['Intermediate']
    
            listHeuristic = list(zip(sizeArray,lineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, stepArray, timeArray, nodecountArray, machineArray, spanningArray,iterArray,TLArray))
            nameColumn = ['Nodes','Lines','inst','totalDemand','ratioDemandRequired','totalRiskArray','minRisk','ratioRisk','Step','Time','Trial','Machine','Spanning','iter','timeLimit']
            heuristicTable = pd.DataFrame(listHeuristic,columns =nameColumn)
            heuristicTable.to_csv(r'result_heuristic/result_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,True,machineName,iteration), index = False)#Check
            
        if bestSolRisk == solRisk:
            same = True
            for l in lineArray:
                if bestSolFunction[l] != solFunction[l]:
                    same = False
                    break
    
            if same == True:
                nLocal += 1
                if RMSD * min(nLocal/20,1) > random.random():
                    relZ = copy.deepcopy(halfZ)  
                    nLocal = 0
                    reset = True
                
            if same == False:
                nLocal = 0
                
        if reset == False:
            alpha = 1 / (1 + math.exp(4 * RMSD))
            for l in lineArray:
                relZ[l] = (1 - alpha) * relZ[l]
                
            for l in bestSolEdges:
                relZ[l] += alpha
                
        toc = time.time()

    sizeArray += [len(nodeArray)]
    lineSizeArray += [len(lineArray)] ###
    instArray += [inst]
    ratioDemandRequiredArray += [ratioDemandRequired] 
    totalDemandArray += [totalDemand]
    totalRiskArray += [totalRisk]
    resultRiskArray += [bestSolRisk]
    ratioRiskArray += [bestSolRisk/totalRisk]
    timeArray += [toc - tic]
    nodecountArray += [bestTrial]
    #securityArray += [security]
    machineArray += [machineName]
    spanningArray += ['spanning']
    iterArray += [iteration] ###
    TLArray += [timeLimit] ###
    stepArray += ['Final']

    listHeuristic = list(zip(sizeArray,lineSizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, stepArray, timeArray, nodecountArray, machineArray, spanningArray,iterArray,TLArray))
    nameColumn = ['Nodes','Lines','inst','totalDemand','ratioDemandRequired','totalRiskArray','minRisk','ratioRisk','Step','Time','Trial','Machine','Spanning','iter','timeLimit']
    heuristicTable = pd.DataFrame(listHeuristic,columns =nameColumn)
    heuristicTable.to_csv(r'result_heuristic/result_size%s_inst%s_span%s_machine%s_iter%s.csv'%(len(nodeArray),inst,True,machineName,iteration), index = False)#Check
        
    return bestSolEdges, bestSolFunction, bestSolRisk, bestTrial, bestTime, iteration, timeLimit


def multigraph(nodes,lines,lineArray):
    key = {}
    G = nx.MultiGraph()
    lineFunction = {}
    for l in lineArray:
        [source_l] = lines.loc[lines['Line']==l,'Source']
        [target_l] = lines.loc[lines['Line']==l,'Target']
        key[l] = G.add_edge(source_l,target_l)
        lineFunction[source_l,target_l,key[l]] = l
        lineFunction[target_l,source_l,key[l]] = l
    return G, key, lineFunction


def profile(lines,nodes):
    totalDemand = 0.0
    totalGeneration = 0.0
    nodeArray = []
    dBound = {}
    minGen = {}
    maxGen = {}
    for i in nodes['Node']:
        nodeArray += [i]
        [load_i] = nodes.loc[nodes['Node']==i,'Load']
        totalDemand += load_i
        dBound[i] = load_i
        [min_gen_i] = nodes.loc[nodes['Node']==i,'min_generation']
        [max_gen_i] = nodes.loc[nodes['Node']==i,'max_generation']
        minGen[i] = min_gen_i
        maxGen[i] = max_gen_i
        totalGeneration += max_gen_i
        
    totalDemand = min(totalDemand, totalGeneration)
        
    lineArray = []
    theSwitchable = []
    theContingent = []
    theMonitored = []
    totalRisk = 0
    riskFunction = {}
    sourceNode = {}
    targetNode = {}
    switchableFunction = {}
    contingentFunction = {}
    monitoredFunction = {}
    for l in lines['Line']:
        lineArray += [l]
        [source_l] = lines.loc[lines['Line']==l,'Source']
        [target_l] = lines.loc[lines['Line']==l,'Target']
        sourceNode[l] = source_l
        targetNode[l] = target_l
        [switch_l] = lines.loc[lines['Line']==l,'Switchable']
        switchableFunction[l] = switch_l
        [contingency_l] = lines.loc[lines['Line']==l,'Vulnerable?']
        contingentFunction[l] = contingency_l
        [capacity_l] = lines.loc[lines['Line']==l,'Normal Flow Limit (MW)']
        monitoredFunction[l] = 0
        [risk_l] = lines.loc[lines['Line']==l,'Risk']
        riskFunction[l] = risk_l
    
        totalRisk += risk_l    
        
        if switch_l == 1:
            theSwitchable += [l]
    
        if contingency_l == 1:
            theContingent += [l]
            
        if capacity_l < 99999:
            theMonitored += [l]
            monitoredFunction[l] = 1
            
    return dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction


def rounding(edgeWeight,lineArray,G,riskFunction):
    sortedLineArray = sorted(lineArray, key = lambda l: edgeWeight[l], reverse = True)
    
    solFunction = {}
    simpleG = nx.Graph(G)
    tree = {}
    for l in lineArray:
        tree[l] = 0
        simpleG[sourceNode[l]][targetNode[l]]['weight'] =  - 1
        solFunction[l] = 0
    for l in lineArray:
        if simpleG[sourceNode[l]][targetNode[l]]['weight'] < edgeWeight[l]:
            simpleG[sourceNode[l]][targetNode[l]]['weight'] = edgeWeight[l]
            simpleG[sourceNode[l]][targetNode[l]]['origin'] = l
    
    T = nx.maximum_spanning_tree(simpleG)
    
    solRisk = 0
    covered = {}
    notCovered = []
    treeEdges = []
    solEdges = []
    for (u,v) in T.edges():
        l = simpleG[u][v]['origin']
        tree[l] = 1
        if contingentFunction[l] == 1:
            covered[l] = 0
            notCovered += [l]
        if contingentFunction[l] == 0:
            covered[l] = 1
        treeEdges += [l]
        solEdges += [l]
        solRisk += riskFunction[l]
        solFunction[l] = 1
        T[u][v]['origin'] = l
    
    coTree = []
    coTreeRemaining = []
    for l in sortedLineArray:
        u = sourceNode[l]
        v = targetNode[l]
        k = key[l]
    
        if len(notCovered) == 0:
            break        
            
        if tree[l] == 0:
            if T.degree(u) % 2 == 1 and T.degree(v) % 2 == 1: 
                coverAnyNew = False
                sp = nx.shortest_path(T,source=sourceNode[l],target=targetNode[l])
                for i in range(1,len(sp)):
                    l_i = T[sp[i-1]][sp[i]]['origin']
                    if covered[l_i] == 0:
                        notCovered.remove(l_i)
                        covered[l_i] = 1
                        coverAnyNew = True
        
                if coverAnyNew == True:
                    coTree += [l]
                    solEdges += [l]
                    solRisk += riskFunction[l]
                    solFunction[l] = 1
            else:
                coTreeRemaining += [l]
    
    for l in coTreeRemaining:
        u = sourceNode[l]
        v = targetNode[l]
        k = key[l]
    
        if len(notCovered) == 0:
            break                
        else:
            coverAnyNew = False
            sp = nx.shortest_path(T,source=sourceNode[l],target=targetNode[l])
            for i in range(1,len(sp)):
                l_i = T[sp[i-1]][sp[i]]['origin']
                if covered[l_i] == 0:
                    notCovered.remove(l_i)
                    covered[l_i] = 1
                    coverAnyNew = True
    
            if coverAnyNew == True:
                coTree += [l]
                solEdges += [l]
                solRisk += riskFunction[l]
                solFunction[l] = 1
    
    return solEdges, solFunction, solRisk    


### code starts here

# machine, date
machineName = socket.gethostname()
print(machineName)
print(datetime.datetime.now())

# network
#(numNodes,folderName) = (30,'matpower/case30')
#(numNodes,folderName) = (118,'matpower/case118')
#(numNodes,folderName) = (300,'matpower/case300')
#(numNodes,folderName) = (1354,'matpower/case1354pegase')

#timeLimit = 30 # seconds

for (numNodes,folderName,timeLimit) in [(1354,'matpower/case1354pegase',7200)]:#[(30,'matpower/case%s'%(30),30),(118,'matpower/case%s'%(118),100),(300,'matpower/case%s'%(300),300)]:#[(30,'matpower/case%s'%(30),30),(118,'matpower/case%s'%(118),100),(300,'matpower/case%s'%(300),300),(1354,'matpower/case1354pegase',1000)]:


    machineArray = []
    methArray = []
    sizeArray = []
    lineSizeArray = []
    instArray = []
    ratioDemandRequiredArray = [] 
    totalDemandArray = []
    totalRiskArray = []
    resultRiskArray = []
    ratioRiskArray = []
    solEdgeArray = []    
    trialArray = []
    timeArray = []
    spanningArray = []
    iterArray = []
    TLArray = []    
    
    
    # instance
    inst = 1
    #iteration = 0
    ratioDemandRequired = 0.9
    slackNode = 1 # root node
    
    for inst in range(1,10+1):
    
        # input data
        lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
        nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))
        
        
        # profile data
        profiling = profile(lines,nodes)        
        dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = profiling
        demandRequired = totalDemand * ratioDemandRequired
        G, key, lineFunction = multigraph(nodes,lines,lineArray)  
        print(datetime.datetime.now())
        print('timeLimit =',timeLimit)
        
        if __name__ == '__main__':
            numCores = mp.cpu_count()
            p = mp.Pool(numCores)
        
            multiArgs = []  
            for iteration in range(numCores):
                multiArgs += [(lineArray,nodeArray,inst,G,timeLimit,profiling,machineName,iteration)]  
        
            results = p.map(ARR2, multiArgs)
        
            grandRisk = totalRisk + 1    
            grandTime = timeLimit * 10
            
            for bestPackage in results:
                bestSolEdges, bestSolFunction, bestSolRisk, bestTrial, bestTime, bestIteration, bestTimeLimit = bestPackage
        
                if grandRisk == bestSolRisk and grandTime > bestTime:
                    grandRisk = bestSolRisk
                    grandSolEdges = copy.deepcopy(bestSolEdges)
                    grandSolFunction = copy.deepcopy(bestSolFunction)
                    grandTrial = bestTrial            
                    grandTime = bestTime
                    grandIteration = bestIteration
                    grandTimeLimit = bestTimeLimit
                
                if grandRisk > bestSolRisk:
                    grandRisk = bestSolRisk
                    grandSolEdges = copy.deepcopy(bestSolEdges)
                    grandSolFunction = copy.deepcopy(bestSolFunction)
                    grandTrial = bestTrial            
                    grandTime = bestTime
                    grandIteration = bestIteration
                    grandTimeLimit = bestTimeLimit
        
            print()
            print(grandRisk,len(grandSolEdges),grandTrial,grandTime,grandIteration,grandTimeLimit)
            print(datetime.datetime.now())
        
        
            machineArray += [machineName]
            methArray += ['ARR']
            sizeArray += [len(nodeArray)]
            lineSizeArray += [len(lineArray)]
            instArray += [inst]
            ratioDemandRequiredArray += [ratioDemandRequired] 
            totalDemandArray += [totalDemand]
            totalRiskArray += [totalRisk]
            resultRiskArray += [grandRisk]
            ratioRiskArray += [grandRisk/totalRisk]
            solEdgeArray += [len(grandSolEdges)]    
            trialArray += [grandTrial]
            timeArray += [grandTime]
            spanningArray += ['spanning']
            iterArray += [grandIteration]
            TLArray += [grandTimeLimit]    
        
            listGrand = list(zip(machineArray,methArray,sizeArray,lineSizeArray,instArray,ratioDemandRequiredArray,totalDemandArray,totalRiskArray,resultRiskArray,ratioRiskArray,solEdgeArray,trialArray,timeArray,spanningArray,iterArray,TLArray))
            nameColumn = ['Machine','Method','Nodes','Lines','inst','ratioDemandRequired','totalDemand','totalRiskArray','minRisk','ratioRisk','solEdges','Trial','Time','Spanning','iter','timeLimit']
            grandTable = pd.DataFrame(listGrand,columns =nameColumn)
            grandTable.to_csv(r'summary/summary_size%s_span%s_machine%s.csv'%(len(nodeArray),True,machineName), index = False)#Check
                
            
            lArray = []
            solArray = []
            for l in lineArray:
                lArray += [l]
                solArray += [grandSolFunction[l]]
            loptTable = pd.DataFrame(list(zip(lArray,solArray)),columns = ['Line','Reference'])
            loptTable.to_csv(r'lopt/lopt_size%s_inst%s_span%s.csv'%(len(nodeArray),inst,True), index = False)#Check
                
                
            
            
            
            
            
            
            
            
            
            
            
        
            
            
        
        
    









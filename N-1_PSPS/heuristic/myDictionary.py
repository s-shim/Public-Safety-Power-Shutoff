import pandas as pd
import random
import copy
import networkx as nx


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
       

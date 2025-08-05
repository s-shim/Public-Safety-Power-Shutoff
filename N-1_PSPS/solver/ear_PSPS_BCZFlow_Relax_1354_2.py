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

for (numNodes,folderName) in [(1354,'matpower/case%spegase'%(1354))]:#[(30,'matpower/case%s'%(30)),(118,'matpower/case%s'%(118)),(300,'matpower/case%s'%(300))]:
    for inst in [2]:#range(1,50+1):
        print(machineName)
        print(datetime.datetime.now())
        
        lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
        nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))
                
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
        
        
        augNeighbors = {}
        for i in nodeArray:
            augNeighbors[i] = []    
        for l in augLineArray:
            augNeighbors[sourceNode[l]] += [l]
            augNeighbors[targetNode[l]] += [l]
            
        
        model = Model('N-1 PSPS')
        #model.setParam('LogFile', 'result%s/grblog/grblog_%s_N-1_PSPS_Span%s_BCZ_Enhanced_%s_inst%s.txt'%(len(nodeArray),machineName,spanning,len(nodeArray),inst))    
        #model.setParam('TimeLimit', TL)  # Set time limit to 60 seconds
        
        model.setParam('MIPFocus',3)   
        model.setParam('Cuts',3)    
        
        x_vars = []
        x_names = []
        for i in nodeArray:
            x_vars += [(i)]
            x_names += ['X[%s]'%i]
        # X = model.addVars(x_vars, vtype = GRB.BINARY, name = x_names)
        X = model.addVars(x_vars, vtype = GRB.CONTINUOUS, name = x_names)
        for i in nodeArray:
            LHS = [(1,X[i])]
            model.addConstr(LinExpr(LHS)<=1, name='boudnX[%s]'%i)
        
        dummy_vars = [(1)]
        dummy_names = ['dummy[%s]'%1]
        dummy = model.addVars(dummy_vars, vtype = GRB.BINARY, name = dummy_names)
        LHS = [(1,dummy[1]),(-1,X[1])]
        model.addConstr(LinExpr(LHS)==0, name='dummy[1] = X[1]')
        
            
        z_vars = []
        z_names = []
        for l in augLineArray:
            z_vars += [(l)]
            z_names += ['Z[%s]'%l]
        # Z = model.addVars(z_vars, vtype = GRB.BINARY, name = z_names)
        Z = model.addVars(z_vars, vtype = GRB.CONTINUOUS, name = z_names)
        for l in augLineArray:
            LHS = [(1,Z[l])]
            model.addConstr(LinExpr(LHS)<=1, name='boudnZ[%s]'%l)
        
        
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
        
        
        # add constraints
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
        model.addConstr(LinExpr(LHS)==1, name='reference node') # if slackNode = 1
        
        for l in augLineArray:
            LHS_source = [(1,Z[l]),(-1,X[sourceNode[l]])]
            LHS_target = [(1,Z[l]),(-1,X[targetNode[l]])]
        
            model.addConstr(LinExpr(LHS_source)<=0, name='source bound (%s)'%l)
            model.addConstr(LinExpr(LHS_target)<=0, name='target bound (%s)'%l)
        
            if switchableFunction[l] == 0:
                print('non-switchable =',l,sourceNode[l],targetNode[l])
                LHS_edge1 = [(-1,Z[l]),(1,X[sourceNode[l]])]
                LHS_edge2 = [(-1,Z[l]),(1,X[targetNode[l]])]
                model.addConstr(LinExpr(LHS_edge1)==0, name='edge inequality1 (%s)'%l)
                model.addConstr(LinExpr(LHS_edge2)==0, name='edge inequality2 (%s)'%l)
        
        
        # flow variables and constraints to connect the nodes to 1
        f_vars = []
        f_names = []
        for l in augLineArray:
            f_vars += [(l,1)]
            f_names += ['F[%s,%s]'%(l,1)]
            f_vars += [(l,-1)]
            f_names += ['F[%s,%s]'%(l,-1)]
        F = model.addVars(f_vars, vtype = GRB.CONTINUOUS, name = f_names)
        
        LHS1 = []
        LHS1_in = []
        for l in augLineArray:
            if sourceNode[l] == 1:
                LHS1 += [(1,F[l,1])]
                LHS1_in += [(1,F[l,-1])]
            if targetNode[l] == 1:
                LHS1 += [(1,F[l,-1])]
                LHS1_in += [(1,F[l,1])]
        for i in nodeArray:
            if i != 1:
                LHS1 += [(-1,X[i])]
        model.addConstr(LinExpr(LHS1)==0, name='ourflow from 1')
        model.addConstr(LinExpr(LHS1_in)==0, name='inflow into 1')
        
        for i in nodeArray:
            if i != 1:
                LHS = [(-1,X[i])]
                for l in augLineArray:
                    if sourceNode[l] == i:
                        LHS += [(-1,F[l,1]),(1,F[l,-1])]
                    if targetNode[l] == i:
                        LHS += [(1,F[l,1]),(-1,F[l,-1])]
                model.addConstr(LinExpr(LHS)==0, name='balance at %s'%i)
                
        for l in augLineArray:
            LHS = [(-(len(G.nodes)-1),Z[l])]
            LHS += [(1,F[l,1]),(1,F[l,-1])]
            model.addConstr(LinExpr(LHS)<=0, name='flow capacity at %s'%l)
        
        
        model._bestTotalRisk = totalRisk
        model._bestRelTotalRisk = 0
        def security(model, where):
            if where == GRB.Callback.MIPSOL:
        
                # make a list of edges selected in the solution
                xVal = {}
                zVal = {}
                for i in nodeArray:
                    xVal[i] = model.cbGetSolution(X[i])
                    
                for l in augLineArray:
                    simpleG[sourceNode[l]][targetNode[l]]['capacity'] = 0
        
                for l in augLineArray:
                    zVal[l] = model.cbGetSolution(Z[l])
                    simpleG[sourceNode[l]][targetNode[l]]['capacity'] += zVal[l]
        
                for i in nodeArray:
                    for j in nodeArray:
                        if i < j:                    
                            cut_value, partition = nx.minimum_cut(simpleG, i, j, capacity = 'capacity')
                            reachable, non_reachable = partition
                            
                            reachableNode = {}
                            for k in reachable:
                                reachableNode[k] = 1
                            for k in non_reachable:
                                reachableNode[k] = -1
                                
                            if cut_value < 2 * xVal[i] + 2 * xVal[j] - 2 - 0.0001:
                                LHS = [(-2,X[i]),(-2,X[j])]
                                for l in augLineArray:
                                    if reachableNode[sourceNode[l]] * reachableNode[targetNode[l]] == -1:
                                        LHS += [(1,Z[l])]
                                model.cbLazy(LinExpr(LHS)>=-2)            
        
        
        
        # objective function
        objTerms = []
        for l in lineArray:
            objTerms += [(riskFunction[l],Z[l])]
        
        # optimize                                
        model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
        model.update()
        model.Params.lazyConstraints = 1
        model.optimize(security)
        
        
        varName = ['relax']
        varVal = [model.getObjective().getValue()]
        varName += ['time']
        varVal += [model.Runtime]
        for v in model.getVars():
            varName += [v.varname]
            varVal += [v.x]
        relaxSolution = pd.DataFrame(list(zip(varName, varVal)),columns =['varName', 'varVal'])
        relaxSolution.to_csv(r'relax/relax_%s_inst%s.csv'%(len(G.nodes),inst), index = False)#Check
        
        
        
        
        
        

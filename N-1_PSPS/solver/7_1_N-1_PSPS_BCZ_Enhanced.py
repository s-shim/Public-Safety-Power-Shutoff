from gurobipy import *
import pandas as pd
import random
import copy
import networkx as nx
import myDictionary as md
import socket
import datetime

machineName = socket.gethostname()
for spanning in [False]:#[True, False]:
    for (numNodes,folderName,TL) in [(118,'matpower/case%s'%(118),600)]:#[(118,'matpower/case%s'%(118),60),(300,'matpower/case%s'%(300),900)]:
        print(machineName)
        print(datetime.datetime.now())
        
        #numNodes = 1354
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
        
        
        for inst in range(1,10+1):#range(1,50+1):
            
            
            lines = pd.read_csv('../../%s/lines_%s_inst%s.csv'%(folderName,numNodes,inst))
            nodes = pd.read_csv('../../%s/nodes_%s.csv'%(folderName,numNodes))
                    
            dBound, minGen, maxGen, nodeArray, riskFunction, sourceNode, targetNode, lineArray, theSwitchable, theContingent,theMonitored, totalDemand, totalRisk, switchableFunction, contingentFunction, monitoredFunction = md.profile(lines,nodes)
            G, key, lineFunction = md.multigraph(nodes, lines, lineArray)
            simpleG = nx.Graph(G)
            print('total risk =',totalRisk)
            
            ratioDemandRequired = 0.9
            demandRequired = totalDemand * ratioDemandRequired
            
                    
            model = Model('N-1 PSPS')
            model.setParam('LogFile', 'result%s/grblog/grblog_%s_N-1_PSPS_Span%s_BCZ_Enhanced_%s_inst%s.txt'%(len(nodeArray),machineName,spanning,len(nodeArray),inst))    
            model.setParam('TimeLimit', TL)  # Set time limit to 60 seconds

            model.setParam('MIPFocus',3)   
            model.setParam('Cuts',3)    

# =============================================================================
#             model.setParam('Cuts',0)    
#             model.setParam('GomoryPasses',150)     
#             model.setParam('MIRCuts',2)    
#             model.setParam('FlowCoverCuts',2)    
#             model.setParam('ZeroHalfCuts',2)    
# =============================================================================


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
            
            for l in lineArray:
                Z[l].Start = 1
            
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
            
            
            def security(model, where):
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
                            if cut_value < xVal[i] - 0.0001:
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
                        if cut_valueC < zVal[c] - 0.0001:
                            LHSC = [(-1,Z[c])]
                            for l in lineArray:
                                if l != c:
                                    if sourceNode[l] in reachableC and targetNode[l] in non_reachableC:
                                        LHSC += [(1,Z[l])]
                                    if sourceNode[l] in non_reachableC and targetNode[l] in reachableC:
                                        LHSC += [(1,Z[l])]
                            model.cbLazy(LinExpr(LHSC)>=0)     
                            

                if where == GRB.callback.MIPNODE:
                    status = model.cbGet(GRB.callback.MIPNODE_STATUS)                
                    if status == GRB.OPTIMAL:
                        # make a list of edges selected in the solution
                        xVal = {}
                        zVal = {}
                        for i in nodeArray:
                            xVal[i] = model.cbGetNodeRel(X[i])
                
                        for l in lineArray:
                            simpleG[sourceNode[l]][targetNode[l]]['capacity'] = 0
                
                        for l in lineArray:
                            zVal[l] = model.cbGetNodeRel(Z[l])
                            simpleG[sourceNode[l]][targetNode[l]]['capacity'] += zVal[l]
                            
                         
                        tolGrand = 0.0001
                        violation = 0
                        for i in nodeArray:
                            if i != 1:
                                cut_value, partition = nx.minimum_cut(simpleG, 1, i, capacity = 'capacity')
                                reachable, non_reachable = partition
                                if cut_value < xVal[i] - tolGrand:
                                    if violation < - cut_value + xVal[i]: 
                                        violation = - cut_value + xVal[i]
                                    LHS = [(-1,X[i])]
                                    for l in lineArray:
                                        if sourceNode[l] in reachable and targetNode[l] in non_reachable:
                                            LHS += [(1,Z[l])]
                                        if sourceNode[l] in non_reachable and targetNode[l] in reachable:
                                            LHS += [(1,Z[l])]
                                    model.cbCut(LinExpr(LHS)>=0)
                                    
                                    
                        for c in theContingent:
                            simpleGC = copy.deepcopy(simpleG)
                            simpleGC[sourceNode[c]][targetNode[c]]['capacity'] = simpleGC[sourceNode[c]][targetNode[c]]['capacity'] - zVal[c]
                            cut_valueC, partitionC = nx.minimum_cut(simpleGC, sourceNode[c], targetNode[c], capacity = 'capacity')
                            reachableC, non_reachableC = partitionC
                            if cut_valueC < zVal[c] - tolGrand:
                                if violation < - cut_valueC + zVal[c]:
                                    violation = - cut_valueC + zVal[c]
                                LHSC = [(-1,Z[c])]
                                for l in lineArray:
                                    if l != c:
                                        if sourceNode[l] in reachableC and targetNode[l] in non_reachableC:
                                            LHSC += [(1,Z[l])]
                                        if sourceNode[l] in non_reachableC and targetNode[l] in reachableC:
                                            LHSC += [(1,Z[l])]
                                model.cbCut(LinExpr(LHSC)>=0)

                            if cut_valueC < xVal[sourceNode[c]] + xVal[targetNode[c]] - 1 - tolGrand:
                                if violation < - cut_valueC + xVal[sourceNode[c]] + xVal[targetNode[c]] - 1:
                                    violation = - cut_valueC + xVal[sourceNode[c]] + xVal[targetNode[c]] - 1
                                LHSC = [(-1,X[sourceNode[c]]),(-1,X[targetNode[c]])]
                                for l in lineArray:
                                    if l != c:
                                        if sourceNode[l] in reachableC and targetNode[l] in non_reachableC:
                                            LHSC += [(1,Z[l])]
                                        if sourceNode[l] in non_reachableC and targetNode[l] in reachableC:
                                            LHSC += [(1,Z[l])]
                                model.cbCut(LinExpr(LHSC)>=-1)


                        print('violation =',violation)
                        if violation < 0.0001:
                            for l in lineArray:
                                tol = 0.0001
                                if - 2 * zVal[l] - 1 + xVal[sourceNode[l]] + xVal[targetNode[l]] > 0.0001:
                                #if violation + 0.85 < - 2 * zVal[l] - 1 + xVal[sourceNode[l]] + xVal[targetNode[l]]:
                                    for c in theContingent:
                                        if c != l:
                                            simpleGC = copy.deepcopy(simpleG)
                                            simpleGC[sourceNode[c]][targetNode[c]]['capacity'] = simpleGC[sourceNode[c]][targetNode[c]]['capacity'] - zVal[c]
                                            simpleGC[sourceNode[l]][targetNode[l]]['capacity'] = simpleGC[sourceNode[l]][targetNode[l]]['capacity'] - zVal[l]
                                            cut_valueC, partitionC = nx.minimum_cut(simpleGC, sourceNode[l], targetNode[l], capacity = 'capacity')
                                            reachableC, non_reachableC = partitionC
                                            if  0.0001 < - cut_valueC - 2 * zVal[l] - 1 + xVal[sourceNode[l]] + xVal[targetNode[l]]:
                                                    
                                                # print('coparallel violation =',- 2 * zVal[l] - 1 + xVal[sourceNode[l]] + xVal[targetNode[l]] - cut_valueC)
                                                LHSC = [(2,Z[l]),(-1,X[sourceNode[l]]),(-1,X[targetNode[l]])]
                                                for l1 in lineArray:
                                                    if l1 != c and l1 != l:
                                                        if sourceNode[l1] in reachableC and targetNode[l1] in non_reachableC:
                                                            LHSC += [(1,Z[l1])]
                                                        if sourceNode[l1] in non_reachableC and targetNode[l1] in reachableC:
                                                            LHSC += [(1,Z[l1])]
                                                model.cbCut(LinExpr(LHSC)>=-1)
                                                

            
            if spanning == True:
                for i in nodeArray:
                    LHS = [(1,X[i])]
                    model.addConstr(LinExpr(LHS)==1, name='Fix X (%s)'%i)
            
            
            
                
            objTerms = []
            for l in lineArray:
                objTerms += [(riskFunction[l],Z[l])]
            
                                            
            model.setObjective(LinExpr(objTerms), GRB.MINIMIZE)
            model.setParam('PreCrush', 1)    
            model.update()
            model.Params.lazyConstraints = 1
            model.optimize(security)
            best_bound = model.ObjBound
            #security = -1
            print('### best bound =',best_bound)
            
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
            optSolution.to_csv(r'result%s/opt/opt_N-1_PSPS_Span%s_BCZ_Enhanced_%s_inst%s.csv'%(len(nodeArray),spanning,len(nodeArray),inst), index = False)#Check
            
                        
            sizeArray += [len(nodeArray)]
            instArray += [inst]
            ratioDemandRequiredArray += [ratioDemandRequired]
            totalDemandArray += [totalDemand]
            totalRiskArray += [totalRisk]
            resultRiskArray += [best_bound]
            ratioRiskArray += [best_bound/totalRisk]
            timeArray += [model.Runtime]
            nodecountArray += [model.NodeCount]
            securityArray += [-1]
            machineArray += [machineName]
            spanningArray += [spanning]
            
            
            resultTable = pd.DataFrame(list(zip(sizeArray,instArray, totalDemandArray, ratioDemandRequiredArray, totalRiskArray, resultRiskArray, ratioRiskArray, timeArray, nodecountArray, securityArray, machineArray, spanningArray)),columns =['Nodes', 'inst', 'totalDemand','ratioDemandRequired', 'totalRiskArray', 'riskBd', 'ratioBd', 'Time', 'B&B Nodes', 'Security', 'Machine', 'Spanning'])
            resultTable.to_csv(r'result%s/result/result_%s_N-1_PSPS_Span%s_BCZ_Enhanced_%s.csv'%(len(nodeArray),machineName,spanning,len(nodeArray)), index = False)#Check
            
            

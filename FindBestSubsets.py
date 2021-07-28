#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:37:48 2019

@author: elinewby
"""

import networkx as nx
from FVS_python import FVS as FVS
import pandas as pd
import itertools
import numpy as np
import scipy as sc
import numpy.linalg as la
import math
import scipy.stats as sp

def findBestFVSSubsets(netName, numNodes, intersections = [['OutDegree','Distance','PRINCE','ModPRINCE','CheiRank','Cycles','SCC'],['PRINCE', 'ModPRINCE','CheiRank'],['ModPRINCE','CheiRank']]):
    net = nx.read_graphml(netName)
    revNet = nx.reverse(net)    
    metricDict = {0:'OutDegree', 1:'Distance', 2:'PRINCE', 3:'ModPRINCE', 4:'CheiRank', 5:'Cycles', 6:'SCC'}
    #functions
    def combos(nodes, num):
        combs = []
        combList = list(itertools.combinations(nodes, num))
        for comb in combList:
            combs.append(sorted(comb))
        return combs
    
    def averageValues(net, arr):
        nodes = list(net.nodes())
        combs = combos(nodes, numNodes)
        avgs = np.zeros(len(combs))
        combNames = []
        for i in range(len(combs)):
            comb = ''
            for j in range(len(combs[i])):
                comb += combs[i][j] + ', '
            combNames.append(comb[:-2])
            values = arr[i]
            j = 0
            while(j < len(values)):
                if(values[j] == -1):
                    values = np.delete(values,j)
                    j -= 1
                j += 1
            avgs[i] = np.sum(abs(values))/len(nodes)
        
        return pd.DataFrame({'Intervention':combNames,'Value':avgs})
    
    def distance(net, numPerts):
        nodes = list(net.nodes())
        combs = combos(nodes, numNodes)
        edges = list(net.edges(data = True))
        adj = np.asarray(nx.adjacency_matrix(net).todense()).astype(np.double)
        signedAdj = np.copy(adj)
        for edge in edges:
            source = edge[0]
            sink = edge[1]
            sign = edge[2]['data']
            if(sign == -1):
                signedAdj[nodes.index(source),nodes.index(sink)] = -1
        dist = np.zeros((len(combs),len(nodes)))
        for i in range(len(combs)):
            for j in range(len(nodes)):
                minDists = []
                for node in combs[i]:
                    try:
                        minDists.append(nx.shortest_path_length(net, node, nodes[j]))
                    except:
                        minDists.append(1E16)
                minDist = min(minDists)
                dist[i][j] = 1/(1+minDist)
        return averageValues(net, dist)
    
    def PRINCEPropagation(net, isModified):
        nodes = list(net.nodes())
        combs = combos(nodes, numNodes)
        edges = list(net.edges(data = True))
        adj = np.asarray(nx.adjacency_matrix(net).todense()).astype(np.double)
        signedAdj = np.copy(adj)
        for edge in edges:
            source = edge[0]
            sink = edge[1]
            sign = edge[2]['data']
            if(sign == -1):
                signedAdj[nodes.index(source),nodes.index(sink)] = -1
        inDegrees = np.zeros(len(nodes))
        outDegrees = np.zeros(len(nodes))
        for i in range(len(nodes)):
            if(net.in_degree(nodes[i]) > 0):
                if(isModified):
                    inDegrees[i] = 1/net.in_degree(nodes[i])
                else:
                    inDegrees[i] = 1/(net.in_degree(nodes[i])**(1/2))
                    #inDegrees[i] = 1/(net.in_degree(nodes[i]))
            if(net.out_degree(nodes[i]) > 0):
                outDegrees[i] = 1/(net.out_degree(nodes[i])**(1/2))
        D1 = np.diag(outDegrees)
        D2 = np.diag(inDegrees)
        if(isModified):
            W = np.matmul(signedAdj, D2)
        else:
            W = np.matmul(D1,np.matmul(signedAdj,D2))
        a = 0.9
        numInputs = 2**(numNodes-1)
        i = 0
        j = 0
        negativeInputs = []
        while(i < numInputs):
            neg = combos(np.linspace(0,numNodes-1,numNodes), j)
            for n in neg:
                negativeInputs.append(n)
                i += 1
                if(i >= numInputs):
                    break  
            j += 1
        averages = np.zeros((len(combs),numInputs))
        for i in range(numInputs):
            Y = np.zeros((len(nodes),len(combs)))
            for c in combs:
                for node in c:
                    if(c.index(node) in negativeInputs[i]):
                        Y[nodes.index(node)][combs.index(c)] = -1
                    else:
                        Y[nodes.index(node)][combs.index(c)] = 1
            S = (1-a)*np.matmul(la.inv(np.identity(len(nodes))-a*W.transpose()),Y)
            combNames = []
            for j in range(len(combs)):
                comb = ''
                for k in range(len(combs[j])):
                    comb += combs[j][k] + ', '
                combNames.append(comb[:-2])
                for k in range(len(nodes)):
                    averages[j][i] += abs(S[k][j])/len(nodes)
        #maxOfAvgs = np.zeros((len(combs),1))
        maxOfAvgs = np.zeros(len(combs))
        for i in range(len(averages)):
            maxOfAvgs[i] = max(averages[i])
        return pd.DataFrame({'Intervention':combNames,'Value':maxOfAvgs})
    
    def pageRank(net):
        nodes = list(net.nodes())
        combs = combos(nodes, numNodes)
        combNames = []
        pR = nx.pagerank(net)
        result = np.zeros(len(combs))
        for i,comb in enumerate(combs):
            combName = ''
            notP = 1
            for c in comb:
                combName += c + ', '
                notP*=(1-pR[c])
            P = 1-notP
            result[i] = np.round(P, decimals = 3)
            combNames.append(combName[:-2])
        return pd.DataFrame({'Intervention':combNames,'Value':result})
    
    
    def nearestNeighbor(net):
        nodes = list(net.nodes())
        combs = combos(nodes, numNodes)
        outDegrees = []
        combNames = []
        for comb in combs:
            combName = ''
            outDegree = 0
            for c in comb:
                combName += c + ', '
                outDegree += net.out_degree(c)
            combNames.append(combName[:-2])
            outDegrees.append(outDegree)
        return pd.DataFrame({'Intervention':combNames,'Value':outDegrees})
    
    def findCycles():
        nodes  = df_cycles['node']
        cycleSize = df_cycles['cycles']
        combNames = []
        cycles = []
        for i in range(len(nodes)):
            combNames.append(nodes[i])
            cycles.append(cycleSize[i])
        return pd.DataFrame({'Intervention':combNames,'Value':cycles})
    
    def findSCC():
        nodes  = df_scc['node']
        sccSize = df_scc['scc']
        combNames = []
        scc = []
        for i in range(len(nodes)):
            combNames.append(nodes[i])
            scc.append(sccSize[i])
        return pd.DataFrame({'Intervention':combNames,'Value':scc})
    
    #Find FVS subsets and calculate FVS-specific metrics
    nodes_num=[]
    nodes_cycle=[]
    nodes_comb=[]
    nodes_combnum=[]
    nodes_scc=[]
    
    itrs=100
    print('Finding FVS Subsets:')
    for i in range(itrs):
        print(i,"%");
        G_FVS=FVS.FVS(net, randomseed = i)
        for n in range(numNodes,numNodes+1):
            addedCombs = []
            for comb in itertools.combinations(G_FVS,n):
                comb_ord=tuple(sorted(list(comb)))
                addedCombs.append(comb_ord)
                FVS2=[x for x in G_FVS if x not in comb_ord]
                H=net.copy()
                H.remove_nodes_from(FVS2)
                if(comb_ord not in nodes_comb):
                    nodes_comb.append(comb_ord)
                    nodes_num.append(1)
                    nodes_combnum.append(n)
                    cycles = nx.simple_cycles(H)
                    numPos = 0
                    for cycle in cycles:
                        negENum = 0
                        for j in range(len(cycle)):
                            if(H.get_edge_data(cycle[j],cycle[(j+1)%len(cycle)])['data'] == -1):
                                negENum += 1
                        if(negENum % 2 == 0):
                            numPos += 1
                    nodes_cycle.append(numPos)
                    nodes_scc.append(len(max(nx.strongly_connected_components(H), key=len)))
                else:
                    cycles = nx.simple_cycles(H)
                    numPos = 0
                    for cycle in cycles:
                        negENum = 0
                        for j in range(len(cycle)):
                            if(H.get_edge_data(cycle[j],cycle[(j+1)%len(cycle)])['data'] == -1):
                                negENum += 1
                        if(negENum % 2 == 0):
                            numPos += 1
                    index=nodes_comb.index(comb_ord)
                    nodes_cycle[index] = nodes_cycle[index] + numPos
                    nodes_num[index]=nodes_num[index]+1
                    nodes_scc[index]=nodes_scc[index]+len(max(nx.strongly_connected_components(H), key=len))
            for comb in itertools.combinations(net, n):
                comb_ord = tuple(sorted(list(comb)))
                if(comb_ord not in addedCombs):
                    if(comb_ord not in nodes_comb):
                        nodes_comb.append(comb_ord)
                        nodes_combnum.append(n)
                        nodes_cycle.append(0)
                        nodes_num.append(0)
                        nodes_scc.append(0)
    
    print("100%")
    nodes_num=[x/itrs for x in nodes_num]
    nodes_cycle=[x/itrs for x in nodes_cycle]
    nodes_cycle2=[nodes_cycle[i]/nodes_num[i] if (nodes_num[i]!=0) else nodes_cycle[i] for i in range(len(nodes_cycle))]
    df_cycles=pd.DataFrame({"cycles": nodes_cycle2,"num": nodes_num,"combnum": nodes_combnum,"node":[", ".join(x) for x in nodes_comb]})
    df_cycles = df_cycles.sort_values(by = ['combnum','cycles','node'],ascending=False)
    
    nodes_scc=[x/itrs for x in nodes_scc]
    nodes_scc2=[nodes_scc[i]/nodes_num[i] if (nodes_num[i]!=0) else nodes_scc[i] for i in range(len(nodes_scc))]
    df_scc=pd.DataFrame({"scc": nodes_scc2, "num": nodes_num,"combnum": nodes_combnum,"node":[", ".join(x) for x in nodes_comb]})
    df_scc = df_scc.sort_values(by=['combnum','scc'],ascending=False)
    
    combs = list(df_cycles['node'])
    inFVS = list(df_cycles['num'])
    FVSsubsets = []
    FVScombos = []
    for i in range(len(combs)):
        if(inFVS[i] > 0):
            FVSsubsets.append(combs[i])
            FVScombos.append(combs[i].split(', '))
    print()
    print('Finding Topological Metric Values:')
    #Create DataFrames for each Topological Metric
    metricValsDict = {}
    metricValsDict['OutDegree'] = nearestNeighbor(net)
    print('Out-Degree Done')
    metricValsDict['Distance'] = distance(net, False)
    print('Distance Done')
    metricValsDict['PRINCE'] = PRINCEPropagation(net, False)
    print('PRINCE Done')
    metricValsDict['ModPRINCE'] = PRINCEPropagation(net, True)
    print('Modified PRINCE Done')
    metricValsDict['CheiRank'] = pageRank(revNet)
    print('CheiRank Done')
    metricValsDict['Cycles'] = findCycles()
    print('Cycles Done')
    metricValsDict['SCC'] = findSCC()
    print('SCC Done')
    
    
    #Calculate Intersection Percentiles
    resultDict = {}
    topologyDict = {}
    FVSOnlyTopologyDict = {}
    
    for metricNum in range(7):
        df2 = metricValsDict[metricDict[metricNum]]
        perts2 = list(df2['Intervention'].values)
        topoVals = list(df2['Value'].values)
    
        for i in range(len(perts2)):
            if(perts2[i] in FVSsubsets):
                if(perts2[i] in FVSOnlyTopologyDict):
                    FVSOnlyTopologyDict[perts2[i]].append(topoVals[i])
                else:
                    FVSOnlyTopologyDict[perts2[i]] = [topoVals[i]]
                resultDict[perts2[i]] = []
            if(perts2[i] in  topologyDict):
                topologyDict[perts2[i]].append(topoVals[i])
            else:
                topologyDict[perts2[i]] = [topoVals[i]]
    
    for metrics in intersections:
        print()
        print('Calculating Percentile Cutoff values for Intersection '+str(intersections.index(metrics)+1))
        goodPerts = []
        for threshold in np.linspace(0.01,1,100):
            print(str(int(100*(threshold-0.01)))+'%')
            goodPertsOld = goodPerts.copy()
            pertDict = {}
            pertDict2 = {}
            for metricNum in range(7):
                reverseTopologyDict = {}
                for k, v in FVSOnlyTopologyDict.items():
                    if(v[metricNum] in reverseTopologyDict):
                        reverseTopologyDict[v[metricNum]].append(k)
                    else:
                        reverseTopologyDict[v[metricNum]] = [k]
                A = []
                for k in FVSsubsets:
                    A.append(topologyDict[k][metricNum])
                topFVSPertsThresh = math.ceil(len(A)*threshold)
                
                topFVSPerts = []
                highestValsFVS = list(dict.fromkeys(sorted(A, reverse = True)[:topFVSPertsThresh]))
                for h in highestValsFVS:
                    for pert in reverseTopologyDict[h]:
                        topFVSPerts.append(pert)
                    if(len(topFVSPerts) >= topFVSPertsThresh):
                        break
                if(metricDict[metricNum] in metrics):
                    for t in topFVSPerts:
                        if(t in pertDict):
                            pertDict[t]+=1
                            pertDict2[t].append(metricDict[metricNum])
                        else:
                            pertDict[t]=1
                            pertDict2[t] = [metricDict[metricNum]]
            for k,v in pertDict.items():
                if(v == len(metrics)):
                    goodPerts.append(k)
            for i in range(len(goodPerts)):
                if(goodPerts[i] not in goodPertsOld):
                    resultDict[goodPerts[i]].append(np.round((1-threshold)*100))
        print('100%')
    resultDict = {k:v for k,v in sorted(resultDict.items(), key = lambda item:item[1])}
    FVSsubsets = list(resultDict.keys())
    
    return pd.DataFrame(topologyDict, index = metricDict.values()).transpose(), pd.DataFrame(resultDict, index = [', '.join(intersection) for intersection in intersections]).transpose(), FVSsubsets
topoVals, interVals, subsets = findBestFVSSubsets('TLGL.graphml', 1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:28:15 2021

@author: elinewby
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 19:46:31 2020

@author: elinewby
"""
import biolqm
import numpy as np
import pandas as pd
import itertools

def calculateToandAwayControl(netName, subsets, wtSims = 10, subsetSims = 10, seed = 0, verbose = False, write = False, writePath = 'ControlVals.xlsx'):
    
    def writePerturbation(subsets):
        perturbations = []
        for subset in subsets:
            nodes = [n.strip() for n in subset.split(',')]
            for i in range(2**len(nodes)):
                form = "{0:0"+str(len(nodes))+"b}"
                pertVals = list(form.format(i))
                perturbation = [nodes[j]+'%'+str(pertVals[j]) for j in range(len(nodes))]
                if(not (perturbation in perturbations) and len(perturbation) == len(nodes)):
                    perturbations.append(perturbation)
        return perturbations

    def willWrite(path):
        pd.set_option('display.max_colwidth', -1)
        ResultFile = open(path, "w")
        ResultFile.write(attDF.transpose().to_string())
        ResultFile.write("\n\n\n")
        ResultFile.write(resultHamming.transpose().to_string())
        ResultFile.close()
    
    lqm = biolqm.load(netName)
    nodes = [n.toString() for n in lqm.getComponents()]
    if(len(biolqm.fixpoints(lqm)) < len(biolqm.trapspaces(lqm))):
        trap = biolqm.trapspace(lqm)
        isComplex = True
    else:
        trap = biolqm.fixpoints(lqm)
        isComplex = False
    attDF = pd.DataFrame([s for s in trap])
    atts = attDF.values.astype(float)
    
    nodesResultOrder = list(attDF.columns)
    nonFixed = {i:[] for i in range(len(attDF))}
    
    if(isComplex):
        for i in range(len(atts)):
            for j in range(len(atts[i])):
                if(atts[i][j] > 1):
                    nonFixed[i].append(nodesResultOrder[j])
                    atts[i][j] = 0.5
    
    attTypes = {}
    
    attTypes['Type 0'] = [0]
    N = len(nodes)
    numAtts = len(atts)
    numTypes = 1
    for i in range(1,len(atts)):
        typeFound = False
        for j in range(numTypes):
            attType = 'Type ' + str(j)
            minDistance = 1
            for k in range(len(attTypes[attType])):
                aNum = attTypes[attType][k]
                distance = np.sum(abs(atts[aNum] - atts[i]))/len(nodes)
                if(distance < minDistance):
                    minDistance = distance
            if(minDistance <= 0.15 and not typeFound):
                attTypes[attType].append(i)
                typeFound = True
        if(not typeFound):
            attTypes['Type ' + str(numTypes)] = [i]
            numTypes += 1
    
    if(numTypes == 1):
        attTypes = {}
        for i in range(len(atts)):
            attTypes['Type ' + str(i)] = [i]
    
    combos = list(itertools.combinations(attTypes.keys(),2))
    diffs = np.zeros(len(combos))
    for i in range(len(combos)):
        type1 = combos[i][0]
        type2 = combos[i][1]
        att1 = atts[attTypes[type1][0]]
        att2 = atts[attTypes[type2][0]]
        diffs[i] = np.sum(abs(att1-att2))
    
    
    
    #Test Subsets
    pertNodes2 = [[]]
    pertNodes2.extend(writePerturbation(subsets))
    resultTableHamming = {}
    resultDict = {}
    
    np.random.seed(seed)
    pertAttTargets = {}
    goodTo = {}
    goodAway = {}
    for subset in subsets:
        goodTo[subset] = False
        goodAway[subset] = False
    
    for perts in pertNodes2:
        if(perts == []):
            numSims = wtSims
        else:
            numSims = subsetSims
        result = np.zeros(numAtts)
        
        isBad = False
        perturbation = ""
        subset = ""
        for a in perts:
            node = a.split('%')[0]
            subset += node+', '
            perturbation += a
        subset = subset[:-2]
        pertAttTargets[perturbation] = []
        for i in range(numAtts):
            isTo = True
            for p in perts:
                node, val = p.split('%')
                val = int(val)
                if(attDF[node][i] != val):
                    isTo = False
            pertAttTargets[perturbation].append(isTo)
            
        if(sum(pertAttTargets[perturbation]) == numAtts):
            isBad = True
        if(sum(pertAttTargets[perturbation]) < numAtts):
            goodAway[subset] = True
            if(sum(pertAttTargets[perturbation]) > 0):
                goodTo[subset] = True
        if(len(perts) == 0):
            isBad = False
        
        if(isBad):
            if(verbose):
                print()
                print(perts,"Bad Intervention")
            resultTableHamming[perturbation] = []
            names = []
            for i in range(numAtts):
                resultTableHamming[perturbation].append(np.round(result[i], decimals = 3))
                if(i == numAtts):
                    names.append("    Middle")
                else:
                    names.append("    Attractor " + str(i))
            resultDict[perturbation] = [0,0]
            continue
        
        L = len(perts)
        pertNodeIndex = []
        if(verbose):
            print("\nPerturbing ", perts,"\n") 
        pert = lqm
        if(int(L)>=1):
            for n in perts:
                pert = biolqm.perturbation(pert, n)
                pertNode = n.split('%')[0]
                pertNodeIndex.append(nodesResultOrder.index(pertNode))
        pertFps = biolqm.fixpoints(pert)
        pertAtts2 = pd.DataFrame([s for s in pertFps])
        pertAtts = pertAtts2.values
        pertIsComplex = False
        
        if(len(biolqm.fixpoints(pert)) < len(biolqm.trapspaces(pert))):
            pertTrap = biolqm.trapspace(pert)
            pertIsComplex = True
        else:
            pertTrap = biolqm.fixpoints(pert)
            pertIsComplex = False
    
        pertAtts2 = pd.DataFrame([s for s in pertTrap])
        pertAtts = pertAtts2.values
        
        pertNonFixed = {i:[] for i in range(len(pertAtts))}
        if(pertIsComplex):
            for i in range(len(pertAtts)):
                for j in range(len(atts[0])):
                    if(pertAtts[i][j] > 1):
                        pertNonFixed[i].append(nodesResultOrder[j])
                        
        for i in range(numSims):
            if(verbose):
                if(numSims < 100):
                    print("{0:2.0f}".format(i/numSims*100),"%")
                else:
                    if(i%(numSims/100)==0):
                        print("{0:2.0f}".format(i/numSims*100),"%")
                    
            initial1 = np.random.randint(0,2**(N-10))
            initial2 = np.random.randint(0,2**10)
            form = "{0:0"+str(N-10)+"b}"
            initString1 = form.format(initial1)
            initString2 = "{0:010b}".format(initial2)
            initString = initString1+initString2
            
            isDone = False
            steps = 0
            while(not isDone and steps < 20):
                hamming = np.zeros(numAtts)
                hammingPerts = np.zeros(len(pertAtts))
                walk = biolqm.random(pert, "-i "+initString+" -m 50 -s "+str(seed))
                Data = pd.DataFrame([s for s in walk])
                S = Data.values[-1]
                steps += 1
                for j in range(len(hamming)):
                    for k in range(len(S)):
                        if(k not in pertNodeIndex and nodesResultOrder[k] not in nonFixed[j]):
                            hamming[j] += abs(atts[j][k]-S[k])
                if(pertIsComplex):
                    for j in range(len(hammingPerts)):
                        for k in range(len(S)):
                            if(nodesResultOrder[k] not in pertNonFixed[j]):
                                hammingPerts[j] += abs(pertAtts[j][k]-S[k])
                else:
                    for j in range(len(hammingPerts)):
                        hammingPerts[j] = np.sum(abs(pertAtts[j]-S))
                
                initString = ""
                for node in nodes:
                    initString += str(Data[node].values[-1])
    
                if(min(hammingPerts) == 0.0):
                    isDone = True  
    
            hamming = hamming/(len(nodes)-(len(pertNodeIndex)))
            minValue = min(hamming)
            
            bestAtt = np.where(hamming == minValue)[0][0]
            
            closeThreshold = min(diffs)/len(nodes)/2
            
            if(minValue <= closeThreshold):
                result[bestAtt] += 1-minValue/closeThreshold
        
        if(verbose):
            print ("100%")
        resultTableHamming[perturbation] = []
        names = []
        for i in range(len(hamming)):
            resultTableHamming[perturbation].append(np.round(result[i]/numSims, decimals = 3))
            names.append("    Attractor " + str(i))
            
        WTTB, WTNTB, TB, NTB = 0,0,0,0
        for i in range(numAtts):
            if(pertAttTargets[perturbation][i]):
                WTTB += resultTableHamming[''][i]
                TB += resultTableHamming[perturbation][i]
            else:
                WTNTB += resultTableHamming[''][i]
                NTB += resultTableHamming[perturbation][i]
                
        toControl, awayControl = 0,0
        if(sum(pertAttTargets[perturbation]) < numAtts and sum(pertAttTargets[perturbation]) > 0):
            if(TB >= WTTB):
                toControl = (TB-WTTB)/(1-WTTB)
            else:
                toControl = (TB-WTTB)/WTTB
            if(NTB >= WTNTB):
                awayControl = (WTNTB-NTB)/(1-WTNTB)
            else:
                awayControl = (WTNTB-NTB)/WTNTB
        
        elif(sum(pertAttTargets[perturbation]) < numAtts):
            awayControl = WTNTB - NTB
        
        resultDict[perturbation] = [np.round(toControl, decimals = 2), np.round(awayControl, decimals = 2)]
    
    resultHamming = pd.DataFrame(resultTableHamming, index = names).transpose()
    #pertResults = pd.DataFrame(resultDict, index = ['To Control','Away Control']).transpose()
    
    combResultDict = {}
    for perts in pertNodes2:
        perturbation = ''
        subset = ''
        for a in perts:
            node = a.split('%')[0]
            perturbation += a
            subset += node+', '
        subset = subset[:-2]
        if(subset in combResultDict):
            if(resultDict[perturbation][0] > combResultDict[subset][0]):
                combResultDict[subset][0] = resultDict[perturbation][0]
            if(resultDict[perturbation][1] > combResultDict[subset][1]):
                combResultDict[subset][1] = resultDict[perturbation][1]
        else:
            combResultDict[subset] = resultDict[perturbation]
            
    results = pd.DataFrame(combResultDict, index = ['To Control','Away Control']).transpose().drop('')
    if(write):
        writer = pd.ExcelWriter(writePath, engine = 'xlsxwriter')
        results.to_excel(writer)
        pd.DataFrame(goodAway.values(), columns = ['Partially Informative']).to_excel(writer, index = False, startcol=3)
        pd.DataFrame(goodTo.values(), columns = ['Fully Informative']).to_excel(writer, index = False, startcol=4)
        writer.save()
    return results, goodTo, goodAway

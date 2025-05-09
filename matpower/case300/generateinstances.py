import pandas as pd
import random
import copy

numNodes = 300
lines = pd.read_csv('lines_%s.csv'%numNodes)
numInstances = 50

for inst in range(11,numInstances+1):
    instance = copy.deepcopy(lines)
    randoms = []
    theSwitchable = []
    for l in lines['Line']:
        randoms += [int(random.random()*10+1)]
        theSwitchable += [1]
    instance['Switchable'] = theSwitchable
    instance['Risk'] = randoms

    instance.to_csv('lines_%s_inst%s.csv'%(numNodes,inst), index=False)
        

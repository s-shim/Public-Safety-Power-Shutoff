# Python program to read
# json file

import json
import pandas as pd

# Opening JSON file
f = open('2017-02-01.json')

# returns JSON object as 
# a dictionary
data = json.load(f)

nodes = pd.read_csv('nodes_%s.csv'%len(data['Buses']))
nodeOfBus = {}
busOfNode = {}
for i in nodes['Node']:
    [bus_i] = nodes.loc[nodes['Node']==i,'Bus']
    nodeOfBus[bus_i] = i
    busOfNode[i] = bus_i

lines = []
sources = []
targets = []
sourcesN = []
targetsN = []
susceptances = []
# Iterating through the json
# list
for l in data["Transmission lines"]:
    lines += [l[1:]]
    sources += [data["Transmission lines"][l]["Source bus"]]
    targets += [data["Transmission lines"][l]["Target bus"]]
    sourcesN += [nodeOfBus[data["Transmission lines"][l]["Source bus"]]]
    targetsN += [nodeOfBus[data["Transmission lines"][l]["Target bus"]]]
    susceptances += [data["Transmission lines"][l]["Susceptance (S)"]]

list_of_tuples = list(zip(lines,sources,targets,sourcesN,targetsN,susceptances))
df = pd.DataFrame(list_of_tuples, columns=['Line','SourceBus','TargetBus','Source','Target','Susceptance'])
df.to_csv('lines_SHS_%s.csv'%len(data["Buses"]), index=False)

# Closing file
f.close()

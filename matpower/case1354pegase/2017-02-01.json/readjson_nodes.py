# Python program to read
# json file

import json
import pandas as pd

# Opening JSON file
f = open('2017-02-01.json')

# returns JSON object as 
# a dictionary
data = json.load(f)

buses = []
busesN = []
loads = []
# Iterating through the json
# list
noBus = 0
for i in data['Buses']:
    noBus += 1
    buses += [i]
    busesN += [noBus]
    if isinstance(data['Buses'][i]['Load (MW)'], float) == True:
        loads += [data['Buses'][i]['Load (MW)']]
    if isinstance(data['Buses'][i]['Load (MW)'], list) == True:
        loads += [sum(data['Buses'][i]['Load (MW)'])/len(data['Buses'][i]['Load (MW)'])]

min_generations = []
max_generations = []
# Iterating through the json
# list
for j in buses:
    min_generation = 0.0
    max_generation = 0.0
    for i in data['Generators']:
        if data['Generators'][i]['Bus'] == j:
            min_generation += data['Generators'][i]['Production cost curve (MW)'][0]
            max_generation += data['Generators'][i]['Production cost curve (MW)'][-1]
    min_generations += [min_generation]
    max_generations += [max_generation]

list_of_tuples = list(zip(busesN,buses,loads,min_generations,max_generations))
df = pd.DataFrame(list_of_tuples, columns=['Node','Bus','Load','min_generation','max_generation'])
df.to_csv('nodes_%s.csv'%len(buses), index=False)

# Closing file
f.close()

# -*- coding: utf-8 -*-
"""
This script is used for getting the negative (or decoy) datasets.

Written by Dong-Sheng Cao

Date: 2013.8.31
"""

import string
import random
import math
import networkx as nx
import cPickle

path = "~/MEF/3-generate-negative-set/"

interaction = []
targetlist = []
druglist = []
f = file(path + "association1.txt", 'rb')
for i in f:
    temp = i.split('\t')
    temp0 = string.strip(temp[0])
    temp1 = string.strip(temp[1])
    interaction.append((temp0, temp1))
    if temp1 not in druglist:
        druglist.append(temp1)
    if temp0 not in targetlist:
        targetlist.append(temp0)
f.close()

G = nx.Graph()
G.add_edges_from(interaction)
# print G.number_of_nodes()
# print G.number_of_edges()
# print G.nodes()
# print G.neighbors("hsa:8856")

interactiondict = {}
for i in targetlist:
    temp = G.neighbors(i)
    interactiondict.update({i: temp})

print len(interaction)
print len(druglist)
print len(targetlist)

index = True
decoy = list()
drugnum = len(druglist)
while index:
    for i in targetlist:
        tindex = int(math.floor(random.random() * drugnum))
        drug1 = druglist[tindex]
        if drug1 not in interactiondict[i]:
            decoy.append([i, drug1])
        if len(decoy) >= len(interaction):
            index = False
            break

print len(decoy)

f = file(path + 'decoy/decoy10.txt', 'wb')
for i in decoy:
    f.write(i[0] + '\t' + i[1] + '\n')
f.close()

decoy1 = []
for i in decoy:
    decoy1.append([int(i[0]), int(i[1])])

f = file(path + 'decoy/decoy10', 'wb')
cPickle.dump(decoy1, f)
f.close()

# -*- coding: utf-8 -*-
"""
This script is used for getting the MedDRA-based ADR similarity measures.

Created on Mon Aug 19 16:06:09 2013

@author: Dong-Sheng Cao
"""

import cPickle
import MedDRA
import string
import scipy

def FindCommonNumber(list1, list2):
    res = 0
    for i in list1:
        if i in list2:
            res = res + 1
    return res

f = file("~/adr/data/adrass", 'r')
adrass = cPickle.load(f)
f.close()

f = file("~/adr/meddra/llt", 'r')
llt = cPickle.load(f)
f.close()

f = file("~/adr/meddra/pt-llt", 'r')
ptllt = cPickle.load(f)
f.close()

f = file("~/adr/meddra/hlt-llt", 'r')
hltllt = cPickle.load(f)
f.close()

f = file("~/adr/meddra/hlgt-llt", 'r')
hlgtllt = cPickle.load(f)
f.close()

f = file("~/adr/meddra/soc-llt", 'r')
socllt = cPickle.load(f)
f.close()

meddraid = []
for i in adrass:
    temp = string.lower(i[1])
    for j in llt:
        if temp == string.lower(llt[j]):
            meddraid.append(j)

##########################################################

soc = []
for i in meddraid:
    temp = MedDRA.GetHighLevelTerm(i, socllt)
    soc.append(temp)

Nsize = len(adrass)
SOCSimilarity = scipy.zeros((Nsize, Nsize))

for i in range(Nsize):
    temp1 = soc[i]
    for j in range(i + 1):
        temp2 = soc[j]
        c = FindCommonNumber(temp1, temp2)
        a = len(temp1)
        b = len(temp2)
        SOCSimilarity[i, j] = (c + 0.0) / (a + b - c)
        SOCSimilarity[j, i] = (c + 0.0) / (a + b - c)

##########################################################

hlgt = []
for i in meddraid:
    temp = MedDRA.GetHighLevelTerm(i, hlgtllt)
    hlgt.append(temp)

Nsize = len(adrass)
HLGTSimilarity = scipy.zeros((Nsize, Nsize))

for i in range(Nsize):
    temp1 = hlgt[i]
    for j in range(i + 1):
        temp2 = hlgt[j]
        c = FindCommonNumber(temp1, temp2)
        a = len(temp1)
        b = len(temp2)
        HLGTSimilarity[i, j] = (c + 0.0) / (a + b - c)
        HLGTSimilarity[j, i] = (c + 0.0) / (a + b - c)

##########################################################

hlt = []
for i in meddraid:
    temp = MedDRA.GetHighLevelTerm(i, hltllt)
    hlt.append(temp)

Nsize = len(adrass)
HLTSimilarity = scipy.zeros((Nsize, Nsize))

for i in range(Nsize):
    temp1 = hlt[i]
    for j in range(i + 1):
        temp2 = hlt[j]
        c = FindCommonNumber(temp1, temp2)
        a = len(temp1)
        b = len(temp2)
        HLTSimilarity[i, j] = (c + 0.0) / (a + b - c)
        HLTSimilarity[j, i] = (c + 0.0) / (a + b - c)

MEDDRASimilarity = (SOCSimilarity + HLGTSimilarity + HLTSimilarity) / 3.0

##########################################################

f = file("~/adr/data/result/SOCSimilarity", 'w')
cPickle.dump(SOCSimilarity, f)
f.close()

f = file("~/adr/data/result/HLGTSimilarity", 'w')
cPickle.dump(HLGTSimilarity, f)
f.close()

f = file("~/adr/data/result/HLTSimilarity", 'w')
cPickle.dump(HLTSimilarity, f)
f.close()

f = file("~/adr/data/result/MEDDRASimilarity", 'w')
cPickle.dump(MEDDRASimilarity, f)
f.close()

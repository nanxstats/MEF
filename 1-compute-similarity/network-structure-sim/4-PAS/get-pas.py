# -*- coding: utf-8 -*-
"""
This script is used for calculating the PAS scores.

Created on Fri Sep 27 15:03:14 2013

@author: Dong-Sheng Cao
"""

import scipy
import cPickle
import copy
import networkx as nx

path="~/MEF/3-generate-negative-set/"

############################################################################

f=file(path+"association1",'r')
passociation=cPickle.load(f)
f.close()

f=file(path+"decoy/decoy1",'r')
nassociation=cPickle.load(f)
f.close()

tassociation=copy.deepcopy(passociation)
tassociation.extend(nassociation)

G=nx.Graph()
G.add_edges_from(passociation)

networkfeature=scipy.zeros((49606,1))
for index,i in enumerate(tassociation):
    if scipy.mod(index,100)==0:
        print index
    drug=i[0]
    adr=i[1]
    drugnei=G.neighbors(drug)
    adrnei=G.neighbors(adr)
    degreeprod=round(G.degree(drug)*G.degree(adr),5)
#    degreesum=round(G.degree(drug)+G.degree(adr),4)
#    degreeratio=round(G.degree(drug)/(G.degree(adr)+0.0),4)
#    degreediff=round(abs(G.degree(drug)-G.degree(adr)),4)
    networkfeature[index]=degreeprod

################################################################

f=file(path+"pas",'w')
cPickle.dump(networkfeature,f)
f.close()

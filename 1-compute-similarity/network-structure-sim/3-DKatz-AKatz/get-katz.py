# -*- coding: utf-8 -*-
"""
This script is used for computing the Katz index based on the
drug-ADR interaction network.

Created on Fri Oct 25 17:42:54 2013

@author: Dong-Sheng Cao
"""

##################################################################

import networkx as nx
import scipy
import cPickle
from scipy import linalg

# beta=0.00175

f=file("~/MEF/1-compute-similarity/network-structure-sim/3-DKatz-AKatz/association1",'r')
passociation=cPickle.load(f)
f.close()

passociation=[['d'+str(i[0]),'a'+str(i[1])] for i in passociation]

##################################################################

G = nx.Graph()
G.add_edges_from(passociation)

drugset,adrset=nx.algorithms.bipartite.sets(G)
drugG=nx.algorithms.bipartite.projected_graph(G,drugset)
adrG=nx.algorithms.bipartite.projected_graph(G,adrset)

adjmat=nx.adjacency_matrix(G)

##################################################################

identity=scipy.diag(scipy.repeat(1,1563))

a=linalg.eig(adjmat)
eigdrug=scipy.amax(a[0])

# drugpath=nx.katz_centrality(drugG,eigdrug)
beta1=0.7*1/eigdrug.real

katz=linalg.inv(identity-beta1*adjmat)-identity

drugsim2=scipy.zeros((746,746))
adrsim2=scipy.zeros((817,817))

for index1,i in enumerate(G.nodes()):
    for index2,j in enumerate(G.nodes()):
        if i[0]=="d" and j[0]=='d':
            ii=int(i[1:])
            jj=int(j[1:])
            drugsim2[ii-1,jj-1]=katz[index1,index2]
        if i[0]=="a" and j[0]=='a':
            ii=int(i[1:])
            jj=int(j[1:])
            adrsim2[ii-1,jj-1]=katz[index1,index2]

##################################################################

path1="~/MEF/result/"

f=file(path1+"drugkatz1",'w')
cPickle.dump(drugsim2,f)
f.close()

f=file(path1+"adrkatz1",'w')
cPickle.dump(adrsim2,f)
f.close()

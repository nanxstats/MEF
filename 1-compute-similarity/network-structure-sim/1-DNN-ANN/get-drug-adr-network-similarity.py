# -*- coding: utf-8 -*-
"""
This script is used for computing the network neighbor-based
similarity based on the underlying drug-ADR interaction network.

Created on Wed Jul 31 10:30:01 2013

@author: Dong-Sheng Cao
"""

import scipy
import cPickle

############################################################

f=file("~/adr/data/result/AdjacentMatrix",'r')
AdjacentMatrix=cPickle.load(f)
f.close()

(Nsize,Msize)=AdjacentMatrix.shape

############################################################

DrugNetworkSimilarity=scipy.zeros((Nsize,Nsize))
for i in range(Nsize):
    temp1=AdjacentMatrix[i,:]
    for j in range(i+1):
        temp2=AdjacentMatrix[j,:]
        c=sum(temp1*temp2)
        a=sum(temp1)
        b=sum(temp2)
        sim=(c+0.0)/(a+b-c)
        DrugNetworkSimilarity[i,j]=sim
        DrugNetworkSimilarity[j,i]=sim

f=file("~/adr/data/result/DrugNetworkSimilarity",'w')
cPickle.dump(DrugNetworkSimilarity,f)
f.close()

############################################################

ADRNetworkSimilarity=scipy.zeros((Msize,Msize))
for i in range(Msize):
    temp1=AdjacentMatrix[:,i]
    for j in range(i+1):
        temp2=AdjacentMatrix[:,j]
        c=sum(temp1*temp2)
        a=sum(temp1)
        b=sum(temp2)
        sim=(c+0.0)/(a+b-c)
        ADRNetworkSimilarity[i,j]=sim
        ADRNetworkSimilarity[j,i]=sim

f=file("~/adr/data/result/ADRNetworkSimilarity",'w')
cPickle.dump(ADRNetworkSimilarity,f)
f.close()

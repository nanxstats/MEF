# -*- coding: utf-8 -*-
"""
This script is used to compute the ADR coexist similarity based on
the adjacent matrix from drug-ADR associations.

Created on Wed Jul 31 11:11:34 2013

@author: Dong-Sheng Cao
"""

import scipy
import cPickle
from scipy import linalg

f = file("~/adr/data/result/AdjacentMatrix", 'r')
AdjacentMatrix = cPickle.load(f)
f.close()

(Nsize, Msize) = AdjacentMatrix.shape

ADRCoexist = scipy.zeros((Msize, Msize))
for i in range(Msize):
    temp1 = AdjacentMatrix[:, i]
    for j in range(Msize):
        temp2 = AdjacentMatrix[:, j]
        a = sum(temp1 * temp2)
        sim = a / sum(temp1)
        ADRCoexist[i, j] = sim

# ADRCoexistSimilarity=scipy.zeros((Msize,Msize))

ADRCoexistSimilarity = scipy.dot(ADRCoexist, ADRCoexist.T)

P1 = scipy.diag(ADRCoexistSimilarity)
P1 = scipy.diag(P1)

pp = linalg.inv(scipy.sqrt(P1))
temp = scipy.dot(pp, ADRCoexistSimilarity.T)
ADRCoexistSimilarity = scipy.dot(temp, pp.T)

f = file("~/adr/data/result/ADRCoexistSimilarity", 'w')
cPickle.dump(ADRCoexistSimilarity, f)
f.close()

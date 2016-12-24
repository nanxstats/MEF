# -*- coding: utf-8 -*-
"""
This script is used for transforming similarity measures to classification features
based on the neighborhood-based collaborative filtering recommendation methods.

Created on Fri Sep 27 10:16:55 2013

@author: Dong-Sheng Cao
"""

import scipy
import cPickle

##################################################################################

def GetDrugMEFFeature(Association, DrugSimilarityMatrix, AdjacentMatrix, k=1):
    """
	The function is used for generating drug-based recommendation scores
	as classification features.

	Input:
		Association: the drug-ADR association list
		DrugSimilarityMatrix: 746 x 746 drug similarity matrix from one kind of evidence
		AdjacentMatrix: 746 x 817 adjacent matrix for the drug-ADR network
    """
    res = []
    for aa, i in enumerate(Association):
        drugid = i[0]
        targetid = i[1]
        temp0 = DrugSimilarityMatrix[drugid, :]
        index2 = scipy.where(AdjacentMatrix[:, targetid] == 1)[0]
        index2 = list(index2)
        if drugid in index2:
            index2.remove(drugid)
        if len(index2) >= 1:
            kk = scipy.minimum(k, len(index2))
            temp = temp0[index2]
            index = temp.argsort()
            cc = len(index) - kk
            index1 = index[cc:]
            simvalue = temp[index1]
            res.append(scipy.mean(simvalue))
        elif len(index2) == 0:
            res.append(0)

    return res

##################################################################################

def GetTargetMEFFeature(Association, TargetSimilarityMatrix, AdjacentMatrix, k=1):
    """
	The function is used for generating ADR-based recommendation scores
	as classification features.

	Input:
		Association: the drug-ADR association list
		DrugSimilarityMatrix: 817 x 817 ADR similarity matrix from one kind of evidence
		AdjacentMatrix: 746 x 817 adjacent matrix for the drug-ADR network
    """
    res = []
    for i in Association:
        drugid = i[0]
        targetid = i[1]
        temp0 = TargetSimilarityMatrix[targetid, :]

        index2 = scipy.where(AdjacentMatrix[drugid, :] == 1)[0]
        index2 = list(index2)

        if targetid in index2:
            index2.remove(targetid)

        if len(index2) >= 1:
            kk = scipy.minimum(k, len(index2))
            temp = temp0[index2]
            index = temp.argsort()
            cc = len(index) - kk
            index1 = index[cc:]
            simvalue = temp[index1]
            res.append(scipy.mean(simvalue))
        elif len(index2) == 0:
            res.append(0)

    return res

##################################################################################

if "__main__" == __name__:
    path = "~/MEF/4-multiple-evidence-fusion/"
    ############################################################################

    f = file(path + "association1", 'r')
    passociation = cPickle.load(f)
    f.close()

    f = file(path + "decoy/decoy1", 'r')
    nassociation = cPickle.load(f)
    f.close()

    tassociation = copy.deepcopy(passociation)
    tassociation.extend(nassociation)

    path1 = "~/MEF/result/"
    ############################################################################
    # molecular structural
    f = file(path1 + "adrkatz", 'r')
    DrugECFPSimilarity = cPickle.load(f)
    f.close()

    tt = []
    for i in tassociation:
        tt.append([i[0] - 1, i[1] - 1])
    label1 = scipy.repeat(1, 24803)
    label2 = scipy.repeat(0, 24803)
    label = scipy.append(label1, label2)

    #################################################################3
    f = file(path1 + "AdjacentMatrix", 'r')
    AdjacentMatrix = cPickle.load(f)
    f.close()
    #################################################################3
    from sklearn import metrics

    res = GetTargetMEFFeature(tt, DrugECFPSimilarity, AdjacentMatrix, k=3)

    fpr, tpr, thresholds = metrics.roc_curve(label, res)
    print metrics.auc(fpr, tpr)

# -*- coding: utf-8 -*-
"""
This script is used for computing the structure-based similarity by using
ECFP4 and topological fingerprints. In our study, we only used the ECFP4
fingerprints to compute the fingerprints.

Created on Wed Jul 31 09:45:01 2013

@author: Dong-Sheng Cao
"""

import string
import scipy
import cPickle
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem

##########################################################################

drugsmi = []
f = file("~/MEF/1-compute-similarity/node-attribute-sim/Drug-Drug-Sim/1-ECFP/drugkegg", 'r')
drugkegg = cPickle.load(f)
f.close()

drugsmi = [i[4] for i in drugkegg]

drugsmi[220] = "OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)O[C@]([H])(C)CNC(=O)CC[C@]1(C)[C@@H](CC(=O)N)[C@@]2([H])N([Co]C#N)\C1=C(C)/C1=N/C(=C\C3=N\C(=C(C)/C4=N[C@]2(C)[C@@](C)(CC(=O)N)[C@@H]4CCC(=O)N)\[C@@](C)(CC(=O)N)[C@@H]3CCC(=O)N)/C(C)(C)[C@@H]1CCC(=O)N)N1C=NC2=CC(C)=C(C)C=C12"
drugsmi[396] = '[H]O[Co+]N1\C2=C(C)/C3=N/C(=C\C4=N\C(=C(C)/C5=N[C@@](C)([C@@]1([H])[C@H](CC(=O)N)[C@@]2(C)CCC(=O)NC[C@@H](C)OP(=O)([O-])O[C@H]1[C@@H](O)[C@H](O[C@@H]1CO)N1C=NC2=CC(C)=C(C)C=C12)[C@@](C)(CC(N)=O)[C@@H]5CCC(=O)N)\[C@@](C)(CC(=O)N)[C@@H]4CCC(=O)N)/C(C)(C)[C@@H]3CCC(=O)N'
molblock = []
for index, i in enumerate(drugsmi):
    print index
    mol = Chem.MolFromSmiles(i)
    molblock.append(mol)

##########################################################################

fps = [FingerprintMols.FingerprintMol(x) for x in molblock]

# fps=[]
# for index,i in enumerate(molblock):
#    print index+1
#    if index+1==221 or index+1==397:
#        continue
#    else:
#        FingerprintMols.FingerprintMol(i)

print len(fps)

DNum = len(fps)
TopoSimilarity = scipy.zeros((DNum, DNum))
for i in range(DNum):
    for j in range(i + 1):
        sim = DataStructs.FingerprintSimilarity(fps[i], fps[j])
        TopoSimilarity[i, j] = sim
        TopoSimilarity[j, i] = sim

f = file("~/adr/data/result/DrugTopoSmilarity", 'w')
cPickle.dump(TopoSimilarity, f)
f.close()

##########################################################################

fps = [AllChem.GetMorganFingerprint(x, 2) for x in molblock]

print len(fps)

DNum = len(fps)
ECFPSimilarity = scipy.zeros((DNum, DNum))
for i in range(DNum):
    for j in range(i + 1):
        sim = DataStructs.DiceSimilarity(fps[i], fps[j])
        ECFPSimilarity[i, j] = sim
        ECFPSimilarity[j, i] = sim

f = file("~/adr/data/result/DrugECFPSimilarity", 'w')
cPickle.dump(ECFPSimilarity, f)
f.close()

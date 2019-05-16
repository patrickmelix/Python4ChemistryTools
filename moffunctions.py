#Some functions to help analysing DUT-8

import numpy as np
import sys
from ase import Atoms, neighborlist
import os
import pickle

def buildNL(mol, path='./', radii=None, save=True):
   #create nl
   if radii is None:
       radii = {}
       radii[ 'H'] = 0.30
       radii[ 'C'] = 0.77
       radii[ 'N'] = 0.70
       radii[ 'O'] = 0.66
       radii['Ni'] = 1.24
   nAtoms = len(mol)

   if (not os.path.isfile(os.path.join(path, 'neighborList.pickle'))) or (not save):
       #create a list of cutoffs
       cutOff = []
       for j in range(0,nAtoms):
          cutOff.append(radii[mol[j].symbol])
       #initiate neighborlist
       neighborList = neighborlist.NeighborList(cutOff,self_interaction=False,bothways=True)
       neighborList.update(mol)
       if save:
           with open(os.path.join(path, 'neighborList.pickle'),'wb') as f:
                pickle.dump(neighborList,f)
   elif save:
       with open(os.path.join(path, 'neighborList.pickle'),'rb') as f:
            neighborList = pickle.load(f)
   print("Bond Map created")
   return neighborList

def getBenzNiAngleList(ana, imI):
    """Make a list of all occurences"""
    if isinstance(imI, int):
        imI = [imI]
    dihedralList = []
    for i in imI:
        dihedralList.append([])
        nAtoms = len(ana.images[i])
        molecule = ana.images[i]
        bondList = ana.all_bonds[i]
        relDihedrals = ana.get_dihedrals('Ni', 'Ni', 'O', 'C')[i]
        for dihed in relDihedrals:
            nextC = [ idx for idx in bondList[dihed[-1]] if molecule[idx].symbol == 'C']
            assert len(nextC) == 1
            attachedCs = [ idx for idx in bondList[nextC[0]] if (molecule[idx].symbol == 'C') and (idx != dihed[-1])]
            assert len(attachedCs) == 2
            dists = [ molecule.get_distance(dihed[0], idx, mic=True) for idx in attachedCs ]
            if dists[0] < dists[1]:
                dihedralList[-1].append((dihed[0], dihed[1], attachedCs[1], attachedCs[0]))
            else:
                dihedralList[-1].append((dihed[0], dihed[1], attachedCs[0], attachedCs[1]))

    return dihedralList

def getAlphaList(ana, imI):
    """Make a list of all occurences"""
    from itertools import combinations
    allIdx = ana._get_symbol_idxs(imI, 'Ni')
    nAtoms = len(ana.images[imI])
    molecule = ana.images[imI]
    bondList = ana.all_bonds[imI]
    dihedrals = ana.get_dihedrals('Ni', 'O', 'C', 'C')[0]
    alphaList = []
    for niIdx in allIdx:
        relDihedrals = [ d for d in dihedrals if d[0] == niIdx ]
        assert len(relDihedrals) == 4
        s = set([d[-1] for d in relDihedrals])
        if s not in alphaList:
            alphaList.append(s)

    r = []
    for l in alphaList:
        r.extend(list(combinations(list(l), 2)))

    return [r]


def checkBonds(mol, bondList, cluster=False):
    for iAtom, atom in enumerate(mol):
        bondedIdx = bondList[iAtom]
        bondedSym = [ mol[idx].symbol for idx in bondedIdx ]
        if atom.symbol == 'C':
            check = [ False, False, False, False ]
            if (set(bondedSym) == set(['C','H'])) and (len(bondedIdx) == 3):
                check[0] = True
            elif (set(bondedSym) == set(['C','H','N'])) and (len(bondedIdx) == 4):
                check[1] = True
            elif (set(bondedSym) == set(['C'])) and (len(bondedIdx) == 3):
                check[2] = True
            elif (set(bondedSym) == set(['C','O'])) and (len(bondedIdx) == 3):
                check[3] = True
            if not any(check):
                raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
        elif atom.symbol == 'H':
            if ((set(bondedSym) != set(['C'])) and (set(bondedSym) != set(['N']))) or (len(bondedIdx) != 1):
                raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
        elif atom.symbol == 'N':
            if cluster:
                if ((set(bondedSym) != set(['C','Ni'])) and (set(bondedSym) != set(['H','Ni'])) and (set(bondedSym) != set(['C']))) or (not len(bondedIdx) in [4,3]):
                    raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
            else:
                if ((set(bondedSym) != set(['C','Ni'])) and (set(bondedSym) != set(['H','Ni']))) or (len(bondedIdx) != 4):
                    raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
        elif atom.symbol == 'O':
            if (set(bondedSym) != set(['C','Ni'])) or (len(bondedIdx) != 2):
                raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
        elif atom.symbol == 'Ni':
            if (set(bondedSym) != set(['Ni','O','N'])) or (len(bondedIdx) != 6):
                raise RuntimeError("Atom {:} ({:}) bonded to {:}, these are {:}".format(iAtom, atom.symbol,str(bondedIdx),str(bondedSym)))
    return True

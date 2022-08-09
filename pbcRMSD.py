#!/usr/bin/env python3
#
# Script to revert atom wrapping in a trajectory.
# by Patrick Melix
# 2022/05/06
#
# You can import the module and then call .main() or use it as a script
from ase.io import read
from ase.geometry import find_mic
import os, glob
import numpy as np

def main(inList):
    mols = []
    for iFile,inFile in enumerate(inList):
        if not os.path.isfile(inFile):
            raise ValueError('File {:} does not exist'.format(inFile))

        print("Reading {:}".format(inFile))
        tmp = read(inFile, index=slice(0,None))
        if isinstance(tmp, list):
            mols.extend(tmp)
        else:
            mols.append(tmp)
        if iFile == 0:
            refCell = mols[0].get_cell()
            refAtoms = mols[0]

    print("Total number of frames: {}".format(len(mols)))
    for i,mol in enumerate(mols):
        if i == 0:
            continue
        if not np.allclose(mol.get_cell(), mols[i-1].get_cell()):
            print("Cell of frame {:} is not the same as the cell before! Results can be useful if cell changes are small enough.".format(i))
        assert len(mol) == len(refAtoms), "Number of atoms in frame {:} is not the same as the reference atoms!".format(i)
        for iAtom in range(len(mol)):
            micVec, _ = find_mic(mol[iAtom].position - refAtoms[iAtom].position, refCell)
            micVec -= mol[iAtom].position - refAtoms[iAtom].position
            if np.linalg.norm(micVec) > 1e-3:
                print("Wrapping atom {:} in frame {:} back.".format(iAtom, i))
                mol[iAtom].position += micVec
        distances = mol.positions - mols[i-1].positions
        distances = np.linalg.norm(distances, axis=1)
        assert len(distances) == len(mol)
        distances *= distances
        rmsd = np.sqrt(np.sum(distances)/len(mol))
        print("RMSD from frame {:} to frame {:}: {:}".format(i-1, i, rmsd))

    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Print RMSD considering PBC.')
    parser.add_argument('files', type=str, nargs='+', help='input file(s)')
    args = parser.parse_args()
    main(args.files)



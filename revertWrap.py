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

def main(inList,outFile):
    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

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
        assert np.allclose(mol.get_cell(), refCell), "Cell of frame {:} is not the same as the reference cell!".format(i)
        assert len(mol) == len(refAtoms), "Number of atoms in frame {:} is not the same as the reference cell!".format(i)
        for iAtom in range(len(mol)):
            micVec, _ = find_mic(mol[iAtom].position - refAtoms[iAtom].position, refCell)
            micVec -= mol[iAtom].position - refAtoms[iAtom].position
            if np.linalg.norm(micVec) > 1e-3:
                print("Wrapping atom {:} in frame {:} back.".format(iAtom, i))
                mol[iAtom].position += micVec



    with open(outFile,'w') as f:
        for mol in mols:
            mol.write(f,format='extxyz')

    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Remove atom wrapping from trajectory.')
    parser.add_argument('files', type=str, nargs='+', help='input file(s)')
    parser.add_argument('-output', type=str, help='output file')
    args = parser.parse_args()
    main(args.files, args.output)



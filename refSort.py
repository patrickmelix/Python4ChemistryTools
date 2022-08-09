#!/usr/bin/env python3
from ase import io, Atoms, Atom
import numpy as np

def main(refFile, inFile, outFile, delta=0.1):
    refMol = io.read(refFile)

    sortMol = io.read(inFile)

    sortList = []
    for iAtom in range(len(refMol)):
        tmpAtoms = sortMol.copy()
        tmpAtoms.extend(refMol[iAtom])
        distances = tmpAtoms.get_distances(len(tmpAtoms)-1, [i for i in range(len(sortMol))], mic=True)
        #print(distances)
        tmp = np.where(distances < delta)
        if len(tmp) > 1:
            print("More than one possible atom found for idx {}".format(iAtom))
            sortList.append(float('nan'))
            continue
        minIdx = np.argmin(distances)
        if distances[minIdx] > delta:
            print("Did not found atom that is closer than {}A for idx {}, element {}".format(delta,iAtom,refMol[iAtom].symbol))
            sortList.append(float('nan'))
        else:
            #print("Unsorted index {} closest to ref index {}".format(minIdx, iAtom))
            sortList.append(minIdx)

    sortedMol = []
    notFound = [idx for idx in range(len(sortMol)) if idx not in sortList]
    #print(sortList)
    for iAtom in sortList:
        if np.isnan(iAtom):
            continue
        else:
            sortedMol.append(Atom(symbol=sortMol[iAtom].symbol, position=sortMol[iAtom].position))

    final = Atoms(symbols=[a.symbol for a in sortedMol], positions=[a.position for a in sortedMol], cell=sortMol.cell, pbc=sortMol.pbc)
    for idx in notFound:
        final.extend(Atom(symbol=sortMol[idx].symbol, position=sortMol[idx].position))
    assert len(final) == len(sortMol)
    final.write(outFile)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Try to sort an ASE structure file based on a reference file.')
    parser.add_argument('-ref', type=str, help='Reference Structure')
    parser.add_argument('-input', type=str, help='Input File Name')
    parser.add_argument('-out', type=str, help='Output File Name')
    parser.add_argument('--dist', type=float, help='Maximum allowed atomic distance.', default=0.1)
    args = parser.parse_args()
    main(args.ref, args.input, args.out, args.dist)

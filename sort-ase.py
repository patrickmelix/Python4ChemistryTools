#!/usr/bin/env python3
from ase import io
import numpy as np

def main(input,output):
    mol = io.read(input)
    newMol = mol.copy()
    del newMol[:]
    allElements = list(set(mol.symbols))
    allElements.sort()
    for element in set(mol.symbols):
        atoms = mol[np.where(mol.symbols == element)[0]]
        newMol.extend(atoms)
    newMol.write(output)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Sort Atoms by Element')
    parser.add_argument('input', type=str, help='path to file')
    parser.add_argument('output', type=str, help='path to file')
    args = parser.parse_args()
    main(args.input, args.output)

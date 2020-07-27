#!/usr/bin/env python3
# creates an ASE xyz file from a cp2k restart file
import os
from ase import io

def main(inFile,outFile):
    if not os.path.isfile(inFile):
        raise ValueError('File {:} does not exist'.format(str(inFile)))

    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    cell = []
    atoms = []

    with open(inFile, 'r') as f:
        for line in f:
            if '&CELL' in line:
                while True:
                    l = f.readline().strip()
                    if l.startswith(('A', 'B', 'C')):
                        cell.append([ s for s in l.split()[1:] ])
                    if '&END' in l:
                        break
            if '&COORD' in line:
                while True:
                    l = f.readline().strip()
                    if '&END' in l:
                        break
                    atoms.append(l)

    with open(outFile,'w') as f:
        f.write("{:}\n".format(len(atoms)))
        f.write("\n")
        for atom in atoms:
            f.write("{}\n".format(atom))
        for i, vec in enumerate(cell):
            f.write("VEC{} {} {} {}\n".format(i+1, *vec))

    #test
    mol = io.read(outFile)
    mol.write(outFile)





if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert Gaussian output to ASE-extxyz trajectory')
    parser.add_argument('input', type=str, help='input restart file')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    main(args.input,args.output)



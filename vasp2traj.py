#!/usr/bin/python3
#
# Script to convert VASP output to ASE-extxyz trajectory
# by Patrick Melix
# 2018/03/13
#
# You can import the module and then call .main() or use it as a script
from ase import io
import os

def main(inFile,outFile):
    if not os.path.isfile(inFile):
        raise ValueError('File {:} does not exist'.format(str(inFile)))

    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    mol = io.read(inFile, format='vasp-out', index=slice(0,None))
    for frame in mol:
        frame.wrap(center=(0.0,0.0,0.0))
        frame.write(outFile,append=True)
    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert VASP output to ASE-extxyz trajectory')
    parser.add_argument('input', type=str, help='input xyz file')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    main(args.input,args.output)



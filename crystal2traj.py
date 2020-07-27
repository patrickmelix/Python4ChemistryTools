#!/usr/bin/env python3
#
# Script to convert Crystal GO to ASE-extxyz trajectory
# by Patrick Melix
# 2018/04/16
#
# You can import the module and then call .main() or use it as a script
from ase import io
import os, glob

def main(inStr,outFile):
    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    files = glob.glob(os.path.join(inStr,'*'))
    files.sort()

    for inFile in files:
        if not os.path.isfile(inFile):
            raise ValueError('File {:} does not exist'.format(str(inFile)))

        print(inFile)
        mol = io.read(inFile, format='crystal')
        mol.write(outFile,append=True)
    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert Crystal GO to ASE-extxyz trajectory')
    parser.add_argument('input', type=str, help='optstory dir')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    main(args.input,args.output)



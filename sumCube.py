#!/usr/bin/python3
#
# Script to sum cube files
# by Patrick Melix
# 2018/06/07
#
# You can import the module and then call .main() or use it as a script
from ase.io.cube import read_cube, write_cube
import os, glob
import numpy as np

def main(inList,outFile, subtract=False):
    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    cubeData = []

    for inFile in inList:
        if not os.path.isfile(inFile):
            raise ValueError('File {:} does not exist'.format(inFile))

        print(inFile)
        with open(inFile) as f:
            dct = read_cube(f)
        cubeData.append(dct['data'].copy())
        atoms = dct['atoms']
        origin = dct['origin']

    if subtract:
        sumData = cubeData[0]
        for i in range(1, len(cubeData)):
            sumData = np.subtract(sumData, cubeData[i])
    else:
        sumData = np.sum(cubeData, axis=0)
    with open(outFile,'w') as f:
        write_cube(f, atoms, data=sumData, origin=origin)

    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Sum up cube files')
    parser.add_argument('--subtract', help='Subtract not Sum', action='store_true')
    parser.add_argument('-file', type=str, nargs='+', default=[], help='optstory dir')
    parser.add_argument('-output', type=str, help='output file')
    args = parser.parse_args()
    main(args.file,args.output, args.subtract)



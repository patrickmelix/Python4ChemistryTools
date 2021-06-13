#!/usr/bin/env python3
#
# Script to generate an ext-xyz trajectory from individual ext-xyz files with varying atom numbers using ASE.
# Appens 'X' atoms at origin to obtain frames with equal lengths
# by Patrick Melix
# 2020/06/08
#

from ase import io, Atom
import os

def main(inList, outFile='traj.xyz', outFormat='extxyz'):
    #if output exists mv to .bak
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    traj = []

    for inFile in inList:
        if not os.path.isfile(inFile):
            raise ValueError('File {:} does not exist'.format(inFile))

        print(inFile)
        traj.append(io.read(inFile))

    maxLen = max([len(frame) for frame in traj])

    for i in range(len(traj)):
        if len(traj[i]) < maxLen:
            for j in range(maxLen-len(traj[i])):
                traj[i].append(Atom('X'))

    with open(outFile,'w') as f:
        for frame in traj:
            frame.write(f, format=outFormat)
    return


#########################
# Functions
########################




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Combine different lengths of XYZ')
    parser.add_argument('--outformat', help='Output ASE Format', default='extxyz')
    parser.add_argument('--outfile', help='Output File', default='traj.xyz')
    parser.add_argument('-files', type=str, nargs='+', default=[], help='All the XYZ Files')
    args = parser.parse_args()
    main(args.files, args.outfile, args.outformat)



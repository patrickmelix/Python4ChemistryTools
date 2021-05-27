#!/usr/bin/env python3
#
# Script to convert SCM AMV trajectory to xyz format trajectory
# by Patrick Melix
# 2020/10
#
# You can import the module and then call .main() or use it as a script
import sys, os, glob
from ase import io

def main(argv):
    inFile = argv[0]
    outFile = argv[1]

    data = []
    nAtoms = None
    with open(inFile) as f:
        for line in f.readlines():
            if line.strip() == '':
                continue
            if 'Geometry' in line:
                data.append([])
            data[-1].append(line)

    with open(outFile,'w') as f:
        for frame in data:
            f.write(str(len(frame)-1)+'\n')
            for i, entry in enumerate(frame):
                if i == 0:
                    f.write("'"+entry.strip()+"'\n")
                else:
                    f.write(entry)


if __name__ == "__main__":
    if '-h' in sys.argv[1]:
        print("Usage: amv2xyz.py <infile> <outfile>")
        sys.exit(0)
    main(sys.argv[1:])



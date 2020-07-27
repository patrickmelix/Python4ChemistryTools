#!/usr/bin/env python3
#
# Script to convert ASE readable structure to ADF readable xyz format
# by Patrick Melix
# 2018/02/13
#
# You can import the module and then call .main() or use it as a script
import sys, os, glob
from ase import io

def main(argv):
    inFiles = []
    for arg in argv:
        tmp = glob.glob(arg)
        for f in tmp:
            if not os.path.isfile(f):
                raise ValueError('File {:} does not exist'.format(os.path.basename(f)))
            inFiles.append(f)

    for inFile in inFiles:
        outFile = os.path.splitext(os.path.basename(inFile))[0] + '.xyz'
        inExt = os.path.splitext(os.path.basename(inFile))[1]
        if inExt == '.xyz':
            outFile = outFile[:-4] + '_adf' + outFile[-4:]
        #if output exists mv to .bak
        if os.path.isfile(outFile):
            print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
            os.rename(outFile, outFile+'.bak')

        inMol = io.read(inFile,index=slice(0,None))

        if type(inMol) != type([]):
            inMol = [inMol]

        for frame in inMol:
            frame.write(outFile,vec_cell=True,append=True)



if __name__ == "__main__":
    if '-h' in sys.argv[1]:
        print("Usage: ase2adf.py <file1> <file2> <file*>")
        sys.exit(0)
    main(sys.argv[1:])



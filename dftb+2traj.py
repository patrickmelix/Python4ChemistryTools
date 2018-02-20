#!/usr/bin/python3
#
# Script to convert ASE readable structure to ADF readable xyz format
# by Patrick Melix
# 2018/02/13
#
# You can import the module and then call .main() or use it as a script
import os
from ase import io
from ase.io.dftb import read_dftb_lattice

def main():
    if not os.path.isfile('geo_end.xyz'):
        raise ValueError('File geo_end.xyz does not exist')
    if not os.path.isfile('md.out'):
        raise ValueError('File md.out does not exist')

    #if output exists mv to .bak
    outFile = 'traj.xyz'
    if os.path.isfile(outFile):
        print('ATTENTION: {:} exists, moving to *.bak'.format(outFile))
        os.rename(outFile, outFile+'.bak')

    mol = io.read('geo_end.xyz', index=slice(0,None))
    read_dftb_lattice(images=mol)
    for frame in mol:
        frame.wrap()
        frame.write(outFile,append=True)



if __name__ == "__main__":
    main()



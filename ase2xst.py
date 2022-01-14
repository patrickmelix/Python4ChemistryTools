#!/usr/bin/env python3
#
# Script to print cell information from an ASE compatible file into an XST file
# by Patrick Melix
# 2021/11/22
#
from ase import io, Atoms
import os

def main(inFile, outFile):
    if not os.path.isfile(inFile):
        raise ValueError('File {:} does not exist'.format(str(inFile)))

    traj = io.read(inFile, index=slice(0,None))
    if isinstance(traj, Atoms):
        traj = [traj]
    assert isinstance(traj, list), "Given file does not contain a trajectory!"

    print("Length of Trajectory: {} Frames".format(len(traj)))
    with open(outFile, 'w') as out:
        out.write("# NAMD extended system trajectory file\n#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z\n")
        for i,frame in enumerate(traj):
            cell = [ x for vec in frame.cell for x in vec ]
            out.write(("{:5} "+"{:10f}"*9+"\n").format(i, *cell))
    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert ASE Compatible Input into XST file.')
    parser.add_argument('input', type=str, help='input file')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    main(args.input, args.output)



#!/usr/bin/env python3
#
# Script to print cell information from an ase trajectory
# by Patrick Melix
# 2021/07/15
#
from ase import io
import os

def main(inFile):
    if not os.path.isfile(inFile):
        raise ValueError('File {:} does not exist'.format(str(inFile)))

    traj = io.read(inFile, index=slice(0,None))
    assert isinstance(traj, list), "Given file does not contain a trajectory!"

    print("Length of Trajectory: {} Frames".format(len(traj)))
    print(("{:10}"*10).format("Volume","A1","A2","A3","B1","B2","B3","C1","C2","C3"))
    for frame in traj:
        cell = [ x for vec in frame.cell for x in vec ]
        v = frame.get_volume()
        print(("{:10f}"*10).format(v, *cell))
    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Print Vell Information from ASE Compatible Input')
    parser.add_argument('input', type=str, help='input file')
    args = parser.parse_args()
    main(args.input)



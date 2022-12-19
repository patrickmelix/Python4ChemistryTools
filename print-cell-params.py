#!/usr/bin/env python3
#
# Script to print cell parameters from an ase trajectory
# by Patrick Melix
# 2022/10/27
#
from ase import io
import os


def main(inFile):
    if not os.path.isfile(inFile):
        raise ValueError('File {:} does not exist'.format(str(inFile)))

    traj = io.read(inFile, index=slice(0, None))
    assert isinstance(traj, list), "Given file does not contain a trajectory!"

    print("Length of Trajectory: {} Frames".format(len(traj)))
    print(("{:10}"*6).format("A", "B", "C", "alpha", "beta", "gamma"))
    for frame in traj:
        cell = [x for x in frame.cell.cellpar()]
        print(("{:10f}"*6).format(*cell))
    return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Print Vell Information from ASE Compatible Input')
    parser.add_argument('input', type=str, help='input file')
    args = parser.parse_args()
    main(args.input)



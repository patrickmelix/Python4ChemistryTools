#!/usr/bin/env python3
#
# Script to calculate the integral of one or multiple cube file
# by Patrick Melix
# 2022/04/04
#
# You can import the module and then call .main() or use it as a script
import re
from ase.io.cube import read_cube
import os, glob
import numpy as np

def main(inFiles, return_value=False, volume=False):
    if return_value:
        values = []
    for iFile,inFile in enumerate(inFiles):
        if not os.path.isfile(inFile):
            raise ValueError('File {:} does not exist'.format(inFile))
        print("Reading {}".format(inFile))
        with open(inFile) as f:
            dct = read_cube(f)
            integral = np.sum(np.abs(dct['data']))
            if volume:
                integral /= dct['atoms'].get_volume()
            if return_value:
                values.append(integral)
            else:
                print("Integral of {} is {}".format(inFile, integral))

    if return_value:
        if len(values) == 1:
            return values[0]
        else:
            return values
    else:
        return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Sum up cube files')
    parser.add_argument('files', type=str, nargs='+', help='input files')
    parser.add_argument('--volume', help='Devide the Data by the Cell Volume', action='store_true')
    args = parser.parse_args()
    main(args.files, return_value=False, volume=args.volume)

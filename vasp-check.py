#!/usr/bin/env python3
#
# Script to assert proper occupations in VASP calculation using ASE
# by Patrick Melix
# 2021/06/15
#
from ase.calculators.vasp.vasp import Vasp
import os
import numpy as np

def main(path):
    assert os.path.isdir(path), "Given path is not a directory"
    calc = Vasp(directory=path)
    xml = calc._read_xml()
    if xml.get_spin_polarized():
        spins = [0, 1]
        electrons = 1.0
    else:
        spins = [0]
        electrons = 2.0

    nkpoints = len(xml.get_ibz_k_points())
    for s in spins:
        for i in range(nkpoints):
            occ = xml.get_occupation_numbers(i, s)
            if occ is None:
                print("No occupations found in vasprun.xml for kpoint #{} and spin {}!".format(i,s))
                return
            test = np.where(np.logical_or(occ == electrons, occ == 0.0), 1, 0)
            if not test.all():
                print("Bad Occupation found")
                return
    print("Seems like there are no bad occupations (only last step is checked).")

    if not calc.read_convergence():
        print("Either SCF or GO did not converge!")
    else:
        print("No convergence issues found (only last step).")
    return



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Check VASP run for proper occupations')
    parser.add_argument('path', type=str, help='path to VASP files', default='./')
    args = parser.parse_args()
    main(args.path)



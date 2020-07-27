#!/usr/bin/env python
#
# Read Output of a CP2K BSSE Calculation and Calculate the Interaction Energy.
# Only for 2 Component System A-B with 5 Subcalculations in one Output.

import os, shutil
from scm.plams import *

def main(name, assign):
    try:
        if isinstance(config['default_jobmanager'], JobManager):
            print('PLAMS already loaded')
        else:
            raise
    except:
        if os.path.isdir('tmp.plams'):
            shutil.rmtree('tmp.plams')
        init(folder='tmp.plams')


    job = Cp2kJob.load_external(name)
    assert job.results.check_scf()

    struct = toASE(Cp2kSettings2Mol(job.settings))

    molInfo = job.results.grep_output(pattern="MOLECULE KIND INFORMATION", options="-A 9")
    nAtoms = molInfo[9::11]
    nAtoms = [ int(s.split()[-1]) for s in nAtoms ]
    nKinds = molInfo[8::11]
    nKinds = [ int(s.split()[-1]) for s in nKinds ]
    #first entry is the generic non-bsse input
    del nKinds[0], nAtoms[0]

    electronInfo = job.results.grep_output(pattern="Number of electrons:", options="")
    nElectrons = [ int(s.split()[-1]) for s in electronInfo ]
    #spin A and B
    nElectrons = [ (nElectrons[i],nElectrons[i+1]) for i in range(0,len(nElectrons),2) ]

    energies = job.results.grep_output(pattern="SCF run converged in", options="-A 15")[15::17]
    energies = [ float(s.split()[-1]) for s in energies ]

    print("Total Number of Atoms: {}".format(len(struct)))
    print("Electrons: {}".format(nElectrons))
    print("Kinds of Atoms: {}".format(nKinds))
    print("Number of Atoms: {}".format(nAtoms))
    print("Energies: {}".format(energies))

    print(assign)


    interaction = energies[assign["AB"]]-energies[assign["A"]]-energies[assign["B"]]
    print("Interaction Energy NOT Corrected: {}".format(interaction))
    interaction_corr = energies[assign["AB"]]-energies[assign["AB_ghost"]]-energies[assign["A_ghostB"]]
    print("Interaction Energy Corrected: {}".format(interaction_corr))
    bsse = energies[assign["AB_ghost"]]-energies[assign["A"]]+energies[assign["A_ghostB"]]-energies[assign["B"]]
    print("BSSE: {}".format(bsse))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Analyze CP2K BSSE Calculation')
    parser.add_argument('input', type=str,  help='Directory of the CP2K BSSE Calculation')
    parser.add_argument('--assign', type=int, nargs=5, help='Order of Energies. Default: 4 0 1 2 3', default=[4,0,1,2,3])
    args = parser.parse_args()

    assign = {}
    for i,key in enumerate(["AB", "A", "B", "AB_ghost", "A_ghostB"]):
         assign[key] = args.assign[i]

    main(args.input, assign)

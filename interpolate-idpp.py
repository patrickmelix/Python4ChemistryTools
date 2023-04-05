#!/usr/bin/env python3

from ase.neb import NEB, idpp_interpolate
#from ase.calculators.emt import EMT
#from ase.calculators.lj import LennardJones
from ase import io
from xtb.ase.calculator import XTB



def main(startFile, finalFile, outputFile, nImages=5):
    # Optimise molecule.
    initial = io.read(startFile)

    # Create final state.
    final = io.read(finalFile)

    # Generate blank images.
    images = [initial]

    for i in range(nImages):
        images.append(initial.copy())

    for image in images:
        #image.calc = LennardJones()
        image.calc = XTB(method='GFNFF', max_iterations=1) #GFNFF, GFN0-xTB

    images.append(final)

    # Run IDPP interpolation.
    neb = NEB(images)
    #neb.interpolate(mic=True) #'idpp',
    idpp_interpolate(neb, mic=True)

    with open(outputFile,'w') as f:
        for i,frame in enumerate(images):
                frame.write(f,format='extxyz')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Interpolate between two structures using IDPP')
    parser.add_argument('start', type=str, help='start structure')
    parser.add_argument('final', type=str, help='final structure')
    parser.add_argument('output', type=str, help='output structure')
    parser.add_argument('--nimages', type=int, help='Number of images', default='5')
    args = parser.parse_args()
    main(args.start, args.final, args.output, nImages=args.nimages)

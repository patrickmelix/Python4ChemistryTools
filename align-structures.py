#!/usr/bin/env python3

from ase import io
from ase.build import minimize_rotation_and_translation




def main(inputFile, refFile, outputFile):
    # align this molecule.
    initial = io.read(inputFile)

    # Create reference state.
    ref = io.read(refFile)

    new = initial.copy()
    minimize_rotation_and_translation(ref, new)

    new.write(outputFile)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Align two structures')
    parser.add_argument('input', type=str, help='input structure')
    parser.add_argument('ref', type=str, help='reference structure')
    parser.add_argument('output', type=str, help='output structure')
    args = parser.parse_args()
    main(args.input, args.ref, args.output)

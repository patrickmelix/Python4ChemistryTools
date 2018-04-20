#!/usr/bin/python3
#
# Script to generate a supercell trajectory from an ase trajectory
# by Patrick Melix
# 2018/04/20
#
# You can import the module and set the input variables by hand and then call .main() or use the built-in STDInput or command line Arguments
# Vars to be set: trajFileIn, trajFileOut, superCell[3]

import sys, os
import numpy as np
from ase import io, Atoms

#variables
global trajFileIn
global trajFileOut
global superCell
superCell = []
trajFileIn = ''
trajFileOut = ''


def main():
    print('Welcome.')

    #global variables
    global superCell
    global trajFileIn
    global trajFileOut

    #get Input from stdin
    if len(sys.argv) is 1:
        trajFileIn, trajFileOut, superCell = getSTDInput()

    #if we have the right amount of arguments
    elif len(sys.argv) is 6:
        print('Provided information as arguments.')
        trajFileIn, trajFileOut, superCell = getArgumentsInput(sys.argv)

    #otherwise something is wrong
    else:
        print('You gave some arguments, but the number of them is wrong. Please provide "<InputFile> <OutFile> <n1> <n2> <n3>"')
        sys.exit(1)

    #check input
    checkInput(trajFileIn, superCell)

    #open files
    inMol = io.read(trajFileIn,index=slice(0,None))

    #if output exists mv to .bak
    if os.path.isfile(trajFileOut):
        print('ATTENTION: Output file exists, moving to *.bak')
        os.rename(trajFileOut, trajFileOut+'.bak')

    #open output for writing
    outFile = open(trajFileOut, 'w')

    #convert input superCell to integer tuple
    try:
        superCell = tuple([int(x) for x in superCell])
    except:
        print('ERROR: SuperCell Vector contains non-integer values!')
        sys.exit(1)

    #check input superCell vector
    if any(x for x in superCell) < 1:
        print('ERROR: SuperCell Tuple contains values smaller than 1!')
        sys.exit(1)

    #calculate number of unit cells
    superCellSize = np.prod(np.array(superCell), dtype='int')
    print('You will end up with '+str(superCellSize)+' unit cells in each frame.')

    #iterate over frames
    nFrames = len(inMol)
    i = 0
    print('Adding Coordinates:')
    for frame in inMol:
        i += 1
        progress(i/nFrames)
        frame *= superCell

    print('Writing to file:')
    i = 0
    for frame in inMol:
        i += 1
        progress(i/nFrames)
        io.extxyz.write_xyz(outFile,frame)

    outFile.close()


#########################
# Functions
########################

#get input from stdin
def getSTDInput():
    #read input trajectory path
    trajFileIn = input('xyz-Trajectory input file: ')
    trajFileOut = input('xyz-Trajectory output file: ')
    #vectors
    supercell = []
    print('Enter now the Supercell vector as integers, where every element gives the number of unit cells in this direction (format: <x> <y> <z>)')
    superCell = input('Vector: ').split()
    return trajFileIn, trajFileOut, superCell

#get input from arguments
def getArgumentsInput(args):
    superCell = []
    trajFileIn = args[1]
    trajFileOut = args[2]
    superCell = args[3:6]
    return trajFileIn, trajFileOut, superCell

#check input for errors
def checkInput(fileIn, superCell):
    #check path
    if not os.path.isfile(fileIn):
        print('ERROR: Input file does not exist?!')
    if len(superCell) is not 3:
        print('ERROR: Length of supercell vector must be three.')
        sys.exit(1)

#A progress bar
def progress(progress):
    length = 10
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress >= 1:
        progress = 1
        status = "Finished!\r\n"
    block = int(round(length*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "*"*block + "-"*(length-block), int(progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()



if __name__ == "__main__":
    main()



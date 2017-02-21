#!/usr/bin/python3
#
# Script to generate a supercell trajectory from an xyz trajectory (e.g. CP2K MD output)
# by Patrick Melix
# 2017/02/20
#
# You can import the module and set the input variables by hand and then call .main() or use the built-in STDInput or command line Arguments
# Vars to be set: trajFileIn, trajFileOut, vecIn[3*[3]], superCell[3]

import sys, os
import numpy as np

#variables
global vecIn
global trajFileIn
global trajFileOut
vecIn = []
superCell = []
trajFileIn = ''
trajFileOut = ''


def main():
    print('Welcome.')
    
    #global variables
    global vecIn
    global superCell
    global trajFileIn
    global trajFileOut
    
    #get arguments or ask for them
    if len(vecIn) is not 0:
        print('You set the input Variables somewhere else, I will try to use them.')
    
    #get Input from stdin
    elif len(sys.argv) is 1:
        trajFileIn, trajFileOut, vecIn, superCell = getSTDInput()
    
    #if we have the right amount of arguments
    elif len(sys.argv) is 15:
        print('Provided information as arguments.')
        trajFileIn, trajFileOut, vecIn, superCell = getArgumentsInput(sys.argv)
    
    #otherwise something is wrong
    else:
        print('You gave some arguments, but the number of them is wrong. Please provide "<InputFile> <a1> <a2> <a3> <b1> <b2> <b3> <c1> <c2> <c3>"')
        sys.exit(1)

    #check input
    checkInput(trajFileIn, vecIn)
    
    
    #transform input vectors to numpy matrix
    matrixIn = np.asmatrix(vecIn, dtype='float')

    #open files
    inFile, outFile = openFiles(trajFileIn, trajFileOut)

    #calculate the inverse
    try:
        matrixInv = np.linalg.inv(matrixIn)
    except:
        print('Cannot calculate the inverse of your unit cell matrix, maybe you misstyped?')
        sys.exit(1)

    #iterate over frames and save info as fractional coords
    finished = False
    atomVec = []
    frames = [] #save all atoms of this frame in fractional coords
    lineN = 0
    #get number of lines
    nLines = rawpycount(trajFileIn)
    print('Number of lines in trajectory: '+str(nLines))
    print('Reading input trajectory and converting to fractional coordinates:')
    while not finished:
        #read first line of input
        nAtoms = inFile.readline().strip()
        #if is empty last frame has been processed
        if nAtoms == '':
            print('EOF reached.')
            finished = True
            break
        else:
            nAtoms = int(nAtoms)
        #line of last line to be read
        lineN += nAtoms+2

        #update progress bar
        progress(lineN/nLines)

        #does the file end before that line?
        if lineN > nLines:
            print('ERROR: Input file ends earlier than expected!')
            sys.exit(1)

        #read comment line
        comment = inFile.readline().rstrip('\n')

        thisFrame = []
        thisFrame.append(nAtoms)
        thisFrame.append(comment)
        for iAtom in range(0, nAtoms):
            line = inFile.readline().split()
            atomVec = np.array(line[1:], dtype='float')
            element = line[0]
            atomVec = np.dot(matrixInv, atomVec).tolist()[0]

            #check if fractional coords are between 0 and 1
            #disabled since unit cell wrapping is not obligatory
            #if any(x > 1.0 for x in atomVec) or any(x < 0.0 for x in atomVec):
            #    print('ERROR: Your cell vectors do not cover all coordinates (at least one atom is outside your unit cell)')
            #    print(len(frames))
            #    print(element)
            #    print(atomVec)
            #    sys.exit(1)
            thisFrame.append([element] + atomVec)
       
        #add this frame to all frames
        frames.append(thisFrame)

        #update progress bar
        progress(lineN/nLines)

        #eof reached?
        if lineN is nLines:
            print('EOF reached.')
            finished = True

    #convert input superCell to integer vector
    try:
        superCell = [int(x) for x in superCell]
    except:
        print('ERROR: SuperCell Vector contains non-integer values!')
        sys.exit(1)

    #check input superCell vector
    if any(x for x in superCell) < 1:
        print('ERROR: SuperCell Vector contains values smaller than 1!')
        sys.exit(1)
    
    #convert to numpy array
    superCell = np.array(superCell)

    #calculate number of unit cells
    superCellSize = np.prod(superCell, dtype='int')
    print('You will end up with '+str(superCellSize)+' unit cells in each frame.')

    #iterate over frames
    nFrames = len(frames)
    i = 0
    print('Adding Coordinates:')
    for frame in frames:
        i += 1
        progress(i/nFrames)
        #iterate over the three dimensions
        for dim in range(0,3):
            #iterate over all atoms
            for atomVec in frame[2:]:
                #iterate over the unitcells we need to add in this dimension
                for cellN in range(1,superCell[dim]):
                    tmpVec = atomVec[1:]
                    tmpVec[dim] += cellN
                    frame.append([atomVec[0]] + tmpVec)
    
    #now all frames should have all fractional coords
    #convert the fractional coord to real space
    print('Converting to real space and writing to file:')
    i = 0
    for frame in frames:
        i += 1
        progress(i/nFrames)
        nAtoms = len(frame)-2
        outFile.write(str(nAtoms)+'\n')
        outFile.write(str(frame[1])+'\n')
        for atomVec in frame[2:]:
            atomVec = [atomVec[0]] + np.dot(matrixIn, atomVec[1:]).tolist()[0]
            string = '   '.join(str(e) for e in atomVec)
            outFile.write(string+'\n')

    inFile.close()
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
    vecIn = []
    print('Please enter now your three unit-cell vectorsi (format: <a1> <a2> <a3>):')
    vecIn.append(input('A: ').split())
    vecIn.append(input('B: ').split())
    vecIn.append(input('C: ').split())
    supercell = []
    print('Enter now the Supercell vector as integers, where every element gives the number of unit cells in this direction (format: <x> <y> <z>)')
    superCell = input('Vector: ').split()
    return trajFileIn, trajFileOut, vecIn, superCell

#get input from arguments
def getArgumentsInput(args):
    vecIn = []
    superCell = []
    trajFileIn = args[1]
    trajFileOut = args[2]
    vecIn[0:2] = [args[3:6], args[6:9], args[9:12]]
    superCell = args[12:15]
    return trajFileIn, trajFileOut, vecIn, superCell

#check input for errors
def checkInput(fileIn, vec):
    #check path
    if not os.path.isfile(fileIn):
        print('ERROR: Input file does not exist?!')
        sys.exit(1)
    #check Vectors
    for vector in vec:
        if len(vector) is not 3:
            print('ERROR: Each vector has to have three dimensions!')
            sys.exit(1)
    #check superCell vector
    if len(superCell) is not 3:
        print('ERROR: Length of supercell vector must be three.')
        sys.exit(1)


#open input and output file
def openFiles(fileIn, fileOut):
    #if output exists mv to .bak
    if os.path.isfile(fileOut):
        print('ATTENTION: Output file exists, moving to *.bak')
        os.rename(fileOut, fileOut+'.bak')

    #open input for reading
    inFile = open(fileIn, 'r')

    #open output for writing
    outFile = open(fileOut, 'w')

    return inFile, outFile

#fast way of getting number of lines
def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)
def rawpycount(filename):
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    nLines = sum( buf.count(b'\n') for buf in f_gen )
    f.close()
    return nLines

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



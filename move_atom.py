#!/usr/bin/env python3
#
# moves a chosen atom from a ase-compatible input around and saves all the coordinate files as xyz
#
# by Patrick Melix
# 2017/03/17
#
# Usage:
# a) Call this script, answer the questions, done.
# b) Call this script with follogin arguments: <deltaX> <deltaY> <deltaZ> <nX> <nY> <nZ> <atomNumber> <input.xyz>
# c) Load the module and set the four global variables and call .main()
#

import sys, os
from ase import Atoms, io
import numpy as np

#global variables
global delta
global nSteps
global atomNumber
global xyzFile
delta = []
nSteps = []
atomNumber = 0
xyzFile = ''


def main():
   print('Hello!')

   #global Vars
   global delta
   global nSteps
   global atomNumber
   global xyzFile

   #input from variable
   if len(delta) is not 0:
      print('You provided input through global variables, trying to use them')

   #input from argument
   elif len(sys.argv) is 9:
      print('Provided input through arguments, trying to use them')
      delta, nSteps, atomNumber, xyzFile = getArgumentsInput(sys.argv)

   #read input from stdin
   elif len(sys.argv) is 1:
      delta, nSteps, atomNumber, xyzFile = getArgsSTDInput()

   #ERROR
   else:
      print('Input could not be gathered, error!')
      sys.exit(1)

   #read xyz
   atoms = io.read(xyzFile)

   #get atom of interest
   atomIdx = int(atomNumber)-1
   atomOfInterest = atoms[atomIdx]

   #loop over x,y,z
   xyz = ['x', 'y', 'z']
   for direct in range(3):
      i = 0
      vec = np.array([0, 0, 0], dtype='float')
      vec[direct] = float(delta[direct])
      newAtomsPlus = atoms.copy()
      newAtomsMinus = atoms.copy()

      #add delta
      while i < int(nSteps[direct]):
         i += 1
         newAtomsPlus[atomIdx].position += vec
         newAtomsMinus[atomIdx].position -= vec
         filenamePlus = xyz[direct] + '+' + str(i) + '.xyz'
         filenameMinus = xyz[direct] + '-' + str(i) + '.xyz'

         #write files
         io.write(filenamePlus, newAtomsPlus)
         io.write(filenameMinus, newAtomsMinus)

   print("Finished!")




#########################
# Functions
#########################

def getArgumentsInput(args):
   delta = args[1:4]
   nSteps = args[4:7]
   atomNumber = args[7]
   xyzFile = args[8]
   return delta, nSteps, atomNumber, xyzFile

def getArgsSTDInput():
   delta = input('Give me three delta values for x, y and z: ')
   nSteps = input('Give me three numbers for the number of steps in each direction for each x, y and z: ')
   atomNumber = input('Which atom should be moved (number in the xyz: ')
   xyzFile = input('Give me the name of the xyz file: ')
   return delta, nSteps, atomNumber, xyzFile



if __name__ == "__main__":
   main()

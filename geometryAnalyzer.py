import numpy as np
from ase import io, Atoms

def convertNames(inFile, outFile):
    xyz = open(outFile,'w')
    with open(inFile) as f:
        for line in f:
            split = line.strip().split()
            if len(split) < 4:
                xyz.write(line)
                continue
            name = split[0]
            split[1:] = [float(x) for x in split[1:]]
            if name[0].lower() == 'c':
                split[0] = 'C'
            elif name[0].lower() == 'o':
                split[0] = 'O'
            elif name[0].lower() == 'h':
                split[0] = 'H'
            elif name[0:2].lower() == 'ni':
                split[0] = 'Ni'
            elif name[0].lower() == 'n':
                split[0] = 'N'
            xyz.write(("{:10} "+"{:20.6f} "*3+"\n").format(*split))
    xyz.close()

def getBonds(A,B,inMol,bondList):
    #get all bonds A-B
    if not isinstance(inMol, list):
        mols = [ inMol ]
    else:
        mols = inMol
    allBonds = []
    for molecule in mols:
        nAtoms = len(molecule)
        bonds = []
        allIdx = []
        for i in range(0,nAtoms):
            if molecule[i].symbol == A:
                allIdx.append(i)
        for iIdx in allIdx:
            try:
                ibonds = bondList[str(iIdx)]
            except:
                continue
            for bonded in ibonds:
                if not molecule[bonded].symbol == B:
                    continue
                bonds.append([iIdx, bonded])
        #delete duplicates if A=B
        if A == B:
            for bond in bonds:
                del bonds[bonds.index(list(reversed(bond)))]
        allBonds.extend(bonds)
    return allBonds

def getAngles(A,B,C,inMol,bondList):
    #get all angles B-A-C
    if not isinstance(inMol, list):
        mols = [ inMol ]
    else:
        mols = inMol
    allAngles = []
    for molecule in mols:
        nAtoms = len(molecule)
        angles = []
        allIdx = []
        for i in range(0,nAtoms):
            if molecule[i].symbol == A:
                allIdx.append(i)
        for iIdx in allIdx:
            try:
                bonds = bondList[str(iIdx)]
            except:
                continue
            for bonded in bonds:
                if not molecule[bonded].symbol == B:
                    continue
                for j in bonds:
                    if j == bonded:
                        continue
                    if molecule[j].symbol == C:
                        angles.append([bonded, iIdx, j])
        #delete duplicates if B=C
        if B == C:
            for angle in angles:
                del angles[angles.index(list(reversed(angle)))]
        allAngles.extend(angles)
    return allAngles


def getDihedrals(A,B,C,D,molecule,bondList):
    """Make a list of all Dihedrals"""
    dihedralList = []
    allIdx = []
    nAtoms = len(molecule)
    for i in range(0,nAtoms):
        if molecule[i].symbol == A:
            allIdx.append(i)
    for idx in allIdx:#A
        try:
            ibonds = bondList[str(idx)]
        except:
            continue
        for j in bondList[str(idx)]:
            if not molecule[j].symbol == B:
                continue
            for k in bondList[str(j)]:
                if not molecule[k].symbol == C:
                    continue
                if idx == k:
                    continue
                for l in bondList[str(k)]:
                    if not molecule[l].symbol == D:
                        continue
                    if (not l == k) and (not l == idx) and (not l == j):
                        dihedralList.append([idx,j,k,l])
    return dihedralList

def bond(v1,v2):
    """ Returns the length of the vector.  """
    return np.linalg.norm(np.subtract(v1,v2))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle(v1, v2):
    """ Angle in Degree"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))/np.pi*180

def dihedral(vec1,vec2,vec3):
    """ Dihedral in Degree"""
    v1_u = unit_vector(vec1)
    v2_u = unit_vector(vec2)
    v3_u = unit_vector(vec3)
    #get the two normal vectors standing on the planes
    v1v2 = np.cross(v1_u,v2_u)
    v2v3 = np.cross(v2_u,v3_u)
    #angle between them is the dihedral
    return angle(v1v2,v2v3)

def getPBCVector(staticVec, vec, box, cut=5.0):
    #find new pbcVec using PBC so that pbcVec-staticVec is less then 5A away
    #test 6 most propable directions first
    pbcVec = np.subtract(vec,staticVec)
    for i in range(0,3):
        for j in [-1,1]:
            newVec = np.add(vec,box[i]*j)
            newVec = np.subtract(newVec,staticVec)
            if np.linalg.norm(newVec) < cut:
                return newVec
    #if not yet exited, perhaps it is one of the boxes on the edges
    #there are eight of them
    for dim in range(0,3):
        dims = list(range(0,3))
        dims.remove(dim)
        for i in [-1,1]:
            for j in [-1,1]:
                translate = np.add(box[dims[0]]*i,box[dims[1]]*j)
                newVec = np.add(vec,translate)
                newVec = np.subtract(newVec,staticVec)
                if np.linalg.norm(newVec) < cut:
                    return newVec
    #check the corner-connected boxes
    for i in [-1,1]:
        for j in [-1,1]:
            for k in [-1,1]:
                translate = np.add(box[0]*i,box[1]*j)
                translate = np.add(translate,box[2]*k)
                newVec = np.add(vec,translate)
                newVec = np.subtract(newVec,staticVec)
                if np.linalg.norm(newVec) < cut:
                    return newVec
    #if there is no result yet something is wrong
    raise ValueError("No matching PBC point found!")
    
    

def getBondValues(inMol,bondLists):
    if not isinstance(inMol, list):
        mols = [ inMol ]
    else:
        mols = inMol
    bonds = {}
    for name in bondLists:
        bonds[name] = []
    for molecule in mols:
        for name, bondList in bondLists.items():
            for item in bondList:
                bonds[name].append(molecule.get_distance(item[0],item[1],mic=True))
    return bonds

def getAngleValues(inMol,angleLists):
    if not isinstance(inMol, list):
        mols = [ inMol ]
    else:
        mols = inMol
    angles = {}
    for name in angleLists:
        angles[name] = []
    for molecule in mols:
        for name, angleList in angleLists.items():
            for item in angleList:
                angles[name].append(molecule.get_angle(item[0],item[1],item[2],mic=True))
    return angles


def getDihedralValues(inMol, dihedralLists):
    if not isinstance(inMol, list):
        mols = [ inMol ]
    else:
        mols = inMol
    dihedrals = {}
    for name in dihedralLists:
        dihedrals[name] = []
    for molecule in mols:
        for name, dihedralList in dihedralLists.items():
            for item in dihedralList:
                dihedrals[name].append(molecule.get_dihedral(item[0],item[1],item[2],item[3],mic=True))
    return dihedrals


#!/usr/bin/python3
#
# Script to convert a DL_POLY HISTORY file to a XYZ-trajectory without loosing information
# Cell Vectors are put in the comment line in an ASE compatible form
# by Patrick Melix
# 2017/07/08
#
# You can import the module and use the functions or call it as a script



def main():
    print('Welcome.')
    atomNames = {}
    atomNames['ca'] = 'C'
    atomNames['cb'] = 'C'
    atomNames['cc'] = 'C'
    atomNames['cd'] = 'C'
    atomNames['ce'] = 'C'
    atomNames['c1'] = 'C'
    atomNames['cn'] = 'C'
    atomNames['co'] = 'C'
    atomNames['o1'] = 'O'
    atomNames['oc'] = 'O'
    atomNames['n1'] = 'N'
    atomNames['ns'] = 'N'
    atomNames['ni'] = 'Ni'
    atomNames['cu'] = 'Cu'
    atomNames['h1'] = 'H'
    atomNames['ha'] = 'H'
    atomNames['ho'] = 'H'
    atomNames['hn'] = 'H'
    convertHistory2XYZ('HISTORY','traj.xyz',atomNames)
    print('Done!')
    
def convertHistory2XYZ(inFile, outFile, atomNames):
    bunchsize = 1000000     # Experiment with different sizes 1000000
    bunch = []
    print('Processing...')
    with open(inFile, "r", bunchsize) as r, open(outFile, "w", bunchsize) as w:
        #skip header line
        next(r)
        line = r.readline().strip().split()
        levcfg,imcon,nAtoms = [int(x) for x in line[0:3]]
        print('levcfg: '+str(levcfg))
        print('imcon: '+str(imcon))
        print('nAtoms: '+str(nAtoms))
        print('Frame: 0...', end='\r')
        nItems = (levcfg+1)*3
        if imcon == 6:
            latticeBool = 'pbc="T T F"'
        elif imcon > 0:
            latticeBool = 'pbc="T T T"'
        properties = 'Properties=species:S:1:pos:R:3 '
        if levcfg == 1:
            properties = 'Properties=species:S:1:pos:R:3:vel:R:3 '
        elif levcfg == 2:
            properties = 'Properties=species:S:1:pos:R:3:vel:R:3:forces:R:3 '
        iLine = 0
        iFrame = 0
        for line in r:
            split = line.split()
            if split[0].lower() == 'timestep':
                bunch.append(str(nAtoms)+'\n')
                iFrame += 1
                print('Frame: '+str(iFrame)+'...', end='\r')
                tmp = []
                #timestep value
                #tmp.append(split[0])
                #tmp.append(split[1])
                if imcon > 0:
                    #three vectors
                    tmp.append('Lattice="')
                    for i in range(0,3):
                        tmp.append(("{:} "*3).format(*[float(x) for x in r.readline().strip().split()]))
                    tmp.append('" '+properties+latticeBool)
                tmp.append('\n')
                bunch.append(''.join(tmp))
            else:
                tmp = []
                for i in range(0,levcfg+1):
                    tmp.extend([float(x) for x in r.readline().strip().split()])
                name = atomNames[split[0].lower()]
                bunch.append(("{:4}"+("{:16f}")*nItems+'\n').format(name,*tmp))
        
            if len(bunch) == bunchsize:
                w.writelines(bunch)
                iLine += bunchsize
                bunch = []
        w.writelines(bunch)
            

if __name__ == "__main__":
    main()



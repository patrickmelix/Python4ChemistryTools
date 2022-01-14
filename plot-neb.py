#!/usr/bin/env python3
#
# Script to plot VASP+TST NEB calculation results
# by Patrick Melix
# 2022/01
#
# You can import the module and then call .main() or use it as a script
import argparse
import numpy as np
from matplotlib.ticker import MaxNLocator
exec(open("/home/patrickm/git/Python4ChemistryTools/mpl-settings.py").read())


def plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, forces, filename, lw=3, s=0, highlight=None):
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Reaction Coordinate')# [Å]
    plt.ylabel(r'$\Delta E$ [eV]')
    #plt.ylim([-10**exp,10**exp])
    #plt.yscale('symlog')
    #plt.gca().yaxis.grid(True)
    plt.plot(reactionCoord, energySpline, color='black', ls=':', label='Cubic Spline', lw=lw)
    plt.scatter(reactionCoordImageAxis, energies, marker='P', color='red', s=(msbig+s)**2, label='NEB Energy')
    dScale = 0.02
    maxX = max(reactionCoordImageAxis)
    delta = dScale*maxX
    for i,x in enumerate(reactionCoordImageAxis):
        tangentX = [x-delta, x+delta]
        tangentY = [energies[i]+(delta*forces[i]),energies[i]-(delta*forces[i])] #invert sign of forces from neb output
        if i == 0:
            label = 'NEB Force'
        else:
            label = None
        plt.plot(tangentX, tangentY, color='green', ls='-', lw=lw, label=label)
    if highlight is not None:
        plt.scatter(reactionCoordImageAxis[highlight], energies[highlight], marker='o', s=(msbig+s+30)**2, facecolors='none', edgecolors='orange', lw=lw+2, clip_on=False)
    #plt.xticks(x, printDirs[:], rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    #plt.show()
    plt.close()


def main(filename='NEB.png', presentation=False, highlight=None, plot_all=False):
    spline = np.loadtxt('spline.dat')
    print("Spline loaded.")
    nebData = np.loadtxt('neb.dat')
    print("Energy and forces loaded.")
    nImages = len(nebData)
    print("{:} data points found.".format(nImages))

    reactionCoord = [ s[1] for s in spline ]
    xAxis = [ s[0] for s in spline ]
    images = [ d[0] for d in nebData ]
    reactionCoordImageAxis = [ d[1] for d in nebData ]
    #reactionCoordImageAxis = [ reactionCoord[xAxis.index(n)] for n in images ]
    energies = [ d[2] for d in nebData ]
    forces = [ d[3] for d in nebData ]
    energySpline = [ s[2] for s in spline ]

    if presentation:
        lw = 5
        s = 3
        plt.rcParams.update({'font.size': 22})
        plt.rcParams.update({'legend.fontsize': 22})
    else:
        lw = 3
        s = 0
    plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, forces, filename, lw=lw, s=s, highlight=highlight)

    if plot_all:
        #plot the main image and then one with every point highlighted
        filename = filename.split('.')
        filename[-2] += "-{:02d}"
        filename = ".".join(filename)
        for i in range(nImages):
            plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, forces, filename.format(i), lw=lw, s=s, highlight=i)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot VASP+TST NEB results')
    parser.add_argument('--file', help='Plot Filename', default='NEB.png')
    parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')
    parser.add_argument('--highlight', help='Circle Point N', type=int, default=None)
    parser.add_argument('--plotall', help='Create main plot and each highlighted plot.', action='store_true')
    args = parser.parse_args()
    main(args.file, args.presentation, args.highlight, args.plotall)

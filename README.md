# Python4ChemistryTools
My personal Python toolbox for theoretical chemistry. Some Bash scripts also sneaked in.

## Trajectories and File Conversion
- ase2adf.py: Convert any type of ASE readable input file (or files, including asterisk expressions) given as arguments to ADF style XYZ files.
- ase2xst.py: Convert any type of ASE readable input file given as argument to XST file (e.g. to read cell info in VMD).
- aseSupercell.py: Read ASE compatible input and save a supercell xyz-Trajectory file.
- cp2krestart2xyz: Extract single xyz file from CP2K restart/input files and save them as an ASE ext-xyz file. 
- crystal2traj.py: Convert Crystal geometry optimization output to a ASE compatible ext-xyz trajectory file.
- dftb+2traj.py: Extract trajectory including periodicity from DFTB+ MD simulation.
- gaussian2traj.py: Convert Gaussian geometry output to ASE compatible ext-xyz trajectory file.
- print-traj-cell.py: Prints volume and cell vectors from an ASE compatible trajectory.
- supercell.py: Generate Supercell from xyz-Trajectory file, can be used as module or standalone.
- vasp2traj.py: Convert VASP geometry optimization output to ASE compatible ext-xyz trajectory file.

## Structure Manipulation and Analysis
- move_atom.py: Move a single atom from an ASE-compatible input around in space and save all resulting coordinates as xyz. Can be used as a script or module.
- geometryAnalyzer: Some functions for geometric analysis of ASE-format molecules and trajectories. Copy functions into your code or import the file.

## Other
- cp2k-bsse.py: Analyze output of a CP2K BSSE calculation with two components using PLAMS.
- dothemath.py: Script to process all math expressions in a text file. Useful for many things, e.g. processing manually written force field definitions.
- moffunctions.py: Some functions for analyzing MOFs, pretty specific to DUT-8.
- mpl-settings.py: My favorite Matplotlib settings.
- plams_defaults: Options for PLAMS to work on clusters.
- sumCube.py: Make the sum of an arbitrary number of cube files. Also supports subtraction.
- vasp-check.py: Assert proper occupations and SCF+GO convergence in VASP using ASE.

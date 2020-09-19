# Python4ChemistryTools
My personal Python toolbox for theoretical chemistry.

- supercell.py: Generate Supercell from xyz-Trajectory file, can be used as module or standalone.
- aseSupercell.py: Read ASE compatible input and save a supercell xyz-Trajectory file.
- move_atom.py: Move a single atom from an ASE-compatible input around in space and save all resulting coordinates as xyz. Can be used as a script or module.
- geometryAnalyzer: Some functions for geometric analysis of ASE-format molecules and trajectories. Copy functions into your code or import the file.
- plams_defaults: Options for PLAMS to work on clusters.
- ase2adf.py: Convert any type of ASE readable input file (or files, including asterisk expressions) given as arguments to ADF style XYZ files.
- dftb+2traj.py: Extract trajectory including periodicity from DFTB+ MD simulation.
- cp2krestart2xyz: Extract single xyz file from CP2K restart/input files and save them as an ASE ext-xyz file. 
- cp2k-bsse.py: Analyze output of a CP2K BSSE calculation with two components using PLAMS.
- crystal2traj.py: Convert Crystal geometry optimization output to a ASE compatible ext-xyz trajectory file.
- gaussian2traj.py: Convert Gaussian geometry output to ASE compatible ext-xyz trajectory file.
- vasp2traj.py: Convert VASP geometry optimization output to ASE compatible ext-xyz trajectory file.
- sumCube.py: Make the sum of an arbitrary number of cube files. Also supports subtraction.
- moffunctions.py: Some functions for analyzing MOFs, pretty specific to DUT-8.
- mpl-settings.py: My favorite Matplotlib settings.
- dothemath.py: Script to process all math expressions in a text file. Useful for many things, e.g. processing manually written force field definitions.

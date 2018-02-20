# Python4ChemistryTools
My personal Python toolbox for theoretical chemistry.

- supercell.py: Generate Supercell from xyz-Trajectory file, can be used as module or standalone.
- move_atom.py: Move a single atom from an ASE-compatible input around in space and save all resulting coordinates as xyz. Can be used as a script or module.
- geometryAnalyzer: Some functions for geometric analysis of ASE-format molecules and trajectories. Copy functions into your code or import the file.
- plams_defaults: Options for PLAMS to work on clusters.
- ase2adf.py: Convert any type of ASE readable input file (or files, including asterisk expressions) given as arguments to ADF style XYZ files.
- dftb+2traj.py: Extract trajectory including periodicity from DFTB+ MD simulation.

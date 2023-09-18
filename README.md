# mdanalysis_protein_cation_pi_bond_detection
A custom class using MDAnalysis to detect cation pi bonds for a protein in solvent.

## List of Files and Descriptions
- **cation_pi_analysis_demo.ipynb** : Jupyter notebook which demonstrates the usage of the function 'CationPiBondAnalysis' defined in cation_pi_analysis.py
- **cation_pi_analysis.py** : Python script containing the custom MDAnalysis class 'CationPiBondAnalysis.'
- **demo_trajectory/ionized.pdb** : Protein data bank file (pdb) containing atom information for protein, water, and ions. See [Introduction to Protein Data Bank Format](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) for more details on pdb file format.
- **demo_trajectory/ionized.psf** : Protein structure file (psf) containing atom and connectivity information for protein, water, and ions See [PSF Files](https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html) for details on psf file format.
- **demo_trajectory/trajectory.dcd** : Molecular dynamics simulation trajectory of strutre in ionized.pdb.
- **mda_cation_pi_env.yml** : Conda YAML file (yml) which can be used to generate a conda environment with MDAnalysis and other packages required to run 'cation_pi_analysis_demo.ipynb' and 'cation_pi_analysis.py.' See [Creating an environment from an environment.yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) in the conda documentation.
- 
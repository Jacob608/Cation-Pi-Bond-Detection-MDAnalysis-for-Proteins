# Cation Pi Bond Detection with MDAnalysis for Proteins
A custom class using MDAnalysis to detect cation pi bonds for a protein in solvent.

## List of Files and Descriptions
- **cation_pi_analysis_demo.ipynb** : Jupyter notebook which demonstrates the usage of the function 'CationPiBondAnalysis' defined in cation_pi_analysis.py
- **cation_pi_analysis.py** : Python script containing the custom MDAnalysis class 'CationPiBondAnalysis.'
- **demo_trajectory/ionized.pdb** : Protein data bank file (pdb) containing atom information for protein, water, and ions. See [Introduction to Protein Data Bank Format](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) for more details on pdb file format.
- **demo_trajectory/ionized.psf** : Protein structure file (psf) containing atom and connectivity information for protein, water, and ions See [PSF Files](https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html) for details on psf file format.
- **demo_trajectory/trajectory.dcd** : Molecular dynamics simulation trajectory of strutre in ionized.pdb.
- **mda_cation_pi_env.yml** : Conda YAML file (yml) which can be used to generate a conda environment with MDAnalysis and other packages required to run 'cation_pi_analysis_demo.ipynb' and 'cation_pi_analysis.py.' See [Creating an environment from an environment.yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) in the conda documentation.

## Instructions

1. Create a conda virtual environment with the appropriate packages by following the instructions in [Creating an environment from an environment.yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
2. [Activate the environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) and follow step 3 at [this link](https://saturncloud.io/blog/how-to-add-a-python-3-kernel-to-jupyter-ipython/) to install the IPython Kernel package.
3. Register the envrionment as a kernel for use in Jupyter Notebook by following step 4 at [this link](https://saturncloud.io/blog/how-to-add-a-python-3-kernel-to-jupyter-ipython/).
4. Open **cation_pi_analysis_demo.ipynb** using Jupyter Notebook and run the script.

## Additional Information

For details regarding distance cutoff and output of class 'CationPiBondAnalysis,' please see the comments in **cation_pi_analysis.py**.

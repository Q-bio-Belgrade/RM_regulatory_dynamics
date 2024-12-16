# RM Regulatory Dynamics

This repository contains the Mathematica and Python code used to reproduce the results presented in the paper:

**Nonlinear Regulatory Dynamics of Bacterial Restriction-Modification Systems Modulates Horizontal Gene Transfer Susceptibility**

*Authors: Magdalena Djordjevic, Lidija Zivkovic, Hong-Yu OU, Marko Djordjevic*  
*Under review in: Nucleic Acids Research, 2024*  

## Repository Overview

The repository is organized to ensure full reproducibility of the results and figures presented in the paper. It includes Mathematica notebooks and Python simulation scripts structured as follows:

- `mathematica_notebooks/`: Contains Mathematica notebooks for symbolic calculations and analytical derivations.
- `python_simulations/`: Contains Python scripts for numerical simulations and data analysis.

---

## Development Environment

To ensure reproducibility, the code was developed and tested in the following environment:

- **Mathematica version**: 11.2  
- **Python version**: 3.9.19  
- **Conda version**: 24.5.0  

Python dependencies are managed through a conda environment.

---

## How to Use

### **Mathematica Notebooks**
- Located in the `mathematica_notebooks/` folder.
- Each notebook corresponds to a specific figure in the paper. For example, `figure1.nb` reproduces Figure 1.
- Notebooks include code, parameters, and results together for clarity.

### **Python Simulations**
- Python code, data, and results are in the `python_simulations/` folder.
- Scripts are located in the `simulation_scripts/` subfolder. Ensure scripts are executable:
  ```bash
  chmod u+x <script_name>
  ```
- Scripts should be run from the parent folder (`python_simulations/`):
  ```bash
  ./simulation_scripts/<script_name>
  ```
- Inputs and outputs are handled via the Python `argparse` module. To view options and usage, use the `-h` flag:
  ```bash
  ./simulation_scripts/simulation_CV_MR_regulated.py -h
  ```
- Default parameter values are set to reproduce paper results. Outputs are saved to the `simulation_results/` folder.

---

## Parallelization
Python simulation scripts are parallelized:
- Default: Uses all available CPUs (`-1`).
- To specify the number of CPUs, use the `-w` flag:
  ```bash
  ./simulation_scripts/simulation_CV_MR_regulated.py -w <number_of_workers>
  ```

---

## Reproducing Figures

### **Figures in the Main Text**
- **Figure 3**: `mathematica_notebooks/figure3.nb`
- **Figure 4**: `mathematica_notebooks/figure4.nb`
- **Figure 5**: `mathematica_notebooks/figure5.nb`
- **Figure 6**: `mathematica_notebooks/figure6_and_S2.nb`
- **Figure 7**: `mathematica_notebooks/figure7.nb`
- **Figure 8**: `mathematica_notebooks/figure8.nb`
- **Figure 9**: Run the following scripts in order:
  1. `figure9_simulation_CV_MR_constitutive.py`
  2. `figure9_simulation_CV_MR_regulated.py`
  3. `figure9_plotting.py`
  - Figure is saved in the `simulation_results/` folder.
- **Figure 10**: `mathematica_notebooks/figure10.nb`
### **Supplementary Figures**
- **Figure S1**:
  1. Run `figureS1_param_inference.py`
  2. Run `figureS1_R_dynamics_simulation.py`
  - Figure is saved in the `simulation_results/` folder.
- **Figure S2**: `mathematica_notebooks/figure6_and_S2.nb`
---

## License
This repository is licensed under MIT license. For more details, see the `LICENSE` file.

---

## Citation
If you use this code in your work, please cite the paper above.

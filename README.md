# Stratum Order-of-Addition Designs

This repository contains all necessary documentation required to replicate the results in:

**Stratum order-of-addition designs**  
by Liushan Zhou, Ze Liu, Min-Qian Liu, and Guanzhou Chen

## Description
 
The repository includes computer programs for:

- **Figure 1**: Position combinations in the first two columns of a COA(10,5,1)
- **Figure 2**: The pairwise plots of the columns of designs X and X' in Example 6
- **Table S5**: Evaluating OofA designs under various space-filling criteria
- **Table S7**: D-efficiencies of OofA designs with m components and N runs under various models
- **Table S8**: MSPEs under different models in the job scheduling problem
- **Table S9**: MSPEs under different models in the job scheduling problem
- **Table S10** The means and standard deviations of (M,S)-efficiencies under various models

## Details

### Figures 1 and 2

- `Fig1and2.sh` 
  Creates a "Rplots.pdf" file which includes **Figures 1 and 2**.


### Table S5

- `Table_S5_SOofA.R`  
  Computes the six space-filling criteria for the **SOofA** design in **Table S5**.

- `Table_S5_SX.R`  
  Computes the six space-filling criteria for the **SX** design in **Table S5**.

### Table S7

- `Table_S7_SX.R`  
  Computes the D-efficiencies of the **SX** design under the **FP**, **QP**, and **SCP** models in **Table S7**.

- `Table_S7_FPQP_SOofA.R`  
  Computes the mean and standard deviation of the D-efficiencies for the **SOofA** design under the **FP** and **QP** models in **Table S7**.

- `Table_S7_SCP_SOofA.R`  
  Computes the mean and standard deviation of the D-efficiencies for the **SOofA** design under the **SCP** model in **Table S7**.

### Table S8

- `TabS8.sh` 
  Generates the data in **Table S8**, with u = 1, 2, 4, 8, 16.

### Table S9

- `TabS9.sh` 
  Generates the data in **Table S9**, with N_train = 504, 72, 36, 24, 18, 9.

### Table S10

- `TabS10.sh` 
  Produces the numerical results of **Table S10**.
  

## Requirements

- The R code may require additional R packages available in the local R environment.
- The C++ code depends on the Eigen library
	(https://eigen.tuxfamily.org/index.php?title=Main_Page),
	with version = 3.3.9 (higher version may or may not work!).
	For convenience, a copy of Eigen-3.3.9 is provided in `program_2.tar.gz`.


## License

Please see the `LICENSE` file in this repository.  
If you plan to make the repository public, make sure the license text and copyright holder information are finalized before upload.

## Citation

If you use this code in academic work, please cite the corresponding article and its supplementary material.

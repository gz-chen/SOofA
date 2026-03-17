# Code for Table S5 and Table S7 in “Stratum order-of-addition designs”

This repository contains the R scripts used to reproduce the results reported in **Table S5** and **Table S7** of the supplementary material for the article:

**Stratum order-of-addition designs**  
Liushan Zhou, Ze Liu, Min-Qian Liu, and Guanzhou Chen

## Description

These are computer programs for producing selected supplementary results from the article above.  
The current repository covers the scripts for:

- **Table S5**: evaluating OofA designs under various space-filling criteria
- **Table S7**: D-efficiencies of OofA designs with different designs and models

## Contents

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

## Requirements

- **R**
- The scripts may require additional R packages available in the local R environment.

## How to run

Run the corresponding script in R or RStudio. For example:

```r
source("Table_S5_SOofA.R")
```

or

```r
source("Table_S7_SX.R")
```

## Notes

- The **code content has not been modified**; only the filenames were standardized so that they match the actual supplementary table numbers.
- This repository is intended to organize and archive the scripts corresponding to the published supplementary results.
- Please check your local R setup before execution.

## Output

These scripts are intended to reproduce the numerical results reported in:

- **Table S5**
- **Table S7**

## License

Please see the `LICENSE` file in this repository.  
If you plan to make the repository public, make sure the license text and copyright holder information are finalized before upload.

## Citation

If you use this code in academic work, please cite the corresponding article and its supplementary material.

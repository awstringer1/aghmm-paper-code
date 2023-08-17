# Approximate Marginal Likelihood Inference in Mixed Models for Grouped Data
Replication code for paper: *Approximate Marginal Likelihood Inference in Mixed Models for Grouped Data* by Alex Stringer.
All code was run on a 2021 M1 MacBook Pro with 10 cores and 64Gb of RAM. I did not test it on other hardware,
but I don't make use of any speciifc features of the Apple silicon either.

## Package

All scripts require the `aghqmm` package, available from GitHub [here](https://github.com/awstringer1/aghqmm).
You should download this repository to a directory on your machine.

## Data

The data are obtained as follows:

- **Smoking**: obtained from the online supplementary material to the paper *A note on marginalization of regression parameters from mixed models of binary outcomes* by Donald Hedeker, Stephen H. C. du Toit, Hakan Demirtas, and Robert D. Gibbons (*Biometrics* (2017) **74**(1) pp. 354–361). Find that paper [here](https://onlinelibrary.wiley.com/doi/10.1111/biom.12707). Data set is **SmkStudy.dat**.
- **Toenail**: data set `toenail` in `R` package `mice`.

## Replication

Instructions to replicate all analyses in the paper as follows:

- Section 4.3, simulation 1 and Supplement 3.1:
  - In the script `01-simulations-absolute.R`, look for several `CHANGE` comments on lines 7-21. Change the paths here to whatever you like.
  - In partcular, make sure you tell the script where to find the `aghqmm` package.
  - Run the script `01-simulations-absolute.R`.
  - In the script `01-summarize-simulations-absolute.R`, look for `CHANGE` comments on lines 8 and 14. Change the paths here to whatever you chose in the previous script.
  - Run the script `01-summarize-simulations-absolute.R`.
- Section 4.3, simulation 2 and Supplements 3.2.1 and 3.2.2:
  - In the script `02-simulations-compare-glmmadaptive.R`, look for several `CHANGE` comments on lines 7-21. Change the paths here to whatever you like.
  - In partcular, make sure you tell the script where to find the `aghqmm` package.
  - Run the script `02-simulations-compare-glmmadaptive.R`.
  - In the script `02-summarize-simulations-compare-glmmadaptive.R`, look for `CHANGE` comments on lines 7 and 13. Change the paths here to whatever you chose in the previous script.
  - Run the script `02-summarize-simulations-compare-glmmadaptive.R`.
- Section 5.1 and Supplement 4.1, smoking data:
  - Download the `SmkStudy.dat` file as described below and save it to your machine.
  - In the script `03-smoking.R`, look for `CHANGE` comments on lines 7 and 11. Change the paths here to wherever you saved the `SmkStudy.dat` file.
  - At the `CHANGE` comment at line 18, tell the script where to find the `aghqmm` package.
  - Run the script `03-smoking.R`.
  - Note: there will be warnings; the model does not converge successfully for all methods and numbers of quadrature points. These are part of the results discussed in the paper.
- Section 5.2 and Supplement 4.2, toenail data:
  - In the script `04-toenail.R`, look for `CHANGE` comments on line 8. Change the path here to whatever you like.
  - At the `CHANGE` comment at line 15, tell the script where to find the `aghqmm` package.
  - Run the script `04-toenail.R`.
  - Note: there will be warnings; the model does not converge successfully for all methods and numbers of quadrature points. These are part of the results discussed in the paper.
- Supplement 3.3:



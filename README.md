# Approximate Marginal Likelihood Inference in Mixed Models for Grouped Data
Replication code for paper: *Approximate Marginal Likelihood Inference in Mixed Models for Grouped Data* by Alex Stringer.
The paper is [on arXiv](https://arxiv.org/abs/2310.01589).

All code was run on a 2021 M1 MacBook Pro with 10 cores and 64Gb of RAM. I did not test it on other hardware,
but I don't make use of any specific features of the Apple silicon either.

## Package

All scripts require the `aghqmm` package, available from GitHub [here](https://github.com/awstringer1/aghqmm).
The reproduction scripts below will install this package using `remotes::install_github()`, which requires only that you have
a GitHub Personal Access Token correctly configured in your `R` environment; see [here](https://carpentries.github.io/sandpaper-docs/github-pat.html) for a tutorial.

## Data

The data are obtained as follows:

- **Smoking**: obtained from the online supplementary material to the paper *A note on marginalization of regression parameters from mixed models of binary outcomes* by Donald Hedeker, Stephen H. C. du Toit, Hakan Demirtas, and Robert D. Gibbons (*Biometrics* (2017) **74**(1) pp. 354–361). Find that paper [here](https://onlinelibrary.wiley.com/doi/10.1111/biom.12707). Data set is **SmkStudy.dat**.
- **Toenail**: data set `toenail` in `R` package `mice`.

## Replication

Instructions to replicate all analyses in the paper as follows. All scripts will save results and figures to appropriate subdirectories of `getwd()`, which will be created if necessary.

The number of simulations run by default is much smaller than what was used to produce the results in the paper, to facilitate the testing of the scripts using continuous integration.
To reproduce the actual results in the paper, you should set `numruns` and `numsims` appropriately within each script; the numbers used in the paper are included as comments in each script.
Note that this takes days for each example, as a lot of simulations were used.

- Section 4.3, simulation 1 and Supplement 3.1:
  - Run `01-summarize-simulations-absolute.R`.
- Section 4.3, simulation 2 and Supplements 3.2.1 and 3.2.2:
  - Run `02-simulations-compare-glmmadaptive.R`.
  - Run `06-simulations-compare-lme4.R`.
- Section 5.1 and Supplement 4.1, smoking data:
  - Download the `SmkStudy.dat` file as described below and save it to your machine.
  - In the script `03-smoking.R`, look for `CHANGE` comments on lines 7 and 11. Change the paths here to wherever you saved the `SmkStudy.dat` file.
  - Run `03-smoking.R`.
  - Note: there will be warnings; the model does not converge successfully for all methods and numbers of quadrature points. These are part of the results discussed in the paper.
- Section 5.2 and Supplement 4.2, toenail data:
  - Run `04-toenail.R`.
  - Note: there will be warnings; the model does not converge successfully for all methods and numbers of quadrature points. These are part of the results discussed in the paper.
- Supplement 3.3:
  - Run `05-simulations-scalar.R`.

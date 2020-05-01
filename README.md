## Replication files for "Asymptotic Inference for the Constrained Quantile Regression Process"

First written 2017-10-27, cleaned up 2020-04-30.

This repository contains files to replicate the figures in the paper.
Two directories contain the files used to run each simulation experiment, and
they depend on two helper files located in the top-level directory.

### The experiment directories

Each directory contains a `sim.R` file, the script for the experiments, well as 
a shell script that was used on the SLURM system to execute the command.  The 
file `sim.R` should produce the data file `normal_sim.rda` or 
`oneregressor_treatment_sim.rda`.  The files `plot_normal.R` and 
`plot_treatment.R` each create two p-value plots corresponding to supremum-norm 
or L2-norm statistics.

### Common files

Both simulation experiments depend on the file `const_inf_utils.cpp`, which 
contains C++ code that is used to speed up the numerical work.  They also share 
the common file `plotting_utils.R` to create plots.

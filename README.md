# Overview
Our approach combines the **Generalized Lotka-Volterra (GLV) model** with **Bayesian inference**. In contrast to the classical GLV approach, which infers model parameters using penalized Ridge regression (as described in this [article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003388)), our method leverages [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) for Bayesian regression. CmdStan provides posterior distributions for model parameters, offering a robust measure of parameter uncertainty. We have applied this approach to study the ecological dynamics of the gut microbiome in response to dietary fiber supplementation. You can find our paper [here](https://academic.oup.com/ismej/article/16/8/2040/7474293).

---

# Prerequisites
The guidelines for installing CmdStan can be found [here](https://mc-stan.org/docs/cmdstan-guide/installation.html).

---

# GlvSolver Library
**GlvSolver** is a library of functions for input data processing, Stan file preparation, and parsing Stan output.

- **GlvSolver_v2** (Recommended) – Uses **matrix multiplication**, making it significantly more efficient.  
- **GlvSolver_v1** – May encounter **compilation difficulties** with complex problems (e.g., when handling a large number of taxa).

---

# Usage

## 1. Preparing Input Data
Please use the Jupyter notebooks included in the test examples as a starting point. The first step is to prepare an input table for GlvSolver with the following columns:

| Column Name      | Description |
|-----------------|-------------|
| **SubjectID**   | Unique identifier for subjects |
| **SampleID**    | Unique identifier for samples |
| **Timepoint**   | Time point of measurement |
| **Perturbation_1, ..., Perturbation_N** | External perturbations (e.g., treatments, interventions) |
| **Taxon_1, ..., Taxon_N** | Abundance of microbial taxa |

You may replace perturbation and taxon IDs with more interpretable labels, as long as they still begin with "Perturbation_" and "Taxon_", respectively.

---

## 2. Computing Log-Derivatives
Once the input table is ready, use the function `compute_dlogy_dt()` to calculate log-derivatives (dlogX/dt). We implemented two approaches for gradient computation:  

- Cubic spline fitting
- Second-order central difference

The computed log-derivatives will be appended as additional columns in the input table.

---

## 3. Generating X and Y Matrices
Next, generate the **X and Y matrices** required for GLV regression using `generate_XY_matrices()`. Then, pass these matrices to `write_stan_input_file()`. This function generates Stan data and model files while allowing:  

- Exclusion of specific interactions (by setting coefficients to 0 using the `pairs_to_exclude` argument).  
- Parameter constraints on self-interactions (`neg_self_int=True` constrains diagonal elements of beta matrix to be non-positive).

---

## 4. Running CmdStan
Navigate to the folder containing the generated Stan files and run `./run_cmdstan.sh`. By default, this script generates 3 Markov chains and produces a summary of posterior distributions from the Stan model. To verify model convergence, inspect the R̂ (R-hat) value in the summary file. Convergence is achieved when R̂ ≈ 1.00 (typically < 1.05).

---

## 6. Analyzing Posterior Distributions
To parse and summarize posterior distributions of model parameters, use `parse_stan_output()`. You can set `bci_cutoff` to define a credible threshold. The significance column is determined as follows:  

- If the entire credible interval is above or below zero, the interaction is significant.  
- Otherwise, it is not significant.

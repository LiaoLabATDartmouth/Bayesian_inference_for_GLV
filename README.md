# Overview
Our approach combines the Generalized Lotka-Volterra model (GLV) with Bayesian inference. In contrast to the canonical GLV approach, which infers model parameters using penalized Ridge regression (as described in this [article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003388)), our method leverages [CMDStan](https://mc-stan.org/users/interfaces/cmdstan) to perform Bayesian regression. CMDStan provides posterior distributions for model parameters as a measure of parameter uncertainties.

We have applied this approach to study the ecological dynamics of the gut microbiome in response to dietary fiber supplementation. You can find our paper [here](https://academic.oup.com/ismej/article/16/8/2040/7474293).

# Prerequisites
The guidelines for installing CmdStan can be found [here](https://mc-stan.org/docs/cmdstan-guide/installation.html).

# GlvSolver library
GlvSolver is a library of functions for input data processing, stan files preparation, and generating summarys of stan output. The second version (GlvSolver_v2) uses matrix multiplication and is more efficient. It is highly recommended over the first version (GlvSolver_v1), which may encoutner difficulties in compiling complex problems (e.g., too many taxa).

# Usage
Please use the Jupyter notebooks included in the test examples as a starting point. The first step is to prepare an input table for GlvSolver. This table should include the following columns: SubjectID, SampleID, Timepoint, Perturbation_1, Perturbation_2, …, Perturbation_N, Taxon_1, Taxon_2, …, Taxon_N. You may replace the perturbation and taxon IDs with more interpretable labels, as long as they still begin with ‘Perturbation_’ and ‘Taxon_’, respectively.

Once this table is available, you can pass this table to the function `compute_dlogy_dt` to compute log-derivatives (i.e., dlogX/dt). We implemented two approaches for gradient computation: derivatives of cubic spline fitting and second-order central difference. The computed log-derivatives will be added to the input table as additional columns.

The next step is generate X and Y matrics for GLV regression using the function `generate_XY_matrics`. Both matrics are further passed into the function `write_stan_input_file` to generate data and model stan files. In this function, we provide flexible to exclude certain interaction parameters by setting their associated coefficients to 0. These interactions can be set by passing a list of tuples as pairs_toi_exclude. We also allow parameter constrains on self-interactions. If
`neg_self_int` is set to True, the self-interactions (diagnoanls of beta) will be constrained to negative (i.e., force intraspecies competition).

Using the genertated stan data and model files, you can go to the folder that contain these files and run `run_cmdstan.sh`. By defaulkt, it generates 3 markov chains and generat esummary of the chains. plEAS check R2_Hat in the last column of the summary file. A converence requires R2_Hatr very close to 1 (e..g,  < 1.05).

To generate statistics of posterior distributiosn of these parameters, we can run the function `parse_stan_output`. You can set bci_cutoff as a way to specific the credible interaction. The signaicne column is determined when the entire credible interveal is above or below zero. Otherwise, it is not significnat.

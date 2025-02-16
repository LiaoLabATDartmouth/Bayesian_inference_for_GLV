# Overview
Our approach combines the Generalized Lotka-Volterra model (GLV) with Bayesian inference. In contrast to the canonical GLV approach, which infers model parameters using penalized Ridge regression (as described in this [article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003388)), our method leverages [CMDStan](https://mc-stan.org/users/interfaces/cmdstan) to perform Bayesian regression. CMDStan provides posterior distributions for model parameters, which allows inference of parameter uncertainties.

__We have applied this approach to study the ecological dynamics of the gut microbiome in response to dietary fiber supplementation. You can find our paper [here](https://academic.oup.com/ismej/article/16/8/2040/7474293).__

__For a demonstration of our algorithm using synthetic data, please refer to the file demonstration.ipynb.__

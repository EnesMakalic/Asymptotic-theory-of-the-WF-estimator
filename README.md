# Wallace-Freeman Estimator Bias Simulations

This repository contains source code for reproducing Monte Carlo simulations and bias plots associated with the paper: "Asymptotic theory and first-order bias of the Wallace-Freeman estimator".

The code compares the empirical bias of the maximum likelihood estimator (MLE) and the Wallace-Freeman estimator for the Weibull shape parameter.

## Overview
The main simulation function is:
```matlab
plot_weibull_bias(niter)
```

where niter is the number of Monte Carlo replications used for each value of the Weibull shape parameter.

## Main Function

```matlab
[fig1, fig2] = plot_weibull_bias(niter)
```

Optional name-value arguments include:

```matlab
plot_weibull_bias(niter, ...
    k_range=[0.1, 10], ...
    n_k_values=20, ...
    sample_size=100, ...
    lambda=1, ...
    show_stderr=false)
```

## Arguments

- niter: number of Monte Carlo replications for each value of the shape parameter.
- k_range: range of true Weibull shape values to evaluate. Default is [0.1, 10].
- n_k_values: number of shape values in the grid. Default is 20.
- sample_size: sample size for each Monte Carlo replication. Default is 100.
- lambda: fixed Weibull scale parameter. Default is 1.
- show_stderr: whether to display approximate standard-error bands. Default is false.

## Example
To reproduce the default simulation with 100000 Monte Carlo replications:

```matlab
[fig1, fig2] = plot_weibull_bias(1e5, sample_size=20);
```

This generates two figures:
1. Empirical versus analytical shape-parameter bias for the MLE and Wallace-Freeman estimator.
2. Direct comparison of empirical MLE and Wallace-Freeman shape-parameter bias.

## Requirements
The code was written for MATLAB and uses standard MATLAB statistical routines, including Weibull random-number generation and maximum likelihood fitting.

Required MATLAB functionality includes:

- wblrnd for generating Weibull random variables;
- wblfit for maximum likelihood estimation of Weibull parameters;
- parfor for parallel Monte Carlo loops, if parallel execution is enabled.

The Wallace-Freeman estimator is computed using the helper routine:

```matlab
estimate_weibull_mml87
```

## Citation

If you use this code, please cite the associated paper:

```bibtex
@misc{Makalic2026b,
  title = {Asymptotic Theory and Bias Correction for the {{Wallace--Freeman}} Estimator},
  author = {Makalic, Enes and Schmidt, Daniel F.},
  year = 2026,
  month = apr,
  howpublished = {https://arxiv.org/abs/2604.01568v1}
}

```

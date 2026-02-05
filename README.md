# Bayesian Unidimensional Scaling (BUDS)

**BUDS** is an R package based on `rstan` for recovering the latent ordering of observations. It performs Bayesian Unidimensional Scaling by mapping high-dimensional objects to 1D latent coordinates using pairwise dissimilarities.

## Overview

BUDS models pairwise dissimilarities between multivariate and high-dimensional objects to retrieve their latent 1D coordinate and ordering along an underlying continuum (gradient or trajectory). A Bayesian model is fitted to the pairwise dissimilarities between datapoints to recover the latent 1D coordinates for the datapoints along an underlying data trajectory. These 1D coordinates can be used to obtain the hidden order of the observations.

## Key Features

- **Latent 1D Mapping**: Fits a Bayesian model to pairwise dissimilarities to recover 1D coordinates.
- **Hidden Order Recovery**: Identifies the hidden order of observations based on recovered coordinates.
- **Stan-based Inference**: Leverages `rstan` for robust Bayesian inference, supporting both MCMC and Variational Bayes (VB) methods.
- **Visualization**: Provides utilities for visualizing the recovered 1D coordinates and the underlying data trajectory.

## Installation

You can install the package from source using `devtools`:

```r
# install.packages("devtools")
devtools::install(".")
```

## Basic Usage

The main function in the package is `fit_buds()`. It requires a square matrix of pairwise dissimilarities as input.

```r
library(buds)

# Assume D is a square matrix of pairwise dissimilarities
# fit <- fit_buds(D, method = "vb")

# To visualize the 1D mapping:
# plot_buds_1D(fit)
```

For more detailed information, please refer to the package documentation.

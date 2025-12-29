# DAPS: Deterministic Approximation of Posteriors

## Overview

`DAPS` is an R package implementing the **DAP-S** (Deterministic Approximation of Posteriors) algorithm for Bayesian genetic fine-mapping.
DAP-S combines the computational efficiency of SuSiE's variational approximation with the accuracy and interpretability of Bayesian Variable Selection Regression (BVSR) model using the DAP framework.

The method provides a **fast, scalable, and well-calibrated** approach for identifying causal genetic variants and constructing **interpretable signal clusters**, which enable flexible construction of credible sets.

## Installation
```r
# install the package via devtools
install.packages("devtools")
devtools::install_github("wnbo9/DAPS")
# load the package
library(DAPS)
```

## Vignettes
See the [vignettes](https://wnbo9.github.io/DAPS/articles/index.html) for a quick start and worked examples.


## BLAS acceleration
`DAPS` relies on linear algebra operations. When linked against an optimized BLAS implementation-such as **OpenBLAS**, **Intel MKL**, and **Apple's Accelerate-framework**-`DAPS` can be substantially faster than when using the reference (Netlib) BLAS shipped with base R. Most modern R installations already use an optimized BLAS, but performance can vary across systems. You can check which BLAS library your R session is using with:
```r
sessionInfo()
```
A simple benchmark to verify that an optimized BLAS in is use:
```r
m <- 5000
A <- matrix(rnorm(m^2), m, m)
system.time(A %*% A)
```
As a rough guideline, an elapsed time on the order of ~1 second indicates that your R installation is likely using an optimized BLAS.
# DAPS

## Deterministic Approximation of Posteriors

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

## Dependencies
`DAPS` requires R (>= 4.0.0), susieR (>= 0.15.0), Rcpp (>= 1.1.0), RcppArmadillo, Matrix, ggplot2 and patchwork.

We recommend running `DAPS` with an optimized BLAS library, such as **OpenBLAS**, **Intel MKL**, or **Apple's Accelerate framework**, which can be substantially faster than the reference BLAS shipped with R. Check your BLAS:
```r
sessionInfo()   # see the "BLAS:" and "LAPACK:" lines
```
A quick benchmark: an elapsed time of about 1 second or less suggests an optimized BLAS is in use.
```r
m <- 5000
A <- matrix(rnorm(m^2), m, m)
system.time(A %*% A)
```
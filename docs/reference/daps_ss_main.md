# DAPS_SS main function in Rcpp

DAPS_SS main function in Rcpp

## Usage

``` r
daps_ss_main(
  XtX,
  Xty,
  yty,
  n,
  prior,
  proposal,
  proposal_thresh,
  grid,
  twas_weight,
  min_abs_corr,
  sigma2,
  sigma2_alpha
)
```

## Arguments

- XtX:

  Cross-product of genotype matrix (p x p)

- Xty:

  Cross-product of genotype and phenotype vector (p x 1)

- yty:

  Cross-product of phenotype vector (scalar)

- n:

  Sample size (scalar)

- prior:

  Prior inclusion probabilities (p x 1)

- proposal:

  Proposal matrix (p x L) from SuSiE

- proposal_thresh:

  Proposal threshold (scalar)

- grid:

  Grid of scaled prior variance values (k x 1)

- twas_weight:

  Logical indicating whether to use TWAS weights

- min_abs_corr:

  Minimum absolute correlation for TWAS weights (scalar)

- sigma2:

  Residual variance estimate from SuSiE (scalar; used for inf or ash)

- sigma2_alpha:

  Estimated variance of small effect sizes (scalar; used for inf or ash)

## Value

A list containing models and input data

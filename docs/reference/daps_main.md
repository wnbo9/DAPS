# DAPS main function in Rcpp

DAPS main function in Rcpp

## Usage

``` r
daps_main(
  X,
  y,
  n,
  prior,
  proposal,
  proposal_thresh,
  grid,
  twas_weight,
  min_abs_corr
)
```

## Arguments

- X:

  Genotype matrix (n x p)

- y:

  Phenotype vector (n x 1)

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

## Value

A list containing models and input data

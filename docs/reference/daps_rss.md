# DAP-S fine-mapping with regression summary statistics

DAP-S fine-mapping using regression summary statistics to recover
sufficient statistics.

## Usage

``` r
daps_rss(
  z = NULL,
  R,
  n = NULL,
  bhat = NULL,
  shat = NULL,
  var_y = NULL,
  estimate_residual_variance = FALSE,
  ...
)
```

## Arguments

- z:

  A p-vector of single-SNP z-scores.

- R:

  A p x p LD matrix among the SNPs.

- n:

  Sample size.

- bhat:

  A p-vector of single-SNP effect size estimates.

- shat:

  A p-vector of standard errors for bhat. bhat and shat together can be
  replaced by z.

- var_y:

  Phenotypic variance, defined as y'y/(n-1). When not provided, the
  algorithm works on standardized data.

- estimate_residual_variance:

  Estimate residual variance in SuSiE. Default is FALSE. If the in
  sample LD matrix is used, we recommend setting this to TRUE.

- ...:

  Other parameters to be passed to `daps_suff_stat`.

## Value

A list of DAP-S fine-mapping results.

## Examples

``` r
z <- c(6,7)
R <- matrix(1,2,2)
rst <- daps_rss(z=z, R=R)
rst$pip
#>      [,1]
#> [1,]  NaN
#> [2,]  NaN
```

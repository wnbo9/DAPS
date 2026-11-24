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
  prior = NULL,
  null_wt = NULL,
  estimate_residual_variance = FALSE,
  proposal_thresh = 1e-06,
  grid = c(0.04, 0.16, 0.64),
  unmappable_effects = c("none", "inf"),
  twas_weight = FALSE,
  coverage = NULL,
  min_abs_corr = 0.5,
  dedup = TRUE,
  verbose = FALSE,
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

- prior:

  Marginal SNP-level prior weights, given as a vector of length p.
  Default is NULL which corresponds to uninformative weights.

- null_wt:

  Null weight, defined as the prior probability of no causal variants in
  the region. Default is calculated as `prod(1 - prior)`. The prior and
  null_wt are used to compute the prior weights and null weight for
  SuSiE, where the prior weights are conditional on the non-null.

- estimate_residual_variance:

  Estimate residual variance in SuSiE. Default is FALSE. If the in
  sample LD matrix is used, we recommend setting this to TRUE.

- proposal_thresh:

  Threshold for proposal acceptance. Default value is 1e-6.

- grid:

  Grid of values for scaled effect size prior variance. Default is
  c(0.04, 0.16, 0.64).

- unmappable_effects:

  Method for handling unmapped effects. Options are "none", "inf".
  Default is "none".

- twas_weight:

  Compute TWAS weights. Default is FALSE.

- coverage:

  Coverage level for credible set (e.g. 0.95). Default is NULL, in which
  the algorithm returns signal clusters that are not constrained by a
  fixed coverage level.

- min_abs_corr:

  Minimum absolute correlation allowed in a signal cluster or credible
  set. Default is 0.5, which corresponds to a squared correlation of
  0.25.

- dedup:

  Whether to remove duplicate effects (effects with the same top 5
  SNPs). Default is TRUE.

- verbose:

  Whether to print progress messages. Default is FALSE.

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
#> [1] 0.0741417 0.9256085
```

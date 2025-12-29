# DAP-S fine-mapping with individual-level data

DAP-S fine-mapping with individual-level genotype and phenotype data as
input.

## Usage

``` r
daps(
  X,
  y,
  L = min(10, ncol(X)),
  standardize = FALSE,
  prior = NULL,
  null_wt = NULL,
  proposal_thresh = 1e-06,
  grid = c(0.04, 0.16, 0.64),
  twas_weight = FALSE,
  coverage = NULL,
  min_abs_corr = 0.5,
  ...
)
```

## Arguments

- X:

  Genotype matrix (n x p) where n is the number of individuals and p is
  the number of SNPs.

- y:

  Phenotype vector of length n.

- L:

  Maximum number of causal SNPs allowed in the model. Default is min(10,
  p).

- standardize:

  Standardize each column of X. Default is FALSE.

- prior:

  SNP-level prior weights, given as a vector of length p. Default is
  NULL which corresponds to uninformative weights.

- null_wt:

  Null weight used in SuSiE, defined as the prior probability of no
  causal variants in the region. Default is calculated as
  `prod(1 - prior)`.

- proposal_thresh:

  Threshold for proposal acceptance. Default value is 1e-6.

- grid:

  Grid of values for scaled effect size prior variance. Default is
  c(0.04, 0.16, 0.64).

- twas_weight:

  Compute TWAS weights. Default is FALSE.

- coverage:

  Coverage leveel for credible set (e.g. 0.95). Default is NULL, in
  which the algorithm returns signal clusters that are not constrained
  by a fixed coverage level.

- min_abs_corr:

  Minimum absolute correlation allowed in a signal cluster or credible
  set. Ddefault is 0.5, which corresponds to a squared correlation of
  0.25.

- ...:

  Additional arguments passed to `susie`.

## Value

A list of DAP-S fine-mapping results.

## Examples

``` r
set.seed(1234)
n = 1000
p = 1000
beta = rep(0, p)
beta[c(1, 200, 500, 800)] = 1
X = matrix(rnorm(n*p), nrow = n, ncol = p)
X = scale(X,center = TRUE, scale = FALSE)
y = X %*% beta + rnorm(n)
rst = daps(X, y, L = 5)
daps_plot(rst)
```

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

- effect_pip:

  A p\*L matrix of SNP posterior probabilities in each signal

- pip:

  A p-vector of marginal SNP-level posterior inclusion probabilities

- prior:

  A p-vector of SNP prior probability of being causal

- models:

  A dataframe of high-probability models identified by DAP-S, together
  with their posterior probabilities

- twas_weights:

  A p-vector of TWAS weights for each SNP

- info:

  A list of information including signal clusters and stored model
  details

- sets:

  A list of signal clusters or credible sets identified by DAP-S with
  the specified coverage level

- variants:

  A variant-level dataframe for eQTL/enloc-style outputs

## Examples

``` r
set.seed(1234)
n <- 1000
p <- 1000
beta <- rep(0, p)
beta[c(1, 200, 500, 800)] <- 1
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X <- scale(X, center = TRUE, scale = TRUE)
y <- X %*% beta + rnorm(n)
rst <- daps(X, y, L = 5)
rst$sets
#> $sets
#> $sets$S1
#> [1] 800
#> 
#> $sets$S2
#> [1] 200
#> 
#> $sets$S3
#> [1] 500
#> 
#> $sets$S4
#> [1] 1
#> 
#> 
#> $purity
#>    min.abs.corr mean.abs.corr median.abs.corr
#> S1            1             1               1
#> S2            1             1               1
#> S3            1             1               1
#> S4            1             1               1
#> 
#> $set_index
#> [1] 1 2 3 4
#> 
#> $coverage
#> [1] 1 1 1 1
#> 
#> $min_abs_corr
#> [1] 0.5
#> 
#> $requested_coverage
#> NULL
#> 
daps_plot(rst)
```

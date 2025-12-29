# Getting started with DAPS

## DAPS with Individual-level Data

``` r
# Load the DAPS package
library(DAPS) 

# Run DAP-S fine-mapping
set.seed(1234)
n = 1000
p = 1000
beta = rep(0, p)
beta[1:4] = 1
X = matrix(rnorm(n*p), nrow = n, ncol = p)
y = X %*% beta + rnorm(n)
X = scale(X,center = TRUE,scale = FALSE)
y = scale(y,center = TRUE,scale = FALSE)
rst = daps(X, y, L = 5)
rst$sets
#> $sets
#> $sets$S1
#> [1] 1
#> 
#> $sets$S2
#> [1] 2
#> 
#> $sets$S3
#> [1] 3
#> 
#> $sets$S4
#> [1] 4
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
```

## DAPS with Sufficient Statistics

DAPS is also designed to also work efficiently using **sufficient
statistics**, which includes
$\mathbf{X}^{T}\mathbf{X},\;\mathbf{X}^{T}\mathbf{y},\;\mathbf{y}^{T}\mathbf{y},\; n$.

``` r
XtX = crossprod(X)
Xty = as.numeric(crossprod(X, y))
yty = sum(y^2)
rst = daps_ss(XtX, Xty, yty, n, L = 5)
rst$sets
#> $sets
#> $sets$S1
#> [1] 1
#> 
#> $sets$S2
#> [1] 2
#> 
#> $sets$S3
#> [1] 3
#> 
#> $sets$S4
#> [1] 4
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
```

## DAPS with Regression Summary Statistics

A common scenario in genetics is having only GWAS **regression summary
statistics** (e.g., Z-scores and LD matrices). We can use these to
recover the necessary inputs for DAP-S fine-mapping.

``` r
rss = susieR::univariate_regression(X, y)
z = rss$betahat / rss$sebetahat
R = cor(X)
rst <- daps_rss(z=z, R=R, n=n)
rst$sets
#> $sets
#> $sets$S1
#> [1] 1
#> 
#> $sets$S2
#> [1] 2
#> 
#> $sets$S3
#> [1] 3
#> 
#> $sets$S4
#> [1] 4
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
```

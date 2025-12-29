# DAP-S variant-level table for eQTL/enloc-style outputs

Construct a variant-level data frame listing variants inside each signal
cluster (coverage = NULL) or credible set (coverage given).

## Usage

``` r
get_enloc(output, coverage = NULL)
```

## Arguments

- output:

  A list of DAP-S fine-mapping results.

- coverage:

  Desired coverage for credible sets. If NULL, use full signal clusters.

## Value

A data.frame with columns: SC_ID, SNP_ID, PIP.

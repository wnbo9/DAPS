# DAP-S variant-level table for eQTL/enloc-style outputs

Construct a variant-level data frame listing variants inside each signal
cluster or credible set (coverage given).

## Usage

``` r
get_enloc(output)
```

## Arguments

- output:

  A list of DAP-S fine-mapping results.

## Value

A data.frame with columns: SC_ID, SNP_ID, PIP.

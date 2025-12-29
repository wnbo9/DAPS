#' DAPS: Deterministic Approximation of Posteriors for Bayesian fine-mapping
#'
#' DAPS implements the DAP-S algorithm for Bayesian genetic fine-mapping,
#' integrating SuSiE's variational approximation with deterministic posterior
#' computation to identify causal variants and construct interpretable
#' signal clusters.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib DAPS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' DAP-S fine-mapping with individual-level data
#'
#' DAP-S fine-mapping with individual-level genotype and phenotype data as
#' input.
#'
#' @param X Genotype matrix (n x p) where n is the number of individuals and
#'  p is the number of SNPs.
#' @param y Phenotype vector of length n.
#' @param L Maximum number of causal SNPs allowed in the model.  Default is
#'  min(10, p).
#' @param standardize Standardize each column of X. Default is FALSE.
#' @param prior SNP-level prior weights, given as a vector of
#'  length p. Default is NULL which corresponds to uninformative weights.
#' @param null_wt Null weight used in SuSiE, defined as the prior
#'  probability of no causal variants in the region. Default is calculated
#'  as \code{prod(1 - prior)}.
#' @param proposal_thresh Threshold for proposal acceptance. Default value
#'  is 1e-6.
#' @param grid Grid of values for scaled effect size prior variance.
#'  Default is c(0.04, 0.16, 0.64).
#' @param twas_weight Compute TWAS weights. Default is FALSE.
#' @param coverage Coverage leveel for credible set (e.g. 0.95). Default is
#'  NULL, in which the algorithm returns signal clusters that are not
#'  constrained by a fixed coverage level.
#' @param min_abs_corr Minimum absolute correlation allowed in a signal
#'  cluster or credible set. Ddefault is 0.5, which corresponds to a
#'  squared correlation of 0.25.
#' @param ... Additional arguments passed to \code{\link{susie}}.
#'
#' @return A list of DAP-S fine-mapping results.
#' @export
#'
#' @examples
#' set.seed(1234)
#' n = 1000
#' p = 1000
#' beta = rep(0, p)
#' beta[c(1, 200, 500, 800)] = 1
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' X = scale(X,center = TRUE, scale = FALSE)
#' y = X %*% beta + rnorm(n)
#' rst = daps(X, y, L = 5)
#' daps_plot(rst)
daps <- function(
  X, y, L = min(10, ncol(X)),
  standardize = FALSE,
  prior = NULL, null_wt = NULL,
  proposal_thresh = 1e-6,
  grid = c(0.04, 0.16, 0.64),
  twas_weight = FALSE,
  coverage = NULL,
  min_abs_corr = 0.5,
  ...
) {

  X <- scale(X, center = TRUE, scale = standardize)
  y <- scale(y, center = TRUE, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(prior)) prior <- rep(1 / p, p)
  if (is.null(null_wt)) null_wt <- prod(1 - prior)
  if (!is.null(coverage) && (coverage <= 0 || coverage > 1))
    stop("Error: coverage must be in (0, 1].\n")

  susie_args <- list(
    X = X, y = y, L = L,
    standardize = standardize,
    prior_weights = prior / (1 - null_wt),
    null_weight = null_wt,
    coverage = 0
  )
  susie_args <- modifyList(susie_args, list(...))
  susie_fit <- do.call(susieR::susie, susie_args)
  proposal <- get_proposal(susie_fit)

  daps_fit <- daps_main(
    X, y, n,
    prior, proposal, proposal_thresh,
    grid, twas_weight, min_abs_corr
  )

  output <- summarize(daps_fit, prior, min_abs_corr)
  output$sets <- get_set(output, coverage)
  output$variants <- get_enloc(output, coverage)

  return(output)
}



#' DAP-S fine-mapping with sufficient statistics
#'
#' DAP-S fine-mapping using sufficient statistics, either provided directly
#' or recovered from regression summary statistics.
#'
#' @param XtX A p x p matrix X'X in which the columns of X are centered
#'  and p is the number of SNPs.
#' @param Xty A p-vector X'y in which y and the columns of X are centered.
#' @param yty A scalar y'y in which y is centered.
#' @param n Sample size.
#' @param L Maximum number of causal SNPs allowed in the model.  Default is
#'  min(10, p).
#' @param standardize Standardize each column of X and y. Default is FALSE.
#' @param prior SNP-level prior weights, given as a vector of
#'  length p. Default is NULL which corresponds to uninformative weights.
#' @param null_wt Null weight used in SuSiE, defined as the prior
#'  probability of no causal variants in the region. Default is calculated
#'  as \code{prod(1 - prior)}.
#' @param estimate_residual_variance Estimate residual variance in SuSiE.
#'  Default is TRUE.
#' @param proposal_thresh Threshold for proposal acceptance. Default value
#'  is 1e-6.
#' @param grid Grid of values for scaled effect size prior variance.
#'  Default is c(0.04, 0.16, 0.64).
#' @param twas_weight Compute TWAS weights. Default is FALSE.
#' @param coverage Coverage level for credible set (e.g. 0.95). Default is
#'  NULL, in which the algorithm returns signal clusters that are not
#'  constrained by a fixed coverage level.
#' @param min_abs_corr Minimum absolute correlation allowed in a signal
#'  cluster or credible set. Default is 0.5, which corresponds to a
#'  squared correlation of 0.25.
#' @param ... Other parameters to be passed to \code{\link{susie_suff_stat}}.
#'
#' @return A list of DAP-S fine-mapping results.
#' @export
#'
#' @examples
#' set.seed(1234)
#' n = 1000
#' p = 1000
#' beta = rep(0, p)
#' beta[c(1, 200, 500, 800)] = 1
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' X = scale(X,center = TRUE,scale = FALSE)
#' y = X %*% beta + rnorm(n)
#' y = scale(y,center = TRUE,scale = FALSE)
#' XtX = crossprod(X)
#' Xty = as.numeric(crossprod(X, y))
#' yty = sum(y^2)
#' rst = daps_ss(XtX, Xty, yty, n, L = 5)
#' daps_plot(rst)
daps_ss <- function(
    XtX, Xty, yty, n,
    L = min(10, ncol(XtX)),
    standardize = FALSE,
    prior = NULL, null_wt = NULL,
    estimate_residual_variance = TRUE,
    proposal_thresh = 1e-6,
    grid = c(0.04, 0.16, 0.64),
    twas_weight = FALSE,
    coverage = NULL,
    min_abs_corr = 0.5,
    ...
) {

  if (standardize) {
    dXtX <- diag(XtX)
    csd <- sqrt(dXtX / (n - 1))
    csd[csd == 0] <- 1
    XtX <- t((1 / csd) * XtX) / csd
    Xty <- Xty / csd
  }
  p <- ncol(XtX)

  if (is.null(prior)) prior <- rep(1 / p, p)
  if (is.null(null_wt)) null_wt <- prod(1 - prior)
  if (!is.null(coverage) && (coverage <= 0 || coverage > 1))
    stop("Error: coverage must be in (0, 1].\n")

  susie_args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
    standardize = standardize,
    prior_weights = prior / (1 - null_wt),
    null_weight = null_wt,
    estimate_residual_variance = estimate_residual_variance,
    coverage = 0
  )
  susie_args <- modifyList(susie_args, list(...))
  susie_fit <- do.call(susieR::susie_ss, susie_args)
  proposal <- get_proposal(susie_fit)

  daps_fit <- daps_ss_main(
    XtX, Xty, yty, n,
    prior, proposal, proposal_thresh,
    grid, twas_weight, min_abs_corr
  )

  output <- summarize(daps_fit, prior, min_abs_corr)
  output$sets <- get_set(output, coverage)
  output$variants <- get_enloc(output, coverage)

  return(output)
}



#' DAP-S fine-mapping with regression summary statistics
#'
#' DAP-S fine-mapping using regression summary statistics to recover
#'  sufficient statistics.
#'
#' @param z A p-vector of single-SNP z-scores.
#' @param R A p x p LD matrix among the SNPs.
#' @param n Sample size.
#' @param bhat A p-vector of single-SNP effect size estimates.
#' @param shat A p-vector of standard errors for bhat. bhat and shat
#'  together can be replaced by z.
#' @param var_y Phenotypic variance, defined as y'y/(n-1). When not
#'  provided, the algorithm works on standardized data.
#' @param estimate_residual_variance Estimate residual variance in SuSiE.
#'  Default is FALSE. If the in sample LD matrix is used, we recommend
#'  setting this to TRUE.
#' @param ... Other parameters to be passed to \code{\link{daps_suff_stat}}.
#'
#' @return A list of DAP-S fine-mapping results.
#' @export
#'
#' @examples
#' z <- c(6,7)
#' R <- matrix(1,2,2)
#' rst <- daps_rss(z=z, R=R)
#' daps_plot(rst)
daps_rss <- function(
  z = NULL, R, n = NULL,
  bhat = NULL, shat = NULL, var_y = NULL,
  estimate_residual_variance = FALSE,
  ...
) {

  if (is.null(z)) {
    if (is.null(bhat) || is.null(shat)) {
      stop("Error: Either z or both bhat and shat must be provided.\n")
    }
    if (length(bhat) != length(shat)) {
      stop("Error: bhat and shat must have the same length.\n")
    }
    z <- bhat / shat
  }
  p <- length(z)
  if (!is.matrix(R) || nrow(R) != p || ncol(R) != p) {
    stop("'R' must be a p x p matrix, where p = length(z).")
  }
  R <- (R + t(R)) / 2

  if (is.null(n)) {
    XtX <- R
    Xty <- z
    yty <- 1
    n <- 2
  } else {
    # when n is provided, we use PVE-adjusted z scores
    adj <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(adj) * z

    if (is.null(shat) && !is.null(var_y)) {
      # var_y, shat and bhat are provided, results on original scale
      XtXdiag <- var_y * adj / (shat^2)
      XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX <- (XtX + t(XtX)) / 2
      Xty <- z * sqrt(adj) * var_y / shat
      yty <- (n - 1) * var_y
    } else {
      # results on standardized X, y scale
      XtX <- (n - 1) * R
      Xty <- sqrt(n - 1) * z
      yty <- n - 1
    }
  }

  daps_args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n,
    estimate_residual_variance = estimate_residual_variance
  )
  daps_args <- modifyList(daps_args, list(...))
  output <- do.call(daps_ss, daps_args)

  return(output)
}
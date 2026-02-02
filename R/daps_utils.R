#' Extract proposal distribution from SuSiE fit
#'
#' Extract the proposal distribution used in DAP-S
#'  from the SuSiE fit.
#'
#' @param susie_fit A SuSiE fit object.
#'
#' @return A matrix of proposal distribution.
#' @noRd
get_proposal <- function(susie_fit) {
  # Extract non-zero effects
  mat <- t(susie_fit$alpha)
  cols <- which(susie_fit$V != 0)
  if (length(cols) == 0) cols <- 1
  mat <- mat[, cols, drop = FALSE]

  # Remove duplicate effects
  top5 <- apply(mat, 2, function(x) sort(order(x, decreasing = TRUE)[1:5]))
  keep <- !duplicated(as.data.frame(t(top5)))

  return(mat[, keep, drop = FALSE])
}

#' Generate DAP-S fine-mapping results output
#'
#' Generate DAP-S fine-mapping results output
#'
#' @param daps_fit A list of DAP-S fine-mapping results.
#' @param prior The prior probabilities used in DAP-S.
#' @param min_abs_corr Minimum absolute correlation among SNPs
#'  in each credible set.
#' @return A list of DAP-S fine-mapping results.
#' @noRd
summarize <- function(daps_fit, prior, min_abs_corr) {

  models <- data.frame(
    model = daps_fit$model_configs,
    log_bf = daps_fit$log_bf,
    log_prior = daps_fit$log_prior,
    log_posterior = daps_fit$log_posterior,
    posterior_prob = daps_fit$posterior_prob
  )

  info <- list(
    signal_clusters = lapply(daps_fit$signal_clusters, function(x) x + 1),
    cpip = daps_fit$cpip,
    min_abs_corr = min_abs_corr,
    LD = daps_fit$LD,
    effect_pip = daps_fit$effect_pip,
    model_index = daps_fit$model_index,
    reg_weights = daps_fit$reg_weights,
    log_nc = daps_fit$log_nc
  )

  return(list(
    pip = daps_fit$pip,
    prior = prior,
    models = models,
    twas_weights = daps_fit$twas_weights,
    info = info
  ))
}

#' DAP-S fine-mapping getting signal clusters or credible sets
#'
#' DAP-S fine-mapping getting signal clusters or credible sets
#'
#' @param output A list of DAP-S fine-mapping results.
#' @param coverage Desired coverage for credible sets.
#'  If coverage is NULL, return signal clusters instead of credible sets.
#' @return A list of signal clusters or credible sets.
#' @export
get_set <- function(output, coverage) {
  sets <- list()
  set_indices <- c()
  set_coverages <- c()

  for (i in seq_along(output$info$signal_clusters)) {
    cluster_snps <- output$info$signal_clusters[[i]]
    cluster_cpip <- output$info$cpip[i]
    if (length(cluster_snps) == 0) next

    if (is.null(coverage)) {
      sets[[length(sets) + 1]] <- cluster_snps
      set_indices <- c(set_indices, i)
      set_coverages <- c(set_coverages, cluster_cpip)

    } else if (coverage > 0) {
      if (cluster_cpip < coverage) next

      ordered_snps <- cluster_snps
      pips <- output$info$effect_pip[ordered_snps, i]
      cumsum_pip <- cumsum(pips)

      k <- which(cumsum_pip >= coverage)[1]
      if (is.na(k)) {
        k <- length(ordered_snps)
      }
      set_snps <- ordered_snps[seq_len(k)]
      actual_coverage <- sum(output$info$effect_pip[set_snps, i])

      sets[[length(sets) + 1]] <- set_snps
      set_indices <- c(set_indices, i)
      set_coverages <- c(set_coverages, actual_coverage)
    }
  }

  if (length(sets) == 0) {
    return(list(
      sets = NULL,
      coverage = NULL,
      min_abs_corr = output$info$min_abs_corr,
      requested_coverage = if (is.null(coverage)) "signal clusters" else coverage
    ))
  }

  names(sets) <- paste0("S", set_indices)

  # Purity table
  purity <- do.call(rbind, lapply(sets, function(set_snps) {
    if (length(set_snps) <= 1) {
      return(c(min.abs.corr = 1, mean.abs.corr = 1, median.abs.corr = 1))
    }

    LD_sub <- output$info$LD[set_snps, set_snps, drop = FALSE]
    upper_tri <- Matrix::triu(LD_sub, k = 1)
    abs_cors <- abs(upper_tri@x)
    if (length(abs_cors) == 0) abs_cors <- 0

    c(min.abs.corr    = min(abs_cors),
      mean.abs.corr   = mean(abs_cors),
      median.abs.corr = stats::median(abs_cors))
  }))
  purity <- as.data.frame(purity)

  return(list(
    sets = sets,
    purity = purity,
    set_index = set_indices,
    coverage = set_coverages,
    min_abs_corr = output$info$min_abs_corr,
    requested_coverage = if (is.null(coverage)) "signal clusters" else coverage
  ))
}


#' DAP-S variant-level table for eQTL/enloc-style outputs
#'
#' Construct a variant-level data frame listing variants inside each
#' signal cluster or credible set (coverage given).
#'
#' @param output A list of DAP-S fine-mapping results.
#' @return A data.frame with columns: SC_ID, SNP_ID, PIP.
#' @export
get_enloc <- function(output) {

  sets <- output$sets$sets
  set_ids <- output$sets$set_index

  if (length(sets) == 0) {
    return(data.frame(SC_ID = integer(), SNP_ID = integer(), PIP = numeric()))
  }

  sc_id <- vector("list", length(sets))
  snp_id <- vector("list", length(sets))
  pip_values <- vector("list", length(sets))

  for (i in seq_along(sets)) {
    current_sc_id <- set_ids[i]
    current_snps <- sets[[i]]

    if (length(current_snps) == 0) next

    current_pips <- output$info$effect_pip[current_snps, current_sc_id]

    sc_id[[i]] <- rep.int(current_sc_id, length(current_snps))
    snp_id[[i]] <- current_snps
    pip_values[[i]] <- current_pips
  }

  data.frame(
    SC_ID = unlist(sc_id),
    SNP_ID = unlist(snp_id),
    PIP = unlist(pip_values),
    row.names = NULL
  )
}

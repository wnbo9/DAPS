#' Extract proposal distribution from SuSiE fit
#'
#' Extract the proposal distribution used in DAP-S
#'  from the SuSiE fit.
#'
#' @param susie_fit A SuSiE fit object.
#' @param dedup Whether to remove duplicate effects
#'  (effects with the same top 5 SNPs). Default is TRUE.
#'
#' @return A matrix of proposal distribution.
#' @noRd
get_proposal <- function(susie_fit, dedup = TRUE) {
  # Extract non-zero effects
  mat <- t(susie_fit$alpha)
  cols <- which(susie_fit$V > 1e-9)
  if (length(cols) == 0) cols <- 1
  mat <- mat[, cols, drop = FALSE]
  if (is.null(susie_fit$null_index) || susie_fit$null_index == 0) {
    mat <- rbind(mat, 0)
  }

  if (dedup) {
    # Remove duplicate effects
    k <- min(5, nrow(mat))
    key <- apply(mat, 2, function(x)
      paste(sort(order(x, decreasing = TRUE)[seq_len(k)]), collapse = ","))
    mat <- mat[, !duplicated(key), drop = FALSE]
  }

  return(mat)
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
summarize <- function(daps_fit, prior, min_abs_corr, snp_names) {

  models <- data.frame(
    model = daps_fit$model_configs,
    log_bf = as.numeric(daps_fit$log_bf),
    log_prior = as.numeric(daps_fit$log_prior),
    log_posterior = as.numeric(daps_fit$log_posterior),
    posterior_prob = as.numeric(daps_fit$posterior_prob)
  )

  effect_pip <- daps_fit$effect_pip
  dimnames(effect_pip) <- list(snp_names, paste0("S", seq_len(ncol(effect_pip))))

  twas_weights <- if (length(daps_fit$twas_weights) == 0) NULL else
    stats::setNames(as.numeric(daps_fit$twas_weights), snp_names)

  info <- list(
    signal_clusters = lapply(daps_fit$signal_clusters, function(x) as.integer(x) + 1L),
    cpip = as.numeric(daps_fit$cpip),
    min_abs_corr = min_abs_corr,
    LD = daps_fit$LD,
    effect_pip = effect_pip,
    model_index = model_index_to_r(daps_fit$model_index, length(prior)),
    reg_weights = reg_weights_to_r(daps_fit$reg_weights),
    log_nc = daps_fit$log_nc
  )

  return(list(
    pip = stats::setNames(as.numeric(daps_fit$pip), snp_names),
    prior = stats::setNames(as.numeric(prior), snp_names),
    models = models,
    twas_weights = twas_weights,
    info = info
  ))
}

# Index conventions for stored models
#   C++ : 0-based SNP index, null slot = p (last row of the proposal)
#   R   : 1-based SNP index, null slot = 0
# Each model keeps its effect-slot layout (length L for sampled models,
# length 1 for the added single-SNP / null models).

#' Convert model index from C++ (0-based, null = p) to R (1-based, null = 0)
#'
#' @param model_index List of integer vectors returned by the C++ code.
#' @param p Number of SNPs.
#' @return List of integer vectors.
#' @noRd
model_index_to_r <- function(model_index, p) {
  lapply(model_index, function(m) {
    m <- as.integer(m)
    out <- m + 1L
    out[m == p] <- 0L
    out
  })
}

#' Convert model index from R (1-based, null = 0) to C++ (0-based, null = p)
#'
#' Use this before passing \code{info$model_index} back to C++
#' (e.g. in a model-update step).
#'
#' @param model_index List of integer vectors in R convention.
#' @param p Number of SNPs.
#' @return List of integer vectors.
#' @noRd
model_index_to_cpp <- function(model_index, p) {
  lapply(model_index, function(m) {
    m <- as.integer(m)
    if (any(m < 0L | m > p)) stop("model_index values must be in 0..p.")
    out <- m - 1L
    out[m == 0L] <- as.integer(p)
    out
  })
}

#' Convert regression weight names from 0-based to 1-based SNP index
#'
#' @param reg_weights List of named numeric vectors (names = 0-based SNP index).
#' @return List of named numeric vectors (names = 1-based SNP index).
#' @noRd
reg_weights_to_r <- function(reg_weights) {
  lapply(reg_weights, function(w) {
    names(w) <- as.character(as.integer(names(w)) + 1L)
    w
  })
}

#' Convert regression weight names from 1-based to 0-based SNP index
#'
#' @param reg_weights List of named numeric vectors (names = 1-based SNP index).
#' @return List of named numeric vectors (names = 0-based SNP index).
#' @noRd
reg_weights_to_cpp <- function(reg_weights) {
  lapply(reg_weights, function(w) {
    names(w) <- as.character(as.integer(names(w)) - 1L)
    w
  })
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

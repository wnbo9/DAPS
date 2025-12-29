#' DAP-S PIP plot with cluster / credible-set annotation
#'
#' DAP-S PIP plot with cluster / credible-set annotation
#'
#' @param res A list of DAP-S fine-mapping results.
#' @param band Fraction of the row gap allowed for lifting between 0 and 1;
#'  default 0.
#'
#' @return A ggplot object (patchwork layout).
#' @export
#' @import ggplot2
#' @import patchwork
daps_plot <- function(res, band = 0) {
  pip <- res$pip
  sc_list <- res$sets$sets
  effect_pip <- res$info$effect_pip[, res$sets$set_index]
  coverage <- res$sets$requested_coverage
  min_abs_corr <- res$sets$min_abs_corr
  CPIP <- res$sets$coverage

  if (!is.numeric(pip)) stop("pip must be a numeric vector.")
  if (!is.list(sc_list)) stop("sc_list must be a list of index vectors.")

  p <- length(pip)
  k <- length(sc_list)

  sc_names <- names(sc_list)
  if (is.null(sc_names) || all(sc_names == "")) {
    sc_names <- paste0("L", seq_len(k))
  }

  df_pip <- data.frame(snp = seq_len(p), pip = pip)

  sc_df <- do.call(
    rbind,
    Map(function(idx, nm) {
      data.frame(snp = idx, cluster = nm, stringsAsFactors = FALSE)
    }, sc_list, sc_names)
  )

  sc_df$cluster <- factor(sc_df$cluster, levels = sc_names)

  row_step <- 1
  sc_df$row_pos <- as.numeric(factor(sc_df$cluster, levels = rev(sc_names))) * row_step

  seg_df <- unique(sc_df[, c("cluster", "row_pos")])
  seg_df$x_start <- 1
  seg_df$x_end <- p

  if (!is.null(CPIP)) {
    if (length(CPIP) != k) stop("Length of CPIP must match length(sc_list).")
    seg_df$CPIP <- CPIP[match(as.character(seg_df$cluster), sc_names)]
  }

  seg_df <- seg_df[order(seg_df$row_pos), ]
  breaks_y <- seg_df$row_pos
  labels_y <- as.character(seg_df$cluster)

  xlab2 <- if (is.null(coverage)) {
    "Signal clusters"
  } else {
    if (!is.numeric(coverage) || coverage <= 0 || coverage > 1)
      stop("coverage must be in (0, 1].")
    paste0(sprintf("%.0f", coverage * 100), "% credible sets")
  }
  if (!is.null(min_abs_corr)) {
    xlab2 <- paste0(xlab2, " with purity ≥ ", min_abs_corr, "")
  }

  sc_df$y <- sc_df$row_pos
  if (!is.null(effect_pip)) {
    if (!is.matrix(effect_pip) || nrow(effect_pip) != p || ncol(effect_pip) != k)
      stop("effect_pip must be a p x length(sc_list) matrix.")
    if (!is.numeric(band) || band < 0 || band > 1)
      stop("band must be in [0, 1].")

    col_idx <- match(as.character(sc_df$cluster), sc_names)
    w <- effect_pip[cbind(sc_df$snp, col_idx)]
    w[!is.finite(w)] <- 0
    w <- pmin(pmax(w, 0), 1)

    eps <- 1e-3
    lift_max <- band * row_step - eps
    sc_df$y <- sc_df$row_pos + lift_max * w
  }

  ring_df <- NULL
  if (nrow(sc_df) > 0) {
    if (!is.null(effect_pip)) {
      mem_df <- data.frame(
        snp = sc_df$snp,
        cluster = as.character(sc_df$cluster),
        stringsAsFactors = FALSE
      )
      mem_df$set_idx <- match(mem_df$cluster, sc_names)
      mem_df$w <- effect_pip[cbind(mem_df$snp, mem_df$set_idx)]
      mem_df$w[!is.finite(mem_df$w)] <- -Inf

      mem_df <- mem_df[order(mem_df$snp, -mem_df$w), ]
      mem_df <- mem_df[!duplicated(mem_df$snp), c("snp", "cluster")]

      ring_df <- merge(mem_df, df_pip, by = "snp", all.x = TRUE)
      ring_df$cluster <- factor(ring_df$cluster, levels = sc_names)
    } else {
      mem_df <- sc_df[, c("snp", "cluster")]
      mem_df$cluster <- as.character(mem_df$cluster)
      mem_df <- mem_df[!duplicated(mem_df$snp), ]
      ring_df <- merge(mem_df, df_pip, by = "snp", all.x = TRUE)
      ring_df$cluster <- factor(ring_df$cluster, levels = sc_names)
    }
  }

  p1 <- ggplot2::ggplot(df_pip, ggplot2::aes(x = snp, y = pip)) +
    ggplot2::geom_point(color = "black", size = 1.3) +
    {
      if (!is.null(ring_df) && nrow(ring_df) > 0) {
        ggplot2::geom_point(
          data = ring_df,
          ggplot2::aes(x = snp, y = pip, color = cluster),
          shape = 21, fill = NA, stroke = 1.0, size = 2.8
        )
      }
    } +
    ggplot2::scale_color_hue(l = 50, c = 100) +
    ggplot2::scale_x_continuous(
      limits = c(1, p),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.2f", x)
    ) +
    ggplot2::labs(x = "SNP index", y = "PIP") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(5, 5, 0, 40)
    )
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg_df,
      ggplot2::aes(x = x_start, xend = x_end,
                   y = row_pos, yend = row_pos, color = cluster),
      alpha = 0.3, linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = sc_df,
      ggplot2::aes(x = snp, y = y, color = cluster),
      size = 2.2
    ) +
    ggplot2::scale_color_hue(l = 50, c = 100) +
    {
      if (!is.null(CPIP)) {
        ggplot2::geom_text(
          data = seg_df,
          ggplot2::aes(x = Inf, y = row_pos,
                       label = sprintf("%.3f", CPIP)),
          hjust = -0.1, vjust = 0.5, size = 3
        )
      }
    } +
    ggplot2::scale_x_continuous(
      limits = c(1, p),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(min(breaks_y) - 0.6 * row_step, max(breaks_y) + 0.6 * row_step),
      breaks = breaks_y,
      labels = labels_y
    ) +
    ggplot2::labs(x = xlab2, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(size = 9),
      axis.ticks.y     = ggplot2::element_blank(),
      axis.line.y      = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(5, 40, 5, 40)
    )

  height_p2 <- max(0.8, k * 0.35)
  p1 / p2 + patchwork::plot_layout(heights = c(3, height_p2))
}

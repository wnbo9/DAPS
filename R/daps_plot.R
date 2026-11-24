#' DAP-S PIP plot with cluster / credible-set annotation
#'
#' Plot SNP-level posterior inclusion probabilities (PIPs) in a locus-plot
#' style, with SNPs colored by signal cluster (or credible set). The lead SNP
#' of each cluster is shown as a diamond and labeled. A track panel below
#' shows the members of each cluster, with its CPIP (or coverage) on the right.
#'
#' @param res A list of DAP-S fine-mapping results.
#' @param band Fraction of the row gap allowed for lifting SNPs in the track
#'  panel by their within-signal PIP, between 0 and 1; default 0.
#'
#' @return A ggplot object (patchwork layout).
#' @export
#' @import ggplot2
#' @import patchwork
daps_plot <- function(res, band = 0) {
  pip <- as.numeric(res$pip)
  snp_names <- names(res$pip)
  sc_list <- res$sets$sets
  effect_pip <- res$info$effect_pip[, res$sets$set_index, drop = FALSE]
  coverage <- res$sets$requested_coverage
  min_abs_corr <- res$sets$min_abs_corr
  CPIP <- res$sets$coverage

  if (!is.numeric(pip)) stop("pip must be a numeric vector.")
  if (!is.numeric(band) || band < 0 || band > 1) stop("band must be in [0, 1].")

  p <- length(pip)
  df_pip <- data.frame(snp = seq_len(p), pip = pip)

  base_theme <- ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      axis.line = ggplot2::element_line(color = "grey30", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "grey30", linewidth = 0.4)
    )
  x_scale <- ggplot2::scale_x_continuous(
    limits = c(1, p),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  )
  y_scale <- ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) sprintf("%.2f", x),
    expand = ggplot2::expansion(mult = c(0.03, 0.08))
  )

  # no signal clusters / credible sets: PIP panel only
  if (length(sc_list) == 0) {
    return(
      ggplot2::ggplot(df_pip, ggplot2::aes(x = snp, y = pip)) +
        ggplot2::geom_point(shape = 16, color = "grey75", size = 1.3) +
        x_scale + y_scale +
        ggplot2::labs(x = "SNP index", y = "PIP") +
        base_theme
    )
  }
  if (!is.list(sc_list)) stop("sc_list must be a list of index vectors.")

  k <- length(sc_list)
  sc_names <- names(sc_list)
  if (is.null(sc_names) || all(sc_names == "")) sc_names <- paste0("L", seq_len(k))

  if (!is.matrix(effect_pip) || nrow(effect_pip) != p || ncol(effect_pip) != k)
    stop("effect_pip must be a p x length(sc_list) matrix.")
  if (!is.null(CPIP) && length(CPIP) != k)
    stop("Length of CPIP must match length(sc_list).")

  # cluster membership with within-signal PIP
  sc_df <- do.call(rbind, Map(function(idx, j) {
    w <- effect_pip[idx, j]
    w[!is.finite(w)] <- 0
    data.frame(snp = idx, cluster = sc_names[j], w = pmin(pmax(w, 0), 1),
               stringsAsFactors = FALSE)
  }, sc_list, seq_len(k)))
  sc_df$cluster <- factor(sc_df$cluster, levels = sc_names)

  # each SNP is colored by the cluster where it has the largest within-signal PIP
  mem_df <- sc_df[order(sc_df$snp, -sc_df$w), ]
  mem_df <- mem_df[!duplicated(mem_df$snp), c("snp", "cluster")]
  mem_df$pip <- pip[mem_df$snp]
  mem_df <- mem_df[order(mem_df$pip), ]  # draw high-PIP SNPs on top

  # lead SNP of each cluster: largest within-signal PIP
  lead_df <- do.call(rbind, lapply(split(sc_df, sc_df$cluster), function(d) {
    d[which.max(d$w), c("snp", "cluster")]
  }))
  lead_df$pip <- pip[lead_df$snp]
  lead_df$label <- if (is.null(snp_names)) as.character(lead_df$snp) else
    snp_names[lead_df$snp]

  bg_df <- df_pip[!df_pip$snp %in% mem_df$snp, ]

  fill_scale <- ggplot2::scale_fill_hue(l = 50, c = 100, drop = FALSE)
  color_scale <- ggplot2::scale_color_hue(l = 50, c = 100, drop = FALSE)

  # top panel: PIP locus plot
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = bg_df, ggplot2::aes(x = snp, y = pip),
      shape = 16, color = "grey75", size = 1.3
    ) +
    ggplot2::geom_point(
      data = mem_df, ggplot2::aes(x = snp, y = pip, fill = cluster),
      shape = 21, color = "black", size = 2.4, stroke = 0.3
    ) +
    ggplot2::geom_point(
      data = lead_df, ggplot2::aes(x = snp, y = pip, fill = cluster),
      shape = 23, color = "black", size = 3.6, stroke = 0.5
    ) +
    ggplot2::geom_text(
      data = lead_df, ggplot2::aes(x = snp, y = pip, label = label),
      vjust = -1.1, size = 3, color = "grey20"
    ) +
    fill_scale + x_scale + y_scale +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = "PIP") +
    base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 45, 0, 5)
    )

  # bottom panel: cluster tracks
  row_pos <- stats::setNames(rev(seq_len(k)), sc_names)
  sc_df$row_pos <- row_pos[as.character(sc_df$cluster)]
  sc_df$y <- sc_df$row_pos + (band - 1e-3 * (band > 0)) * sc_df$w

  seg_df <- data.frame(cluster = factor(sc_names, levels = sc_names),
                       row_pos = row_pos, x_start = 1, x_end = p)
  if (!is.null(CPIP)) seg_df$CPIP <- CPIP

  xlab2 <- if (is.null(coverage) || identical(coverage, "signal clusters")) {
    "Signal clusters"
  } else {
    if (!is.numeric(coverage) || coverage <= 0 || coverage > 1)
      stop("coverage must be in (0, 1].")
    paste0(sprintf("%.0f", coverage * 100), "% credible sets")
  }
  if (!is.null(min_abs_corr)) {
    xlab2 <- paste0(xlab2, " with purity \u2265 ", min_abs_corr)
  }

  p2 <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg_df,
      ggplot2::aes(x = x_start, xend = x_end, y = row_pos, yend = row_pos,
                   color = cluster),
      alpha = 0.3, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = sc_df, ggplot2::aes(x = snp, y = y, fill = cluster),
      shape = 21, color = "black", size = 2, stroke = 0.2
    ) +
    {
      if (!is.null(CPIP)) {
        ggplot2::geom_text(
          data = seg_df,
          ggplot2::aes(x = Inf, y = row_pos,
                       label = sprintf("%.3f", floor(CPIP * 1000 + 1e-5) / 1000)),
          hjust = -0.2, vjust = 0.5, size = 3
        )
      }
    } +
    fill_scale + color_scale + x_scale +
    ggplot2::scale_y_continuous(
      limits = c(0.4, k + 0.6),
      breaks = row_pos,
      labels = sc_names
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "SNP index", y = NULL, subtitle = xlab2) +
    base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey30"),
      plot.margin = ggplot2::margin(10, 45, 5, 5)
    )

  height_p2 <- max(0.8, k * 0.35)
  p1 / p2 + patchwork::plot_layout(heights = c(3, height_p2))
}

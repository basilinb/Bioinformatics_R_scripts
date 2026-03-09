#' Make a contrast scatterplot comparing effect estimates between two contrasts
#'
#' Creates a scatterplot of estimates for two specified contrasts, classifies
#' genes by significance (based on FDR), selects top genes for labeling, and
#' returns the processed data and ggplot object.
#'
#' @param dat Data frame with columns: gene, hgnc_symbol, comparison, estimate, FDR
#' @param comp_x Character. Name of the comparison to use for x-axis (must match values in `comparison` column)
#' @param comp_y Character. Name of the comparison to use for y-axis (must match values in `comparison` column)
#' @param xlab Character. Label for x-axis (defaults to comp_x)
#' @param ylab Character. Label for y-axis (defaults to comp_y)
#' @param group_x Character. Human-friendly group name corresponding to comp_x (used in significance labels)
#' @param group_y Character. Human-friendly group name corresponding to comp_y (used in significance labels)
#' @param fdr_cutoff Numeric. FDR threshold used to call significance (default 0.05)
#' @param top_n Integer. Number of top genes per significance group to label (default 10)
#' @return A list with elements:
#'   - scatter_data: processed data used for plotting
#'   - top_genes: subset used for labels
#'   - plot: ggplot2 scatterplot object
make_contrast_scatterplot <- function(
  dat,
  comp_x,
  comp_y,
  xlab = comp_x,
  ylab = comp_y,
  group_x,
  group_y,
  fdr_cutoff = 0.05,
  top_n = 10 #,
  #limits = c(-4.4, 4.4)
) {
  # --- reshape data ---
  # keep only rows for the two comparisons, select needed columns, and pivot so each row has estimates and FDRs for both contrasts
  scatter_dat <- dat %>%
    subset(comparison %in% c(comp_y, comp_x)) %>%
    select(gene, hgnc_symbol, comparison, estimate, FDR) %>%
    pivot_wider(names_from = comparison, values_from = c(estimate, FDR)) %>%
    mutate(
      # classify each gene by whether it is significant in neither, one, or both contrasts
      significance = case_when(
        !!sym(paste0("FDR_", comp_y)) <= fdr_cutoff &
          !!sym(paste0("FDR_", comp_x)) <= fdr_cutoff ~ "Significant in both",

        !!sym(paste0("FDR_", comp_y)) <= fdr_cutoff &
          !!sym(paste0("FDR_", comp_x)) > fdr_cutoff ~ paste(
          "Significant in",
          group_y,
          "only"
        ),

        !!sym(paste0("FDR_", comp_y)) > fdr_cutoff &
          !!sym(paste0("FDR_", comp_x)) <= fdr_cutoff ~ paste(
          "Significant in",
          group_x,
          "only"
        ),

        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(significance)

  # --- top genes for labeling ---
  # For genes that are significant in at least one contrast, compute the maximum absolute effect
  # across the two contrasts and pick the top N per significance class for labeling
  top_genes <- scatter_dat %>%
    filter(significance != "Not significant") %>%
    mutate(
      max_abs_effect = pmax(
        abs(!!sym(paste0("estimate_", comp_y))),
        abs(!!sym(paste0("estimate_", comp_x)))
      )
    ) %>%
    group_by(significance) %>%
    slice_max(order_by = max_abs_effect, n = top_n) %>%
    ungroup()

  # compute symmetric plotting limits from the maximum absolute estimate across both contrasts
  limits <- max(c(
    abs(scatter_dat[[paste0("estimate_", comp_y)]]),
    abs(scatter_dat[[paste0("estimate_", comp_x)]])
  ))

  # --- plot ---
  # scatter plot of estimates (x = comp_x, y = comp_y), colored by significance,
  # with labeled top genes and reference dashed lines (zero, y=x, y=-x)
  p <- scatter_dat %>%
    ggplot(aes(
      x = !!sym(paste0("estimate_", comp_x)),
      y = !!sym(paste0("estimate_", comp_y))
    )) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    theme_classic_custom +
    labs(
      x = xlab,
      y = ylab
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = hgnc_symbol),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.4,
      max.overlaps = Inf
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
    xlim(c(-limits, limits)) +
    ylim(c(-limits, limits))

  # return processed data and plot for further inspection or saving
  return(list(
    scatter_data = scatter_dat,
    top_genes = top_genes,
    plot = p
  ))
}

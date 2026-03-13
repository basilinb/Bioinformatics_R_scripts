library(dplyr)
library(ggplot2)
library(purrr)
library(rlang)

#' Plot enrichment comparisons across multiple pairwise comparisons
#'
#' df must contain columns: group, FDR, pathway, and `k/K`.
#' comparisons: a named list; each element is length-2 character vector
#'              with the two group names that appear in df$group.
#' fdr_cutoff: include rows with FDR <= fdr_cutoff.
#' Example:
#' comparisons <- list(
#'  "Group 1:T2 High vs Group 2:T2 Low \n(Non-Asthma)" = c(
#'     "Up in T2 High in Non-Asthma",
#'     "Up in T2 Low in Non-Asthma"
#'   ),
#'   "Group 1:T2 High vs Group 2:T2 Low \n(Asthma)" = c(
#'     "Up in T2 High in Asthma",
#'     "Up in T2 Low in Asthma"
#'   ),
#'
#'   "Group 1:Asthma vs Group 2:Non-Asthma \n(T2 High)" = c(
#'     "Up in Asthma in T2 High",
#'     "Up in Non-Asthma in T2 High"
#'   )
#' )
#'
#' plot_enrichment_comparisons(
#'   df = enrichment_Significant_genes_eos1_based_T2_High_vs_Low_and_Asthma_Non_Asthma_Screening_lavage_GO_H_muppits_list$H,
#'   comparisons = comparisons,
#'   fdr_cutoff = 0.2
#' )

plot_enrichment_comparisons <- function(df, comparisons, fdr_cutoff = 0.05) {
  # ---- validate inputs ----------------------------------------------------
  req_cols <- c("group", "FDR", "pathway", "k/K")
  miss <- setdiff(req_cols, colnames(df))
  if (length(miss)) {
    stop("df is missing required columns: ", paste(miss, collapse = ", "))
  }

  if (!is.list(comparisons) || is.null(names(comparisons))) {
    stop(
      "`comparisons` must be a named list (each element a length-2 character vector)."
    )
  }

  # ---- build plotting dataframe -------------------------------------------
  # For each named comparison, keep significant rows for each of the two groups
  plot_df <- imap_dfr(comparisons, function(groups, comp_name) {
    if (length(groups) < 2) {
      stop("Comparison '", comp_name, "' must contain two group names.")
    }
    g1 <- groups[1]
    g2 <- groups[2]

    df %>%
      filter(FDR <= fdr_cutoff & (group == g1 | group == g2)) %>%
      # label which comparison this row belongs to and the "direction" using the exact group name
      mutate(
        comparison = comp_name,
        direction = group
      ) %>%
      # keep only rows for the two groups (defensive; already filtered)
      filter(group %in% c(g1, g2))
  })

  if (nrow(plot_df) == 0) {
    stop(
      "No rows passed the FDR cutoff; increase fdr_cutoff or check group names."
    )
  }

  # create a numeric size metric from the `k/K` column
  plot_df <- plot_df %>%
    mutate(size_metric = as.numeric(`k/K`))

  # preserve the order of comparisons as provided and keep pathway factor stable
  plot_df <- plot_df %>%
    mutate(
      comparison = factor(comparison, levels = unique(names(comparisons))),
      pathway = factor(pathway, levels = unique(pathway))
    )

  # ---- build and return ggplot object -------------------------------------
  p <- ggplot(
    plot_df,
    aes(x = comparison, y = pathway, size = size_metric, shape = direction)
  ) +
    geom_point(
      aes(color = direction, fill = direction),
      alpha = 0.6,
      stroke = 0.5
    ) +
    # use filled shapes so color+fill both work well
    scale_shape_manual(
      values = set_names(c(21, 22, 23, 24, 25), nm = unique(plot_df$direction))
    ) +
    scale_size(range = c(3, 9), guide = "legend") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      x = "Comparison",
      y = "Pathway",
      size = "k/K",
      shape = "Group (direction)",
      color = "Group (direction)",
      fill = "Group (direction)"
    )

  p
}

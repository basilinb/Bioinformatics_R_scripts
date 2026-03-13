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
  plot_df <- lapply(names(comparisons), function(comp_name) {
    groups <- comparisons[[comp_name]]

    g1 <- groups[1]
    g2 <- groups[2]

    g1_df <- df %>%
      filter(group == g1, FDR <= fdr_cutoff) %>%
      mutate(
        comparison = comp_name,
        direction = "Up in Group 1"
      )

    g2_df <- df %>%
      filter(group == g2, FDR <= fdr_cutoff) %>%
      mutate(
        comparison = comp_name,
        direction = "Down in Group 1"
      )

    bind_rows(g1_df, g2_df)
  }) %>%
    bind_rows()

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
    geom_point(aes(color = direction, fill = direction), alpha = 0.5) +
    scale_shape_manual(
      values = c(
        "Up in Group 1" = 24, # upward triangle
        "Down in Group 1" = 25 # downward triangle
      )
    ) +
    scale_color_manual(
      values = c(
        "Up in Group 1" = "darkred",
        "Down in Group 1" = "darkblue"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Up in Group 1" = "darkred",
        "Down in Group 1" = "darkblue"
      )
    ) +
    scale_size(range = c(3, 10)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Comparison",
      y = "Pathway",
      size = "k/K",
      shape = "Direction"
    )

  p
}

# Running GSEA geneset Plot
# Author: Basilin Benson, bbenson@benaroyaresearch.org
#
#

# © Basilin Benson 2024
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# ----------------------------------------
# GSEA plot for indivdual genesets
# ----------------------------------------
# Makes the GSEA line plot
#
# Inputs:
#   - FC_list: named list of genes and estimates from model.
#   - geneset_db: geneset to use, "H","C2","C5"
#   - geneset_name: geneset of interest e.g. "HALLMARK_INTERFERON_RESPONSE"
#   - right_label: label on left side bar eg."SE"
#   - left_label: label on right side bar eg."NSE"
#   - highlight_gene: Highlight one or more genes (character vector)
#   - right_label_col: color on left side bar eg."SE"
#   - left_label_col: color on right side bar eg."NSE"
#
# Outputs:
#   - A named list with:
#     - plot: GSEA line plot for geneset (patchwork)
#     - results: dataframe of pathway genes with columns: gene, rank, stat, leading_edge
#
# Example:
# out <- GSEAgenesetPlot(FC.ls_Severe_vs_Non_baseline$`Severe vs Non Severe `,
#                        geneset_db = "H",
#                        geneset_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#                        left_label = "NSE", right_label = "SE",
#                        highlight_gene = c("CXCL10", "ISG15"), left_label_col = "darkred",
#                        right_label_col = "black")
# out$plot
# out$results                          # full dataframe
# out$results %>% filter(leading_edge) # leading edge genes only
#
# Custom database example:
# out <- GSEAgenesetPlot(FC.ls_Severe_vs_Non_baseline$`Severe vs Non Severe `,
#                        custom_db = customDatabase,
#                        geneset_name = "custom_geneset",
#                        left_label = "NSE", right_label = "SE",
#                        highlight_gene = "CXCL10", left_label_col = "darkred",
#                        right_label_col = "black")
#
# This is what a custom database should look like:
# cluster_gs
#     gs_name       gene_symbol
#     Cluster_1      TSPAN6
#     Cluster_2        DPM1
#     Cluster_3       SCYL3

GSEAgenesetPlot <- function(
  FC_list,
  geneset_db = "H",
  geneset_name,
  custom_db = NULL,
  right_label = "Upregulated",
  left_label = "Downregulated",
  highlight_gene = NULL,
  right_label_col = "red",
  left_label_col = "blue"
) {
  require(ggplot2)
  require(msigdbr)
  require(fgsea)
  require(dplyr)
  require(stringr)
  require(patchwork)

  # Load gene set database
  if (!is.null(custom_db)) {
    if (!all(c("gs_name", "gene_symbol") %in% colnames(custom_db))) {
      stop("Custom database must have columns: gs_name and gene_symbol")
    }
    temp_gmt <- custom_db %>%
      split(.$gs_name) %>%
      lapply(function(x) x$gene_symbol)
  } else {
    temp_gmt <- msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = geneset_db
    ) %>%
      split(.$gs_name) %>%
      lapply(function(x) x$gene_symbol)
  }

  if (!(geneset_name %in% names(temp_gmt))) {
    stop("Gene set not found in the database. Check spelling or availability.")
  }

  # Prepare ranked list
  ranks <- sort(FC_list, decreasing = TRUE)

  # Run enrichment plot
  plt_enrich <- fgsea::plotEnrichment(temp_gmt[[geneset_name]], ranks) +
    labs(title = str_to_title(gsub("_", " ", geneset_name))) +
    xlab("Rank") +
    ylab("Enrichment Score")

  # If genes are highlighted, add annotations with vertical staggering
  if (!is.null(highlight_gene)) {
    valid_genes <- highlight_gene[highlight_gene %in% names(ranks)]
    if (length(valid_genes) > 0) {
      y_base <- 0.04
      y_step <- 0.035
      for (i in seq_along(valid_genes)) {
        gene <- valid_genes[i]
        gene_rank <- which(names(ranks) == gene)
        y_end <- y_base + (i - 1) * y_step
        y_text <- y_end + 0.019
        plt_enrich <- plt_enrich +
          annotate(
            "segment",
            x = gene_rank,
            xend = gene_rank + 60,
            y = 0,
            yend = y_end,
            colour = "black"
          ) +
          annotate(
            "text",
            x = gene_rank + 200,
            y = y_text,
            label = gene,
            fontface = "italic",
            size = 4.5
          )
      }
    }
  }

  # Get enrichment data used for the plot
  fgsea_enrich_data <- fgsea::plotEnrichmentData(
    pathway = temp_gmt[[geneset_name]],
    stats = ranks
  )

  rank_stats <- fgsea_enrich_data$stats

  mean_stat <- mean(rank_stats$stat)
  sd_stat <- sd(rank_stats$stat)

  # Build pathway gene dataframe with ranks, stats, and leading edge flag
  pathway_genes <- temp_gmt[[geneset_name]]
  pathway_genes_in_ranks <- pathway_genes[pathway_genes %in% names(ranks)]
  gene_ranks <- match(pathway_genes_in_ranks, names(ranks))

  # Determine leading edge from the enrichment curve
  # The curve data includes ES at each rank; leading edge genes are those
  # before the peak (positive ES) or after the trough (negative ES)
  curve_data <- fgsea_enrich_data$curve
  if (max(curve_data$ES) > abs(min(curve_data$ES))) {
    peak_rank <- curve_data$rank[which.max(curve_data$ES)]
    leading_edge_flag <- gene_ranks <= peak_rank
  } else {
    trough_rank <- curve_data$rank[which.min(curve_data$ES)]
    leading_edge_flag <- gene_ranks >= trough_rank
  }

  results_df <- data.frame(
    gene = pathway_genes_in_ranks,
    rank = gene_ranks,
    stat = ranks[pathway_genes_in_ranks],
    leading_edge = leading_edge_flag,
    row.names = NULL
  ) %>%
    arrange(rank)

  # Build bottom panels
  m <- rank_stats %>%
    ggplot(aes(x = rank, y = stat)) +
    geom_bar(stat = "identity", color = "black") +
    theme_classic() +
    geom_vline(
      xintercept = median(rank_stats$rank),
      linetype = 2,
      color = "grey25"
    ) +
    xlab("Rank") +
    ylab("Statistic")

  # Calculate range for positioning labels
  xpos <- quantile(rank_stats$rank, probs = seq(0, 1, 0.10))

  n <- rank_stats %>%
    mutate(z_score = (stat - mean_stat) / sd_stat) %>%
    ggplot(aes(x = rank, y = 1, fill = stat)) +
    geom_tile() +
    scale_fill_gradient2(
      low = right_label_col,
      mid = "white",
      high = left_label_col
    ) +
    ylab("") +
    theme_classic() +
    theme(legend.position = "none") +
    geom_text(
      aes(x = xpos[[2]], y = 1, label = left_label),
      color = left_label_col,
      fontface = "bold",
      size = 5,
      check_overlap = TRUE
    ) +
    geom_text(
      aes(x = xpos[[10]], y = 1, label = right_label),
      color = right_label_col,
      fontface = "bold",
      size = 5,
      check_overlap = TRUE
    ) +
    xlab("Rank")

  # Layout for combined plot
  des <- "
    111
    111
    111
    111
    111
    222
    333
    333
  "

  combined_plot <- plt_enrich +
    n +
    m +
    plot_layout(nrow = 3, design = des, axes = "collect")

  return(list(
    plot = combined_plot,
    results = results_df
  ))
}

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
#   - highlight_gene: Highlight a gene
#   - right_label_col: color on left side bar eg."SE"
#   - left_label_col: color on right side bar eg."NSE"
#
# Outputs:
#   - GSEA line plot for geneset
#
# Example:
# plt_gsea_ifn_alpha_test <- run_gsea_analysis(FC.ls_Severe_vs_Non_baseline$`Severe vs Non Severe `, 
#                                             geneset_db = "H", 
#                                             geneset_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#                                             left_label = "NSE",right_label = "SE",
#                                             highlight_gene = "CXCL10",left_label_col = "darkred",
#                                             right_label_col = "black")

# Print or save plot
#print(plt_gsea_ifn_alpha)



GSEAgenesetPlot <- function(FC_list, geneset_db = "H", geneset_name, 
                              right_label = "Upregulated", left_label = "Downregulated", 
                              highlight_gene = NULL,right_label_col = "red",left_label_col="blue") {
  
  # Load the gene set database using msigdbr
  temp_gmt <- msigdbr(species = "Homo sapiens", category = geneset_db) %>%
    split(.$gs_name ) %>% 
    lapply(function(x) x$gene_symbol)
  
  if (!(geneset_name %in% names(temp_gmt))) {
    stop("Gene set not found in the database. Check spelling or availability.")
  }
  
  # Prepare ranked list
  ranks <- sort(FC_list, decreasing = TRUE)
  
  # Run enrichment plot
  plt_enrich <- plotEnrichment(temp_gmt[[geneset_name]], ranks) + 
    labs(title = str_to_title(gsub("_"," ",geneset_name))) +xlab("Rank") + ylab("Enrichment Score")
  
  # If a gene is highlighted, find its rank and add annotation
  if (!is.null(highlight_gene) && highlight_gene %in% names(ranks)) {
    gene_rank <- which(names(ranks) == highlight_gene)
    gene_stat <- ranks[highlight_gene]
    
    plt_enrich <- plt_enrich +
      annotate("segment", x = gene_rank, xend = gene_rank + 60, 
               y = 0, yend = 0.04, colour = "black") +
      annotate("text", x = gene_rank + 200, y = 0.059, 
               label = highlight_gene, fontface = "italic", size = 4.5)
  }
  
  # Run enrichment statistics
  fgsea_res <- fgsea::plotEnrichmentData(pathway = temp_gmt[[geneset_name]], stats = ranks)
  
  # Custom visualization
  mean_stat <- mean(fgsea_res$stats$stat)
  sd_stat <- sd(fgsea_res$stats$stat)
  
  m <- fgsea_res$stats %>%
    ggplot(aes(x = rank, y = stat)) + 
    geom_bar(stat = "identity", color = "black") + 
    theme_classic() + 
    geom_vline(xintercept = median(fgsea_res$stats$rank), linetype = 2, color = "grey25") +xlab("Rank")+ ylab("Statistic")
  
  # Calculate range for positioning labels
  xpos<-quantile(fgsea_res$stats$rank,probs = seq(0, 1, 0.10))
  
  n <- fgsea_res$stats %>%
    mutate(z_score = (stat - mean_stat) / sd_stat) %>%
    ggplot(aes(x = rank, y = 1, fill = stat)) + 
    geom_tile() + 
    scale_fill_gradient2(low = right_label_col, mid = "white", high = left_label_col) + 
    ylab("") + theme_classic() + theme(legend.position = "none") +
    geom_text(aes(x = xpos[[2]], y = 1, label = left_label),
              color = left_label_col, fontface = "bold", size = 5, check_overlap = TRUE) +
    geom_text(aes(x = xpos[[10]], y = 1, label = right_label),
              color = right_label_col, fontface = "bold", size = 5, check_overlap = TRUE) +xlab("Rank")
  
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
  
  combined_plot <- plt_enrich + n + m + plot_layout(nrow = 3, design = des,axes = "collect")
  
  return(combined_plot)
}

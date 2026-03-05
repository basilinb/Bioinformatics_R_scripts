volcano_funcxn <- function(
  limma_output,
  coeff = NULL,
  fdr = 0.01,
  n_top_gene = 10,
  add_genes_int = NULL
) {
  for (package in c('tidyverse', 'ggplot2', 'ggrepel')) {
    if (!require(package, character.only = T, quietly = T)) {
      install.packages(package)
      library(package, character.only = T)
    }
  }

  if (is.null(coeff)) {
    coeff <- limma_output$volcano_comp %>% unique()
  } else {
    coeff <- coeff
  }

  options(ggrepel.max.overlaps = Inf)
  # Read in file generated from extract pvalue function for model of interest
  p_val <- limma_output
  vol_plt_list <- list()
  # If you want to run function for all coefficients

  # running for loop for each given coefficient
  for (i in coeff) {
    pval_coef <- subset(p_val, volcano_comp == i)
    # getting all hgnc symbol for up regulated genes
    up_reg <- head(
      subset(pval_coef, estimate > 0 & FDR < fdr) %>% arrange(FDR),
      n_top_gene
    ) %>%
      rename(delabel = hgnc_symbol) %>%
      select(gene, delabel)
    # getting all hgnc symbol for down regulated genes
    down_reg <- head(
      subset(pval_coef, estimate < 0 & FDR < fdr) %>% arrange(FDR),
      n_top_gene
    ) %>%
      rename(delabel = hgnc_symbol) %>%
      select(gene, delabel)
    add_gene_reg <- head(
      subset(pval_coef, hgnc_symbol %in% add_genes_int) %>% arrange(FDR),
      n_top_gene
    ) %>%
      rename(delabel = hgnc_symbol) %>%
      select(gene, delabel)
    # making the dataframe to plot
    pval_plot <- left_join(pval_coef, down_reg, by = "gene") %>%
      left_join(., up_reg, by = "gene") %>%
      left_join(., add_gene_reg, by = "gene") %>%
      mutate(delabel_final = coalesce(delabel.x, delabel.y, delabel)) %>%
      mutate(FC.group = ifelse(estimate > 0, "up", "down")) %>%
      mutate(
        FC.group = case_when(
          FDR <= fdr ~ FC.group,
          FDR > fdr & !is.na(delabel_final) ~ "Gene of Interest(NS)",
          FDR > fdr ~ "Not significant"
        )
      ) %>%
      #mutate(FC.group = factor(FC.group,levels = c("Not significant","Gene of Interest","down","up"))) %>%
      arrange(desc(FC.group)) %>%
      mutate(FDR = ifelse(FDR <= 0, FDR + 1e-320, FDR))
    cols <- c(
      "up" = "red",
      "down" = "blue",
      "Gene of Interest(NS)" = "steelblue4",
      "Not significant" = "grey"
    )
    # running ggplot to make volcano plot
    plt_vol <- ggplot(
      data = pval_plot,
      aes(x = estimate, y = -log10(FDR), col = FC.group)
    ) +
      geom_point(size = 3, alpha = 0.7) +
      theme_classic() +
      geom_text_repel(aes(label = delabel_final), color = "black") +
      scale_color_manual(values = cols) +
      geom_vline(xintercept = 0, col = "black", linetype = 2) + #xlim(c(-6.1,6.7)) + ylim(c(0,14.9)) +
      geom_hline(yintercept = -log10(fdr), col = "black", linetype = 2) +
      ggtitle(i)
    vol_plt_list[[i]] <- plt_vol
  }
  return(vol_plt_list)
}

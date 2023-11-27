volcano_funcxn <- function(limma_output, coeff = "all",fdr = 0.01, n_top_gene = 10) {
  for (package in c('tidyverse', 'ggplot2', 'ggrepel')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  
  
  options(ggrepel.max.overlaps = Inf)
  # Read in file generated from extract pvalue function for model of interest 
  p_val <- limma_output
  vol_plt_list <- list()
  # If you want to run function for all coefficients


    # running for loop for each given coefficient
    for (i in coeff) {
      pval_coef <- subset(p_val, variable == i)
      # getting all hgnc symbol for up regulated genes
      up_reg <-  head(subset(pval_coef, estimate > 0 & FDR < fdr) %>% arrange(FDR), n_top_gene) %>% 
        rename(delabel = gene) %>%  
        select(ensgene, delabel)
      # getting all hgnc symbol for down regulated genes
      down_reg <-  head(subset(pval_coef, estimate < 0 & FDR < fdr) %>% arrange(FDR), n_top_gene) %>% 
        rename(delabel = gene) %>%  
        select(ensgene, delabel)
      # making the dataframe to plot
      pval_plot <- left_join(pval_coef,down_reg, by="ensgene") %>% 
        left_join(.,up_reg,by="ensgene") %>% 
        mutate(delabel = coalesce(delabel.x,delabel.y)) %>% mutate(FC.group = ifelse(estimate>0,"up","down")) %>% 
        mutate(FC.group = case_when(FDR<=fdr ~ FC.group,
                                    FDR>fdr ~ "Not significant"))
      cols <- c("up"="red","down"="blue","Not significant"="grey")
      # running ggplot to make volcano plot
      plt_vol<-ggplot(data=pval_plot, aes(x=estimate, y=-log10(FDR), col=FC.group)) +
        geom_point() + 
        theme_classic() +
        geom_text_repel(aes(label=delabel),color="black") +
        scale_color_manual(values=cols) +
        geom_vline(xintercept=0, col="black",linetype=2) +
        geom_hline(yintercept=-log10(fdr), col="black",linetype=2) 
      vol_plt_list[[i]] <- plt_vol
    }
  return(vol_plt_list)
  }


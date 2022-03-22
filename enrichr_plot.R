"EnrichR Plots function

#################

Basilin Benson
University of Washington, basilinb@uw.edu
Copyright (C) 2021 Basilin Benson
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
    file_name = Thsi is the file generated from extract_pval function


OPTIONAL
  coeff = By default the function will print all coefficients from the linear model. can also specify specific coefficients e.g. c('coeff1', 'coeff2')
  fdr = You can set the FDR correction by default it is 0.2
  n_top_gene = number of up regulated and down regulated genes to label with hgnc symbol and color. default is 10
   
Example
  volcano_funcxn(file_name='results/test.pval_w_hgnc.csv',
             coeff='coeff1',
             fdr='0.05', n_top_gene=15,mod_gene='gene',
             out_dir='./')
"

# All the libraries
enrichr_plot <- function(enrich_file = NULL, type = "gene",dbs = c("MSigDB_Hallmark_2020"),Base_name = "Gene") {
  for (package in c('tidyverse', 'ggplot2', 'ggrepel','enrichR')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  plot_list <- list()
  
  for (d in dbs) {
    file_e <- enrich_file[[d]] 
    max_odds <- max(file_e[["Odds.Ratio"]])
    plt <- file_e %>% head(20) %>% 
      arrange(Adjusted.P.value) %>% 
      mutate(FDR = ifelse(Adjusted.P.value<=0.001,"<=0.001",ifelse(Adjusted.P.value>=0.001 & Adjusted.P.value<=0.01,"<=0.01",
                                                                   ifelse(Adjusted.P.value>=0.01 & Adjusted.P.value<=0.05,"<=0.05",ifelse(Adjusted.P.value>=0.05 & Adjusted.P.value<=0.1,"<=0.1",                                     ifelse(Adjusted.P.value>=0.1 & Adjusted.P.value<=0.2,"<=0.2","NS")))))) %>% 
      separate(Overlap,c("Gene_count","Size_Term"),remove = FALSE) %>% transform(Gene_count = as.numeric(Gene_count)) %>% 
      mutate(Term = fct_reorder(Term,Adjusted.P.value,.desc = TRUE)) %>%
      ggplot( aes(x=Term, y=Odds.Ratio,color=FDR)) +
      geom_segment( aes(xend=Term, yend=0)) + geom_label(aes(label=str_wrap(Overlap,1), y=0), size=3,color="black",hjust = 1) +
      geom_point( size=4) + ylim(c(-1,max_odds+0.1)) +
      scale_color_manual(values = c("<=0.001" = "red4","<=0.01"="red1","<=0.05"="orange3","<=0.1"="orange","<=0.2"="khaki2","NS"="grey")) +
      coord_flip() +
      theme_bw() +
      xlab("Enriched Terms") + ylab("Odds Ratio") + ggtitle(paste("Enrichr plot for Gene in",Base_name,"in",d,sep = " ")) 
    plot_list[[paste("plt",Base_name,d, sep="_")]] <- plt
  }
  return(plot_list)
}
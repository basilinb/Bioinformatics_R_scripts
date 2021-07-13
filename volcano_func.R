"Volcano PLot Function

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
             fdr='0.05', n_top_gene=15,
             out_dir='./')
"
# All the libraries
volcano_funcxn <- function(file_name, coeff = "all",fdr = 0.2, n_top_gene = 10, out_dir = "./") {
  for (package in c('tidyverse', 'ggplot2', 'ggrepel')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  options(ggrepel.max.overlaps = Inf)
  # Read in file generated from extract pvalue function for model of interest 
  p_val <- read_csv(file_name) %>% arrange(adj.P.Val)
  # If you want to run function for all coefficients
  if(coeff == "all") {
    coeff <- unique(p_val$group)
  }
  else {
  }
  # running for loop for each given coefficient
  for (i in coeff) {
    pval_coef <- subset(p_val, group == i)
    # getting all hgnc symbol for up regulated genes
    up_reg <-  head(subset(pval_coef, FC.group == "up" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>% 
              rename(delabel = hgnc_symbol) %>%  
              select(geneName, delabel)
    # getting all hgnc symbol for down regulated genes
    down_reg <-  head(subset(pval_coef, FC.group == "down" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>% 
                rename(delabel = hgnc_symbol) %>%  
                select(geneName, delabel)
    # making the dataframe to plot
    pval_plot <- left_join(pval_coef,down_reg, by="geneName") %>% 
                left_join(.,up_reg,by="geneName") %>% 
                mutate(delabel = coalesce(delabel.x,delabel.y)) %>% 
                mutate(FC.group = case_when(is.na(delabel) == FALSE ~ FC.group,
                                  is.na(delabel) == TRUE ~ "No"))
    cols <- c("up"="red","down"="blue","No"="black")
    # running ggplot to make volcano plot
    plt_vol<-ggplot(data=pval_plot, aes(x=logFC, y=-log10(adj.P.Val), col=FC.group, label=delabel)) +
      geom_point() + 
      theme_minimal() +
      geom_text_repel() +
      scale_color_manual(values=cols) +
      geom_vline(xintercept=0, col="black") +
      geom_hline(yintercept= 0, col="black")
    ggsave(paste(out_dir,i,"_volcano_plot.pdf"))
  }
}

"Make GSEA plots

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu
Copyright (C) 2020 Kim Dill-Mcfarland
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
  gene_list = named list object with named numeric vectors of gene symbols and logFC
  gmt_file = character giving full path the gene ontology GMT file
  or
  gmt_ls = list object with gene ontology data

OPTIONAL
  nperm = number of permutations for P-value calculations. Default is 1000
  name = character to add to output file names. Default NULL
  outdir = Output directory. Default 'results/GSEA/'
  plot = logical if should produce a plot of enrichment sources. Default FALSE
  plot.fdr = Significance cutoff for terms to include in plot. Default 0.05
  plot.groups = List of groups in which to apply FDR cutoff
  plotdir = Output directory for plot. Default the same as outdir except
            in 'figs' instead of 'results'
"

#################


GSEA_plot <- function(gsea_result,plot.fdr=0.05, plot.groups=NULL){
  #### Setup ####
  require(tidyverse)
  require(ggplot2)
  #### Plot ####
    #Filter to significant results
    if(is.null(plot.groups)){
      to.plot <- gsea_result %>% 
        filter(fgsea.FDR <= plot.fdr)
    } else {
      #First group
      terms.to.plot <- gsea_result %>% 
        filter(fgsea.FDR <= plot.fdr & group %in% plot.groups[[1]]) %>% 
        distinct(pathway) %>% unlist(use.names=FALSE)
      #Remaining groups
      for(i in c(2:length(plot.groups))){
        terms.temp <- gsea_result %>% 
          filter(fgsea.FDR <= plot.fdr & group %in% plot.groups[[i]]) %>% 
          distinct(pathway) %>% unlist(use.names=FALSE)
        
        terms.to.plot <- intersect(terms.to.plot, terms.temp)
      }
      to.plot <- gsea_result %>% 
        filter(pathway %in% terms.to.plot)
    }
    
    
    if(nrow(to.plot > 0)){
      plot.dat <- gsea_result %>% 
        #Significant terms
        filter(pathway %in% to.plot$pathway) %>% 
        #color by significance
        mutate(Significance = ifelse(fgsea.FDR <= 0.01, "FDR < 0.01",
                                     ifelse(fgsea.FDR <= 0.05, "FDR < 0.05",
                                            ifelse(fgsea.FDR <= 0.1, "FDR < 0.1",
                                                   ifelse(fgsea.FDR <= 0.2, "FDR < 0.2",
                                                          "NS"))))) 
      
      #Enrichment score limits
      plot.lim <- max(abs(plot.dat$fgsea.NES))+0.2
      
      plot <- plot.dat %>%  
        ggplot(aes(reorder(pathway, fgsea.NES), fgsea.NES)) +
        geom_segment(aes(reorder(pathway, fgsea.NES), 
                         xend=pathway, y=0, yend=fgsea.NES)) +
        geom_point(size=3, aes(fill = Significance),
                   shape=21, stroke=1) +
        geom_hline(yintercept = 0) +
        
        scale_fill_manual(values=c("FDR < 0.01"="#a50026", 
                                   "FDR < 0.05"="#f46d43",
                                   "FDR < 0.1"="#fdae61",
                                   "FDR < 0.2"="#ffffbf",
                                   "NS"="grey")) +
        lims(y=c(-plot.lim,plot.lim)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="GSEA", fill = "FGSEA significance") + 
        facet_grid( ~ group) +
        theme_bw()
    } else{
      message("No significant terms plotted.")
    }
  return(plot)
}
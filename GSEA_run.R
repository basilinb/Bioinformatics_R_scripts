"Run GSEA with fgsea

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
            
Available Geneset and how to Input in group and subgroup
Group     Sub-group
C1        
C2        CGP
C2        CP
C2        CP:BIOCARTA
C2        CP:KEGG
C2        CP:PID
C2        CP:REACTOME
C2        CP:WIKIPATHWAYS
C5        GO:BP
C5        GO:CC
C5        GO:MF
C5        HPO
C6
H
"

#################

#' Title
#'
#' @param gene_list 
#' @param gmt_ls 
#' @param nperm 
#' @param group 
#' @param subgroup 
#' @param name 
#' @param writeout 
#' @param outdir 
#'
#' @return
#' @export
#'
#' @examples



GSEA_run <- function(gene_list, gmt_ls=NULL, nperm=1000,group = NULL, subgroup = NULL, S = "human"){
  #### Setup ####
  for (package in c('tidyverse','fgsea','msigdbr')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }

  
  #Blank list to hold results
  all.results <- list()
  
  #### Data ####
  #Load gene ontology
  if(!is.null(group) & !is.null(subgroup)){
    #myGO <- fgsea::gmtPathways(gmt_file)
    gene_sets = msigdbr(species = "human", category = group, subcategory = subgroup)
    myGO <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  } else if(!is.null(group) & is.null(subgroup)){
    gene_sets = msigdbr(species = "human", category = group)
    myGO <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  }else if(!is.null(gmt_ls)){
    myGO <- gmt_ls
  } else {
    stop("Please provide gene set data as file for list object.")
  }
  
  #### Loop ####
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
    message(genes)
    #Extract 1 gene list
    genes.temp <- gene_list[[genes]]
    #Order by fold change
    genes.temp <- sort(genes.temp, decreasing = TRUE)
    
    #### FGSEA ####
    #Set score type based on fold change
    if(min(genes.temp) < 0 & max(genes.temp) > 0){
      scoreType <- "std"
    } else if(max(genes.temp) <= 0){
      scoreType <- "neg"
    } else if(min(genes.temp) >= 0){
      scoreType <- "pos"
    } else{
      stop("Could not determine score type from fold changes.")
    }
    
    #Run GSEA with fgsea
    fg.result <- fgsea::fgseaSimple(pathways = myGO, 
                                    stats = genes.temp,
                                    nperm=nperm,
                                    #eps=0,
                                    scoreType=scoreType) %>% 
      as.data.frame()
    
    
    
    
    #### Combine results ####
    gsea.result <- fg.result %>% 
      dplyr::select(pathway, pval, padj, ES, NES, size, leadingEdge) %>% 
      mutate(fgsea.FC = ifelse(NES < 0, "down","up")) %>% 
      dplyr::rename(fgsea.pval = pval, fgsea.FDR = padj, fgsea.NES = NES, 
                    fgsea.ES=ES, fgsea.size=size, fgsea.leadingEdge=leadingEdge) %>% 
      #format leading edge list
      unnest(cols = c(fgsea.leadingEdge)) %>% 
      group_by(pathway,fgsea.pval,fgsea.FDR,fgsea.ES,
               fgsea.NES,fgsea.size,fgsea.FC) %>% 
      dplyr::summarise(fgsea.leadingEdge=paste(unique(fgsea.leadingEdge), collapse=";"),
                       .groups="drop")
    
    #### Save ####
    all.results[[genes]] <- gsea.result
  }
  
  #### Format output ####
  #Unlist results into 1 df
  

  
  if(group!= "c6_GSEA.result"){
    all.results.df <- do.call(rbind.data.frame, all.results) %>% 
      rownames_to_column("group") %>% 
      mutate(group = gsub("[.][0-9]{0,4}","",group)) %>% mutate(pathway = sub("[A-Z]*_","",pathway)) %>% mutate(pathway = gsub("_"," ",pathway))
  } else {
    all.results.df <- do.call(rbind.data.frame, all.results) %>% 
      rownames_to_column("group") %>% 
      mutate(group = gsub("[.][0-9]{0,4}","",group))
  }
  
  return(all.results.df)
}
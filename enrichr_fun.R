"EnrichR Function 

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
enrichr_fun <- function(gene_mod_list = NULL, type = "gene",gene_id = "HGNC",dbs = c("MSigDB_Hallmark_2020"),mod_of_interest = "all") {
  for (package in c('tidyverse','enrichR','biomaRt')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  #Download Ensembl gene list to get HGNC symbols
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  all_genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol'), mart = ensembl) %>% filter(gene_biotype == "protein_coding")
  setEnrichrSite("Enrichr")
  enrich_list <- list()
  gene_list <- list()
  
  if (type!="gene"){
    if(mod_of_interest == "all"){
      mod_of_interest = gene_mod_list$mods %>% pull(module) %>% unique() %>% sort()
    }
  } else {
    if (gene_id == "HGNC") {
      gene_list <- gene_mod_list
    }else{
      gene_list <- all_genes %>% filter(ensembl_gene_id %in% gene_mod_list) %>% pull(hgnc_symbol)
    }
  }
  
  for (d in dbs){
  if (type == "gene"){
    enriched <- enrichr(gene_list,d)
    enrich_list <-  append(enrich_list,enriched)
  }
    else {
      for (m in mod_of_interest) {
        gene_list <- gene_mod_list$mods %>% subset(module == m) %>% pull(hgnc_symbol) 
        enriched <- enrichr(gene_list,d)
        enrich_list[[paste("module", m,sep="_")]] <- enriched
      }
    }
  }

return(enrich_list)
}
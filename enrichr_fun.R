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
    gene_mod_list = Gene list or Module object


OPTIONAL
  type = By default the function will run for gene list. Need to specify mod to run on modules
  gene_id = HGNC symbol or Ensembl ID. Default is HGNC
  dbs = Geneset database to run enrichr on. Default is MSigDB_Hallmark_2020
  mod_of_interest = If running function on module object. Which Modules do you want to run enrichr on, by default will run on all modules.
   
Example
  enrichr_run(gene_mod_list = NULL, 
  type = gene,gene_id = HGNC,dbs = c(MSigDB_Hallmark_2020),mod_of_interest = all)
"



enrichr_run <- function(gene_mod_list = NULL, type = "gene",gene_id = "HGNC",dbs = c("MSigDB_Hallmark_2020"),mod_of_interest = "all") {
  # This loop loads all the packages needed for function, if package is not installed it will automatically do it
  for (package in c('tidyverse','enrichR','biomaRt')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  #Download Ensembl gene list to get HGNC symbols
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  all_genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol'), mart = ensembl) %>% filter(gene_biotype == "protein_coding")
  #Setting the website for Enrichr to run
  setEnrichrSite("Enrichr")
  #Defining two empty lists 
  enrich_list <- list()
  gene_list <- list()
  
  #Checking if user wants to run on gene list or modules and formatting the data
  if (type!="gene"){
    if(mod_of_interest == "all"){
      mod_of_interest = gene_mod_list$mods %>% pull(module) %>% unique() %>% sort()
    }
  } else {
    # Check to see if its HGNC symbol or Ensembl ID and formatting the data 
    if (gene_id == "HGNC") {
      gene_list <- gene_mod_list
    }else{
      gene_list <- all_genes %>% filter(ensembl_gene_id %in% gene_mod_list) %>% pull(hgnc_symbol)
    }
  }
  #Loop through the databases provided by user
  for (d in dbs){
  if (type == "gene"){
    #Running Enrichr
    enriched <- enrichr(gene_list,d)
    # Append to our main list we want to return, So if we provide multiple databases this will append them to the main list. 
    enrich_list <-  append(enrich_list,enriched)
  }
    else {
      #Loop through the modules and genes within the modules
      for (m in mod_of_interest) {
        #Getting gene list from the module
        gene_list <- gene_mod_list$mods %>% subset(module == m) %>% pull(hgnc_symbol)
        #Running Enrichr 
        enriched <- enrichr(gene_list,d)
        # Append to our main list we want to return, So if we provide multiple databases this will append them to the main list.
        enrich_list[[paste("module", m,sep="_")]] <- enriched
      }
    }
  }

return(enrich_list)
}
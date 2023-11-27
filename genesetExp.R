genesetExp <- function(
    data_voom = NULL,custom_dat=NULL, genesetDb_cat="H",genesetDb_sub=NULL,genesetName=NULL,custom_geneset = NULL,custom_geneset_name="custom_geneset"){
  #currently custom geneset needs to be hgnc symbol to work
  ########## LOAD PACKAGES ############# 
  set.seed(2828)
  # Data manipulation and figures
  #data_voom = dat.voom_subset_control
  require(tidyverse, quietly = TRUE)
  require(msigdbr,quietly = TRUE)
  geneset_list <- list()
  genesetName_temp <- cusom_geneset_temp <- NULL
  geneset_median_list <-list()
  
  
  db.format <- msigdbr::msigdbr("human", category=genesetDb_cat,subcategory = genesetDb_sub)
  
  if(length(genesetName) > 0 & length(custom_geneset) > 0 ){
    cusom_geneset_temp <- data_voom$genes %>% subset(hgnc_symbol %in% custom_geneset) %>% pull(geneName)
    genesetName_temp <- genesetName
    geneset_list[[custom_geneset_name]] <- cusom_geneset_temp
  } else if (length(genesetName) > 0) {
    genesetName_temp <- genesetName
  } else if (length(custom_geneset) > 0) {
    cusom_geneset_temp <- data_voom$genes %>% subset(hgnc_symbol %in% custom_geneset) %>% pull(geneName)
    geneset_list[[custom_geneset_name]] <- cusom_geneset_temp
  } else{
    stop("Please provide a geneset name or a custom geneset")
  }
  
  for(i in genesetName_temp) {
    genes_in_geneset <- db.format %>% subset(gs_name == i) %>% 
      pull(human_ensembl_gene) %>% unique()
    geneset_list[[i]] <- genes_in_geneset
  }
  

  for(x in 1:length(geneset_list)){
  intersect_of_genes_btw_data_and_geneset <- intersect(data_voom$genes$geneName,geneset_list[[x]])
  data_voom_temp <- data_voom[intersect_of_genes_btw_data_and_geneset,]
  geneset_median_by_lib <- data_voom_temp$E %>% as.data.frame() %>% select(starts_with("lib")) %>% 
    summarise(across(everything(),median)) %>% 
    #summarise_all(funs(median)) %>% 
    t() %>% as.data.frame()  %>% rownames_to_column("libid") %>% 
    rename(median_expression = V1)
  geneset_mean_by_lib <- data_voom_temp$E %>% as.data.frame() %>% select(starts_with("lib")) %>% 
    summarise(across(everything(),mean)) %>% 
    #summarise_all(funs(median)) %>% 
    t() %>% as.data.frame()  %>% rownames_to_column("libid") %>% 
    rename(mean_expression = V1) %>% left_join(geneset_median_by_lib) %>% left_join(data_voom_temp$targets)
  geneset_median_list[[names(geneset_list)[x]]] <- geneset_mean_by_lib
  }  
 return(geneset_median_list) 
}



vizGO <- function(gsea_result,
        variable_of_interest="all" , FDR_cutoff=0.01,rrvgo_threshold=0.7
){
  #gsea_result = enrich_module_C5
  ########## LOAD PACKAGES ############# 
  set.seed(2828)
  temp_result_list <- list()
  all_result <- list()
  #variable_of_interest <- c("module_01","module_02")
  # Data manipulation and figures
  require(tidyverse, quietly = TRUE)
  # mediation analysis
  require(rrvgo, quietly = TRUE)
  # Progress tracking
  require(msigdbr, quietly = TRUE)
  
  db.format <- msigdbr::msigdbr("human", "C5") %>% select(gs_name,gs_exact_source) %>% rename(pathway = gs_name)

  
  if(str_detect(variable_of_interest,"all") %>% unique()) {
    var_temp <- gsea_result %>% pull(group) %>% unique()
  }else{
  var_temp <- variable_of_interest
}
  
  for(g_of_interest in var_temp){
    print(paste("running rrvgo for",g_of_interest,sep = " "))
    gsea_dat_temp <- gsea_result %>% subset(group == g_of_interest & FDR <= FDR_cutoff)
    if(length(gsea_dat_temp$group) > 1){
    gsea_dat_temp_id <- left_join(gsea_dat_temp,db.format)
    simMatrix_temp <- calculateSimMatrix(gsea_dat_temp_id$gs_exact_source,
                                          orgdb="org.Hs.eg.db",
                                          ont="BP",
                                          method="Rel")
    scores_temp <- setNames(-log10(gsea_dat_temp_id$qvalue), gsea_dat_temp_id$gs_exact_source)
    reducedTerms_temp <- reduceSimMatrix(simMatrix_temp,
                                         scores_temp,
                                    threshold=rrvgo_threshold,
                                    orgdb="org.Hs.eg.db")
    temp_result_list$simMatrix <- simMatrix_temp
    temp_result_list$reducedTerms <- reducedTerms_temp
    all_result[[paste("rrvgo",g_of_interest,sep = "_")]] <- temp_result_list
    }
    else{
      
      print(paste("Module",g_of_interest,"has no significant pathways at FDR",FDR_cutoff,"\n try increasing the FDR",sep = " "))
      all_result[["Failed_modules"]] <- append(all_result[["Failed_modules"]],g_of_interest)
    }
  }
  
  return(all_result)
}



# Running GAM loop 
# Author: Basilin Benson, bbenson@benaroyaresearch.org
# Started Oct 2023
#

# © Basilin Benson 2023
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#


# DESCRIPTION:
# Contains function run GAM looping through either modules or gene
#
# 
# Outputs:  list of lists and Dataframes:
#           Fixed effects = Gives the parametric effects
#           Smooth effects = Gives dataframe of smooth effects
#           Gene model summary = The GAM model summaries
#           Model_object = The GAM Model Object
#           Failed Gene List = If a gene failes while running GAM then it is added to this list
            

#---------------------------
# 
#--------------------------

###############################
#   GAM_loop   #
###############################

# REQUIRED

# voom_object  = (voom object) voom object containign RNASeq dataset in which module coherence is to be tested. Use this for gene level GAM
# module_dat_w_meta = (Dataframe) Module expression data frame with metatdata, columns with module expression where the rows libids and also has the metadata joined. Use one or the other not both
# gam_formula = GAM formula

# OPTIONAL
# geneSet (character string) = name of the column in module_gene_sets which defines modules. Defaults to "geneSet"
# module_set (character string) = name the study from which modules were built. This is used simply for labeling of outputs. Default = "STUDY"
# sample_set (character string) = name of the study from which the data come. This is used simply for labeling of outputs. Default = "STUDY"
# remove_sets (vector) = a vector of character strings naming modules which you want removed from the analysis (eg. "0").
# return_plot (logical) - logical indicating whether plot should be printed when function runs. Defaults = TRUE.
# R_cutoff (numeric) = a vector or single numeric value at which to draw red cutoff lines in correlation (R) plot. Default = 0.3
# P_cutoff (numeric) = a vector or single numeric value at which to draw red cutoff lines in significance (P) plot. Default = 0.01

# EXAMPLE USAGE
#For Gene Level
#Gene_gam <- GAM_loop(voom_object  = dat.voom_combined_w_batch2_asthma_no_treatment,
#                     gam_formula = '~treatment+s(time,k=5,bs="cr")+s(time,k=5,bs="cr",by = Infection)+ s(donor,bs="re")')

#For Module Level
#Moudle_gam <- GAM_loop(module_dat_w_meta = module_dat_asthma_no_treatment_w_meta,
#                     gam_formula = '~treatment+s(time,k=5,bs="cr")+s(time,k=5,bs="cr",by = Infection)+ s(donor,bs="re")')

########### DEFINE INPUTS ###############


GAM_loop <- function(voom_object=NULL,module_dat_w_meta = NULL,gam_formula=NULL) {
  require(tidyverse, quietly = TRUE)
  
  require(mgcv, quietly = TRUE)
  
  
  if (!is.null(module_dat_w_meta)){
    dat.temp <- module_dat_w_meta
    gene_name <- dat.temp %>% select(starts_with("module")) %>% colnames()
  }else{
    voom_obj <- voom_object
    voom_obj$targets <- voom_obj$targets %>% mutate(donorId = as.factor(donorId))
    dat.temp <- voom_obj$E %>% t() %>% as.data.frame() %>% rownames_to_column("libid") %>% 
      left_join(voom_obj$targets)
    gene_name <- voom_obj$E %>% rownames() %>% unique()
  }
  
  
  n_iter <- gene_name %>% unique() %>% length()
  
  pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  x = 1
  
  gmt_dat_s_table <- data.frame()
  gmt_dat_p_table <- data.frame()
  gmt_aic_bic_table <-data.frame() #voom_obj$genes %>% select(geneName,hgnc_symbol)
  gmt_list_summary <- list()
  gmt_list <- list()
  gmt_list_final <- list()
  failed_gene_list <- list()
  gam_form <- gam_formula #'~exacerbation+s(time,k=5,bs="cr")+s(time,k=5,bs="cr",by = exacerbation) + s(donorId,bs="re")'
  
  
  for(ens_gene in gene_name){
    formula_gam <- as.formula(paste(ens_gene,gam_form))
    tryCatch({gmt <- mgcv::gam(formula_gam,data=dat.temp)
    gmt_summ<- summary(gmt)
    gmt_p_table <- gmt_summ$p.table %>% as.data.frame() %>% rownames_to_column("variable") %>% mutate(geneName = ens_gene)
    gmt_s_table <- gmt_summ$s.table %>% as.data.frame() %>% rownames_to_column("variable") %>% mutate(geneName = ens_gene)
    gmt_dat_p_table <- rbind(gmt_dat_p_table,gmt_p_table)
    gmt_dat_s_table <- rbind(gmt_dat_s_table,gmt_s_table)
    #gmt_aic_bic_table <- rbind(gmt_aic_bic_table,
    #gene_name %>% as.data.frame() %>% rownames_to_column("geneName") %>% 
    #subset(geneName == ens_gene) %>%
    #mutate(AIC_metric = AIC(gmt)) %>% mutate(BIC_metric = BIC(gmt)))
    gmt_list_summary[[ens_gene]] <- gmt_summ
    gmt_list[[ens_gene]] <- gmt},
    error = function(e) {an.error.occured <<- ens_gene})
    failed_gene_list <- append(failed_gene_list,an.error.occured)
    #print(paste("Running gene",ens_gene,":",x,"of",n_iter))
    
    setTxtProgressBar(pb, x)
    x=x+1
  }
  close(pb)
  gmt_dat_p_table_corrected <- gmt_dat_p_table %>% group_by(variable) %>% dplyr::mutate(FDR=stats::p.adjust(`Pr(>|t|)`, method="BH")) %>% arrange(FDR)
  gmt_dat_s_table_corrected <- gmt_dat_s_table %>% group_by(variable) %>% dplyr::mutate(FDR=stats::p.adjust(`p-value`, method="BH")) %>% arrange(FDR)
  gmt_list_final[["Fixed effects"]] <- gmt_dat_p_table_corrected
  gmt_list_final[["Smooth effects"]] <- gmt_dat_s_table_corrected
  gmt_list_final[["Gene model summary"]] <- gmt_list_summary
  gmt_list_final[["Model_object"]] <- gmt_list
  gmt_list_final[["Failed Gene List"]] <- failed_gene_list
  return(gmt_list_final)
}
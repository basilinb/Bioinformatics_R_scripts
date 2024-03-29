causalFDR2 <- function(
    mediation_output, 
    method_fdr = "BH"){
  #mediation_output <- mediation.analysis_tns_w_pre_int
  ########## LOAD PACKAGES ############# 
  set.seed(2828)
  
  # Data manipulation and figures
  require(tidyverse, quietly = TRUE)
  require(stringi, quietly = TRUE)
  require(stringr, quietly = TRUE)
  
  ########## RUN FUNCTION ############# 
  
  # Assign specified contrasts 
  mediation.contrasts <-mediation_output$contrasts
  
  #-----------------------------------------------------
  # GET CONTRAST NAMES IN DF FORMAT AND LIST OUTPUT DFS
  #-----------------------------------------------------
  
  # Set up lists to store elements
  contrast_names<-list()
  contrast_output_dfs<-list()
  contrast_output_dfs_mediators<- list()
  output_dfs_contrasts<-list()
  
  # Pull formatted contrast names and the output dataframes they correspond to
  for(i in 1:length(mediation.contrasts)){
    # formatted names
    contrast_names[i]<-paste("_T-", mediation.contrasts[[i]][1], "_C-", mediation.contrasts[[i]][2], sep = "")
    # corresponding output dataframes
    contrast_output_dfs[[i]]<-names(mediation_output$summary)[which(grepl(contrast_names[i], names(mediation_output$summary)))]
    # mediators for each output
    contrast_output_dfs_mediators[[i]]<-
      unlist(lapply(contrast_output_dfs[[i]], function(x){str_split(x, "_T-post_C-pre")[[1]][1]}))
    # assign output names to mediators
    names(contrast_output_dfs_mediators[[i]])<-contrast_output_dfs[[i]]
    
    # List contrast for each df list (should be indentical)
    output_dfs_contrasts[[i]]<-rep(contrast_names[i], length(contrast_output_dfs[[i]]))
    # Assign names
    names(output_dfs_contrasts[[i]])<-contrast_output_dfs[[i]]
    
  }
  
  # Concanenate lists
  contrast_output_dfs_mediators_all<-unlist(contrast_output_dfs_mediators)
  output_dfs_contrasts_all<-unlist(output_dfs_contrasts)
  
  ##----------------------------------------##
  ##     RESHAPE AND AGGREGATE OUPUT      
  ##----------------------------------------##
  
  ## convert summaries to wide format
  summaries_transformed<-
    lapply(mediation_output$summary, 
           
           function(x){
             # Pull the dataframe name 
             df.name <- deparse(substitute(x)) 
             # Isolate contrast name
             contrast_name<- unlist(lapply(strsplit(df.name, '`', fixed = TRUE), '[', 2))
             # Isolate mediator name
             mediator_name<- paste(unlist(str_split(contrast_name, "_"))[1],unlist(str_split(contrast_name, "_"))[2],sep = "_")
             # Coerce df to widest format
             x%>%    
               as.data.frame()%>%
               rownames_to_column("Effect")%>%
               mutate(mediator = mediator_name)%>% # This doesn't currently work
               reshape2::melt(measure.vars = c("Estimate", "95% CI Lower", "95% CI Upper", "p-value"))%>%
               unite("Effect_Measure", c(Effect, variable))%>%
               spread(key = "Effect_Measure", value = "value")
           })
  
  # Assign gene names into each wide summary
  for(i in names(summaries_transformed)){
    summaries_transformed[[i]]$mediator<-contrast_output_dfs_mediators_all[[i]]
  }
  
  # Assign contrast names into each wide summary
  for(i in names(summaries_transformed)){
    summaries_transformed[[i]]$contrast <- str_sub(output_dfs_contrasts_all[[i]], 2)
  }
  
  
  # Bind summaries
  summaries_bound<-
    summaries_transformed%>%
    bind_rows()%>%
    dplyr::select(mediator, contrast, everything())
  
  
  ##----------------------------------------------##
  ##             ADJUST P VALUES
  ##----------------------------------------------##
  
  summaries_bound_adjusted<-
    summaries_bound%>%
    group_by(contrast)%>%
    mutate(`ACME (average)_p-value.adjusted` = p.adjust(`ACME (average)_p-value`, method = method_fdr))%>%
    mutate(`ACME (control)_p-value.adjusted` = p.adjust(`ACME (control)_p-value`, method = method_fdr))%>%
    mutate(`ACME (treated)_p-value.adjusted` = p.adjust(`ACME (treated)_p-value`, method = method_fdr))%>%
    mutate(`ADE (average)_p-value.adjusted` = p.adjust(`ADE (average)_p-value`, method = method_fdr))%>%
    mutate(`ADE (control)_p-value.adjusted` = p.adjust(`ADE (control)_p-value`, method = method_fdr))%>%
    mutate(`ADE (treated)_p-value.adjusted` = p.adjust(`ADE (treated)_p-value`, method = method_fdr))%>%
    mutate(`Prop. Mediated (average)_p-value.adjusted` = p.adjust(`Prop. Mediated (control)_p-value`,
                                                                  method = method_fdr))%>%
    mutate(`Prop. Mediated (control)_p-value.adjusted` = p.adjust(`Prop. Mediated (control)_p-value`, 
                                                                  method = method_fdr))%>%
    mutate(`Prop. Mediated (treated)_p-value.adjusted` = p.adjust(`Prop. Mediated (treated)_p-value`, 
                                                                  method = method_fdr))%>%
    mutate(`Total Effect_p-value.adjusted` = p.adjust(`Total Effect_p-value`, method = method_fdr))%>%
    
    # reorder dataframe
    dplyr::select(
      'mediator', 'contrast', 
      
      'ACME (average)_95% CI Lower', 'ACME (average)_95% CI Upper', 'ACME (average)_Estimate', 
      
      'ACME (average)_p-value', 'ACME (average)_p-value.adjusted', 
      
      'ACME (control)_95% CI Lower', 'ACME (control)_95% CI Upper', 'ACME (control)_Estimate', 
      
      'ACME (control)_p-value', 'ACME (control)_p-value.adjusted', 
      
      'ACME (treated)_95% CI Lower', 'ACME (treated)_95% CI Upper', 'ACME (treated)_Estimate', 
      
      'ACME (treated)_p-value', 'ACME (treated)_p-value.adjusted', 
      
      'ADE (average)_95% CI Lower', 'ADE (average)_95% CI Upper', 'ADE (average)_Estimate', 
      
      'ADE (average)_p-value', 'ADE (average)_p-value.adjusted', 
      
      'ADE (control)_95% CI Lower', 'ADE (control)_95% CI Upper', 'ADE (control)_Estimate', 
      
      'ADE (control)_p-value', 'ADE (control)_p-value.adjusted', 
      
      'ADE (treated)_95% CI Lower', 'ADE (treated)_95% CI Upper', 'ADE (treated)_Estimate', 
      
      'ADE (treated)_p-value', 'ADE (treated)_p-value.adjusted', 
      
      'Prop. Mediated (average)_95% CI Lower', 'Prop. Mediated (average)_95% CI Upper', 'Prop. Mediated (average)_Estimate', 
      
      'Prop. Mediated (average)_p-value', 'Prop. Mediated (average)_p-value.adjusted', 
      
      'Prop. Mediated (control)_95% CI Lower', 'Prop. Mediated (control)_95% CI Upper', 'Prop. Mediated (control)_Estimate', 
      
      'Prop. Mediated (control)_p-value', 'Prop. Mediated (control)_p-value.adjusted', 
      
      'Prop. Mediated (treated)_95% CI Lower', 'Prop. Mediated (treated)_95% CI Upper', 'Prop. Mediated (treated)_Estimate', 
      
      'Prop. Mediated (treated)_p-value', 'Prop. Mediated (treated)_p-value.adjusted', 
      
      'Total Effect_95% CI Lower', 'Total Effect_95% CI Upper', 'Total Effect_Estimate', 
      
      'Total Effect_p-value', 'Total Effect_p-value.adjusted'
      
    )
  
  
  
  ##---------------------------
  ##       RETURN OUTPUT
  ##---------------------------
  
  summaries_bound_adjusted
  
}


##################################
## 2) FDRcorrectedTable()   #####
##################################

# DESCRIPTION:
# Function receives dataframe output from mediationFDRcorrect() and a mediator of interest, and outputs reformatted (publication friendly) tables. 

# REQUIRED
# FDR_corrected_output = (list) An output from mediationFDRcorrect()
# mediator = (character) A mediator present in "FDR_corrected_output", for which a reformatted table is desired. 

FDRcorrectedTable = function(FDR_corrected_output, mediator_of_interest){
  
  # Subset dataframe to mediator of interest
  df_mediator<- FDR_corrected_output%>%dplyr::filter(mediator == mediator_of_interest)
  
  # Create list to save tables
  tables_formatted<-list()
  
  for(i in unique(df_mediator)$contrast){
    tables_formatted[[i]]<-df_mediator%>%dplyr::filter(contrast == i)%>%
      reshape2::melt(measure.vars = colnames(df_mediator)[-c(1,2)])%>%
      tidyr::separate(variable, c("variable", "measure"),sep = "_", remove=TRUE)%>%
      spread(key = "measure", value = "value")
  }
  
  ##---------------------------
  ##       RETURN OUTPUT
  ##---------------------------
  
  tables_formatted
  
}
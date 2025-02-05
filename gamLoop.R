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
GAMLoop <- function(voomObject = NULL, moduleDatWMeta = NULL, gamFormula = NULL) {
  require(tidyverse, quietly = TRUE)
  require(mgcv, quietly = TRUE)
  
  # Prepare input data depending on which object is provided
  if (!is.null(moduleDatWMeta)) {
    datTemp <- moduleDatWMeta
    geneName <- datTemp %>% select(starts_with("module")) %>% colnames()
  } else {
    voomObj <- voomObject
    voomObj$targets <- voomObj$targets %>% mutate(donorId = as.factor(donorId))
    
    # Transform expression data and merge with metadata
    datTemp <- voomObj$E %>% t() %>% as.data.frame() %>% rownames_to_column("libid") %>% 
      left_join(voomObj$targets)
    
    # Extract unique gene names
    geneName <- voomObj$E %>% rownames() %>% unique()
  }
  
  # Determine number of iterations (genes to process)
  nIter <- length(unique(geneName))
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 1, max = nIter, style = 3, width = 50, char = "=")
  
  # Initialize empty result containers
  gmtDatSTable <- data.frame()
  gmtDatPTable <- data.frame()
  gmtAicBicTable <- data.frame()
  gmtListSummary <- list()
  gmtList <- list()
  gmtListFinal <- list()
  failedGeneList <- list()
  
  # Define GAM formula
  gamForm <- gamFormula 
  
  x <- 1  # Iteration counter
  
  for (ensGene in geneName) {
    formulaGAM <- as.formula(paste(ensGene, gamForm))
    
    # Try fitting GAM model and catching any errors
    tryCatch({
      gmt <- mgcv::gam(formulaGAM, data = datTemp)
      gmtSumm <- summary(gmt)
      
      # Extract fixed effects p-values
      gmtPTable <- gmtSumm$p.table %>% as.data.frame() %>% rownames_to_column("variable") %>% mutate(geneName = ensGene)
      
      # Extract smooth effects p-values
      gmtSTable <- gmtSumm$s.table %>% as.data.frame() %>% rownames_to_column("variable") %>% mutate(geneName = ensGene)
      
      # Append results to main dataframes
      gmtDatPTable <- rbind(gmtDatPTable, gmtPTable)
      gmtDatSTable <- rbind(gmtDatSTable, gmtSTable)
      
      # Store model summary and object in lists
      gmtListSummary[[ensGene]] <- gmtSumm
      gmtList[[ensGene]] <- gmt
    }, error = function(e) {
      # Capture failed genes
      failedGeneList <- append(failedGeneList, ensGene)
    })
    
    # Update progress bar
    setTxtProgressBar(pb, x)
    x <- x + 1
  }
  
  # Close progress bar
  close(pb)
  
  # Apply multiple testing correction (FDR) to fixed effects
  gmtDatPTableCorrected <- gmtDatPTable %>% group_by(variable) %>% 
    dplyr::mutate(FDR = stats::p.adjust(`Pr(>|t|)`, method = "BH")) %>% 
    arrange(FDR)
  
  # Apply multiple testing correction (FDR) to smooth effects
  gmtDatSTableCorrected <- gmtDatSTable %>% group_by(variable) %>% 
    dplyr::mutate(FDR = stats::p.adjust(`p-value`, method = "BH")) %>% 
    arrange(FDR)
  
  # Store results in final output list
  gmtListFinal[["FixedEffects"]] <- gmtDatPTableCorrected
  gmtListFinal[["SmoothEffects"]] <- gmtDatSTableCorrected
  gmtListFinal[["GeneModelSummary"]] <- gmtListSummary
  gmtListFinal[["ModelObject"]] <- gmtList
  gmtListFinal[["FailedGeneList"]] <- failedGeneList
  
  return(gmtListFinal)
}

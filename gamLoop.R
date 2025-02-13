# Running GAM loop 
# Author: Basilin Benson, bbenson@benaroyaresearch.org
# 
#

# © Basilin Benson 2024
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#


# DESCRIPTION:
# Contains function run GAM looping through either modules or gene

# ----------------------------------------
# Generalized Additive Model (GAM) Loop
# ----------------------------------------
# This function fits a GAM model to each gene/module in the dataset
# and extracts key statistical outputs for further analysis.
#
# Inputs:
#   - voomObject: (optional) A voom-transformed object containing gene expression data.
#   - moduleDatwMeta: (optional) A data frame with gene/module data and metadata.
#   - gamFormula: A string specifying the GAM formula (e.g., "~ s(Time, bs='cs')").
#
# Outputs:
#   - `Fixed Effects`: Data frame with fixed-effect p-values and adjusted FDR.
#   - `Smooth Effects`: Data frame with smooth-effect p-values and adjusted FDR.
#   - `Gene Model Summary`: List of GAM model summaries for each gene.
#   - `Model Object`: List of fitted GAM models.
#   - `Failed Gene List`: Vector of genes that failed model fitting.
#   - `ExecutionTime`: Time taken to run the function.


GAMLoop <- function(voomObject = NULL, moduleDatwMeta = NULL, gamFormula = NULL) {
  # Load required packages
  require(tidyverse, quietly = TRUE)  # For data wrangling
  require(mgcv, quietly = TRUE)       # For GAM modeling
  require(data.table, quietly = TRUE) # For efficient list-to-dataframe conversion
  
  startTime <- Sys.time()  # Start time tracking
  
  # Determine which dataset to use: module-level or gene-level
  if (!is.null(moduleDatwMeta)) {
    # If a module data matrix is provided, use it
    datTemp <- moduleDatwMeta
    geneNames <- colnames(datTemp %>% select(starts_with("module")))  # Extract module names
  } else {
    # Otherwise, process gene expression data from voomObject
    voomObject$targets <- voomObject$targets %>% mutate(donorId = as.factor(donorId))  # Ensure donor ID is a factor
    datTemp <- voomObject$E %>% t() %>% as.data.frame() %>%
      rownames_to_column("libid") %>% left_join(voomObject$targets, by = "libid")  # Merge metadata
    geneNames <- rownames(voomObject$E)  # Get gene names from the expression matrix
  }
  
  # Initialize a progress bar to track computation progress
  nIter <- length(geneNames)
  pb <- txtProgressBar(min = 1, max = nIter, style = 3, width = 50, char = "=")
  
  # Initialize preallocated memory lists for outputs
  gmtPTableList <- vector("list", nIter)  # Stores p-values for fixed effects
  gmtSTableList <- vector("list", nIter)  # Stores p-values for smooth terms
  gmtSummaryList <- vector("list", nIter) # Stores full model summaries
  gmtModelList <- vector("list", nIter)   # Stores full GAM model objects
  failedGenes <- c()       # Track genes that fail during model fitting
  processedGenes <- c()    # Track genes that are successfully processed
  
  gamForm <- gamFormula  # Store the GAM formula passed by the user
  
  # Loop through each gene and fit a GAM model
  for (i in seq_along(geneNames)) {
    ensGene <- geneNames[i]  # Get the gene/module name
    # Construct the GAM formula dynamically
    formulaGam <- as.formula(paste(ensGene, gamForm))  
    
    tryCatch({
      # Fit the GAM model
      gmt <- mgcv::gam(formulaGam, data = datTemp)
      gmtSummary <- summary(gmt)  # Get model summary
      
      # Extract and store p-values for fixed effects
      gmtPTableList[[ensGene]] <- gmtSummary$p.table %>%
        as.data.frame() %>% rownames_to_column("variable") %>%
        mutate(geneName = ensGene)
      
      # Extract and store p-values for smooth terms
      gmtSTableList[[ensGene]] <- gmtSummary$s.table %>%
        as.data.frame() %>% rownames_to_column("variable") %>%
        mutate(geneName = ensGene)
      
      # Store the full model summary and model object
      gmtSummaryList[[i]] <- gmtSummary
      gmtModelList[[i]] <- gmt  
      
      # Add successfully processed gene to the tracking list
      processedGenes <- c(processedGenes, ensGene)  
      
    }, error = function(e) {
      # If an error occurs, add the gene to the failed list
      failedGenes <<- append(failedGenes, ensGene)  
      return(NULL)
    })
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)  # Close progress bar
  
  
  # Remove NULL entries from lists before binding
  gmtSummaryList <- gmtSummaryList[!sapply(gmtSummaryList, is.null)]
  gmtModelList <- gmtModelList[!sapply(gmtModelList, is.null)]
  
  # Convert stored p-value lists into data frames
  gmtDatPTable <- data.table::rbindlist(gmtPTableList, fill = TRUE)
  gmtDatSTable <- data.table::rbindlist(gmtSTableList, fill = TRUE)
  
  # Adjust p-values for multiple comparisons using the Benjamini-Hochberg (BH) method
  gmtDatPTable <- gmtDatPTable %>%
    group_by(variable) %>%
    mutate(FDR = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
    arrange(FDR)
  
  gmtDatSTable <- gmtDatSTable %>%
    group_by(variable) %>%
    mutate(FDR = p.adjust(`p-value`, method = "BH")) %>%
    arrange(FDR)
  
  # Assign names only to successfully processed genes
  names(gmtSummaryList) <- processedGenes
  names(gmtModelList) <- processedGenes
  
  # Calculate total execution time
  endTime <- Sys.time()
  executionTime <- endTime - startTime
  
  # Return a structured list containing all results
  return(list(
    `FixedEffects` = gmtDatPTable,   # Processed fixed effects results
    `SmoothEffects` = gmtDatSTable,  # Processed smooth effects results
    `GeneModelSummary` = gmtSummaryList,  # Full model summaries
    `ModelObject` = gmtModelList,    # Fitted model objects
    `FailedGeneList` = failedGenes,  # List of genes that failed to fit
    `ExecutionTime` = executionTime  # Total time taken
  ))
}



runLmerLoop <- function(expressionDatawMeta, lmerformula=NULL,nameCol = NULL,expCol=NULL) {
  require(lmerTest, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(tibble, quietly = TRUE)
  
  # Error handling: Ensure formula and colName are provided
  if (is.null(formula) && is.null(nameCol) && is.null(expCol)) {
    stop("Error: Please provide both an lmer formula and a gene/module column name.")
  }
  if (is.null(formula)) {
    stop("Error: Please provide an lmer formula.")
  }
  if (is.null(nameCol)) {
    stop("Error: Please provide a gene/module column name.")
  }
  
  if (is.null(expCol)) {
    stop("Error: Please provide a gene/module expression column name.")
  }

  
  # Extract unique module names
  varNames <- unique(expressionDatawMeta[[nameCol]])
  
  # Preallocate list to store results
  resultsList <- vector("list", length(varNames))
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 1, max = length(varNames), style = 3, width = 50, char = "=")
  
  # Loop through each module and fit the LMER model
  for (i in seq_along(varNames)) {
    varName <- varNames[i]
    
    # Define formula dynamically
    formula <- as.formula(paste0(expCol, lmerformula))
    
    # Fit the mixed effects model
    dat_sub <- expressionDatawMeta %>% subset(get(nameCol) == varName)
    modelSummary <- summary(lmerTest::lmer(formula, data = dat_sub))
    
    # Extract coefficient table and format output
    resultsList[[i]] <- modelSummary$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      mutate(gene = varName)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)  # Close progress bar after loop finishes
  
  # Combine all results into a single data frame
  resultsDf <- bind_rows(resultsList)
  #Run multiple testing correction
  resultsDf <- resultsDf %>% rename("pval" = "Pr(>|t|)","estimate"="Estimate") %>%
    group_by(variable) %>%
    mutate(FDR = stats::p.adjust(pval, method = "BH")) %>%
    arrange(FDR)
  
  return(resultsDf)
}

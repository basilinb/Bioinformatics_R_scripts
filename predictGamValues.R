# Running GAM Prediction
# Author: Basilin Benson, bbenson@benaroyaresearch.org
# 
#

# © Basilin Benson 2024
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# ----------------------------------------
# Generalized Additive Model (GAM) Prediction
# ----------------------------------------
# This function predicts values from fitted GAM models, excluding specified terms 
# (e.g., random effects) and returning the output in a wide matrix format.
#
# Inputs:
#   - moduleGAM: A list containing fitted GAM models (output from `GAMLoop`).
#   - excludeTerms: Character vector specifying terms to exclude from prediction 
#                   (default: `c("s(donorId)")`).
#
# Outputs:
#   - A matrix where:
#       - Rows represent genes/modules.
#       - Columns correspond to treatment-infection-time combinations.
#       - Values are predicted GAM fits.
#
# Details:
#   - Uses `tidymv::predict_gam()` to compute predictions while excluding 
#     specified terms (e.g., random effects).
#   - If prediction fails for a gene, it is skipped, and a message is displayed.
#   - The final output is transformed into a wide-format matrix for easy downstream use.
predictGamValues <- function(moduleGAM, excludeTerms = c("s(donorId)"),
                             treat = NULL, time = NULL) {
  require(tidymv, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  
  # Initialize list to store predicted values
  predictedList <- vector("list", length(moduleGAM$ModelObject))
  
  # Loop through gene models and predict values
  for (i in seq_along(predictedList)) {
    gene <- names(moduleGAM$ModelObject)[i]
    model <- moduleGAM$ModelObject[[gene]]
    #print(gene)  # Print progress
    
    # Predict and handle potential errors
    tryCatch({
      predictedList[[i]] <- tidymv::predict_gam(
        model, exclude_terms = excludeTerms, values = list(donorId = NULL)
      ) %>%
        mutate(geneName = gene)
    }, error = function(e) {
      message(paste("Prediction failed for gene:", gene))
    })
  }
  
  # Combine results into a single data frame
  predictedValues <- data.table::rbindlist(predictedList)
  
  # Convert to wide format
  predictedValuesWide <- predictedValues %>%
    select(geneName, treat, time, fit) %>%
    pivot_wider(names_from = c(treat, time), values_from = fit) %>%
    column_to_rownames("geneName") %>%
    as.matrix()
  
  return(predictedValuesWide)
}

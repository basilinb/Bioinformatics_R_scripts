# ----------------------------------------
# Generalized Additive Model (GAM) Loop
# ----------------------------------------
# Author: Basilin Benson  <bbenson@benaroyaresearch.org>
# License: GNU General Public License (GPL v3 or later)
# © Basilin Benson 2024
#
# This software is licensed under the GNU General Public License. You may modify
# and/or redistribute it under GPL version 3 or later. For license details, visit:
# https://www.gnu.org/licenses/gpl-3.0.html
#
# DESCRIPTION:
# This function loops over gene-level or module-level expression data and fits
# a Generalized Additive Model (GAM) to each gene/module using the provided formula.
# It supports flexible output options including fixed effects, smooth effects,
# model summaries, and the fitted model objects themselves. Failed model fits
# are tracked and reported.
#
# This function is useful for modeling longitudinal gene expression or
# module-level summaries in time-course or repeated-measures designs.
#
# PARAMETERS:
#   - voomObject:
#       A voom-transformed object (from limma) containing the expression matrix (E)
#       and associated metadata in `targets`. Used when moduleDatwMeta is not provided.
#
#   - moduleDatwMeta:
#       A data frame containing module-level expression data and metadata.
#       If provided, the function will use this instead of voomObject.
#
#   - gamFormula:
#       A string containing the formula to be passed to mgcv::gam.
#       Example: "~ s(Time, bs = 'cs') + donorId"
#
#   - returnFixedEffects:
#       Logical. If TRUE, returns a data frame of p-values and FDR for fixed effects.
#
#   - returnSmoothEffects:
#       Logical. If TRUE, returns a data frame of p-values and FDR for smooth terms.
#
#   - returnGeneModelSummary:
#       Logical. If TRUE, returns a named list of GAM summary objects for each gene/module.
#
#   - returnModelObject:
#       Logical. If TRUE, returns a named list of fitted GAM model objects.
#
# RETURNS:
# A named list containing:
#   - FixedEffects:
#       Data frame of fixed-effect terms and FDR values (if returnFixedEffects = TRUE).
#
#   - SmoothEffects:
#       Data frame of smooth-effect terms and FDR values (if returnSmoothEffects = TRUE).
#
#   - GeneModelSummary:
#       List of GAM summary objects for each gene/module (if returnGeneModelSummary = TRUE).
#
#   - ModelObject:
#       List of fitted GAM model objects (if returnModelObject = TRUE).
#
#   - FailedGeneList:
#       Character vector of genes or modules that failed model fitting.
#
#   - ExecutionTime:
#       The total time taken to complete the GAM fitting loop.
#
# EXAMPLE USAGE:
#   gamResults <- GAMLoop(
#     voomObject = v,
#     gamFormula = "~ s(Timepoint, bs = 'cs') + Treatment + donorId",
#     returnFixedEffects = TRUE,
#     returnSmoothEffects = TRUE,
#     returnGeneModelSummary = TRUE,
#     returnModelObject = TRUE
#   )

GAMLoop <- function(
  voomObject = NULL,
  moduleDatwMeta = NULL,
  gamFormula = NULL,
  returnFixedEffects = TRUE,
  returnSmoothEffects = TRUE,
  returnGeneModelSummary = TRUE,
  returnModelObject = TRUE
) {
  require(tidyverse, quietly = TRUE)
  require(mgcv, quietly = TRUE)
  require(data.table, quietly = TRUE)

  startTime <- Sys.time()

  if (!is.null(moduleDatwMeta)) {
    datTemp <- moduleDatwMeta
    geneNames <- colnames(datTemp %>% select(starts_with("module")))
  } else {
    voomObject$targets <- voomObject$targets %>%
      mutate(donorId = as.factor(donorId))
    datTemp <- voomObject$E %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("libid") %>%
      left_join(voomObject$targets, by = "libid")
    geneNames <- rownames(voomObject$E)
  }

  nIter <- length(geneNames)
  pb <- txtProgressBar(min = 1, max = nIter, style = 3, width = 50, char = "=")

  gmtPTableList <- if (returnFixedEffects) vector("list", nIter) else NULL
  gmtSTableList <- if (returnSmoothEffects) vector("list", nIter) else NULL
  gmtSummaryList <- if (returnGeneModelSummary) vector("list", nIter) else NULL
  gmtModelList <- if (returnModelObject) vector("list", nIter) else NULL
  failedGenes <- c()

  print("Starting GAM loop:")

  for (i in seq_along(geneNames)) {
    ensGene <- geneNames[i]
    formulaGam <- as.formula(paste(ensGene, gamFormula))

    tryCatch(
      {
        gmt <- mgcv::gam(formulaGam, data = datTemp)

        # ── Compute summary once, reuse for both effect tables ──────────────
        if (
          returnFixedEffects || returnSmoothEffects || returnGeneModelSummary
        ) {
          gmtSummary <- summary(gmt)
        }

        if (returnFixedEffects) {
          gmtPTableList[[i]] <- gmtSummary$p.table %>%
            as.data.frame() %>%
            rownames_to_column("variable") %>%
            mutate(geneName = ensGene)
        }

        if (returnSmoothEffects) {
          gmtSTableList[[i]] <- gmtSummary$s.table %>%
            as.data.frame() %>%
            rownames_to_column("variable") %>%
            mutate(geneName = ensGene)
        }

        if (returnGeneModelSummary) {
          gmtSummaryList[[i]] <- gmtSummary
        }

        if (returnModelObject) {
          gmtModelList[[i]] <- gmt
        }
      },
      error = function(e) {
        failedGenes <<- append(failedGenes, ensGene)
      }
    )

    setTxtProgressBar(pb, i)
  }

  close(pb)

  # ── Identify which genes succeeded BEFORE removing NULLs ───────────────────
  # Use the model list if available, otherwise fall back to summary list.
  # This is the source of truth for which indices actually succeeded.
  ref_list <- if (returnModelObject) gmtModelList else gmtSummaryList
  succeededIdx <- which(!sapply(ref_list, is.null))

  # ── Remove NULLs and assign names from the correct indices ─────────────────
  if (returnGeneModelSummary) {
    gmtSummaryList <- gmtSummaryList[succeededIdx]
    names(gmtSummaryList) <- geneNames[succeededIdx]
  }

  if (returnModelObject) {
    gmtModelList <- gmtModelList[succeededIdx]
    names(gmtModelList) <- geneNames[succeededIdx]
  }

  # ── Bind data frames (NULLs in rbindlist are handled via fill = TRUE) ───────
  gmtDatPTable <- if (returnFixedEffects) {
    data.table::rbindlist(gmtPTableList, fill = TRUE)
  } else {
    NULL
  }

  gmtDatSTable <- if (returnSmoothEffects) {
    data.table::rbindlist(gmtSTableList, fill = TRUE)
  } else {
    NULL
  }

  if (returnFixedEffects) {
    gmtDatPTable <- gmtDatPTable %>%
      group_by(variable) %>%
      mutate(FDR = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
      arrange(FDR)
  }

  if (returnSmoothEffects) {
    gmtDatSTable <- gmtDatSTable %>%
      group_by(variable) %>%
      mutate(FDR = p.adjust(`p-value`, method = "BH")) %>%
      arrange(FDR)
  }

  executionTime <- Sys.time() - startTime

  results <- list(
    FailedGeneList = failedGenes,
    ExecutionTime = executionTime
  )
  if (returnFixedEffects) {
    results[["FixedEffects"]] <- gmtDatPTable
  }
  if (returnSmoothEffects) {
    results[["SmoothEffects"]] <- gmtDatSTable
  }
  if (returnGeneModelSummary) {
    results[["GeneModelSummary"]] <- gmtSummaryList
  }
  if (returnModelObject) {
    results[["ModelObject"]] <- gmtModelList
  }

  return(results)
}

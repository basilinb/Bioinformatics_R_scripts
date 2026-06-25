#' Run linear mixed-effects or linear models across genes/modules
#'
#' Fits an lmer (or lm) model independently for each unique gene or module in
#' the input data, extracts coefficient tables, and applies BH FDR correction
#' within each model term.
#'
#' @param expressionDatawMeta A data frame containing gene/module expression
#'   values, metadata covariates, and a grouping column (e.g. gene name).
#' @param lmerformula A character string for the right-hand side of the model
#'   formula, e.g. \code{"~ viral_load + sex + (1|donorId)"}. The expression
#'   column is prepended automatically.
#' @param nameCol Character. Column name in \code{expressionDatawMeta} that
#'   identifies each gene or module to iterate over.
#' @param expCol Character. Column name for the expression/response variable.
#' @param use_lmer Logical. If \code{TRUE} (default), fits
#'   \code{lmerTest::lmer}; if \code{FALSE}, fits \code{stats::lm}.
#' @param keep_models Logical. If \code{TRUE}, returns a list containing both
#'   the results data frame and a named list of fitted model objects. If
#'   \code{FALSE} (default), returns only the results data frame.
#'
#' @return If \code{keep_models = FALSE}, a tibble with columns \code{variable},
#'   \code{estimate}, \code{Std. Error}, \code{df} (lmer only), \code{t value},
#'   \code{pval}, \code{gene}, and \code{FDR}. If \code{keep_models = TRUE}, a
#'   named list with elements \code{resultsDf} and \code{modelObjList}.
#'
#' @details Each gene/module is subset from \code{expressionDatawMeta} and fit
#'   independently. FDR correction is applied per model term (grouped by
#'   \code{variable}) using \code{stats::p.adjust(method = "BH")}. Note that
#'   storing model objects (\code{keep_models = TRUE}) can be very
#'   memory-intensive, as each lmer object retains a copy of its input data.
#'
#' @examples
#' # Basic usage with lmer
#' results <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   use_lmer = TRUE
#' )
#'
#' # Fixed-effects only (lm)
#' results_lm <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ viral_load + sex",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   use_lmer = FALSE
#' )
#'
#' # Keep model objects for diagnostics
#' results_full <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   keep_models = TRUE
#' )
#' results_full$resultsDf
#' results_full$modelObjList[["IFNL1"]]
#'
#' @importFrom lmerTest lmer
#' @importFrom data.table rbindlist
#' @importFrom dplyr rename group_by mutate arrange
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom stats as.formula lm p.adjust
#' @export
runLmerLoop <- function(
  expressionDatawMeta,
  lmerformula = NULL,
  nameCol = NULL,
  expCol = NULL,
  use_lmer = TRUE,
  keep_models = FALSE
) {
  require(lmerTest, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(tibble, quietly = TRUE)
  require(data.table, quietly = TRUE)

  if (is.null(lmerformula) && is.null(nameCol) && is.null(expCol)) {
    stop(
      "Error: Please provide both an lmer formula and a gene/module column name."
    )
  }
  if (is.null(lmerformula)) {
    stop("Error: Please provide an lmer formula.")
  }
  if (is.null(nameCol)) {
    stop("Error: Please provide a gene/module column name.")
  }
  if (is.null(expCol)) {
    stop("Error: Please provide a gene/module expression column name.")
  }

  varNames <- unique(expressionDatawMeta[[nameCol]])

  resultsList <- vector("list", length(varNames))
  if (keep_models) {
    modelObjList <- vector("list", length(varNames))
  }

  pb <- txtProgressBar(
    min = 1,
    max = length(varNames),
    style = 3,
    width = 50,
    char = "="
  )

  for (i in seq_along(varNames)) {
    varName <- varNames[i]
    formula <- as.formula(paste0(expCol, lmerformula))
    dat_sub <- expressionDatawMeta %>% subset(get(nameCol) == varName)

    if (use_lmer) {
      model <- lmerTest::lmer(formula, data = dat_sub)
    } else {
      model <- lm(formula, data = dat_sub)
    }

    modelSummary <- summary(model)

    resultsList[[i]] <- modelSummary$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      mutate(gene = varName)

    if (keep_models) {
      modelObjList[[varName]] <- model
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  resultsDf <- data.table::rbindlist(resultsList) %>%
    as_tibble() %>%
    rename("pval" = "Pr(>|t|)", "estimate" = "Estimate") %>%
    group_by(variable) %>%
    mutate(FDR = stats::p.adjust(pval, method = "BH")) %>%
    arrange(FDR)

  if (keep_models) {
    return(list("resultsDf" = resultsDf, "modelObjList" = modelObjList))
  } else {
    return(resultsDf)
  }
}

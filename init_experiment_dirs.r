#' Initialize experiment subdirectory structure (with checks & warnings)
#'
#' Create a consistent set of experiment-specific folders under a project's
#' plots and results directories. Accepts either the list returned by
#' init_project_dirs() (uses $plots and $results) or a bare path string.
#'
#' If experiment_name is missing or NULL a warning is issued and folders are
#' created directly in the parent plots/results paths. At least one of
#' plot_subdirs or result_subdirs must be non-NULL.
#'
#' @param parent_dir list or character. Either the list from init_project_dirs()
#'   (must contain $plots and $results) or a base path string to use for both
#'   plots and results.
#' @param experiment_name character or NULL. Name for the experiment folder
#'   (e.g. "IL13_dose"). If NULL or missing, a warning is issued and the
#'   experiment-level folders are the parent plots/results paths.
#' @param plot_subdirs character or NULL. Vector of subdirectory names to create
#'   under the experiment plots directory. Must be a non-empty character vector
#'   when provided.
#' @param result_subdirs character or NULL. Vector of subdirectory names to
#'   create under the experiment results directory. Must be a non-empty
#'   character vector when provided.
#' @return invisible named list of created paths (plots, results, and entries
#'   named plots_<subdir> / results_<subdir> when provided).
#' @examples
#' proj <- init_project_dirs("P513-5")
#'
#' exp_dirs <- init_experiment_dirs(
#'   proj,
#'   "IL13_dose",
#'   plot_subdirs  = c("volcano", "heatmap", "enrichment",
#'                     "scatterplot", "boxplot", "GAM"),
#'   result_subdirs = c("GSEA", "Norm")
#' )
#'
#' plot_exp_dirs <- init_experiment_dirs(
#'   proj$plots,
#'   "IL13_dose",
#'   plot_subdirs = c("volcano", "heatmap", "enrichment",
#'                    "scatterplot", "boxplot", "GAM")
#' )
#' @export

init_experiment_dirs <- function(
  parent_dir,
  experiment_name = NULL,
  plot_subdirs = NULL,
  result_subdirs = NULL
) {
  # --- Resolve plots_base / results_base ---------------------------------
  # Accept the list from init_project_dirs() (must have $plots and $results)
  # OR a bare path string (used as the base for whichever side is needed).
  if (is.list(parent_dir) && !is.data.frame(parent_dir)) {
    if (!all(c("plots", "results") %in% names(parent_dir))) {
      stop(
        "When parent_dir is a list it must contain both $plots and $results ",
        "(i.e. the object returned by init_project_dirs())"
      )
    }
    plots_base <- parent_dir$plots
    results_base <- parent_dir$results
  } else if (is.character(parent_dir) && length(parent_dir) == 1) {
    # Only assign the base we actually need to avoid creating phantom dirs
    plots_base <- if (!is.null(plot_subdirs)) parent_dir else NULL
    results_base <- if (!is.null(result_subdirs)) parent_dir else NULL
  } else {
    stop(
      "parent_dir must be either the list returned by init_project_dirs() ",
      "or a single character path string"
    )
  }

  # --- experiment_name ---------------------------------------------------
  # Default is NULL; missing() check is no longer needed
  if (is.null(experiment_name)) {
    warning(
      "No experiment name provided — folders will be created in the ",
      "parent plots/results directories"
    )
    exp_plots <- plots_base
    exp_results <- results_base
  } else {
    if (
      !is.character(experiment_name) ||
        length(experiment_name) != 1 ||
        nchar(experiment_name) == 0
    ) {
      stop("experiment_name must be a single, non-empty string or NULL")
    }
    exp_plots <- if (!is.null(plots_base)) {
      file.path(plots_base, experiment_name)
    } else {
      NULL
    }
    exp_results <- if (!is.null(results_base)) {
      file.path(results_base, experiment_name)
    } else {
      NULL
    }
  }

  # --- Require at least one subdir vector --------------------------------
  if (is.null(plot_subdirs) && is.null(result_subdirs)) {
    stop(
      "Both plot_subdirs and result_subdirs are NULL — ",
      "specify at least one non-NULL subdirectory vector"
    )
  }

  # --- Validate and build sub-paths --------------------------------------
  if (!is.null(plot_subdirs)) {
    if (!is.character(plot_subdirs) || length(plot_subdirs) == 0) {
      stop("plot_subdirs must be a non-empty character vector or NULL")
    }
    plot_sub_paths <- file.path(exp_plots, plot_subdirs)
  } else {
    plot_sub_paths <- character(0)
  }

  if (!is.null(result_subdirs)) {
    if (!is.character(result_subdirs) || length(result_subdirs) == 0) {
      stop("result_subdirs must be a non-empty character vector or NULL")
    }
    result_sub_paths <- file.path(exp_results, result_subdirs)
  } else {
    result_sub_paths <- character(0)
  }

  # --- Create directories ------------------------------------------------
  dirs <- unique(c(exp_plots, exp_results, plot_sub_paths, result_sub_paths))
  dirs <- dirs[!is.null(dirs) & nzchar(dirs)] # drop any NULLs / empty strings

  lapply(dirs, \(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

  # --- Return named list of paths ----------------------------------------
  out <- list()
  if (!is.null(exp_plots)) {
    out$plots <- exp_plots
  }
  if (!is.null(exp_results)) {
    out$results <- exp_results
  }

  if (length(plot_sub_paths) > 0) {
    out <- c(
      out,
      setNames(as.list(plot_sub_paths), paste0("plots_", plot_subdirs))
    )
  }
  if (length(result_sub_paths) > 0) {
    out <- c(
      out,
      setNames(as.list(result_sub_paths), paste0("results_", result_subdirs))
    )
  }

  invisible(out)
}

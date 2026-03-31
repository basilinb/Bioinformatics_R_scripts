#' Initialize project directory structure
#'
#' Create a standard project folder layout.
#' The function creates (if missing) the project root and subfolders:
#' Figs, Results, Data_Clean, and Data_Raw.
#'
#' @param project_number character. Project identifier used as the folder name
#'   under the lab folder (e.g. "P513-5").
#' @param box_base character. Base path to Box storage
#'   (e.g. "~/Library/CloudStorage/Box-Box/BRI_Projects").
#' @return invisible list with absolute paths to created directories:
#'   root, plots, results, data_clean, data_raw.
#' @examples
#' proj <- init_project_dirs(
#'   project_number = "P513-5",
#'   box_base = "~/Library/CloudStorage/Box-Box/BRI_Projects"
#' )
#' @export
init_project_dirs <- function(
  project_number = NULL,
  box_base = NULL
) {
  # Validate box_base
  if (
    is.null(box_base) ||
      !is.character(box_base) ||
      length(box_base) != 1 ||
      !nzchar(box_base)
  ) {
    stop(
      "box_base must be a single non-empty path string ",
      "(e.g. '~/Library/CloudStorage/Box-Box/BRI_Projects')"
    )
  }

  # Validate project_number
  if (
    is.null(project_number) ||
      !is.character(project_number) ||
      length(project_number) != 1 ||
      !nzchar(project_number)
  ) {
    stop(
      "project_number must be a single non-empty string (e.g. 'P100-1')"
    )
  }

  # Build paths
  dir_root <- file.path(box_base, project_number)
  dir_plots <- file.path(dir_root, "Figs")
  dir_results <- file.path(dir_root, "Results")
  dir_data_clean <- file.path(dir_root, "Data_Clean")
  dir_data_raw <- file.path(dir_root, "Data_Raw")

  # Create directories (recursive = TRUE handles dir_root itself)
  dirs <- c(dir_root, dir_plots, dir_results, dir_data_clean, dir_data_raw)
  lapply(dirs, \(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

  # Return paths invisibly so callers can use them
  invisible(list(
    root = dir_root,
    plots = dir_plots,
    results = dir_results,
    data_clean = dir_data_clean,
    data_raw = dir_data_raw
  ))
}

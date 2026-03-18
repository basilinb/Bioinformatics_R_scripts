library(tidyverse)
library(pak)
# Get all installed packages (user library only, excluding base/recommended)
pkgs <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
user_pkgs <- pkgs[is.na(pkgs$Priority), "Package"] # excludes base & recommended

# Read DESCRIPTION metadata for each package
get_pkg_info <- function(pkg) {
  desc <- tryCatch(packageDescription(pkg), error = function(e) NULL)
  if (is.null(desc)) {
    return(NULL)
  }

  data.frame(
    package = pkg,
    remote_type = desc$RemoteType %||% NA_character_,
    remote_user = desc$RemoteUsername %||% NA_character_,
    remote_repo = desc$RemoteRepo %||% NA_character_,
    remote_ref = desc$RemoteRef %||% NA_character_,
    bioc_views = desc$biocViews %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# %||% is the "null coalescing" operator — define it if not already loaded
`%||%` <- function(x, y) if (is.null(x)) y else x

info <- lapply(user_pkgs, get_pkg_info)

# Drop any NULLs (packages where packageDescription() failed entirely)
info <- Filter(Negate(is.null), info)

info_df <- do.call(rbind, info)

info_df$bioc_views[info_df$bioc_views == ""] <- NA_character_


# GitHub/GitLab: has a remote_type we want to install from source
gh_pkgs <- info_df[
  !is.na(info_df$remote_type) &
    info_df$remote_type %in% c("github", "gitlab"),
]

# xgit packages — incomplete metadata, flag for manual review
xgit_pkgs <- info_df[
  !is.na(info_df$remote_type) &
    info_df$remote_type == "xgit",
]

# "standard" remote type = CRAN mirror, treat as CRAN
standard_pkgs <- info_df[
  !is.na(info_df$remote_type) &
    info_df$remote_type == "standard",
]

# Bioconductor: no remote_type (or not github/gitlab) but has bioc_views
bioc_pkgs <- info_df[
  is.na(info_df$remote_type) &
    !is.na(info_df$bioc_views),
]

# CRAN: everything else with no remote_type and no bioc_views
cran_pkgs <- info_df[
  is.na(info_df$remote_type) &
    is.na(info_df$bioc_views),
]

# Sanity check — these should sum to nrow(info_df)
cat("GitHub:  ", nrow(gh_pkgs), "\n")
cat("xgit:    ", nrow(xgit_pkgs), "\n")
cat("standard:", nrow(standard_pkgs), "\n")
cat("Bioc:    ", nrow(bioc_pkgs), "\n")
cat("CRAN:    ", nrow(cran_pkgs), "\n")
cat(
  "Total:   ",
  nrow(gh_pkgs) +
    nrow(xgit_pkgs) +
    nrow(standard_pkgs) +
    nrow(bioc_pkgs) +
    nrow(cran_pkgs),
  "vs",
  nrow(info_df),
  "\n"
)

# These need manual attention — print package names and look them up
print(xgit_pkgs$package)
xgit_pkgs <- xgit_pkgs %>%
  mutate(
    remote_type = "github",
    remote_user = "BenaroyaResearch",
    remote_repo = package,
    remote_ref = "HEAD"
  )


# Merge with github packages (same install method) and add a column for the full GitHub repo spec
gh_pkgs <- bind_rows(gh_pkgs, xgit_pkgs)


#Build the refs
cran_refs <- c(cran_pkgs$package, standard_pkgs$package)
bioc_refs <- paste0("bioc::", bioc_pkgs$package)
gh_refs <- paste0(gh_pkgs$remote_user, "/", gh_pkgs$remote_repo)

# For gh_refs, add @ref if it's not HEAD
gh_refs <- ifelse(
  !is.na(gh_pkgs$remote_ref) & gh_pkgs$remote_ref != "HEAD",
  paste0(gh_refs, "@", gh_pkgs$remote_ref),
  gh_refs
)


all_refs <- c(cran_refs, bioc_refs, gh_refs)


# Save
saveRDS(all_refs, "~/pkg_list.rds")
writeLines(all_refs, "~/pkg_list.txt")

# I wanted to trim the list of packages so I manually did it
all_refs_trimmed <- read_csv("Package_install/packages_updated.csv") %>%
  pull(package) %>%
  unique()
saveRDS(all_refs_trimmed, "Package_install/pkg_list.rds")
writeLines(all_refs_trimmed, "Package_install/pkg_list.txt")
# Reinstall all packages from the list
# Install pak first
install.packages("pak")

# Load your saved list
all_refs <- readRDS("Package_install/pkg_list.rds")

#Remove NAs if any
all_refs <- all_refs[!is.na(all_refs)]

# set repos / Bioc version and check library paths
options(repos = BiocManager::repositories())
options(pkg.bioc_version = BiocManager::version())
.libPaths()
install.packages("pak", repos = "https://cloud.r-project.org")
install.packages("BiocManager", repos = "https://cloud.r-project.org")


# split list into CRAN/Bioc vs GitHub remotes and install in two phases
gh <- grep("/", all_refs, value = TRUE) # GitHub entries (user/repo)
cran_bioc <- setdiff(all_refs, gh) # CRAN/Bioconductor entries
cran_bioc_clean <- gsub("^bioc::", "", cran_bioc)
# install CRAN/Bioc first (no interactive prompts)
pak::pkg_install(cran_bioc_clean)

# then install GitHub remotes (one-by-one or in small batches)
pak::pkg_install(gh)
library(dplyr)
`%notin%` <- Negate(`%in%`)
cran_bioc_clean_n <- cran_bioc_clean[
  cran_bioc_clean %notin% c("duckdb", "tm", "BiocSingular")
]
pak::pkg_install(cran_bioc_clean_n)
pak

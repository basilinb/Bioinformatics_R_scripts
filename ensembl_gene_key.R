#' Reproducible, version-pinned Ensembl protein-coding gene keys
#'
#' Replaces the hand-rolled `if (file.exists("key_pc_mouse.csv"))` blocks that
#' were duplicated across projects. Three properties this provides that the old
#' pattern did not:
#'
#' 1. The Ensembl release is pinned explicitly, so the annotation cannot drift
#'    when Ensembl publishes a new release (~quarterly).
#' 2. The release is part of the cache filename, so bumping it forces a real
#'    re-download instead of silently reusing a stale file, and previous
#'    releases stay on disk so old results remain reproducible.
#' 3. The release is recorded in the data itself, in an `ensembl_release`
#'    column derived from the same variable passed to `version =`, so the
#'    filename, the pin and the data can never disagree.
#'
#' Requires: HoardR, biomaRt, dplyr, curl, jsonlite. Assumes dplyr/magrittr
#' are attached (uses `%>%`), which is true for every calling document here.

# ---- 1. Live release detection ---------------------------------------------

#' Latest Ensembl release, or NA if it cannot be determined
#'
#' Tries three independent sources and takes the first that answers.
#' rest.ensembl.org intermittently returns HTTP 500 -- observed live, the same
#' URL 500'd and then succeeded minutes later -- so a single source is not
#' reliable enough to base a staleness warning on. All three sources were
#' verified to agree on release 116.
#'
#' biomaRt::listEnsemblArchives() is deliberately NOT used: it probes every
#' archive host in turn and times out (observed hanging on
#' jun2026.archive.ensembl.org).
#'
#' The answer is memoised for `cache_hours`, so a typical run makes no network
#' call at all and Ensembl is not hammered.
#'
#' @param timeout numeric seconds allowed per source. Default 10. A value <= 0
#'   means "no time budget", and returns NA immediately without any request --
#'   note libcurl treats CURLOPT_TIMEOUT = 0 as *unlimited*, so this is handled
#'   explicitly rather than passed through
#' @param cache_hours numeric hours a previous answer stays good. Default 24.
#'   Pass 0 to force a fresh lookup (needed when testing the network path).
#' @param dir directory for the memo file. Defaults to the key cache dir, or
#'   tempdir() if that is unset -- the memo is a convenience, not critical state
#' @param sources list of list(url=, parse=) pairs. Defaults to
#'   .ensembl_release_sources; exposed so the fallback behaviour is testable
#'
#' @return integer Ensembl release, or NA_integer_ if every source failed
#'
#' @examples
#' ensembl_live_release()                              # memoised, no network
#' ensembl_live_release(cache_hours = 0)               # force a real lookup
#' ensembl_live_release(timeout = 0, cache_hours = 0)  # NA, no request made
ensembl_live_release <- function(timeout = 10, cache_hours = 24,
                                 dir = getOption("ensembl_key_dir", tempdir()),
                                 sources = .ensembl_release_sources) {
  memo <- file.path(dir, ".ensembl_live_release")
  if (cache_hours > 0 && file.exists(memo) &&
      difftime(Sys.time(), file.mtime(memo), units = "hours") < cache_hours) {
    v <- suppressWarnings(as.integer(readLines(memo, warn = FALSE)[1]))
    if (!is.na(v)) return(v)
  }

  # No time budget means no check. Guarded explicitly because libcurl reads
  # CURLOPT_TIMEOUT = 0 as "wait forever", the opposite of the intent.
  if (!is.finite(timeout) || timeout <= 0) return(NA_integer_)

  ms <- max(1L, as.integer(round(timeout * 1000)))
  fetch <- function(url) {
    h <- curl::new_handle(timeout_ms = ms, connecttimeout_ms = ms,
                          followlocation = TRUE)
    r <- curl::curl_fetch_memory(url, handle = h)
    if (r$status_code != 200L) stop("HTTP ", r$status_code, call. = FALSE)
    rawToChar(r$content)
  }

  for (s in sources) {
    v <- tryCatch(s$parse(fetch(s$url)), error = function(e) NA_integer_)
    if (length(v) == 1L && !is.na(v) && v > 0L) {
      try(writeLines(as.character(v), memo), silent = TRUE)
      return(as.integer(v))
    }
  }
  NA_integer_
}

# Three independent sources, tried in order. rest.ensembl.org intermittently
# returns HTTP 500, so no single source is load-bearing. Assigned after the
# function so it can be swapped in tests: .ensembl_release_sources[[1]]$url <- ...
.ensembl_release_sources <- list(
  # {"releases":[116]}
  list(url   = "https://rest.ensembl.org/info/data/?content-type=application/json",
       parse = function(txt) max(as.integer(jsonlite::fromJSON(txt)$releases))),
  # {"release":116}
  list(url   = "https://rest.ensembl.org/info/software?content-type=application/json",
       parse = function(txt) as.integer(jsonlite::fromJSON(txt)$release)),
  # directory listing containing release-116/
  list(url   = "https://ftp.ensembl.org/pub/",
       parse = function(txt) max(as.integer(sub("release-", "",
                 regmatches(txt, gregexpr("release-[0-9]+", txt))[[1]]))))
)

# ---- 2. Species table and cache location -----------------------------------

# The mouse/human copy-paste, collapsed into data.
.ensembl_species <- list(
  mouse = list(dataset = "mmusculus_gene_ensembl", symbol = "mgi_symbol"),
  human = list(dataset = "hsapiens_gene_ensembl",  symbol = "hgnc_symbol")
)

#' Canonical cache directory for Ensembl gene keys
#'
#' Absolute, so it does not depend on the working directory -- knitting (wd =
#' the document's folder) and running interactively (wd = whichever .Rproj is
#' open) must resolve to the same place.
#'
#' Deliberately not here::here(): a project with two .Rproj files resolves
#' here() differently depending on which one is open, which is precisely the
#' bug this is meant to prevent.
#'
#' Set with options(ensembl_key_dir = "...") or the ENSEMBL_KEY_DIR env var.
#'
#' @return absolute path, as a character scalar
ensembl_key_dir <- function() {
  d <- getOption("ensembl_key_dir", Sys.getenv("ENSEMBL_KEY_DIR", ""))
  if (!nzchar(d)) {
    stop("Set options(ensembl_key_dir = ) to an absolute cache directory ",
         "before requesting a gene key.", call. = FALSE)
  }
  path.expand(d)
}

# ---- 3. The expensive, cached step: the full unfiltered key -----------------

#' Unfiltered BioMart gene key for one species at one pinned Ensembl release
#'
#' Cached as .rds only. It is an intermediate rather than a share artifact, and
#' .rds preserves NA exactly -- readr writes "" as an empty field but NA as the
#' literal string "NA", and reads both back as NA, so a CSV cache would return
#' a different object on a hit than on a miss.
#'
#' The cheap protein-coding filter is deliberately left OUT of the cache, so
#' filter policy can change without re-downloading.
#'
#' @param species "mouse" or "human"
#' @param release integer Ensembl release. Becomes part of the cache filename
#' @param dir cache directory. Default ensembl_key_dir()
#' @param overwrite logical force a re-download. Default FALSE
#'
#' @return data.frame of geneName, <symbol>, gene_biotype, ensembl_release
get_ensembl_full_key <- function(species, release, dir = ensembl_key_dir(),
                                 overwrite = FALSE) {
  species <- match.arg(species, names(.ensembl_species))
  stopifnot(length(release) == 1L, !is.na(release))
  release <- as.integer(release)
  spec <- .ensembl_species[[species]]

  HoardR::cache(
    paths = HoardR::make_path(
      dir      = dir,
      project  = paste0("ensembl", release),   # the release IS the cache key
      analysis = paste0(species, "_full_gene_key"),
      ext      = "rds"
    ),
    {
      mart <- biomaRt::useEnsembl(
        biomart = "genes",
        dataset = spec$dataset,
        version = release            # the ONE place a release is chosen
      )
      message("BioMart host for release ", release, ": ", mart@host)

      biomaRt::getBM(
        attributes = c("ensembl_gene_id", spec$symbol, "transcript_biotype"),
        mart       = mart
      ) %>%
        # rename to match count data
        dplyr::rename(geneName = "ensembl_gene_id",
                      gene_biotype = "transcript_biotype") %>%
        # biomaRt returns "" (never NA) for a gene with no assigned symbol.
        # Normalise once, here, so a cache hit and a fresh download agree.
        dplyr::mutate(dplyr::across(dplyr::all_of(spec$symbol),
                                    ~ dplyr::na_if(trimws(.x), ""))) %>%
        # Provenance DERIVED from the same variable passed to version= above,
        # so the pin, the filename and the data cannot disagree.
        dplyr::mutate(ensembl_release = release)
    },
    overwrite = overwrite
  )
}

# ---- 4. The cheap derived artifact -----------------------------------------

#' Protein-coding gene key for one species at one pinned Ensembl release
#'
#' @param species "mouse" or "human"
#' @param release integer Ensembl release
#' @param dir cache directory. Default ensembl_key_dir()
#' @param drop_missing_symbol logical drop genes with no assigned symbol.
#'   Default FALSE. For left_join consumers dropping is byte-for-byte neutral
#'   (those genes yield NA either way), but it would remove ~592 human / ~94
#'   mouse real protein-coding genes from any `%in%` style gene universe.
#' @param check_release logical compare the pin against the live release and
#'   warn if behind. Default TRUE
#' @param overwrite logical force a re-download. Default FALSE
#'
#' @return data.frame of geneName, <symbol>, gene_biotype, ensembl_release,
#'   one row per protein-coding gene, sorted by geneName
#'
#' @examples
#' options(ensembl_key_dir = "~/gene_annotation")
#' key.pc_mouse <- get_ensembl_pc_key("mouse", release = 115L)
get_ensembl_pc_key <- function(species, release, dir = ensembl_key_dir(),
                               drop_missing_symbol = FALSE,
                               check_release = TRUE, overwrite = FALSE) {
  species <- match.arg(species, names(.ensembl_species))
  spec    <- .ensembl_species[[species]]
  release <- as.integer(release)

  if (isTRUE(check_release)) {
    live <- ensembl_live_release()
    if (is.na(live)) {
      message("Could not reach Ensembl to check for a newer release; ",
              "skipping staleness check. Pinned to ", release, ".")
    } else if (live > release) {
      warning("Pinned Ensembl release ", release, " is behind the live release ",
              live, ". This is intentional and reproducible. To update, set the ",
              "release to ", live, " and re-run: a new cache file is written and ",
              "the old one is kept.", call. = FALSE)
    } else {
      message("Ensembl release check: pinned release ", release, " is current.")
    }
  }

  out <- get_ensembl_full_key(species, release, dir, overwrite = overwrite) %>%
    # Keep protein coding genes only
    dplyr::filter(gene_biotype == "protein_coding") %>%
    # Remove duplicate entries
    dplyr::distinct(geneName, .keep_all = TRUE) %>%
    dplyr::arrange(geneName)

  if (isTRUE(drop_missing_symbol)) {
    n0  <- nrow(out)
    out <- out[!is.na(out[[spec$symbol]]), , drop = FALSE]
    message("Dropped ", n0 - nrow(out), " ", species,
            " protein-coding genes with no ", spec$symbol, ".")
  }
  out
}

# ---- 5. Consumer-side helpers ----------------------------------------------

#' Infer species from a vector of Ensembl gene IDs
#'
#' Inspects ALL ids, not just the first. The common
#' `startsWith(rownames(counts)[1], "ENSMU")` idiom checks only row 1, which
#' misfires on custom-genome datasets where a custom contig such as "SS78_CAR"
#' can be the first row -- a mouse dataset would then be handed the human key.
#'
#' @param gene_ids character vector of gene IDs
#' @return "mouse" or "human"
ensembl_species_from_ids <- function(gene_ids) {
  ids <- stats::na.omit(as.character(gene_ids))
  if (!length(ids)) stop("No gene IDs supplied.", call. = FALSE)
  if (any(startsWith(ids, "ENSMUS"))) return("mouse")
  if (any(startsWith(ids, "ENSG")))   return("human")
  stop("Could not infer species from gene IDs (e.g. '", ids[1], "'). ",
       "Expected ENSMUSG... (mouse) or ENSG... (human).", call. = FALSE)
}

#' Narrow a gene key to exactly geneName + gene_symbol
#'
#' This is the contract that un-keyed `left_join(key)` calls depend on: they
#' natural-join on shared column names, so the key must expose only the join
#' column and the payload. Carrying `ensembl_release` on the full key is safe
#' precisely because consumers narrow through this function first.
#'
#' @param key a key from get_ensembl_pc_key()
#' @return data.frame with columns geneName, gene_symbol
ensembl_gene_symbol_map <- function(key) {
  sym <- intersect(c("mgi_symbol", "hgnc_symbol"), names(key))
  if (length(sym) != 1L) {
    stop("Expected exactly one of mgi_symbol / hgnc_symbol in the key; found: ",
         paste(sym, collapse = ", "), call. = FALSE)
  }
  key %>%
    dplyr::select(dplyr::all_of(c("geneName", sym))) %>%
    dplyr::rename(gene_symbol = dplyr::all_of(sym))
}

#' Pinned-vs-live release table, for printing in a rendered report
#'
#' Prints on every run so that "no warning" is never mistaken for "current" --
#' a failed check reports status "unknown", not "current".
#'
#' @param release the pinned release
#' @return one-row data.frame
ensembl_release_report <- function(release) {
  live <- ensembl_live_release()
  data.frame(
    pinned_release = as.integer(release),
    live_release   = live,
    status = if (is.na(live)) "unknown (could not reach Ensembl)"
             else if (live <= release) "current"
             else paste0("stale by ", live - release, " release(s)"),
    checked_on = as.character(Sys.Date())
  )
}

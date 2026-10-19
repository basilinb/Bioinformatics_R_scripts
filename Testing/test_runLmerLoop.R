# Verification script for runLmerLoop.R
# Run from the repo root:  Rscript Testing/test_runLmerLoop.R
#
# Plain stopifnot() assertions — the repo has no testthat harness. The
# reference is the pre-matrix/parallel version of runLmerLoop from git
# (commit c4e5600), so default calls must reproduce it exactly.

suppressPackageStartupMessages({
  library(dplyr)
  library(lmerTest)
})

oldEnv <- new.env()
eval(
  parse(text = system("git show c4e5600:runLmerLoop.R", intern = TRUE)),
  envir = oldEnv
)
runLmerLoop_old <- oldEnv$runLmerLoop
source("runLmerLoop.R")

pass <- function(label) cat("  PASS  ", label, "\n")
section <- function(label) cat("\n== ", label, " ==\n", sep = "")

# Silence progress bars and expected messages/warnings
quiet <- function(expr) {
  capture.output(res <- suppressWarnings(suppressMessages(expr)))
  res
}

# Expect a warning/message matching a pattern while evaluating expr
expect_condition <- function(expr, pattern, type = c("warning", "message")) {
  type <- match.arg(type)
  seen <- FALSE
  handler <- function(c) {
    if (grepl(pattern, conditionMessage(c))) seen <<- TRUE
    invokeRestart(if (type == "warning") "muffleWarning" else "muffleMessage")
  }
  capture.output(
    res <- if (type == "warning") {
      withCallingHandlers(suppressMessages(expr), warning = handler)
    } else {
      withCallingHandlers(suppressWarnings(expr), message = handler)
    }
  )
  stopifnot(seen)
  res
}

# Compare coefficient tables by gene + variable, ignoring row order
same_coefs <- function(a, b, tol = 1e-8) {
  key <- function(d) {
    d <- as.data.frame(ungroup(d))
    d$gene <- as.character(d$gene)
    d[order(d$gene, d$variable), ]
  }
  a <- key(a)
  b <- key(b)
  rownames(a) <- rownames(b) <- NULL
  isTRUE(all.equal(a, b, tolerance = tol))
}

# ── Simulated data ───────────────────────────────────────────────────────────
set.seed(42)
nGenes <- 150
nDonors <- 20
nTimes <- 3
targets <- data.frame(
  libId = sprintf("lib%03d", seq_len(nDonors * nTimes)),
  donorId = factor(rep(sprintf("d%02d", seq_len(nDonors)), each = nTimes)),
  x = rnorm(nDonors * nTimes),
  grp = factor(rep(c("A", "B", "C"), length.out = nDonors * nTimes)),
  stringsAsFactors = FALSE
)
rownames(targets) <- targets$libId
donorEffect <- matrix(
  rnorm(nGenes * nDonors, sd = 0.7),
  nGenes,
  nDonors
)[, as.integer(targets$donorId)]
E <- 5 +
  outer(rnorm(nGenes, sd = 0.5), targets$x) +
  donorEffect +
  matrix(rnorm(nGenes * nrow(targets)), nGenes, nrow(targets))
rownames(E) <- sprintf("GENE%03d", seq_len(nGenes))
colnames(E) <- targets$libId
voomObj <- list(E = E, targets = targets)

to_long <- function(E, targets) {
  data.frame(
    geneName = rep(rownames(E), times = ncol(E)),
    geneExp = as.vector(E),
    targets[rep(seq_len(ncol(E)), each = nrow(E)), ],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
longDf <- to_long(E, targets)

lmForm <- "~ x + grp"
lmerForm <- "~ x + grp + (1|donorId)"
longArgs <- list(
  expressionDatawMeta = longDf,
  nameCol = "geneName",
  expCol = "geneExp"
)
run_old <- function(...) quiet(do.call(runLmerLoop_old, c(longArgs, list(...))))
run_new <- function(...) quiet(do.call(runLmerLoop, c(longArgs, list(...))))

# ── 0. Backward compatibility ────────────────────────────────────────────────
section("0. Backward compatibility (default arguments)")

oldLm <- run_old(lmerformula = lmForm, use_lmer = FALSE)
stopifnot(identical(run_new(lmerformula = lmForm, use_lmer = FALSE), oldLm))
pass("lm loop output identical to old version")

oldLmer <- run_old(lmerformula = lmerForm)
stopifnot(identical(run_new(lmerformula = lmerForm), oldLmer))
pass("lmer loop output identical to old version")

oldKeep <- run_old(lmerformula = lmerForm, keep_models = TRUE)
newKeep <- run_new(lmerformula = lmerForm, keep_models = TRUE)
stopifnot(
  identical(names(newKeep), c("resultsDf", "modelObjList")),
  identical(newKeep$resultsDf, oldKeep$resultsDf),
  identical(names(newKeep$modelObjList), rownames(E)),
  all(vapply(
    rownames(E),
    function(g) {
      isTRUE(all.equal(
        lme4::fixef(newKeep$modelObjList[[g]]),
        lme4::fixef(oldKeep$modelObjList[[g]])
      ))
    },
    logical(1)
  ))
)
pass("keep_models: same results, models by name match (old list had 2n entries, new has n)")

longFactor <- longDf
longFactor$geneName <- factor(longFactor$geneName)
stopifnot(identical(
  quiet(runLmerLoop(longFactor, lmForm, "geneName", "geneExp", use_lmer = FALSE)),
  quiet(runLmerLoop_old(longFactor, lmForm, "geneName", "geneExp", use_lmer = FALSE))
))
pass("factor nameCol: identical to old version")

longBad <- longDf
longBad$geneExp[longBad$geneName == "GENE010"] <- NA
oldErr <- tryCatch(
  quiet(runLmerLoop_old(longBad, lmForm, "geneName", "geneExp", use_lmer = FALSE)),
  error = function(e) "error"
)
newErr <- tryCatch(
  quiet(runLmerLoop(longBad, lmForm, "geneName", "geneExp", use_lmer = FALSE)),
  error = function(e) conditionMessage(e)
)
stopifnot(identical(oldErr, "error"), grepl("GENE010", newErr))
pass("a failing gene still stops the call (and now names the gene)")

# ── 1-3. Matrix lm ───────────────────────────────────────────────────────────
section("1-3. Matrix lm")

matVoom <- quiet(runLmerLoop(voomObj = voomObj, lmerformula = lmForm, use_lmer = FALSE))
stopifnot(same_coefs(matVoom, oldLm), identical(names(matVoom), names(oldLm)))
pass("1. voom input matches lm gene by gene")

matLong <- run_new(lmerformula = lmForm, use_lmer = FALSE, sampleCol = "libId")
stopifnot(same_coefs(matLong, oldLm), identical(names(matLong), names(oldLm)))
pass("2. long df + sampleCol matches lm gene by gene")

stopifnot(identical(
  quiet(runLmerLoop(
    voomObj = voomObj, lmerformula = lmForm, use_lmer = FALSE,
    genes = rownames(E)[1:20]
  ))$gene %>% unique() %>% sort(),
  rownames(E)[1:20]
))
pass("genes argument subsets the fit")

# ── 4. NA and rank handling ──────────────────────────────────────────────────
section("4. NA and rank handling")

naCovTargets <- targets
naCovTargets$x[c(3, 17)] <- NA
naCovLong <- to_long(E, naCovTargets)
stopifnot(same_coefs(
  quiet(runLmerLoop(
    voomObj = list(E = E, targets = naCovTargets), lmerformula = lmForm,
    use_lmer = FALSE
  )),
  quiet(runLmerLoop_old(naCovLong, lmForm, "geneName", "geneExp", use_lmer = FALSE))
))
pass("NA covariate drops that library, as lm does")

naExpE <- E
naExpE[5, 3] <- NA
naExpE[40, c(1, 9)] <- NA
naExpOld <- quiet(runLmerLoop_old(
  to_long(naExpE, targets), lmForm, "geneName", "geneExp", use_lmer = FALSE
))
naExpNew <- expect_condition(
  runLmerLoop(
    voomObj = list(E = naExpE, targets = targets), lmerformula = lmForm,
    use_lmer = FALSE
  ),
  "2 gene\\(s\\) have missing expression",
  "message"
)
stopifnot(same_coefs(naExpNew, naExpOld))
pass("NA expression: those genes fall back to lm, results match")

naLongDrop <- to_long(E, targets)
naLongDrop <- naLongDrop[!(naLongDrop$geneName == "GENE007" & naLongDrop$libId == "lib004"), ]
stopifnot(same_coefs(
  expect_condition(
    runLmerLoop(naLongDrop, lmForm, "geneName", "geneExp", use_lmer = FALSE, sampleCol = "libId"),
    "1 gene\\(s\\) have missing expression",
    "message"
  ),
  quiet(runLmerLoop_old(naLongDrop, lmForm, "geneName", "geneExp", use_lmer = FALSE))
))
pass("long df with a missing gene/sample row: falls back for that gene")

rankTargets <- targets
rankTargets$x2 <- 2 * rankTargets$x
rankForm <- "~ x + x2 + grp"
rankNew <- expect_condition(
  runLmerLoop(
    voomObj = list(E = E, targets = rankTargets), lmerformula = rankForm,
    use_lmer = FALSE
  ),
  "rank deficient",
  "message"
)
stopifnot(same_coefs(
  rankNew,
  quiet(runLmerLoop_old(to_long(E, rankTargets), rankForm, "geneName", "geneExp", use_lmer = FALSE))
))
pass("rank-deficient design falls back to lm, results match")

badMeta <- longDf
badMeta$x[1] <- badMeta$x[1] + 1
stopifnot(inherits(
  tryCatch(
    quiet(runLmerLoop(badMeta, lmForm, "geneName", "geneExp", use_lmer = FALSE, sampleCol = "libId")),
    error = function(e) e
  ),
  "error"
))
pass("sampleCol with metadata that varies within a sample is an error")

# ── 5. Parallel lmer ─────────────────────────────────────────────────────────
section("5. Parallel lmer")

parLong <- run_new(lmerformula = lmerForm, n_cores = 2)
stopifnot(identical(parLong, oldLmer))
pass("long df, n_cores = 2: identical to old version")

parVoom <- quiet(runLmerLoop(voomObj = voomObj, lmerformula = lmerForm, n_cores = 2))
stopifnot(same_coefs(parVoom, oldLmer, tol = 1e-10), "df" %in% names(parVoom))
pass("voom input, n_cores = 2: matches old version, including df")

# ── 6. keep_models ───────────────────────────────────────────────────────────
section("6. keep_models")

parKeep <- run_new(lmerformula = lmerForm, keep_models = TRUE, n_cores = 2)
stopifnot(
  identical(parKeep$resultsDf, oldKeep$resultsDf),
  identical(names(parKeep$modelObjList), rownames(E)),
  inherits(parKeep$modelObjList[["GENE001"]], "lmerModLmerTest")
)
pass("lmer, n_cores = 2: named model list returned")

lmKeep <- expect_condition(
  runLmerLoop(voomObj = voomObj, lmerformula = lmForm, use_lmer = FALSE, keep_models = TRUE),
  "does not save model objects",
  "warning"
)
stopifnot(
  identical(names(lmKeep), c("resultsDf", "modelObjList")),
  is.null(lmKeep$modelObjList),
  same_coefs(lmKeep$resultsDf, oldLm)
)
pass("matrix lm + keep_models: warning, modelObjList is NULL")

# ── 7. on_error ──────────────────────────────────────────────────────────────
section("7. on_error")

for (cores in c(1, 2)) {
  stopErr <- tryCatch(
    quiet(runLmerLoop(longBad, lmForm, "geneName", "geneExp", use_lmer = FALSE, n_cores = cores)),
    error = function(e) conditionMessage(e)
  )
  stopifnot(is.character(stopErr), grepl("GENE010", stopErr))

  warnRes <- expect_condition(
    runLmerLoop(
      longBad, lmForm, "geneName", "geneExp", use_lmer = FALSE,
      n_cores = cores, on_error = "warn"
    ),
    "1 model\\(s\\) failed",
    "warning"
  )
  stopifnot(
    !"GENE010" %in% warnRes$gene,
    length(unique(warnRes$gene)) == nGenes - 1
  )
}
pass("stop (default) stops; warn skips the failed gene, serial and parallel")

# ── 8. Timing (informational) ────────────────────────────────────────────────
section("8. Timing (informational, not asserted)")

timeIt <- function(expr) unname(system.time(quiet(expr))["elapsed"])

bigGenes <- 2000
bigE <- matrix(rnorm(bigGenes * nrow(targets)), bigGenes, nrow(targets))
rownames(bigE) <- sprintf("G%04d", seq_len(bigGenes))
colnames(bigE) <- targets$libId
bigLong <- to_long(bigE, targets)

cat(sprintf(
  "  lm, %d genes:   old %.2fs | new loop %.2fs | matrix (long) %.2fs | matrix (voom) %.2fs\n",
  bigGenes,
  timeIt(runLmerLoop_old(bigLong, lmForm, "geneName", "geneExp", use_lmer = FALSE)),
  timeIt(runLmerLoop(bigLong, lmForm, "geneName", "geneExp", use_lmer = FALSE)),
  timeIt(runLmerLoop(bigLong, lmForm, "geneName", "geneExp", use_lmer = FALSE, sampleCol = "libId")),
  timeIt(runLmerLoop(voomObj = list(E = bigE, targets = targets), lmerformula = lmForm, use_lmer = FALSE))
))

lmerGenes <- 500
lmerLong <- to_long(bigE[seq_len(lmerGenes), ], targets)
nc <- max(1L, parallel::detectCores() - 1L)
cat(sprintf(
  "  lmer, %d genes: old %.2fs | new 1 core %.2fs | new %d cores %.2fs\n",
  lmerGenes,
  timeIt(runLmerLoop_old(lmerLong, lmerForm, "geneName", "geneExp")),
  timeIt(runLmerLoop(lmerLong, lmerForm, "geneName", "geneExp")),
  nc,
  timeIt(runLmerLoop(lmerLong, lmerForm, "geneName", "geneExp", n_cores = nc))
))

cat("\nAll assertions passed.\n")

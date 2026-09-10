# Verification script for make_color_palette.R
# Run from the repo root:  Rscript Testing/test_make_color_palette.R
#
# Plain stopifnot() assertions — the repo has no testthat harness.

source("make_color_palette.R")

pass <- function(label) cat("  PASS  ", label, "\n")
section <- function(label) cat("\n== ", label, " ==\n", sep = "")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Min CIEDE2000 distance between two sets of colors (0 = identical color present)
min_delta_e <- function(a, b) {
  min(farver::compare_colour(
    farver::decode_colour(a, to = "lab"),
    farver::decode_colour(b, to = "lab"),
    from_space = "lab",
    method = "cie2000"
  ))
}

# Compare ignoring the redundant "FF" alpha suffix paletteer appends
same_colors <- function(a, b) {
  identical(
    unname(.normalize_colors(a)),
    unname(.normalize_colors(b))
  )
}

# ── 1. Auto-selection ────────────────────────────────────────────────────────
section("1. Auto-selection")

# CIEDE2000 scoring must actually produce a score. farver's argument is
# `method`, not `metric`; when this was wrong every score fell through the
# candidate loop's tryCatch to -1 and which.max silently returned the FIRST
# candidate row, so no palette was ever really scored. Guard against that
# regressing.
stopifnot(.score_palette_distinctiveness(c("#000000", "#FFFFFF")) > 50)
stopifnot(.score_palette_distinctiveness(c("#000000", "#010101")) < 5)
stopifnot(.score_palette_distinctiveness(character(0)) == 0)
pass("distinctiveness scoring returns real CIEDE2000 distances")

# Scoring being live means an identical-color palette must lose to a distinct
# one, which is only true if scores differ across candidates.
scores <- vapply(
  c("RColorBrewer::Dark2", "RColorBrewer::Paired", "ggthemes::colorblind"),
  function(p) {
    .score_palette_distinctiveness(
      as.character(paletteer::paletteer_d(p))[1:5]
    )
  },
  numeric(1)
)
stopifnot(length(unique(round(scores, 6))) == 3, all(scores > 0))
pass("candidate palettes receive distinct, positive scores")

# Locked-in outputs for the common cases, so future changes to scoring or to
# the candidate pool show up as a deliberate diff rather than a silent shift.
expected <- list(
  list(
    vars = c("Control", "24hr", "48hr", "96hr", "10Day"),
    args = list(),
    want = c("#4D4D4DFF", "#5DA5DAFF", "#FAA43AFF", "#60BD68FF", "#F17CB0FF")
  ),
  list(vars = c("M", "F"), args = list(), want = c("#000000FF", "#E69F00FF")),
  list(
    vars = LETTERS[1:3],
    args = list(),
    want = c("#000000FF", "#E69F00FF", "#56B4E9FF")
  ),
  list(
    vars = LETTERS[1:8],
    args = list(),
    want = c(
      "#000000FF",
      "#E69F00FF",
      "#56B4E9FF",
      "#009E73FF",
      "#F0E442FF",
      "#0072B2FF",
      "#D55E00FF",
      "#CC79A7FF"
    )
  ),
  list(
    vars = LETTERS[1:10],
    args = list(),
    want = c(
      "#A6CEE3FF",
      "#1F78B4FF",
      "#B2DF8AFF",
      "#33A02CFF",
      "#FB9A99FF",
      "#E31A1CFF",
      "#FDBF6FFF",
      "#FF7F00FF",
      "#CAB2D6FF",
      "#6A3D9AFF"
    )
  ),
  list(
    vars = LETTERS[1:15],
    args = list(),
    want = c(
      "#222222",
      "#F3C300",
      "#875692",
      "#F38400",
      "#A1CAF1",
      "#BE0032",
      "#C2B280",
      "#848482",
      "#008856",
      "#E68FAC",
      "#0067A5",
      "#F99379",
      "#604E97",
      "#F6A600",
      "#B3446C"
    )
  ),
  list(
    vars = LETTERS[1:5],
    args = list(colorblind_safe = FALSE),
    want = c("#E96C0CFF", "#CDEE40FF", "#8BD4F7FF", "#1927DDFF", "#0D0006FF")
  )
)

for (case in expected) {
  got <- quiet(do.call(
    make_color_palette,
    c(list(case$vars, verbose = FALSE), case$args)
  ))
  stopifnot(identical(unname(got), case$want))
  stopifnot(identical(names(got), unique(as.character(case$vars))))
}
pass(sprintf("%d auto-selection cases match the recorded output", length(expected)))

# ── 2. Cross-variable avoidance ──────────────────────────────────────────────
section("2. Cross-variable avoidance via avoid_keys")

tp <- quiet(make_color_palette(
  c("Control", "24hr", "48hr", "96hr"),
  registry_key = "timepoint",
  verbose = FALSE
))
sex <- quiet(make_color_palette(
  c("M", "F"),
  registry_key = "sex",
  avoid_keys = "timepoint",
  verbose = FALSE
))
trt <- quiet(make_color_palette(
  c("Drug", "Placebo", "Untreated"),
  registry_key = "treatment",
  avoid_keys = "all",
  verbose = FALSE
))

all_three <- .normalize_colors(c(tp, sex, trt))
stopifnot(length(unique(all_three)) == length(all_three))
pass("timepoint / sex / treatment share no colors")

cat(sprintf(
  "        min deltaE  timepoint-sex: %.1f  timepoint-treatment: %.1f  sex-treatment: %.1f\n",
  min_delta_e(tp, sex),
  min_delta_e(tp, trt),
  min_delta_e(sex, trt)
))
stopifnot(min_delta_e(tp, sex) > 5, min_delta_e(tp, trt) > 5, min_delta_e(sex, trt) > 5)
pass("all cross-variable pairs are > 5 deltaE apart")

# Without avoid_keys the collision is real — confirms the feature does work
sex_naive <- quiet(make_color_palette(c("M", "F"), verbose = FALSE))
stopifnot(min_delta_e(tp, sex_naive) < 1)
pass("without avoid_keys the same colors are reused (the problem being fixed)")

# avoid_colors, explicit
no_black <- quiet(make_color_palette(
  LETTERS[1:4],
  avoid_colors = c("black", "#E69F00"),
  verbose = FALSE
))
stopifnot(!any(.normalize_colors(no_black) %in% c("#000000", "#E69F00")))
pass("avoid_colors excludes named and hex colors alike")

# A typo'd key warns rather than silently doing nothing
# NB: suppressWarnings() must not wrap the call here — its inner handler would
# muffle the warning before the outer one below ever sees it.
warned <- FALSE
invisible(withCallingHandlers(
  suppressMessages(make_color_palette(
    LETTERS[1:3],
    avoid_keys = "no_such_key",
    verbose = FALSE
  )),
  warning = function(w) {
    if (grepl("not found in the registry", conditionMessage(w))) {
      warned <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
))
stopifnot(warned)
pass("unknown avoid_keys warns")

# ── 3. append_colors ─────────────────────────────────────────────────────────
section("3. append_colors")

# RColorBrewer::Dark2 has 8 colors; ask for 10 levels with 2 appended.
dark2 <- as.character(paletteer::paletteer_d("RColorBrewer::Dark2"))
appended <- quiet(make_color_palette(
  LETTERS[1:10],
  palette = "RColorBrewer::Dark2",
  append_colors = c("black", "grey"),
  verbose = FALSE
))
stopifnot(same_colors(appended[1:8], dark2))
stopifnot(same_colors(appended[9:10], c("black", "grey")))
pass("pinned palette + append_colors covers more levels than the palette alone")

# The same call without append_colors is an error (palette too short) —
# and the error points at append_colors as the remedy.
err <- tryCatch(
  quiet(make_color_palette(
    LETTERS[1:10],
    palette = "RColorBrewer::Dark2",
    verbose = FALSE
  )),
  error = conditionMessage
)
stopifnot(is.character(err), grepl("append_colors", err))
pass("too-short palette errors and suggests append_colors")

# Auto-selection + append: appended colors land on the trailing levels and the
# scored palette avoids them.
auto_app <- quiet(make_color_palette(
  LETTERS[1:5],
  append_colors = "black",
  verbose = FALSE
))
stopifnot(same_colors(auto_app[5], "black"))
stopifnot(!"#000000" %in% .normalize_colors(auto_app[1:4]))
pass("auto-selection reserves the appended color and scores around it")

# ── 4. Custom palette vectors ────────────────────────────────────────────────
section("4. Custom palette vectors")

mine <- quiet(make_color_palette(
  c("Drug", "Placebo", "Untreated"),
  palette = c("#E1BD6D", "black", "grey40"),
  verbose = FALSE
))
stopifnot(identical(
  unname(mine),
  c("#E1BD6D", "#000000", "#666666")
))
stopifnot(identical(names(mine), c("Drug", "Placebo", "Untreated")))
pass("raw color vector is used verbatim and normalized to hex")

bad <- tryCatch(
  quiet(make_color_palette(
    LETTERS[1:2],
    palette = c("#E1BD6D", "notacolor"),
    verbose = FALSE
  )),
  error = conditionMessage
)
stopifnot(is.character(bad), grepl("notacolor", bad))
pass("an invalid color errors with the offending value named")

# A custom palette plus a registry_key persists the assignment
invisible(quiet(make_color_palette(
  c("Drug", "Placebo"),
  palette = c("#E1BD6D", "black"),
  registry_key = "custom demo",
  verbose = FALSE
)))
stopifnot(exists("custom_demo", envir = .color_registry))
pass("custom palette + registry_key stores the assignment (key sanitized)")

# The sanitized key is now reachable from forget_color_registry()
quiet(forget_color_registry("custom demo"))
stopifnot(!exists("custom_demo", envir = .color_registry))
pass("forget_color_registry() sanitizes the key like make_color_palette() does")

# ── 5. registry_overwrite ────────────────────────────────────────────────────
section("5. registry_overwrite")

quiet(forget_color_registry("ow"))
first <- quiet(make_color_palette(
  c("A", "B", "C"),
  palette = c("red", "green", "blue"),
  registry_key = "ow",
  verbose = FALSE
))
redone <- quiet(make_color_palette(
  c("A", "B"),
  palette = c("cyan", "magenta"),
  registry_key = "ow",
  registry_overwrite = TRUE,
  verbose = FALSE
))
cache <- get("ow", envir = .color_registry)

stopifnot(same_colors(redone, c("cyan", "magenta")))
stopifnot(same_colors(cache[c("A", "B")], c("cyan", "magenta")))
stopifnot(same_colors(cache["C"], "blue"))
stopifnot(identical(sort(names(cache)), c("A", "B", "C")))
pass("only the levels passed in are reassigned; C survives untouched")

# Without the flag, the cache still wins
kept <- quiet(make_color_palette(
  c("A", "B"),
  palette = c("yellow", "orange"),
  registry_key = "ow",
  verbose = FALSE
))
stopifnot(same_colors(kept, c("cyan", "magenta")))
pass("registry_overwrite = FALSE leaves cached colors in place")

# Flag without a key warns
warned <- FALSE
invisible(withCallingHandlers(
  suppressMessages(make_color_palette(
    LETTERS[1:2],
    registry_overwrite = TRUE,
    verbose = FALSE
  )),
  warning = function(w) {
    if (grepl("no effect without a registry_key", conditionMessage(w))) {
      warned <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
))
stopifnot(warned)
pass("registry_overwrite without registry_key warns")

# ── 6. Idempotency and incremental growth ────────────────────────────────────
section("6. Idempotency and incremental growth")

quiet(forget_color_registry("grow"))
v1 <- quiet(make_color_palette(
  c("Control", "24hr", "48hr"),
  registry_key = "grow",
  verbose = FALSE
))
v1_again <- quiet(make_color_palette(
  c("Control", "24hr", "48hr"),
  registry_key = "grow",
  verbose = FALSE
))
stopifnot(identical(v1, v1_again))
pass("repeated calls under one key are identical")

v2 <- quiet(make_color_palette(
  c("Control", "24hr", "48hr", "96hr"),
  registry_key = "grow",
  verbose = FALSE
))
stopifnot(identical(v2[names(v1)], v1))
stopifnot(!.normalize_colors(v2[["96hr"]]) %in% .normalize_colors(v1))
cat(sprintf(
  "        new level is %.1f deltaE from its siblings\n",
  min_delta_e(v2[["96hr"]], v1)
))
stopifnot(min_delta_e(v2[["96hr"]], v1) > 5)
pass("a level added later keeps the old colors and gets a distinct new one")

# Factor level order is respected and does not disturb assignment
f <- factor(c("48hr", "Control", "24hr"), levels = c("Control", "24hr", "48hr"))
vf <- quiet(make_color_palette(f, registry_key = "grow", verbose = FALSE))
stopifnot(identical(names(vf), c("Control", "24hr", "48hr")))
stopifnot(identical(vf, v1[names(vf)]))
pass("factor input reuses the same assignment in level order")

# ── 7. Overflow behavior ─────────────────────────────────────────────────────
section("7. Overflow behavior")

k21 <- quiet(make_color_palette(LETTERS[1:21], verbose = FALSE))
stopifnot(length(k21) == 21, !anyNA(k21))
stopifnot(length(unique(.normalize_colors(k21))) == 21)
pass("21 levels fills from pals::kelly with no NA and no duplicates")

cap <- tryCatch(
  quiet(make_color_palette(LETTERS[1:22], verbose = FALSE)),
  error = conditionMessage
)
stopifnot(is.character(cap), grepl("Maximum supported is 21", cap))
pass("22 levels errors with the corrected cap")

# Avoidance that cannot be satisfied warns instead of silently colliding
quiet(forget_color_registry("hog"))
invisible(quiet(make_color_palette(
  LETTERS[1:8],
  palette = "RColorBrewer::Dark2",
  registry_key = "hog",
  verbose = FALSE
)))
warned <- FALSE
squeezed <- withCallingHandlers(
  suppressMessages(make_color_palette(
    c("X", "Y"),
    palette = "RColorBrewer::Dark2",
    avoid_keys = "hog",
    verbose = FALSE
  )),
  warning = function(w) {
    if (grepl("unused by the avoided variables", conditionMessage(w))) {
      warned <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
)
stopifnot(warned, length(squeezed) == 2, !anyNA(squeezed))
pass("an exhausted palette warns and still returns usable colors")

# ── 8. Save / load round trip ────────────────────────────────────────────────
section("8. Registry save / load round trip")

tmp <- tempfile(fileext = ".rds")
before <- as.list(.color_registry)
quiet(save_color_registry(tmp))
rm(list = ls(.color_registry), envir = .color_registry)
stopifnot(length(ls(.color_registry)) == 0)
quiet(load_color_registry(tmp))
stopifnot(identical(before[order(names(before))], as.list(.color_registry)[order(names(as.list(.color_registry)))]))
unlink(tmp)
pass("registry survives a save / load round trip")

cat("\nAll checks passed.\n")

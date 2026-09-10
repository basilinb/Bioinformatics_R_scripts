# make_color_palette (and registry helpers)

**File:** `make_color_palette.R`

## Description

Generates stable, perceptually distinct named color palettes for categorical variables. Automatically selects the best-scoring qualitative palette from `paletteer` (optionally restricted to colorblind-safe options) using CIEDE2000 perceptual distance scoring. Supports a persistent cross-session color registry so the same variable always gets the same color across plots and analyses.

When a document colors several variables at once, `avoid_keys` steers each variable away from the colors already assigned to the others, so no two legends share a color. Palettes can be pinned, extended with colors of your own via `append_colors`, or supplied outright as a character vector.

Also provides ggplot2 drop-in scale wrappers (`scale_color_auto`, `scale_fill_auto`) and registry management helpers.

---

## make_color_palette

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `variables` | character or factor vector | — | The variable values to assign colors to. Unique values are extracted in first-appearance (or factor level) order. Max 21 unique values. |
| `palette` | character | `NULL` | Palette to use instead of auto-selection. Either a `"package::name"` string resolved through paletteer (e.g. `"wesanderson::AsteroidCity1"`), **or** a character vector of your own colors / hex codes (e.g. `c("#E1BD6D", "black", "grey40")`). |
| `append_colors` | character | `NULL` | Colors appended after the palette's own colors, for the trailing levels — e.g. `c("black", "grey")` for control / reference groups. Works with a pinned palette (extends it) and with auto-selection (that many fewer colors are drawn from the scored palette, and the appended colors are themselves avoided while scoring). |
| `colorblind_safe` | logical | `TRUE` | Restrict to colorblind-safe palettes. Falls back to full pool if no safe palette covers the required number of colors. |
| `registry_key` | character | `NULL` | Key for the color registry. Reuses previously assigned colors for known variables; assigns new colors only to new variables. |
| `avoid_keys` | character | `NULL` | Registry key(s) whose colors this variable should steer away from, so several variables in one document do not share colors. Use `"all"` for every other key currently in the registry. |
| `avoid_colors` | character | `NULL` | Explicit colors to steer away from, in addition to `avoid_keys`. Accepts names or hex. |
| `registry_overwrite` | logical | `FALSE` | Reassign the levels passed in this call even if they are already registered. Other levels stored under the key are left untouched. |
| `verbose` | logical | `TRUE` | Print the selected palette name and scoring info. |

### Returns

A named character vector of hex color codes, one per unique variable value.

### How avoidance works

Avoidance is folded into the CIEDE2000 scoring, not applied after the fact — so it influences *which palette is chosen*, not just which colors are taken from it. A palette whose colors sit on top of an avoided color is penalized as hard as one whose colors sit on top of each other. On top of that:

- The colors already stored under this call's own `registry_key` are **always** avoided, so incrementally adding a level cannot duplicate one of its siblings.
- Colors are taken in the palette's own order once avoided ones are dropped, preserving the palette's intended look.
- If a palette cannot supply enough unused colors, you get a `warning()` and a distance-maximizing fallback rather than a silent collision.

> **Tip:** `avoid_keys` is opt-in and needs the other variables to have been registered first. Assign your variables in a fixed order at the top of a document, then reuse the registry keys downstream.

---

## Registry helpers

| Function | Description |
|----------|-------------|
| `save_color_registry(path)` | Save all current registry entries to an RDS file (default: `~/.color_registry.rds`) |
| `load_color_registry(path)` | Load a previously saved registry from disk |
| `view_color_registry()` | Print all current registry entries |
| `forget_color_registry(key)` | Remove a key from the registry entirely (forces re-assignment on next call) |

> **Re-coloring:** to change only *some* levels of a key while leaving the rest alone, pass `registry_overwrite = TRUE` to `make_color_palette` instead of forgetting the whole key. Keys are sanitized the same way in both places, so `"time point"` and `time_point` refer to the same entry.

---

## ggplot2 scale wrappers

### `scale_color_auto(variables, palette, append_colors, colorblind_safe, registry_key, avoid_keys, avoid_colors, registry_overwrite, ...)`
Drop-in replacement for `scale_color_manual`. Passes additional arguments to `ggplot2::scale_color_manual`.

### `scale_fill_auto(variables, palette, append_colors, colorblind_safe, registry_key, avoid_keys, avoid_colors, registry_overwrite, ...)`
Drop-in replacement for `scale_fill_manual`.

---

## Dependencies

`paletteer`, `pals`, `dplyr`, `farver`, `ggplot2`, `RColorBrewer`, `ggthemes`, `ggsci`, `ochRe`, `colorBlindness`

Install all at once:
```r
install.packages(c(
  "paletteer", "pals", "dplyr", "farver", "ggplot2",
  "RColorBrewer", "ggthemes", "ggsci", "colorBlindness"
))
# ochRe is GitHub-only:
# remotes::install_github("hollylkirk/ochRe")
```

The palette-source packages are only needed for the colorblind-safe candidates that live in them; missing ones are skipped during scoring.

---

## Example

```r
source("make_color_palette.R")

# Auto-select best colorblind-safe palette
cols <- make_color_palette(c("Control", "24hr", "48hr", "96hr", "10Day"))

# Use in ggplot
ggplot(df, aes(x = Timepoint, y = Expression, color = Timepoint)) +
  geom_point() +
  scale_color_manual(values = cols)

# Stable registry — same variables always get same colors across sessions
cols <- make_color_palette(df$Timepoint, registry_key = "timepoint")
save_color_registry()  # persist to ~/.color_registry.rds

# Next session
load_color_registry()
cols <- make_color_palette(df$Timepoint, registry_key = "timepoint")
# Returns same colors as before

# ── Several variables in one document, none sharing colors ────────────────
tp_cols  <- make_color_palette(df$TimePoint, registry_key = "timepoint")
sex_cols <- make_color_palette(df$Sex,
                               registry_key = "sex",
                               avoid_keys   = "timepoint")
trt_cols <- make_color_palette(df$Treatment,
                               registry_key = "treatment",
                               avoid_keys   = "all")

# Same thing through the wrappers
ggplot(df, aes(x, y, color = Sex)) +
  geom_point() +
  scale_color_auto(df$Sex, registry_key = "sex", avoid_keys = "timepoint")

# ── Extend a palette with colors of your own ──────────────────────────────
# AsteroidCity1 for the treated groups, black + grey for the controls
cols <- make_color_palette(
  df$Group,
  palette       = "wesanderson::AsteroidCity1",
  append_colors = c("black", "grey")
)

# Appended colors work with auto-selection too — the scored palette is
# chosen to stay clear of them
cols <- make_color_palette(df$Group, append_colors = "black")

# ── Bring your own palette ────────────────────────────────────────────────
cols <- make_color_palette(
  df$Treatment,
  palette      = c("#E1BD6D", "#EABE94", "black", "grey40"),
  registry_key = "treatment"   # the assignment is registered like any other
)

# ── Re-color just some levels, keeping the rest of the key intact ─────────
make_color_palette(
  df$Treatment,
  palette            = "RColorBrewer::Dark2",
  registry_key       = "treatment",
  registry_overwrite = TRUE
)

# Pin a specific palette
cols_pinned <- make_color_palette(
  df$Treatment,
  palette = "RColorBrewer::Dark2"
)
```

---

## Notes

- **Fixed in this version:** the CIEDE2000 scoring passed `metric =` to `farver::compare_colour()`, whose argument is `method =`. Every candidate therefore errored into the loop's `tryCatch` fallback score of `-1`, and `which.max()` silently returned the *first* candidate row — no palette was ever actually scored. Auto-selected palettes differ from previous runs as a result. Colors already stored in a saved registry are untouched.
- The maximum number of levels is **21**, not 22: `pals::kelly(22)` minus its near-white first entry.
- `Testing/test_make_color_palette.R` covers auto-selection, avoidance, `append_colors`, custom vectors, `registry_overwrite`, idempotency, overflow, and the save/load round trip. Run it with `Rscript Testing/test_make_color_palette.R` from the repo root.

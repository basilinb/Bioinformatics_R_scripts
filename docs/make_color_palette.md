# make_color_palette (and registry helpers)

**File:** `make_color_palette.R`

## Description

Generates stable, perceptually distinct named color palettes for categorical variables. Automatically selects the best-scoring qualitative palette from `paletteer` (optionally restricted to colorblind-safe options) using CIEDE2000 perceptual distance scoring. Supports a persistent cross-session color registry so the same variable always gets the same color across plots and analyses.

Also provides ggplot2 drop-in scale wrappers (`scale_color_auto`, `scale_fill_auto`) and registry management helpers.

---

## make_color_palette

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `variables` | character or factor vector | — | The variable values to assign colors to. Unique values are extracted in first-appearance (or factor level) order. Max 22 unique values. |
| `palette` | character | `NULL` | Pin a specific palette in `"package::name"` format (e.g. `"wesanderson::AsteroidCity1"`). Bypasses automatic scoring. |
| `colorblind_safe` | logical | `TRUE` | Restrict to colorblind-safe palettes. Falls back to full pool if no safe palette covers the required number of colors. |
| `registry_key` | character | `NULL` | Key for the color registry. Reuses previously assigned colors for known variables; assigns new colors only to new variables. |
| `verbose` | logical | `TRUE` | Print the selected palette name and scoring info. |

### Returns

A named character vector of hex color codes, one per unique variable value.

---

## Registry helpers

| Function | Description |
|----------|-------------|
| `save_color_registry(path)` | Save all current registry entries to an RDS file (default: `~/.color_registry.rds`) |
| `load_color_registry(path)` | Load a previously saved registry from disk |
| `view_color_registry()` | Print all current registry entries |
| `forget_color_registry(key)` | Remove a specific key from the registry (forces re-assignment on next call) |

---

## ggplot2 scale wrappers

### `scale_color_auto(variables, palette, colorblind_safe, registry_key, ...)`
Drop-in replacement for `scale_color_manual`. Passes additional arguments to `ggplot2::scale_color_manual`.

### `scale_fill_auto(variables, palette, colorblind_safe, registry_key, ...)`
Drop-in replacement for `scale_fill_manual`.

---

## Dependencies

`paletteer`, `pals`, `dplyr`, `colorspace`, `farver`, `ggplot2`, `RColorBrewer`, `ggthemes`, `ggsci`, `ochRe`, `colorBlindness`

Install all at once:
```r
install.packages(c("paletteer", "pals", "dplyr", "colorspace", "farver"))
```

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

# Drop-in ggplot2 wrappers
ggplot(df, aes(x, y, color = Treatment)) +
  geom_point() +
  scale_color_auto(df$Treatment, registry_key = "treatment")

ggplot(df, aes(x, y, fill = Treatment)) +
  geom_col() +
  scale_fill_auto(df$Treatment, registry_key = "treatment")

# Pin a specific palette
cols_pinned <- make_color_palette(
  df$Treatment,
  palette = "RColorBrewer::Dark2"
)
```

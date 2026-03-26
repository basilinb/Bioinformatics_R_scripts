# calculate_module_coherence

**File:** `calculate_module_coherence.R`

## Description

Evaluates the internal coherence (co-expression) of gene modules by computing all pairwise Pearson correlations between genes within each module and testing their significance. Returns correlation data and three `ggplot2` boxplot visualizations. Useful for validating module definitions in a new dataset.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mods` | data.frame | — | Module definition table. Must contain a gene name column and a module assignment column. |
| `mods_title` | character | `"STUDY1"` | Label for the module source used in plot titles. |
| `mods_var` | character | `"module.char"` | Column name in `mods` containing module assignments. |
| `remove_mods` | character/numeric vector | `c(0,"0","00","grey")` | Module labels to exclude (e.g. noise/unassigned modules). |
| `dat` | EList or matrix | — | Expression data. If an `EList` (voom output), `$E` is used automatically. |
| `gene_var` | character | `"geneName"` | Column name in `mods` containing gene names that match `rownames(dat)`. |
| `dat_title` | character | `"STUDY2"` | Label for the expression dataset used in plot titles. |
| `return_plot` | logical | `TRUE` | If `TRUE`, prints the combined faceted plot. |
| `r_cutoff` | numeric | `0.3` | Correlation cutoff shown as a dashed reference line. |
| `p_cutoff` | numeric | `0.01` | P-value cutoff shown as a dashed reference line. |

## Returns

A named list:

| Element | Description |
|---------|-------------|
| `coherence_boxplot_combined` | Faceted boxplot (Correlation + −log10 P side-by-side) |
| `coherence_boxplot_cor` | Boxplot of correlation strength only |
| `coherence_boxplot_p` | Boxplot of −log10(p-value) only |
| `subgene_correlation_df` | data.frame with all gene pairs, correlations, and p-values |
| `GSabsent` | (only if missing genes exist) list of genes in module definitions but absent from expression data |

## Dependencies

`dplyr`, `ggplot2`, `tidyr`, `gtools`, `stats`, `utils`

## Example

```r
source("calculate_module_coherence.R")

# mods: data frame with columns "geneName" and "module.char"
# v: limma voom EList

coh <- calculate_module_coherence(
  mods       = module_df,
  mods_title = "WGCNA_Modules",
  mods_var   = "module.char",
  dat        = v,
  gene_var   = "geneName",
  dat_title  = "ValidationCohort",
  r_cutoff   = 0.3,
  p_cutoff   = 0.01
)

# Access plots
print(coh$coherence_boxplot_combined)
print(coh$coherence_boxplot_cor)

# Access correlation data
head(coh$subgene_correlation_df)
```

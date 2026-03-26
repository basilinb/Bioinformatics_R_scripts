# filter_rare

**File:** `filter_rare.R`

## Description

Filters out low-abundance (rare) genes from an edgeR `DGEList` object. A gene is retained only if its CPM exceeds `min_CPM` in at least `min_sample` samples (or `min_pct`% of samples). Optionally plots mean-variance trends before and after filtering using `limma::voom`.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dat` | DGEList | — | edgeR `DGEList` object (output of `edgeR::DGEList`). |
| `min_CPM` | numeric | — | Minimum counts per million threshold. |
| `gene_var` | character | `"ensembl_gene_id"` | Column name in `dat$genes` that matches the row names of the count matrix. |
| `min_sample` | integer | `NULL` | Minimum number of samples that must meet `min_CPM`. Provide either this **or** `min_pct`. |
| `min_pct` | numeric | `NULL` | Minimum **percentage** of samples (0–100) that must meet `min_CPM`. Provide either this **or** `min_sample`. |
| `plot` | logical | `FALSE` | If `TRUE`, prints mean-variance trend plots (unfiltered vs filtered) via `limma::voom`. |

> **Deprecated parameters:** `min.CPM`, `gene.var`, `min.sample`, `min.pct` (snake_case versions are preferred).

## Returns

A filtered `DGEList` object with rare genes removed from both `$counts` and `$genes`. Prints a message reporting how many genes were removed and the percentage.

## Dependencies

`edgeR`, `limma` (for plot), `ggplot2` (for plot)

## Example

```r
source("filter_rare.R")

# Filter by absolute sample count
dat_filtered <- filter_rare(
  dat        = my_dge,
  min_CPM    = 0.1,
  min_sample = 3,
  gene_var   = "geneName"
)

# Filter by percentage of samples
dat_filtered <- filter_rare(
  dat      = my_dge,
  min_CPM  = 0.1,
  min_pct  = 10,       # at least 10% of samples
  plot     = TRUE,
  gene_var = "geneName"
)
```

> **Tip:** A common starting point for bulk RNA-seq is `min_CPM = 0.1` in `min_pct = 10` (at least 10% of samples). Adjust based on your smallest group size.

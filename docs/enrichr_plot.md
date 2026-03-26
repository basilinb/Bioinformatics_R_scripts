# enrichr_plot

**File:** `enrichr_plot.R`
**Author:** Basilin Benson (UW)
**License:** GPL v3

## Description

Creates lollipop plots from Enrichr results (output of `enrichr_run`). For each specified database, plots the top 20 enriched terms ordered by adjusted p-value, with odds ratios on the x-axis, overlap counts as labels, and points colored by FDR significance tier.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enrich_file` | named list | `NULL` | Output from `enrichr_run`. Each element is a data frame from one Enrichr database. |
| `type` | character | `"gene"` | Passed through; currently not used in plot logic but kept for consistency with `enrichr_run`. |
| `dbs` | character vector | `c("MSigDB_Hallmark_2020")` | Which databases in `enrich_file` to plot. Must match keys in the list. |
| `Base_name` | character | `"Gene"` | Label used in the plot title (e.g. gene list name or module ID). |

## Returns

A named list of `ggplot2` objects, one per database. Keys follow the pattern `"plt_<Base_name>_<db>"`.

Color scale:

| Color | FDR |
|-------|-----|
| Dark red | ≤ 0.001 |
| Red | ≤ 0.01 |
| Orange | ≤ 0.05 |
| Light orange | ≤ 0.1 |
| Yellow | ≤ 0.2 |
| Grey | NS |

## Dependencies

`tidyverse`, `ggplot2`, `ggrepel`, `enrichR`

## Example

```r
source("enrichr_fun.R")
source("enrichr_plot.R")

enrich_res <- enrichr_run(
  gene_mod_list = my_genes,
  dbs           = c("MSigDB_Hallmark_2020")
)

plt_list <- enrichr_plot(
  enrich_file = enrich_res,
  dbs         = c("MSigDB_Hallmark_2020"),
  Base_name   = "MyGeneList"
)

print(plt_list[["plt_MyGeneList_MSigDB_Hallmark_2020"]])
```

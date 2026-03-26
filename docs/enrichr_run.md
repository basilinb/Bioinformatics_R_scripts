# enrichr_run

**File:** `enrichr_fun.R`
**Author:** Basilin Benson (UW)
**License:** GPL v3

## Description

Convenience wrapper for running [Enrichr](https://maayanlab.cloud/Enrichr/) on a gene list or on module gene lists from a module object. Supports HGNC symbols directly or Ensembl IDs (which are converted to HGNC via `biomaRt`). Returns a list of Enrichr result data frames.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_mod_list` | character vector or module object | `NULL` | Gene list (character vector of IDs) or a module object containing `$mods` with `module` and `hgnc_symbol` columns. |
| `type` | character | `"gene"` | `"gene"` to run on a flat gene list; any other value (e.g. `"mod"`) to run on modules. |
| `gene_id` | character | `"HGNC"` | Gene identifier: `"HGNC"` for symbols, anything else triggers Ensembl → HGNC conversion via `biomaRt`. |
| `dbs` | character vector | `c("MSigDB_Hallmark_2020")` | Enrichr database(s) to query. See [Enrichr databases](https://maayanlab.cloud/Enrichr/#stats) for available options. |
| `mod_of_interest` | character | `"all"` | When `type != "gene"`, which modules to run. `"all"` runs all modules; otherwise provide a character vector of module names. |

## Returns

A named list of Enrichr result data frames. When running on a gene list, the list is keyed by database name. When running on modules, it is keyed by `"module_<name>"`.

## Dependencies

`tidyverse`, `enrichR`, `biomaRt`

## Example

```r
source("enrichr_fun.R")
source("enrichr_plot.R")

# Run on a gene list (HGNC symbols)
enrich_res <- enrichr_run(
  gene_mod_list = c("STAT1", "IRF7", "MX1", "IFIT1"),
  type          = "gene",
  gene_id       = "HGNC",
  dbs           = c("MSigDB_Hallmark_2020", "GO_Biological_Process_2023")
)

# Plot results
plt_list <- enrichr_plot(
  enrich_file = enrich_res,
  dbs         = c("MSigDB_Hallmark_2020"),
  Base_name   = "IFN genes"
)
print(plt_list[[1]])

# Run on modules
enrich_mod_res <- enrichr_run(
  gene_mod_list  = module_object,
  type           = "mod",
  dbs            = c("MSigDB_Hallmark_2020"),
  mod_of_interest = c("module_01", "module_03")
)
```

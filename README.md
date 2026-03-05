# Bioinformatics R scripts

A small collection of reusable R functions/scripts I use for common bioinformatics / transcriptomics workflows (RNA-seq), including:

- gene set enrichment (GSEA / fgsea, MSigDB via **msigdbr**)
- over-representation/enrichment workflows (e.g., **clusterProfiler** enricher, **Enrichr**)
- longitudinal / repeated-measures modeling (GAM loops via **mgcv**, LMM loops via **lmerTest**)
- plotting helpers (volcano plots, GSEA plots)
- module / gene-set utilities (module coherence, module overlap / comparison)

> Note: This repository is currently a scripts collection (not an R package). Functions are sourced directly from `.R` files.

---

## Contents (high level)

### Modeling / statistics
- **`GAMLoop.R`** — Run a generalized additive model (GAM) per gene or per module and return tables + model objects.
- **`GAM_loop.R`** — Older GAM loop version (header indicates it should not be used).
- **`runLmerLoop.R`** — Loop over genes/modules and fit an `lmer()` (or `lm()`) model per feature; returns results + model objects.
- **`predictGamValues.R`** — Predict fitted values from GAM objects (e.g., excluding random-effect terms) and return wide matrix output.

### Enrichment / gene sets
- **`BIGprofiler.R`** — Over-representation analysis wrapper around `clusterProfiler::enricher()` using MSigDB collections via `msigdbr` (or a custom TERM2GENE-like data frame).
- **`GSEA_run.R`** — Run GSEA using `fgsea` with MSigDB gene sets via `msigdbr` (or a custom GMT list).
- **`GSEA_plot.R`** — Plot GSEA results (NES by pathway) with significance categories.
- **`GSEAgenesetPlot.R`** — Plot enrichment curve for a single gene set and ranked statistics (optionally highlight a gene).
- **`enrichr_fun.R`** — Convenience wrapper for running Enrichr on genes or on module gene lists.
- **`enrichr_plot.R`** — Plot Enrichr results (odds ratios, overlap labels, FDR bins).
- **`vizGo.R`** — GO term reduction/visualization helper using **rrvgo** (requires GO IDs from MSigDB C5 mapping).

### Expression / filtering utilities
- **`filter_rare.R`** — Filter low-abundance genes in an edgeR `DGEList` by CPM threshold and sample count/percent.
- **`genesetExp.R`** — Summarize expression over a gene set (mean/median per library) using MSigDB or a custom gene set.

### Module utilities
- **`calculate_module_coherence.R`** — Compute within-module gene–gene correlation coherence (and plots).
- **`module_compareModules.r`** — Compare module sets across studies using Sorensen similarity + hypergeometric tests; includes plots.
- **`overlapGroups.R`** — Helper to extract overlap groups from an UpSet-style list input.

### Plotting
- **`Volcano_fxn.R`** — Volcano plot generator from differential expression-like results (labels top hits, optional genes of interest).

---

## Installation / requirements

Most scripts use CRAN/Bioconductor packages. Common dependencies across the repo include:

- `tidyverse`
- `mgcv`
- `data.table`
- `fgsea`
- `msigdbr`
- `clusterProfiler`
- `enrichR`
- `edgeR`, `limma`
- `lmerTest`
- `rrvgo`
- (others appear in individual scripts)

You can install packages as needed. Example (not exhaustive):

```r
install.packages(c("tidyverse", "mgcv", "data.table", "ggrepel", "patchwork"))
# Bioconductor examples:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("edgeR", "limma", "fgsea", "clusterProfiler"))
```

---

## Usage

Because this is not packaged, you typically load functions with `source()`.

```r
source("GAMLoop.R")
source("GSEA_run.R")
source("GSEA_plot.R")
```

### Example: run a GAM per gene (RNA-seq voom object)

```r
source("GAMLoop.R")

# voomObject should look like limma::voom output:
# - voomObject$E expression matrix
# - voomObject$targets metadata with libid + covariates
res <- GAMLoop(
  voomObject = v,
  gamFormula = "~ s(Timepoint, bs = 'cs') + Treatment + donorId",
  returnFixedEffects = TRUE,
  returnSmoothEffects = TRUE,
  returnGeneModelSummary = FALSE,
  returnModelObject = FALSE
)

head(res$FixedEffects)
head(res$SmoothEffects)
```

### Example: run fgsea-based GSEA using MSigDB Hallmark

```r
source("GSEA_run.R")
source("GSEA_plot.R")

# gene_list: named list of named numeric vectors (e.g., logFC named by gene symbol)
 gsea_res <- GSEA_run(gene_list = gene_list, group = "H")
plt <- GSEA_plot(gsea_res, plot.fdr = 0.05)
print(plt)
```

---

## Notes / conventions

- Many scripts call `require()` internally and may attempt to install missing CRAN packages.
- Several functions assume certain column names in inputs (e.g., `donorId`, `libid`, `hgnc_symbol`, etc.). Check the function headers and parameter docs in each file.
- `GAM_loop.R` is labeled as an older version; prefer `GAMLoop.R`.

---

## License

Some files include GPL headers (GPL v3 or later) and some code includes attribution (e.g., GSEA scripts credited to Kim Dill-McFarland, 2020). If you plan to redistribute this repository as a single bundled work, consider adding a top-level `LICENSE` file and ensuring the licensing/attribution across scripts is consistent.

---

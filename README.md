# Bioinformatics R scripts

A collection of reusable R functions for common bioinformatics / transcriptomics workflows (bulk RNA-seq), including:

- Gene set enrichment (GSEA / fgsea, MSigDB via **msigdbr**)
- Over-representation analysis (**clusterProfiler**, **Enrichr**)
- Longitudinal / repeated-measures modeling (GAM via **mgcv**, LMM via **lmerTest**)
- Mediation analysis and causal inference
- Plotting helpers (volcano plots, GSEA plots, enrichment comparisons)
- Module / gene-set utilities (coherence, overlap)
- Color palette management

> **Note:** This is a scripts collection, not an R package. Functions are sourced directly from `.R` files.

Individual function documentation lives in [`docs/`](docs/).

---

## Function reference

### Modeling / statistics

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `GAMLoop` | `GAMLoop.R` | Fit a GAM per gene or module; returns fixed effects, smooth effects, model objects | [docs](docs/GAMLoop.md) |
| `runLmerLoop` | `runLmerLoop.R` | Fit `lmer` (or `lm`) per gene/module; returns coefficient table + model objects | [docs](docs/runLmerLoop.md) |
| `predictGamValues` | `predictGamValues.R` | Predict fitted values from GAM objects, excluding specified terms (e.g. random effects) | [docs](docs/predictGamValues.md) |
| `GAM_loop` | `GAM_loop.R` | *(Deprecated — use `GAMLoop` instead)* Older GAM loop version | — |

### Enrichment / gene sets

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `GSEA_run` | `GSEA_run.R` | Run GSEA with `fgsea` using MSigDB or a custom GMT list | [docs](docs/GSEA_run.md) |
| `GSEA_plot` | `GSEA_plot.R` | Lollipop plot of GSEA results (NES by pathway, colored by FDR) | [docs](docs/GSEA_plot.md) |
| `GSEAgenesetPlot` | `GSEAgenesetPlot.R` | Detailed three-panel GSEA figure for a single gene set | [docs](docs/GSEAgenesetPlot.md) |
| `BIGprofiler` | `BIGprofiler.R` | Over-representation analysis via `clusterProfiler::enricher` with MSigDB or custom DB | [docs](docs/BIGprofiler.md) |
| `enrichr_run` | `enrichr_fun.R` | Run Enrichr on a gene list or module gene lists | [docs](docs/enrichr_run.md) |
| `enrichr_plot` | `enrichr_plot.R` | Lollipop plot of Enrichr results (odds ratio, FDR colors) | [docs](docs/enrichr_plot.md) |
| `vizGO` | `vizGo.R` | Reduce redundant GO terms via `rrvgo` semantic similarity | [docs](docs/vizGO.md) |
| `plot_enrichment_comparisons` | `plot_enrichment_comparisons.R` | Bubble plot of enrichment results across multiple pairwise comparisons | [docs](docs/plot_enrichment_comparisons.md) |

### Expression / filtering utilities

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `filter_rare` | `filter_rare.R` | Filter low-abundance genes from a `DGEList` by CPM threshold | [docs](docs/filter_rare.md) |
| `genesetExp` | `genesetExp.R` | Summarize per-library median/mean expression over a gene set | [docs](docs/genesetExp.md) |

### Mediation / causal inference

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `causalFDR2` | `causalFDR2.R` | FDR-correct mediation analysis output across all effect types | [docs](docs/causalFDR2.md) |
| `FDRcorrectedTable` | `causalFDR2.R` | Reformat FDR-corrected mediation output into publication-ready tables | [docs](docs/causalFDR2.md) |

### Module utilities

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `calculate_module_coherence` | `calculate_module_coherence.R` | Compute within-module gene–gene correlation coherence; returns boxplots + data | [docs](docs/calculate_module_coherence.md) |
| `overlapGroups` | `overlapGroups.R` | Extract overlap groups from an UpSet-style list input | [docs](docs/overlapGroups.md) |

### Plotting

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `volcano_funcxn` | `Volcano_fxn.R` | Volcano plot from limma DE results; labels top hits and genes of interest | [docs](docs/volcano_funcxn.md) |
| `make_contrast_scatterplot` | `make_contrast_scatterplot.R` | Scatterplot comparing effect estimates between two DE contrasts | [docs](docs/make_contrast_scatterplot.md) |

### Utilities

| Function | File | Description | Docs |
|----------|------|-------------|------|
| `make_color_palette` | `make_color_palette.R` | Auto-select colorblind-safe palette with cross-session registry | [docs](docs/make_color_palette.md) |
| `scale_color_auto` / `scale_fill_auto` | `make_color_palette.R` | ggplot2 drop-in color/fill scale wrappers | [docs](docs/make_color_palette.md) |

---

## Installation / requirements

Most scripts use CRAN/Bioconductor packages. Common dependencies:

- `tidyverse`, `ggplot2`, `ggrepel`, `patchwork`
- `mgcv`, `lmerTest`
- `data.table`
- `fgsea`, `msigdbr`
- `clusterProfiler`, `enrichR`
- `edgeR`, `limma`
- `rrvgo`, `org.Hs.eg.db`
- `paletteer`, `pals`, `farver` (for `make_color_palette`)

```r
install.packages(c("tidyverse", "mgcv", "lmerTest", "data.table",
                   "ggrepel", "patchwork", "paletteer", "pals", "farver"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "fgsea", "clusterProfiler",
                       "rrvgo", "org.Hs.eg.db"))
```

See [`Reinstalling_all_R_packages.R`](Reinstalling_all_R_packages.R) for a utility to catalog and reinstall all packages by source.

---

## Usage

Because this is not packaged, load functions with `source()`:

```r
source("GAMLoop.R")
source("GSEA_run.R")
source("GSEA_plot.R")
```

### Example: GAM per gene (RNA-seq voom object)

```r
source("GAMLoop.R")

res <- GAMLoop(
  voomObject = v,
  gamFormula = "~ s(Timepoint, bs = 'cs') + Treatment + s(donorId, bs = 're')",
  returnFixedEffects    = TRUE,
  returnSmoothEffects   = TRUE,
  returnGeneModelSummary = FALSE,
  returnModelObject      = FALSE
)

head(res$FixedEffects)
head(res$SmoothEffects)
```

### Example: GSEA with Hallmark gene sets

```r
source("GSEA_run.R")
source("GSEA_plot.R")

# gene_list: named list of named numeric vectors (logFC named by gene symbol)
gsea_res <- GSEA_run(gene_list = gene_list, group = "H")
plt <- GSEA_plot(gsea_res, plot.fdr = 0.05)
print(plt)
```

---

## Notes / conventions

- Many scripts call `require()` internally and may attempt to install missing CRAN packages.
- Several functions assume specific column names (e.g. `donorId`, `libid`, `hgnc_symbol`, `geneName`). Check the [function docs](docs/) for input requirements.
- `GAM_loop.R` is an older version; prefer `GAMLoop.R`.

---

## License

Some files include GPL headers (GPL v3 or later). Some code includes attribution (e.g. GSEA scripts credited to Kim Dill-McFarland, 2020).

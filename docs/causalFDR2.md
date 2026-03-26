# causalFDR2 & FDRcorrectedTable

**File:** `causalFDR2.R`

---

## causalFDR2

### Description

Transforms the summary output of a mediation analysis (from `mediationGeneSet.R`) from long to wide format and applies FDR correction across all mediators within each contrast. Corrects p-values for ACME (average/control/treated), ADE (average/control/treated), Total Effect, and Proportion Mediated terms.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mediation_output` | list | — | Output object from the mediation analysis workflow. Must contain `$contrasts` (list of contrast definitions) and `$summary` (named list of mediation summary data frames). |
| `method_fdr` | character | `"BH"` | FDR correction method passed to `p.adjust`. Common options: `"BH"` (Benjamini-Hochberg), `"bonferroni"`. |

### Returns

A wide-format data frame with one row per mediator × contrast combination. Columns include estimates, 95% CIs, raw p-values, and BH-adjusted p-values for all mediation effect types.

### Example

```r
source("causalFDR2.R")

fdr_results <- causalFDR2(
  mediation_output = mediation.analysis_results,
  method_fdr       = "BH"
)

head(fdr_results)
```

---

## FDRcorrectedTable

### Description

Takes the output of `causalFDR2` and reformats it into publication-friendly long-format tables for a specific mediator of interest. One table is returned per contrast.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `FDR_corrected_output` | data.frame | — | Output from `causalFDR2`. |
| `mediator_of_interest` | character | — | A single mediator name present in the `mediator` column of `FDR_corrected_output`. |

### Returns

A named list of data frames (one per contrast) in long format with columns: `variable`, `measure` (Estimate / CI Lower / CI Upper / p-value / adjusted p-value), and `value`.

### Example

```r
source("causalFDR2.R")

fdr_results <- causalFDR2(mediation_output = med_res)

# Get formatted table for a specific gene/mediator
tables <- FDRcorrectedTable(
  FDR_corrected_output = fdr_results,
  mediator_of_interest = "STAT1"
)

tables[["T-post_C-pre"]]
```

## Dependencies

`tidyverse`, `stringi`, `stringr`, `reshape2`

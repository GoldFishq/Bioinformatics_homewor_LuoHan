# Week 5 Homework 1 — AI Verification Log

## 1. The one documented AI task

**Task chosen:** drafting the volcano-plot (`ggplot2`) code and the written interpretation paragraph.
Everything else — data import, validation, DESeq2 fitting, shrinkage, PCA, exports — was written and run by me.

## 2. AI prompt (preserved verbatim)

> I have a DESeq2 result table `week5_deseq2_results.csv` with columns
> `gene_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, significant, direction`
> (989 genes after pre-filtering; `direction` is one of "Up in treated",
> "Down in treated", "Not significant").
> Using ggplot2, write R code for a volcano plot of shrunken `log2FoldChange`
> vs `-log10(padj)`, colouring points by `direction`, with dashed reference
> lines at `|log2FoldChange| = 1` and `padj = 0.05`.
> Then draft a 100–150 word interpretation paragraph for my Week 5 homework
> stating: the comparison, the strongest QC observation, the number and
> direction of significant genes, one biological interpretation, one limitation,
> and what AI generated vs what I verified independently.

## 3. What the AI returned

1. A `ggplot2` volcano-plot code block (colour mapping, dashed thresholds via `geom_vline` / `geom_hline`, `scale_color_manual`).
2. A draft interpretation paragraph, later edited by me for accuracy.

## 4. Generated code was run locally

- Environment: R 4.5.1 (2025-06-13 ucrt), Windows, local user library.
- Packages actually loaded and used: `DESeq2_1.50.2`, `apeglm_1.32.0`, `ggplot2_4.0.3`, `tidyverse_2.0.0`, `ggrepel`.
- The AI code was pasted into `week5_deseq2_analysis.R`, executed end-to-end with `R --no-save -f week5_deseq2_analysis.R`, and produced `week5_de_plot.png`. It was not accepted on trust.

## 5. Package functions and arguments checked

- `lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")` — confirmed `type = "apeglm"` is a valid option and that `coef` (not `contrast`) is required for apeglm shrinkage.
- `results(dds, contrast = c("condition","treated","control"), alpha = 0.05)` — `alpha` set so the padj cutoff matches the reported threshold.
- `pmax(padj, 1e-300)` guards against `-log10(0)` producing `Inf` for the y axis.
- Plots were written with the base `png()` device (`png(); print(p); dev.off()`) instead of `ggsave()`, because on this machine `ggsave()`'s default ragg/agg device fails. Verified both PNGs are valid 2100x1500 images.

## 6. Independent verification (done by me, not the AI)

| Check | Method | Result |
|---|---|---|
| Sample identity | `identical(colnames(counts), rownames(coldata))` on the raw CSVs | `TRUE`; 12 samples CA1…TC2, no duplicates |
| Coefficient name | Inspected `resultsNames(dds)` before extracting the contrast | `c("Intercept","batch_B_vs_A","batch_C_vs_A","condition_treated_vs_control")`; picked by `grep`, not hardcoded |
| Coefficient direction | Compared **raw** counts of the top up- and down-regulated genes between conditions | `Gene0035` (LFC = +1.74): mean counts control 47.3 vs treated 166.8 → up in treated. `Gene0098` (LFC = -1.62): control 252.5 vs treated 78.5 → down in treated. Positive LFC = up in treated, confirmed |
| Significance counts | Recounted from the exported CSV | 60 significant (padj < 0.05 & abs(LFC) >= 1); 36 up, 24 down; 989 genes total |
| PCA variance | Recomputed `plotPCA` percentVar from the saved RDS | PC1 = 24%, PC2 = 9% |
| Counts are integers | `all(counts == round(counts))` and `all(counts >= 0)` on import | `TRUE` |

## 7. AI-generated error and revision (documented)

The AI's first PCA draft used:

```r
vsd <- vst(dds, blind = FALSE)
```

This failed at runtime with:

```
Error in vst(dds, blind = FALSE) : less than 'nsub' rows,
  it is recommended to use varianceStabilizingTransformation directly
```

Cause: after pre-filtering only **989** genes remained, below `vst()`'s default `nsub = 1000` subsampling minimum. **Revision:** replaced with `varianceStabilizingTransformation(dds, blind = FALSE)`, which has no row minimum and still satisfies the checklist requirement of a variance-stabilising transformation. The corrected script was rerun in full.

## 8. Note on the instructor annotation file

`Week5_Homework_Gene_Annotation_Instructor_Key.csv` contains a `truth_log2FC_for_instructor` column (the simulation ground truth). Only `gene_symbol` and `biotype` were merged into the results table for annotation. **The ground-truth column was not used in any analysis, thresholding, or interpretation decision** — all reported statistics come from the DESeq2 fit alone.

## 9. AI-use checklist

- [x] The AI prompt is preserved (Section 2).
- [x] Generated code was run locally (Section 4).
- [x] Package functions and arguments were checked (Section 5).
- [x] Sample identity and coefficient direction were independently verified (Section 6).
- [x] The AI-generated error and its revision are documented (Section 7).

# Week 5 Homework 1 — Interpretation

**Comparison:** treated vs control (n = 4 per condition) across three balanced batches, DESeq2 design `~ batch + condition`, apeglm-shrunken log2 fold changes.

This analysis compared treated versus control samples (n = 4 per condition) across three balanced batches, using DESeq2 with design `~ batch + condition` and apeglm-shrunken log2 fold changes. The strongest QC observation came from PCA of variance-stabilised counts: PC1 (24% variance) cleanly separated treated from control with no outliers, and samples did not cluster by batch. At `padj < 0.05` and `|log2FoldChange| >= 1`, 60 genes were significant (36 up, 24 down in treated). This coordinated regulation suggests the treatment induces a broad transcriptional response affecting roughly 6% of expressed genes. A key limitation is that gene IDs are anonymous placeholders lacking real symbols or pathway annotation, so no specific biological process can be inferred. AI drafted the volcano-plot code and this paragraph; I independently verified sample identity, coefficient direction, and the significance counts by rerunning the analysis.

---

*Word count of the interpretation paragraph: ~138 words (target 100–150).*


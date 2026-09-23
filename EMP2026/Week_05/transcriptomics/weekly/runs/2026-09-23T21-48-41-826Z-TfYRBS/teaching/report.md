# Week 5 Transcriptomics — Analysis Report

Student: Luo Han (SUAT24000110)

## Experiment
- Name: RNAseq_LIPUS
- Omics: transcriptomics
- Differential analysis: DESeq2, comparison `DMSO` (reference) vs `T4400` (test), pairwise
- Genes tested: 16757

## Outputs
- `data/RNAseq_LIPUS_assay.csv` — count matrix (features x samples)
- `data/RNAseq_LIPUS_metadata.csv` — sample metadata / colData
- `results/RNAseq_LIPUS_diff_analysis.csv` — DESeq2 DEG table (baseMean, log2FC, pvalue, padj)
- `plots/heatmap_RNAseq_LIPUS.pdf` — expression heatmap
- `plots/volcano_RNAseq_LIPUS.pdf` — volcano plot

Note: regenerated locally from the EMP analysis session for manual GitHub upload.

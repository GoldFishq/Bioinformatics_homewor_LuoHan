# ================================================================
# Week 5 Homework 1 - Bulk RNA-seq differential expression with DESeq2
# Student analysis script
# Design: ~ batch + condition   (treated vs control, 3 balanced batches)
# ================================================================

# ---- locate directories (works with Rscript and RStudio) ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
setwd(script_dir)

# Input files live in the sibling "input" folder (NOT overwritten here).
input_dir <- file.path(dirname(script_dir), "input")

count_file    <- file.path(input_dir, "Week5_Homework_Count_Matrix.csv")
metadata_file <- file.path(input_dir, "Week5_Homework_Sample_Metadata.csv")
anno_file     <- file.path(input_dir, "Week5_Homework_Gene_Annotation_Instructor_Key.csv")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
  library(ggrepel)
})

cat("Working directory:", getwd(), "\n")

# ------------------------------------------------------------------
# 1. Import
# ------------------------------------------------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
print(head(counts[, 1:6]))

# ------------------------------------------------------------------
# 2. Mandatory validation
# ------------------------------------------------------------------
stopifnot("ncol(counts) == nrow(coldata)" =
            ncol(counts) == nrow(coldata))
stopifnot("colnames(counts) match rownames(coldata)" =
            identical(colnames(counts), rownames(coldata)))
stopifnot("no duplicated sample IDs" =
            !any(duplicated(rownames(coldata))))
stopifnot("all counts are non-negative" = all(counts >= 0))
stopifnot("all counts are integers" =
            all(as.matrix(counts) == round(as.matrix(counts))))

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch <- factor(coldata$batch)

cat("\nBatch x condition table:\n")
print(table(coldata$batch, coldata$condition))
cat("\nLibrary size (colSums) summary:\n")
print(summary(colSums(counts)))

# ------------------------------------------------------------------
# 3. Construct DESeq2 object
#    batch is included to remove the balanced batch effect so that
#    the treatment effect is not confounded with it.
# ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ batch + condition
)

# ------------------------------------------------------------------
# 4. Pre-filter
#    Keep genes with at least 10 counts in at least 3 samples.
# ------------------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3
cat("\nGenes before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("Genes after filtering :", nrow(dds), "\n")

# ------------------------------------------------------------------
# 5. Fit model
# ------------------------------------------------------------------
dds <- DESeq(dds)
coef_names <- resultsNames(dds)
cat("\nresultsNames(dds):\n")
print(coef_names)

# DO NOT assume the coefficient name - read it from the model.
target_coef <- grep("condition_treated_vs_control", coef_names, value = TRUE)
if (length(target_coef) != 1) {
  stop("Expected exactly one treated-vs-control coefficient. Inspect resultsNames(dds).")
}
cat("\nUsing coefficient:", target_coef, "\n")

# ------------------------------------------------------------------
# 6. Extract treated vs control, then shrink LFC with apeglm
# ------------------------------------------------------------------
res <- results(
  dds,
  contrast = c("condition", "treated", "control"),
  alpha = 0.05
)

res_shrunk <- lfcShrink(
  dds,
  coef = target_coef,
  type = "apeglm"
)

# ------------------------------------------------------------------
# 7. Build the full result table (all filtered genes, incl. non-sig)
#    Annotate with gene_symbol and biotype ONLY.
#    NOTE: the instructor key also contains truth_log2FC_for_instructor
#    (the simulation ground truth). It is deliberately NOT used in any
#    analysis or significance decision - see week5_AI_verification_log.md.
# ------------------------------------------------------------------
anno <- read.csv(anno_file, check.names = FALSE) |>
  dplyr::select(gene_id, gene_symbol, biotype)

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  left_join(anno, by = "gene_id") |>
  mutate(
    significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(padj)

write.csv(res_df, "week5_deseq2_results.csv", row.names = FALSE)

n_sig <- sum(res_df$significant)
n_up  <- sum(res_df$direction == "Up in treated")
n_down <- sum(res_df$direction == "Down in treated")
cat("\nSignificant genes (padj<0.05 & |LFC|>=1):", n_sig,
    " (Up:", n_up, ", Down:", n_down, ")\n")
print(table(res_df$direction))

# ------------------------------------------------------------------
# 8. PCA (VST transformation, blind = FALSE)
# ------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df,
                aes(x = PC1, y = PC2, color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    title = "Week 5 RNA-seq PCA (VST, blind = FALSE)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition", shape = "Batch"
  ) +
  theme_bw(base_size = 12)

png("week5_pca.png", width = 7, height = 5, units = "in", res = 300)
print(p_pca)
dev.off()
cat("\nPCA saved. PC1 =", percent_var[1], "%, PC2 =", percent_var[2], "%\n")

# ------------------------------------------------------------------
# 9. Volcano plot (shrunken LFC vs -log10(padj))
# ------------------------------------------------------------------
plot_df <- res_df |>
  mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))

p_volcano <- ggplot(plot_df,
                    aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Up in treated" = "#C0392B",
    "Down in treated" = "#2F6DB3",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Differential expression: treated vs control",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p value",
    color = NULL
  ) +
  theme_bw(base_size = 12)

png("week5_de_plot.png", width = 7, height = 5, units = "in", res = 300)
print(p_volcano)
dev.off()
cat("Volcano plot saved.\n")

# ------------------------------------------------------------------
# 10. Reproducibility artifacts
# ------------------------------------------------------------------
saveRDS(dds, "week5_deseq2_object.rds")
capture.output(sessionInfo(), file = "session_info.txt")
cat("\nDESeq2 object and sessionInfo saved.\n")
cat("Done.\n")

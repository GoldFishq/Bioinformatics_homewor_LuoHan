# ============================================================
# Week 6 — 16S differential-abundance + PERMANOVA statistics.
# Uses the EMP-web backend engine for import + taxonomy prep to
# Genus, then computes the "interpret" statistics with vegan /
# wilcoxon (mirroring EMP's wilcox method) + effect sizes.
# ============================================================
REPO    <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web"
BACKEND <- file.path(REPO, "webapp", "backend")
TESTS   <- file.path(REPO, "tests")
OUT     <- "D:/学习资料/大三上/生物信息学/作业/Week6/analysis"
OUT_RES <- file.path(OUT, "outputs", "results")
dir.create(OUT_RES, showWarnings = FALSE, recursive = TRUE)

Sys.setenv(EMP_ROOT        = REPO)
Sys.setenv(EMP_DATA_DIR     = file.path(REPO, ".local_run", "data"))
Sys.setenv(BACKEND_DIR      = BACKEND)
Sys.setenv(EMP_BACKEND_DIR  = BACKEND)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S"), ""), ..., "\n", sep = "")

suppressPackageStartupMessages({
  library(EasyMultiProfiler); library(MultiAssayExperiment)
  library(SummarizedExperiment); library(S4Vectors)
  library(vegan); library(ggplot2)
})

.helpers <- c("storage.R","session.R","utils.R","auth.R","projects.R","plot_theme.R",
              "import.R","analysis.R","viz.R","workflow_registry.R",
              "workflow_metabolomics.R","workflow_metagenomics.R",
              "workflow_transcriptomics.R","workflow_chipseq.R",
              "workflow_microbiome_16s.R","workflow_microbiome_16s_api.R",
              "clinical.R","jobs.R","runall.R")
for (h in .helpers) { f <- file.path(BACKEND, "helpers", h); if (file.exists(f)) source(f) }

EXP <- "Microbiome_16S"; GROUP_VAR <- "Group"

## 1) session + import (species level)
sid <- create_session(); emp_register_session_owner(sid, "local")
imp <- import_omics_files(
  data_file       = file.path(TESTS, "16S_level-7.csv"),
  metadata_file   = file.path(TESTS, "16S_mapping.csv"),
  experiment_name = EXP, data_type = "tax", assay_name = "counts",
  start_level = "Species", tax_sep = ";", session_id = sid, owner_id = "local")
logmsg("imported: ", imp$n_samples, " samples x ", imp$n_features, " features")

## 2) prepare taxonomy to Genus, keep ALL genera (not just top 40)
m16s_prepare_taxonomy_step(sid, EXP, collapse_level = "Genus",
                           keep_top_n = 0L, drop_unassigned = FALSE)
empt <- load_empt(sid, EXP)
ad   <- as.matrix(SummarizedExperiment::assays(empt)[[1]])   # genus x sample
cd   <- as.data.frame(SummarizedExperiment::colData(empt))
logmsg("genus matrix: ", nrow(ad), " genera x ", ncol(ad), " samples")

## 3) keep only samples that have a Group (drop the 2 ungrouped K_XYL samples)
keep <- which(!is.na(cd[[GROUP_VAR]]) & nzchar(as.character(cd[[GROUP_VAR]])))
ad <- ad[, keep, drop = FALSE]; cd <- cd[keep, , drop = FALSE]
cd[[GROUP_VAR]] <- factor(cd[[GROUP_VAR]])
logmsg("grouped samples kept: ", ncol(ad), " | group sizes: ",
      paste(names(table(cd[[GROUP_VAR]])), table(cd[[GROUP_VAR]]), sep="=", collapse=", "))
write.csv(cd, file.path(OUT_RES, "coldata_grouped.csv"), row.names = TRUE)

## 4) relative abundance (compositional closure) per sample
rel <- ad
rel <- sweep(rel, 2, colSums(rel, na.rm = TRUE), `/`)
rel[!is.finite(rel)] <- 0

## 5) alpha diversity (Shannon, Observed, Simpson)
alpha_df <- data.frame(
  sample   = colnames(rel),
  Group    = cd[[GROUP_VAR]],
  shannon  = diversity(t(rel), index = "shannon"),
  observed = specnumber(t(rel)),
  simpson  = diversity(t(rel), index = "simpson"),
  stringsAsFactors = FALSE)
write.csv(alpha_df, file.path(OUT_RES, "alpha_by_sample.csv"), row.names = FALSE)
# per-group summary + Kruskal-Wallis + pairwise wilcox (before vs after)
alpha_summary <- function(metric) {
  v <- alpha_df[[metric]]; g <- alpha_df$Group
  kw <- kruskal.test(v ~ g)$p.value
  grps <- levels(g)
  rows <- data.frame()
  for (lv in grps) rows <- rbind(rows, data.frame(group = lv, mean = mean(v[g == lv]),
                                                  sd = sd(v[g == lv]), median = median(v[g == lv])))
  # disease-specific before vs after
  pw <- data.frame()
  for (dz in c("IBS", "UC")) {
    a <- v[g == paste0(dz, "_before")]; b <- v[g == paste0(dz, "_after")]
    wt <- wilcox.test(a, b, exact = FALSE)
    pw <- rbind(pw, data.frame(contrast = paste0(dz, "_before vs _after"),
                               mean_before = mean(a), mean_after = mean(b),
                               p = wt$p.value))
  }
  list(kw = kw, by_group = rows, paired = pw)
}
ash <- alpha_summary("shannon"); aob <- alpha_summary("observed"); asi <- alpha_summary("simpson")
write.csv(ash$by_group,  file.path(OUT_RES, "alpha_shannon_bygroup.csv"), row.names = FALSE)
write.csv(ash$paired,    file.path(OUT_RES, "alpha_shannon_paired.csv"), row.names = FALSE)
write.csv(aob$by_group,  file.path(OUT_RES, "alpha_observed_bygroup.csv"), row.names = FALSE)
write.csv(aob$paired,    file.path(OUT_RES, "alpha_observed_paired.csv"), row.names = FALSE)
logmsg("Shannon KW p = ", signif(ash$kw, 4),
       " | IBS b>a p = ", signif(ash$paired$p[1], 4),
       " | UC b>a p = ", signif(ash$paired$p[2], 4))

## 6) beta diversity: Bray-Curtis + PCoA + PERMANOVA (+ dispersion check)
bray <- vegdist(t(rel), method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)
pcoa_xy <- as.data.frame(pcoa$points)
names(pcoa_xy) <- c("PCo1", "PCo2")
pcoa_xy$sample <- rownames(pcoa_xy); pcoa_xy$Group <- cd[[GROUP_VAR]]
write.csv(pcoa_xy, file.path(OUT_RES, "pcoa_coords.csv"), row.names = FALSE)

# derive Disease + Time factors
cd$Disease <- sub("_(before|after)", "", cd[[GROUP_VAR]])
cd$Time    <- sub("^(IBS|UC)_", "", cd[[GROUP_VAR]])
perm_group <- adonis2(bray ~ Group, data = cd, permutations = 999)
perm_dz    <- adonis2(bray ~ Disease, data = cd, permutations = 999)
perm_tm    <- adonis2(bray ~ Time, data = cd, permutations = 999)
perm_dz_tm <- adonis2(bray ~ Disease * Time, data = cd, permutations = 999)
capt <- function(o) data.frame(term = rownames(o), R2 = o$R2, F = o$F, p = o$`Pr(>F)`,
                               stringsAsFactors = FALSE)
perm_all <- rbind(data.frame(model = "Group",    capt(perm_group),   stringsAsFactors = FALSE),
                  data.frame(model = "Disease",  capt(perm_dz),      stringsAsFactors = FALSE),
                  data.frame(model = "Time",     capt(perm_tm),      stringsAsFactors = FALSE),
                  data.frame(model = "Disease*Time", capt(perm_dz_tm), stringsAsFactors = FALSE))
write.csv(perm_all, file.path(OUT_RES, "permanova.csv"), row.names = FALSE)
# dispersion homogeneity (betadisper / PERMDISP)
bd <- betadisper(bray, cd[[GROUP_VAR]])
permdisp <- permutest(bd, permutations = 999)
disp_row <- data.frame(test = "betadisper(Group)",
                       F = permdisp$tab$F, p = permdisp$tab$`Pr(>F)`, stringsAsFactors = FALSE)
write.csv(disp_row, file.path(OUT_RES, "dispersion.csv"), row.names = FALSE)
logmsg("PERMANOVA Group R2 = ", signif(perm_group$R2[1], 4), " p = ", signif(perm_group$`Pr(>F)`[1], 4),
       " | Time R2 = ", signif(perm_tm$R2[1], 4), " p = ", signif(perm_tm$`Pr(>F)`[1], 4),
       " | dispersion p = ", signif(permdisp$tab$`Pr(>F)`, 4))

## 7) differential abundance (Wilcoxon on relative abundance), disease-stratified
cliff <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(p = NA, delta = NA, auc = NA))
  wt <- wilcox.test(x, y, exact = FALSE)
  U  <- as.numeric(wt$statistic); auc <- U / (n1 * n2)
  list(p = wt$p.value, delta = 2 * auc - 1, auc = auc)
}
diff_by_disease <- function(dz) {
  g_before <- paste0(dz, "_before"); g_after <- paste0(dz, "_after")
  b <- rel[, cd[[GROUP_VAR]] == g_before, drop = FALSE]
  a <- rel[, cd[[GROUP_VAR]] == g_after, drop = FALSE]
  out <- data.frame()
  for (taxa in rownames(rel)) {
    mb <- as.numeric(b[taxa, ]); ma <- as.numeric(a[taxa, ])
    if (sum(mb > 0, na.rm = TRUE) < 3 && sum(ma > 0, na.rm = TRUE) < 3) next
    cf <- cliff(mb, ma)
    out <- rbind(out, data.frame(
      genus = taxa,
      mean_rel_before = mean(mb, na.rm = TRUE),
      mean_rel_after  = mean(ma, na.rm = TRUE),
      log2fc_rel = log2((mean(ma, na.rm = TRUE) + 1e-6) / (mean(mb, na.rm = TRUE) + 1e-6)),
      cliffs_delta = cf$delta, p = cf$p))
  }
  out$q <- p.adjust(out$p, method = "BH")
  out[order(out$p), ]
}
dib <- diff_by_disease("IBS")
diu <- diff_by_disease("UC")
write.csv(dib, file.path(OUT_RES, "diff_taxa_IBS.csv"), row.names = FALSE)
write.csv(diu, file.path(OUT_RES, "diff_taxa_UC.csv"),  row.names = FALSE)
logmsg("IBS diff taxa (BH q<0.05): ", sum(dib$q < 0.05, na.rm = TRUE),
       " | UC diff taxa (BH q<0.05): ", sum(diu$q < 0.05, na.rm = TRUE))

## 8) top taxa overall (for bar/heatmap context)
top <- data.frame(genus = rownames(rel), mean_rel = rowMeans(rel, na.rm = TRUE))
top <- top[order(-top$mean_rel), ][1:40, ]
write.csv(top, file.path(OUT_RES, "top40_genera.csv"), row.names = FALSE)

## 9) baseline disease difference (IBS_before vs UC_before)
db <- cliff(rel[, cd[[GROUP_VAR]] == "IBS_before", drop = FALSE],
            rel[, cd[[GROUP_VAR]] == "UC_before", drop = FALSE])
logmsg("DONE. results in ", OUT_RES)
saveRDS(list(session_id = sid, experiment = EXP), file.path(OUT, "session_info_stats.rds"))

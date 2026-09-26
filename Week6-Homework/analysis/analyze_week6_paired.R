# ============================================================
# Week 6 — PAIRED 16S analysis (correct design: each patient
# has a before _01 and after _02 sample). EMP-web engine for
# import + Genus prep; paired statistics with vegan / wilcox.
# ============================================================
REPO    <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web"
BACKEND <- file.path(REPO, "webapp", "backend"); TESTS <- file.path(REPO, "tests")
OUT     <- "D:/学习资料/大三上/生物信息学/作业/Week6/analysis"
OUT_RES <- file.path(OUT, "outputs", "results"); dir.create(OUT_RES, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(EMP_ROOT = REPO, EMP_DATA_DIR = file.path(REPO, ".local_run", "data"),
           BACKEND_DIR = BACKEND, EMP_BACKEND_DIR = BACKEND)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                            paste(..., collapse = " "), "\n", sep = "")
suppressPackageStartupMessages({ library(EasyMultiProfiler); library(MultiAssayExperiment)
  library(SummarizedExperiment); library(S4Vectors); library(vegan); library(ggplot2) })
.helpers <- c("storage.R","session.R","utils.R","auth.R","projects.R","plot_theme.R","import.R",
  "analysis.R","viz.R","workflow_registry.R","workflow_metabolomics.R","workflow_metagenomics.R",
  "workflow_transcriptomics.R","workflow_chipseq.R","workflow_microbiome_16s.R",
  "workflow_microbiome_16s_api.R","clinical.R","jobs.R","runall.R")
for (h in .helpers) { f <- file.path(BACKEND, "helpers", h); if (file.exists(f)) source(f) }
EXP <- "Microbiome_16S"; GROUP_VAR <- "Group"
sid <- create_session(); emp_register_session_owner(sid, "local")
imp <- import_omics_files(data_file = file.path(TESTS, "16S_level-7.csv"),
  metadata_file = file.path(TESTS, "16S_mapping.csv"), experiment_name = EXP,
  data_type = "tax", assay_name = "counts", start_level = "Species", tax_sep = ";",
  session_id = sid, owner_id = "local")
m16s_prepare_taxonomy_step(sid, EXP, collapse_level = "Genus", keep_top_n = 0L, drop_unassigned = FALSE)
empt <- load_empt(sid, EXP)
ad <- as.matrix(SummarizedExperiment::assays(empt)[[1]])
cd <- as.data.frame(SummarizedExperiment::colData(empt))
keep <- which(!is.na(cd[[GROUP_VAR]]) & nzchar(as.character(cd[[GROUP_VAR]])))
ad <- ad[, keep, drop = FALSE]; cd <- cd[keep, , drop = FALSE]
cd[[GROUP_VAR]] <- factor(cd[[GROUP_VAR]])
cd$patient <- sub("_\\d+$", "", colnames(ad))
cd$Disease <- sub("_(before|after)", "", cd[[GROUP_VAR]])
cd$Time    <- sub("^(IBS|UC)_", "", cd[[GROUP_VAR]])
cd$Pheno   <- sub("^(IBS|UC)_(before|after)_", "", cd$Group_sub)  # poor / great
logmsg("grouped samples: ", ncol(ad), " patients: ", length(unique(cd$patient)))
rel <- ad; rel <- sweep(rel, 2, colSums(rel, na.rm = TRUE), `/`); rel[!is.finite(rel)] <- 0

## ---- Alpha (paired within disease) ----
alpha_df <- data.frame(sample = colnames(rel), Group = cd[[GROUP_VAR]],
  shannon = diversity(t(rel), "shannon"), observed = specnumber(t(rel)),
  simpson = diversity(t(rel), "simpson"), stringsAsFactors = FALSE)
alpha_paired <- function(metric) {
  v <- alpha_df[[metric]]; g <- alpha_df$Group; pidv <- cd$patient; out <- data.frame()
  for (dz in c("IBS", "UC")) {
    bv <- v[g == paste0(dz, "_before")]; names(bv) <- pidv[g == paste0(dz, "_before")]
    av <- v[g == paste0(dz, "_after")];  names(av) <- pidv[g == paste0(dz, "_after")]
    common <- intersect(names(bv), names(av))
    bv <- bv[common]; av <- av[common]
    if (length(common) >= 3) {
      wt <- wilcox.test(bv, av, paired = TRUE, exact = FALSE)
      out <- rbind(out, data.frame(disease = dz, n_pairs = length(common),
        mean_before = mean(bv), mean_after = mean(av),
        median_delta = median(av - bv), wilcox_p = wt$p.value, stringsAsFactors = FALSE))
    }
  }; out
}
ap_sh <- alpha_paired("shannon"); ap_ob <- alpha_paired("observed"); ap_si <- alpha_paired("simpson")
write.csv(ap_sh, file.path(OUT_RES, "paired_alpha_shannon.csv"), row.names = FALSE)
write.csv(ap_ob, file.path(OUT_RES, "paired_alpha_observed.csv"), row.names = FALSE)
write.csv(ap_si, file.path(OUT_RES, "paired_alpha_simpson.csv"), row.names = FALSE)
logmsg("PAIRED alpha Shannon: IBS p=", signif(ap_sh$wilcox_p[1],4), " (delta ",
      signif(ap_sh$median_delta[1],3), "); UC p=", signif(ap_sh$wilcox_p[2],4),
      " (delta ", signif(ap_sh$median_delta[2],3), ")")

## ---- Beta: PCoA + PERMANOVA (independent, community-level) + paired shift ----
bray <- vegdist(t(rel), "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)
pxy <- as.data.frame(pcoa$points); names(pxy) <- c("PCo1","PCo2")
pxy$sample <- rownames(pxy); pxy$Group <- cd[[GROUP_VAR]]
write.csv(pxy, file.path(OUT_RES, "pcoa_coords.csv"), row.names = FALSE)
perm_g <- adonis2(bray ~ Group, data = cd, permutations = 999)
perm_dz <- adonis2(bray ~ Disease, data = cd, permutations = 999)
perm_tm <- adonis2(bray ~ Time, data = cd, permutations = 999)
perm_dz_tm <- adonis2(bray ~ Disease * Time, data = cd, permutations = 999)
capt <- function(o) data.frame(term = rownames(o), R2 = o$R2, F = o$F, p = o$`Pr(>F)`, stringsAsFactors = FALSE)
perm_all <- rbind(data.frame(model="Group",capt(perm_g),stringsAsFactors=FALSE),
                  data.frame(model="Disease",capt(perm_dz),stringsAsFactors=FALSE),
                  data.frame(model="Time",capt(perm_tm),stringsAsFactors=FALSE),
                  data.frame(model="Disease*Time",capt(perm_dz_tm),stringsAsFactors=FALSE))
write.csv(perm_all, file.path(OUT_RES, "permanova.csv"), row.names = FALSE)
# dispersion homogeneity
bd <- betadisper(bray, cd[[GROUP_VAR]]); pd <- permutest(bd, permutations = 999)
write.csv(data.frame(test="betadisper(Group)", F=pd$tab$F, p=pd$tab$`Pr(>F)`, stringsAsFactors=FALSE),
          file.path(OUT_RES, "dispersion.csv"), row.names = FALSE)
# paired within-patient shift vs population spread
pid_of <- cd$patient; g_of <- cd[[GROUP_VAR]]
# build before/after index per patient
bidx <- setNames(which(g_of == "IBS_before" | g_of == "UC_before"), pid_of[g_of == "IBS_before" | g_of == "UC_before"])
aidx <- setNames(which(g_of == "IBS_after"  | g_of == "UC_after"),  pid_of[g_of == "IBS_after"  | g_of == "UC_after"])
common_p <- intersect(names(bidx), names(aidx))
wd <- sapply(common_p, function(p) as.numeric(bray[bidx[p], aidx[p]]))
mean_within <- mean(wd, na.rm = TRUE); mean_overall <- mean(as.dist(bray), na.rm = TRUE)
write.csv(data.frame(metric="within_patient_Bray", value = mean_within, n_pairs = length(wd),
                     ratio_to_population = mean_within / mean_overall, stringsAsFactors = FALSE),
          file.path(OUT_RES, "paired_beta_shift.csv"), row.names = FALSE)
logmsg("PERMANOVA Group R2=", signif(perm_g$R2[1],4), " p=", signif(perm_g$`Pr(>F)`[1],4),
       " | Time R2=", signif(perm_tm$R2[1],4), " p=", signif(perm_tm$`Pr(>F)`[1],4),
       " | dispersion p=", signif(pd$tab$`Pr(>F)`,4),
       " | within-patient Bray=", signif(mean_within,4), " (", signif(mean_within/mean_overall,3), "x population)")

## ---- Differential abundance (PAIRED Wilcoxon signed-rank) ----
diff_paired <- function(dz) {
  g <- cd$Group; pidv <- cd$patient
  bsel <- colnames(rel)[g == paste0(dz, "_before")]; asel <- colnames(rel)[g == paste0(dz, "_after")]
  bpid <- setNames(pidv[g == paste0(dz, "_before")], bsel); apid <- setNames(pidv[g == paste0(dz, "_after")], asel)
  common <- intersect(bpid, apid)
  bsel <- bsel[bpid %in% common]; asel <- asel[apid %in% common]
  bmat <- rel[, bsel, drop = FALSE]; amat <- rel[, asel, drop = FALSE]
  ob <- order(bpid[bsel]); oa <- order(apid[asel])
  bmat <- bmat[, ob]; amat <- amat[, oa]
  stopifnot(all(bpid[bsel][ob] == apid[asel][oa]))
  out <- data.frame()
  for (taxa in rownames(rel)) {
    x <- as.numeric(bmat[taxa, ]); y <- as.numeric(amat[taxa, ])
    if (sum(x > 0, na.rm = TRUE) < 3 && sum(y > 0, na.rm = TRUE) < 3) next
    wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
    out <- rbind(out, data.frame(genus = taxa,
      mean_rel_before = mean(x, na.rm = TRUE), mean_rel_after = mean(y, na.rm = TRUE),
      log2fc_rel = log2((mean(y, na.rm = TRUE) + 1e-6) / (mean(x, na.rm = TRUE) + 1e-6)),
      wilcox_p = wt$p.value, stringsAsFactors = FALSE))
  }
  out$q <- p.adjust(out$wilcox_p, "BH")
  out[order(out$wilcox_p), ]
}
dib <- diff_paired("IBS"); diu <- diff_paired("UC")
write.csv(dib, file.path(OUT_RES, "diff_paired_IBS.csv"), row.names = FALSE)
write.csv(diu, file.path(OUT_RES, "diff_paired_UC.csv"),  row.names = FALSE)
logmsg("PAIRED diff (BH q<0.05): IBS=", sum(dib$q < 0.05, na.rm = TRUE),
       " UC=", sum(diu$q < 0.05, na.rm = TRUE))

## ---- Phenotype (poor vs great) baseline difference (independent) ----
pheno_alpha <- data.frame()
for (dz in c("IBS","UC")) {
  for (met in c("shannon","observed")) {
    v <- alpha_df[[met]]; ph <- cd$Pheno; gg <- cd$Group
    sel <- gg == paste0(dz, "_before")
    a <- v[sel & ph == "poor"]; b <- v[sel & ph == "great"]
    if (length(a) >= 3 && length(b) >= 3) {
      wt <- wilcox.test(a, b, exact = FALSE)
      pheno_alpha <- rbind(pheno_alpha, data.frame(disease = dz, metric = met,
        pheno = "poor_vs_great(before)", mean_poor = mean(a), mean_great = mean(b),
        wilcox_p = wt$p.value, stringsAsFactors = FALSE))
    }
  }
}
write.csv(pheno_alpha, file.path(OUT_RES, "pheno_alpha.csv"), row.names = FALSE)
logmsg("phenotype baseline alpha poor_vs_great: ", paste(pheno_alpha$wilcox_p, collapse=", "))

## ---- Top taxa ----
top <- data.frame(genus = rownames(rel), mean_rel = rowMeans(rel, na.rm = TRUE))
top <- top[order(-top$mean_rel), ][1:40, ]; write.csv(top, file.path(OUT_RES, "top40_genera.csv"), row.names = FALSE)
logmsg("DONE -> ", OUT_RES)

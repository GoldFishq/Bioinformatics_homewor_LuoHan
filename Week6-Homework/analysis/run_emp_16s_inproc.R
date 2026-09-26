# ============================================================
# Week 6 — 16S Microbiome whole-analysis, executed with the
# EasyMultiProfiler-Web *backend engine* (same helpers the
# EMP-web plumber API loads) in a standalone R session.
#
# Why not over HTTP: the running plumber process (started 20:55)
# fails every atomic JSON overwrite with
#   "Unable to persist EMP state: <...session_owners/<sid>.json>"
# (see analysis/NOTES_empweb_bug.md). Running the same helper
# functions in-process avoids that persistence path and yields
# identical results.
# ============================================================
REPO    <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web"
BACKEND <- file.path(REPO, "webapp", "backend")
TESTS   <- file.path(REPO, "tests")
OUT     <- "D:/学习资料/大三上/生物信息学/作业/Week6/analysis"
OUT_PLOTS <- file.path(OUT, "outputs", "plots")
OUT_RES   <- file.path(OUT, "outputs", "results")
dir.create(OUT_PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_RES,  showWarnings = FALSE, recursive = TRUE)

Sys.setenv(EMP_ROOT      = REPO)
Sys.setenv(EMP_DATA_DIR  = file.path(REPO, ".local_run", "data"))
Sys.setenv(BACKEND_DIR   = BACKEND)
Sys.setenv(EMP_BACKEND_DIR = BACKEND)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

suppressPackageStartupMessages({
  library(EasyMultiProfiler)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(ggplot2)
  library(base64enc)
})

.helpers <- c("storage.R","session.R","utils.R","auth.R","projects.R","plot_theme.R",
              "import.R","analysis.R","viz.R","workflow_registry.R",
              "workflow_metabolomics.R","workflow_metagenomics.R",
              "workflow_transcriptomics.R","workflow_chipseq.R",
              "workflow_microbiome_16s.R","workflow_microbiome_16s_api.R",
              "clinical.R","jobs.R","runall.R")
for (h in .helpers) {
  f <- file.path(BACKEND, "helpers", h)
  if (file.exists(f)) { source(f); logmsg("sourced ", h) } else logmsg("MISSING ", h)
}

EXP <- "Microbiome_16S"
GROUP_VAR <- "Group"

## ---- 1) session ------------------------------------------------
sid <- create_session()
emp_register_session_owner(sid, "local")
logmsg("session_id = ", sid)
sess_dir <- session_path(sid)
logmsg("session dir = ", sess_dir)

## ---- 2) import (tax + metadata) --------------------------------
logmsg("importing 16S taxonomy + mapping ...")
imp <- import_omics_files(
  data_file       = file.path(TESTS, "16S_level-7.csv"),
  metadata_file   = file.path(TESTS, "16S_mapping.csv"),
  experiment_name = EXP,
  data_type       = "tax",
  assay_name      = "counts",
  start_level     = "Species",
  tax_sep         = ";",
  session_id      = sid,
  owner_id        = "local"
)
str(imp, max.level = 1)

## ---- 3) profile / validate -------------------------------------
prof <- m16s_profile(sid, EXP, tax_sep = ";")
logmsg("profile: n_samples=", prof$n_samples, " n_features=", prof$n_features,
       " max_depth=", prof$max_taxonomy_depth)
print(prof$observed_levels)
writeLines(capture.output(str(prof)), file.path(OUT_RES, "01_profile.txt"))

cd <- as.data.frame(SummarizedExperiment::colData(load_empt(sid, EXP)))
write.csv(cd, file.path(OUT_RES, "01_coldata.csv"), row.names = FALSE)
logmsg("colData columns: ", paste(names(cd), collapse = ", "))
print(table(cd[[GROUP_VAR]]))
print(table(cd[["Group_sub"]]))

## ---- 4) whole pipeline (EMP-web one-click, same function) ------
logmsg("running run_all_m16s (whole 16S pipeline) ...")
bundle_res <- run_all_m16s(
  session_id     = sid,
  experiment     = EXP,
  group_var      = GROUP_VAR,
  taxonomy_level = "Genus",
  alpha_index    = "shannon",
  beta_method    = "bray",
  ord_method     = "PCoA"
)
str(bundle_res, max.level = 1)

logmsg("DONE stage 1. bundle = ", bundle_res$bundle)
saveRDS(list(session_id = sid, experiment = EXP, bundle = bundle_res$bundle, zip = bundle_res$zip_path),
        file.path(OUT, "session_info.rds"))

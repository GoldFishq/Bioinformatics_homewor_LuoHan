# ============================================================
# Week 6 Homework — 16S Microbiome analysis driven through
# the EasyMultiProfiler-Web (EMP-web) backend API.
#
# Data : tests/16S_level-7.csv  (species-level taxonomy, EMP format)
#        tests/16S_mapping.csv  (SampleID, Group, Group_sub)
# Backend: http://127.0.0.1:8000  (EMP-web v9.0.4, plumber)
# ============================================================
suppressPackageStartupMessages({
  library(httr2); library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
BASE <- "http://127.0.0.1:8000"
DATA_DIR <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web/tests"
OUT <- "D:/学习资料/大三上/生物信息学/作业/Week6/analysis"
OUT_PLOTS <- file.path(OUT, "outputs", "plots")
OUT_RES  <- file.path(OUT, "outputs", "results")
dir.create(OUT_PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_RES,  showWarnings = FALSE, recursive = TRUE)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

post_json <- function(path, body, timeout = 1800) {
  req <- request(paste0(BASE, path)) |>
    req_body_json(body) |>
    req_timeout(timeout) |>
    req_error(is_error = function(r) FALSE)
  resp <- req_perform(req)
  out <- resp_body_json(resp)
  out
}
get_json <- function(path, timeout = 300) {
  req <- request(paste0(BASE, path)) |> req_timeout(timeout) |>
    req_error(is_error = function(r) FALSE)
  resp_body_json(req_perform(req))
}

## ---------------------------------------------------------
## 1) Create a fresh analysis session
## ---------------------------------------------------------
logmsg("Creating session ...")
s <- post_json("/api/session", list())
sid <- s$session_id
logmsg("session_id = ", sid)
stopifnot(!is.null(sid))

## ---------------------------------------------------------
## 2) Import 16S taxonomy table + sample metadata
## ---------------------------------------------------------
logmsg("Importing 16S data (tax) ...")
req <- request(paste0(BASE, "/api/import")) |>
  req_body_multipart(
    data_file       = curl::form_file(file.path(DATA_DIR, "16S_level-7.csv")),
    metadata_file   = curl::form_file(file.path(DATA_DIR, "16S_mapping.csv")),
    experiment_name = "Microbiome_16S",
    data_type       = "tax",
    assay_name      = "counts",
    start_level     = "Species",
    tax_sep         = ";",
    session_id      = sid
  ) |> req_timeout(1800) |> req_error(is_error = function(r) FALSE)
imp <- resp_body_json(req_perform(req))
logmsg("import result:")
str(imp, max.level = 2)

EXP <- "Microbiome_16S"

## ---------------------------------------------------------
## 3) Profile / validate taxonomy depth
## ---------------------------------------------------------
logmsg("Profiling taxonomy ...")
prof <- post_json("/api/workflows/microbiome_16s/profile",
                  list(session_id = sid, experiment = EXP, tax_sep = ";"))
print(prof)
val <- post_json("/api/workflows/microbiome_16s/validate",
                 list(session_id = sid, experiment = EXP, tax_sep = ";"))
print(val)

## Summary of the imported experiment
summ <- get_json(sprintf("/api/summary/%s/%s", sid, EXP))
logmsg("summary:")
str(summ, max.level = 2)
write_json(summ, file.path(OUT_RES, "01_experiment_summary.json"), pretty = TRUE, auto_unbox = TRUE)

## colData (sample metadata as seen by EMP)
cd <- get_json(sprintf("/api/coldata/%s/%s", sid, EXP))
write_json(cd, file.path(OUT_RES, "01_coldata.json"), pretty = TRUE, auto_unbox = TRUE)

## ---------------------------------------------------------
## 4) Prepare taxonomy: collapse to Genus + keep top 40
##    (mirrors the run_all pipeline step 1)
## ---------------------------------------------------------
logmsg("Prepare taxonomy -> Genus ...")
prep <- post_json("/api/workflows/microbiome_16s/prepare/taxonomy",
                  list(session_id = sid, experiment = EXP, collapse_level = "Genus",
                       min_total_abundance = 0, drop_unassigned = FALSE,
                       keep_top_n = 40L, tax_sep = ";", normalize_method = "none"))
logmsg("prepare result:")
str(prep[setdiff(names(prep), "preview_data")], max.level = 2)

## ---------------------------------------------------------
## 5) One-click whole-analysis pipeline (alpha + beta + top taxa + diff)
## ---------------------------------------------------------
logmsg("Submitting run_all (whole 16S pipeline) ...")
ra <- post_json("/api/workflows/microbiome_16s/run_all",
                list(session_id = sid, experiment = EXP, group_var = "Group",
                     taxonomy_level = "Genus", alpha_index = "shannon",
                     beta_method = "bray", ord_method = "PCoA"))
logmsg("run_all job: ", ra$job_id)
stopifnot(!is.null(ra$job_id))

## poll
repeat {
  Sys.sleep(5)
  st <- get_json(sprintf("/api/jobs/%s", ra$job_id))
  status <- st$job$status
  logmsg("  job status = ", status, "  progress = ", st$job$progress %||% "?", "  ",
         st$job$message %||% "")
  if (!is.null(status) && status %in% c("done", "error", "failed", "cancelled")) break
}
jobres <- get_json(sprintf("/api/jobs/%s/result", ra$job_id))
logmsg("run_all finished:")
str(jobres[setdiff(names(jobres), "data")], max.level = 2)
write_json(jobres[setdiff(names(jobres), "data")],
           file.path(OUT_RES, "02_runall_result.json"), pretty = TRUE, auto_unbox = TRUE)

## ---------------------------------------------------------
## 6) Download the produced bundle (all plots + tables)
## ---------------------------------------------------------
logmsg("Listing bundles ...")
bun <- get_json(sprintf("/api/bundles/%s", sid))
print(bun)
if (length(bun$bundles)) {
  bn <- bun$bundles[[length(bun$bundles)]]$name
  logmsg("Downloading bundle: ", bn)
  req <- request(sprintf("%s/api/bundles/%s/%s", BASE, sid, bn)) |>
    req_timeout(600) |> req_error(is_error = function(r) FALSE)
  zipfile <- file.path(OUT, "outputs", bn)
  req_perform(req, path = zipfile)
  logmsg("bundle saved -> ", zipfile, " (", round(file.size(zipfile)/1024), " KB)")
}

saveRDS(list(session_id = sid, experiment = EXP), file.path(OUT, "session_info.rds"))
logmsg("DONE. session_id = ", sid)

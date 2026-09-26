suppressPackageStartupMessages(library(jsonlite))
d <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web/.local_run/data/projects/session_owners"
atomic <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".emp-json-", tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  jsonlite::write_json(value, tmp, auto_unbox = TRUE, null = "null", pretty = TRUE)
  ok <- isTRUE(file.rename(tmp, path))
  cat("  rename ok:", ok, " tmp existed:", file.exists(tmp), "\n")
  ok
}
f <- file.path(d, "aaa_twocall.json")
cat("call 1 (dest absent):\n"); atomic(f, list(session_id = "aaa", n = 1))
cat("call 2 (dest exists):\n"); atomic(f, list(session_id = "aaa", n = 2))
cat("call 3 (dest exists):\n"); atomic(f, list(session_id = "aaa", n = 3))
cat("final content:", paste(readLines(f), collapse = " "), "\n")
unlink(f)
cat("tempdir:", tempdir(), "\n")
cat("Sys.getlocale(LC_CTYPE):", Sys.getlocale("LC_CTYPE"), "\n")

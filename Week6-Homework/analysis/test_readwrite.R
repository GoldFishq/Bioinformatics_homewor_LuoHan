suppressPackageStartupMessages(library(jsonlite))
d <- "D:/学习资料/大三上/生物信息学/EasyMultiProfiler-Web/.local_run/data/projects/session_owners"
f <- file.path(d, "aaa_readwrite.json")
unlink(f)
# 1st write
dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile(".emp-json-", tmpdir = dirname(f))
write_json(list(v = 1), tmp, auto_unbox = TRUE, pretty = TRUE)
cat("write1 rename:", file.rename(tmp, f), "\n")

# read it back (like emp_assert_session_owner does)
ex <- read_json(f, simplifyVector = FALSE)
cat("read ok, v =", ex$v, "\n")

# 2nd write
tmp2 <- tempfile(".emp-json-", tmpdir = dirname(f))
write_json(list(v = 2), tmp2, auto_unbox = TRUE, pretty = TRUE)
cat("write2 rename:", file.rename(tmp2, f), "\n")

# 3rd: read again then write
ex <- read_json(f, simplifyVector = FALSE)
tmp3 <- tempfile(".emp-json-", tmpdir = dirname(f))
write_json(list(v = 3), tmp3, auto_unbox = TRUE, pretty = TRUE)
cat("write3 rename:", file.rename(tmp3, f), "\n")
cat("content:", readLines(f)[grep("v", readLines(f))], "\n")
unlink(f)

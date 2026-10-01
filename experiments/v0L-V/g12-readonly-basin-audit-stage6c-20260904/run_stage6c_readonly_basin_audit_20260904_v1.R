#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
parse_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  value <- arguments[startsWith(arguments, prefix)]
  if (!length(value)) return("")
  sub(prefix, "", value[[1L]], fixed = TRUE)
}

stage6a_root <- parse_argument("stage6a-root")
stage6b_root <- parse_argument("stage6b-root")
output_root <- parse_argument("output-root")
tools_path <- "stage6c_basin_tools_20260904_v1.R"
if (!file.exists(tools_path)) {
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                          value = TRUE)[[1L]])
  tools_path <- file.path(dirname(script_path), basename(tools_path))
}
source(tools_path, local = FALSE)

if (!nzchar(stage6a_root) || !nzchar(stage6b_root) || !nzchar(output_root)) {
  stop("--stage6a-root, --stage6b-root and --output-root are required.",
       call. = FALSE)
}

started <- Sys.time()
result <- s6c_audit(stage6a_root, stage6b_root, output_root)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
cat("STAGE6C_TRUTH_FREE_AUDIT_PASS endpoints=", nrow(result$endpoints),
    " pairs=", nrow(result$pairs), " elapsed_seconds=",
    sprintf("%.3f", elapsed), "\n", sep = "")

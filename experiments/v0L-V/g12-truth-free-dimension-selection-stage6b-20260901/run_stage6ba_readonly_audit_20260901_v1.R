arguments <- commandArgs(trailingOnly = TRUE)
parse_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  value <- arguments[startsWith(arguments, prefix)]
  if (!length(value)) return("")
  sub(prefix, "", value[[1L]], fixed = TRUE)
}

stage6a_root <- parse_argument("stage6a-root")
output_root <- parse_argument("output-root")
if (!nzchar(stage6a_root) || !nzchar(output_root)) {
  stop("--stage6a-root and --output-root are required.", call. = FALSE)
}

script_argument <- commandArgs(trailingOnly = FALSE)
file_argument <- script_argument[startsWith(script_argument, "--file=")]
script_file <- if (length(file_argument)) {
  sub("--file=", "", file_argument[[1L]], fixed = TRUE)
} else {
  "run_stage6ba_readonly_audit_20260901_v1.R"
}
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
source(file.path(script_root,
                 "stage6b_truth_free_tools_20260901_v1.R"), local = TRUE)

s6b_audit_stage6a(stage6a_root, output_root)
cat("STAGE6B_A_READONLY_AUDIT_PASS output_root=", output_root, "\n", sep = "")

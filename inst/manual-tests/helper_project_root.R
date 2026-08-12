.working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(.working_dir, "DESCRIPTION"))) {
  # Prefer the working directory on Windows: a Chinese --file path can be
  # decoded with the startup code page before R switches to UTF-8.
  .candidate <- .working_dir
} else {
  .test_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(.test_script_arg)) {
    stop("Run this test from the multiFSYNC package root.")
  }
  .candidate <- dirname(normalizePath(
    sub("^--file=", "", .test_script_arg[1]), winslash = "/", mustWork = TRUE
  ))
}
while (!file.exists(file.path(.candidate, "DESCRIPTION"))) {
  .parent <- dirname(.candidate)
  if (identical(.parent, .candidate)) stop("Could not locate the multiFSYNC package root.")
  .candidate <- .parent
}
project_root <- normalizePath(.candidate, winslash = "/", mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
.libPaths(c(normalizePath(args[[1L]], mustWork = TRUE), .libPaths()))
cat("lib_paths=", paste(.libPaths(), collapse = "|"), "\n", sep = "")
library(multiFSYNC)
cat("package_path=", system.file(package = "multiFSYNC"), "\n", sep = "")
cat("package_version=", as.character(packageVersion("multiFSYNC")), "\n", sep = "")
cat("pre_score_exported=",
    "bayesSYNC_multi_pre_score" %in% getNamespaceExports("multiFSYNC"),
    "\n", sep = "")

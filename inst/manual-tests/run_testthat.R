library(splines)
library(parallel)
library(testthat)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
src_dir <- file.path(project_root, "R")
for (f in list.files(src_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

test_dir(file.path(project_root, "tests", "testthat"))

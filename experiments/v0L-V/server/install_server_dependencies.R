# Install only missing CRAN dependencies for the v0L-V rc1 snapshot.
# Target Linux versions were recorded after native server acceptance.

script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run this installer with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]),
                winslash = "/", mustWork = TRUE)
}

if (getRversion() < "4.0.0") {
  stop("multiFSYNC requires R >= 4.0.0; found ", getRversion(), ".")
}

dependency_file <- file.path(
  dirname(script_path()), "SERVER_DEPENDENCIES.csv"
)
dependencies <- utils::read.csv(
  dependency_file, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "package", "role", "windows_reference_version",
  "is_base_or_recommended", "source", "linux_actual_version"
)
if (!identical(names(dependencies), required) ||
    anyDuplicated(dependencies$package) || any(!nzchar(dependencies$package))) {
  stop("SERVER_DEPENDENCIES.csv is malformed.")
}

distribution_packages <- dependencies$package[
  dependencies$source == "R_distribution"
]
missing_distribution <- distribution_packages[!vapply(
  distribution_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_distribution)) {
  stop(
    "Required base/recommended packages are absent from this R installation: ",
    paste(missing_distribution, collapse = ", ")
  )
}

cran_packages <- dependencies$package[dependencies$source == "CRAN"]
missing_cran <- cran_packages[!vapply(
  cran_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_cran)) {
  repositories <- getOption("repos")
  if (!length(repositories) || is.na(repositories[["CRAN"]]) ||
      repositories[["CRAN"]] == "@CRAN@") {
    repositories <- c(CRAN = "https://cloud.r-project.org")
  }
  tryCatch(
    utils::install.packages(missing_cran, repos = repositories,
                            dependencies = TRUE),
    error = function(error) {
      stop("CRAN dependency installation failed: ", conditionMessage(error))
    }
  )
}

still_missing <- cran_packages[!vapply(
  cran_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(still_missing)) {
  stop("Missing CRAN packages after installation: ",
       paste(still_missing, collapse = ", "))
}

installed_version <- vapply(dependencies$package, function(package) {
  if (package == "R") return(as.character(getRversion()))
  if (!requireNamespace(package, quietly = TRUE)) return("not_installed_yet")
  as.character(utils::packageVersion(package))
}, character(1))

result <- dependencies[c("package", "role", "source")]
result$actual_version <- installed_version
utils::write.csv(result, row.names = FALSE)

if (any(result$source %in% c("CRAN", "R_distribution") &
        result$actual_version == "not_installed_yet")) {
  stop("Dependency verification failed.")
}

# Record software versions without machine-specific installation paths.
record_session <- function(destination) {
  info <- capture.output(sessionInfo())
  info <- info[!grepl("^(BLAS:|LAPACK:|time zone:)", info)]
  writeLines(info, destination)
}

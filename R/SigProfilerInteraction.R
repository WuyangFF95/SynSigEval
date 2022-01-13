# SigProfilerInteraction.R

## Turn this into a test
## ICAReadCatalog("example-SP-signatures.txt")

#' Read a file containing exposures attributed by SigProfiler/Python
#'
#' @param file The name of the file to read.
#'
#' @return The corresponding signature matrix in standard internal
#' representation.
#'
#' @importFrom utils read.table
#'
#' @export

ReadSigProfilerExposure <- function(file) {

  x <- utils::read.table(
    file = file, sep = "\t",
    as.is = TRUE, header = TRUE)

  y <- t(x)
  colnames(y) <- y[1,,drop = F]
  y <- y[-1,,drop = F]

  return(y)
}

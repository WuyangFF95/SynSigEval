# SignatureAnalyzerInteraction.R


#' Standardize SignatureAnalyzer signature names
#'
#' For example, change \code{BI_COMPOSITE_SNV_SBS83_P}
#' to \code{BI_COMPOSITE_SBS83_P}
#'
#' This is necessary because
#' for COMPOSITE signatures we rbind coordinated
#' "SNV", "DNP", and "INDEL" signatures.
#'
#' @param sig.names Vector of signature names
#'
#' @return Vector of signatures names with "_SNV" removed.
#'
#' @export

FixSASigNames <- function(sig.names) {
  return(gsub("_SNV_", "_", sig.names, fixed = TRUE))
}


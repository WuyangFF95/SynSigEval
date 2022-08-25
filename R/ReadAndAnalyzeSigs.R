#' @title Assess how well extracted signatures match input signatures
#'
#' @description We assume that in many cases extraction programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
#'
#' @param reference.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment to only those signatures in \code{reference.sigs}
#'  that were actually represented in the exposures.
#'
#' @return See \code{\link[mSigTools]{TP_FP_FN_avg_sim}}
#'
#' @details Generates output files by calling
#' \code{\link[mSigTools]{TP_FP_FN_avg_sim}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           reference.sigs,
           ground.truth.exposures) {

    # Import extracted sigs and reference sigs --------------------------------
    ex.sigs <- ICAMS::ReadCatalog(extracted.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")
    ref.sigs <- ICAMS::ReadCatalog(reference.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")




    # If file ground.truth.exposures exists, ----------------------------------
    # Remove reference signatures with zero ground-truth exposures
    # Also rename extracted signatures.
    if (file.exists(ground.truth.exposures)) {
      # Rows are signatures, columns are samples.
      exposure <- mSigTools::read_exposure(
        ground.truth.exposures, check.names = F)
      if (!is.null(exposure)) {
        # Remove signatures that are have zero exposures
        # or not present in ground-truth exposure matrix
        exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
        # Make sure we do not have any signatures in exposures that
        # are not in ref.sigs.
        stopifnot(
          setequal(setdiff(exposed.sig.names, colnames(ref.sigs)), c()))
        ref.sigs <- ref.sigs[  , exposed.sig.names]
      }
    }


    # Calculate extraction measures from mSigTools::TP_FP_FN_avg_sim ----------
    retval <- mSigTools::TP_FP_FN_avg_sim(ex.sigs, ref.sigs,
                                          similarity.cutoff = 0.9)



    # Concatenate filtered ex.sigs and ref.sigs as-is --------------------------
    retval$ex.sigs <- ex.sigs
    retval$ref.sigs <- ref.sigs



    return(retval)
  }

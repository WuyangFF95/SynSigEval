#' @title Assess how well extracted signatures match input signatures
#'
#' @description We assume that in many cases extraction programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
#'
#' @param ground.truth.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment to only those signatures in \code{ground.truth.sigs}
#'  that were actually represented in the exposures.
#'
#' @return See \code{\link[ICAMSxtra]{TP_FP_FN_avg_sim}}
#'
#' @details Generates output files by calling
#' \code{\link[ICAMSxtra]{TP_FP_FN_avg_sim}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           ground.truth.sigs,
           ground.truth.exposures) {

    # Import ex sigs, gt sigs -------------------------------------------------
    ex.sigs <- ICAMS::ReadCatalog(extracted.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")
    gt.sigs <- ICAMS::ReadCatalog(ground.truth.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")




    # If file ground.truth.exposures exists, ----------------------------------
    # Remove gt sigs with zero ground-truth exposures
    # Also rename ex sigs.
    if (file.exists(ground.truth.exposures)) {
      # Rows are signatures, columns are samples.
      exposure <- mSigAct::ReadExposure(
        ground.truth.exposures, check.names = F)
      if (!is.null(exposure)) {
        # Remove signatures that are have zero exposures
        # or not present in ground-truth exposure matrix
        exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
        # Make sure we do not have any signatures in exposures that
        # are not in gt.sigs.
        stopifnot(
          setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
        gt.sigs <- gt.sigs[  , exposed.sig.names]
      }
    }


    # Calculate extraction measures from ICAMSxtra::TP_FP_FN_avg_sim ----------
    retval <- ICAMSxtra::TP_FP_FN_avg_sim(ex.sigs, gt.sigs,
                                          similarity.cutoff = 0.9)



    # Concatenate filtered gt.sigs and ex.sigs as-is --------------------------
    retval$ex.sigs <- ex.sigs
    retval$gt.sigs <- gt.sigs



    return(retval)
  }

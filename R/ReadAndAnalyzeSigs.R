#' @title Assess how well extracted signatures match input signatures
#'
#' We assume that in many cases extraction programs will be run
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
#' @return See \code{\link[ICAMSxtra]{MatchSigsAndRelabel}}
#'
#' @details Generates output files by calling
#' \code{\link[ICAMSxtra]{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           ground.truth.sigs,
           ground.truth.exposures) {
    ex.sigs <- ICAMS::ReadCatalog(extracted.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")
    # read.extracted.sigs.fn(extracted.sigs)
    gt.sigs <- ICAMS::ReadCatalog(ground.truth.sigs,
                                  region = "unknown",
                                  catalog.type = "counts.signature")
    # read.ground.truth.sigs.fn(ground.truth.sigs)
    exposure <- ICAMSxtra::ReadExposure(
      ground.truth.exposures,check.names = F)
    # Rows are signatures, columns are samples.

    retval <- ICAMSxtra::MatchSigsAndRelabel(ex.sigs, gt.sigs, exposure)

    ## If input gt.sigs is a ICAMS catalog,
    ## Move all the attributes of gt.sigs to retval::gt.sigs.
    ## Otherwise (e.g. gt.sigs is a matrix), do nothing.
    if(!is.null(attr(gt.sigs, "catalog.type"))){

      catalog.type <- attr(gt.sigs, "catalog.type")
      region       <- attr(gt.sigs, "region")
      ref.genome   <- attr(gt.sigs, "ref.genome")
      abundance    <- attr(gt.sigs, "abundance")

      retval$gt.sigs <-
        ICAMS::as.catalog(retval$gt.sigs,
                          catalog.type = catalog.type,
                          region       = region,
                          ref.genome   = ref.genome,
                          abundance    = abundance)
   }


    ## If input ex.sigs is a ICAMS catalog,
    ## Move all the attributes of ex.sigs to retval::ex.sigs.
    ## Otherwise (e.g. ex.sigs is a matrix), do nothing.
    if(!is.null(attr(ex.sigs, "catalog.type"))){

      catalog.type <- attr(ex.sigs, "catalog.type")
      region       <- attr(ex.sigs, "region")
      ref.genome   <- attr(ex.sigs, "ref.genome")
      abundance    <- attr(ex.sigs, "abundance")

      retval$ex.sigs <-
        ICAMS::as.catalog(retval$ex.sigs,
                          catalog.type = catalog.type,
                          region       = region,
                          ref.genome   = ref.genome,
                          abundance    = abundance)
   }

    return(retval)
  }

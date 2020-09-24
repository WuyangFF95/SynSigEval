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


#' @title Assess how well inferred exposures match input exposures
#'
#' We assume that in many cases attribution programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
#'
#' @param ground.truth.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param inferred.exp.path File containing mutation counts (exposures)
#' of synthetic tumors which are inferred to extracted or input signatures.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment of signature exposures to only those signatures in
#'  \code{ground.truth.sigs} that were actually represented in the exposures.
#'
#' @return A \code{\link{data.frame}} recording:
#'
#' \code{Ground.truth.exposure}: sum of ground truth exposures of
#' all tumors to all ground-truth signatures.
#'
#' \code{Inferred.exposure}: sum of inferred exposures of
#' all tumors to all ground-truth signatures.
#' Here, inferred exposure of a tumor to a ground-truth
#' signature equals to the sum of the exposures of this tumor
#' to all extracted signatures which are most similar to
#' a ground-truth signature.
#' If there is no extracted signature resembling an ground-truth
#' signature, the inferred exposure of this ground-truth
#' signature will be \code{0}.
#'
#' \code{Absolute.difference}: sum of absolute difference between
#' ground-truth exposure and inferred exposure of all tumors
#' to all ground-truth signatures.
#'
#'
#' @details Generates output files by calling
#' \code{\link[ICAMSxtra]{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeExposures <-
  function(extracted.sigs,
           ground.truth.sigs,
           inferred.exp.path,
           ground.truth.exposures) {

    ## Bilaterally matching between ground-truth and extracted signatures
    sigMatch <- ReadAndAnalyzeSigs(extracted.sigs,
                                   ground.truth.sigs,
                                   ground.truth.exposures)


    ## Read in ground-truth and inferred exposures in ICAMS format
    gtExposures <- ICAMSxtra::ReadExposure(
      ground.truth.exposures,check.names = F)
    inferredExposures <- ICAMSxtra::ReadExposure(
      inferred.exp.path,check.names = F)

    ## Names of ground-truth signatures
    gtSigsNames <- colnames(sigMatch$gt.sigs)

    ## Initialize an empty data.frame for exposure difference
    exposureDiff <- data.frame(matrix(0,nrow = length(gtSigsNames),ncol = 4))
    rownames(exposureDiff) <- gtSigsNames
    colnames(exposureDiff) <- c("Ground.truth.exposure", ## Sum of all tumor's ground-truth exposure to gtSigsName
                                "Inferred.exposure", ## Sum of all tumor's inferred exposure to gtSigsName
                                "Absolute.difference", ## Sum of absolute difference of two exposure values for each tumor
                                "Manhattan.distance") ## L1-difference betwen ground-truth exposure and inferred exposure.
    ## = Absolute.difference/Ground.truth.Exposure

    ## For each of the ground-truth signature, calculate the absolute difference
    ## between its input (ground-truth) exposure and its inferred exposure.
    ## Inferred exposure of a input signature equals to the sum of
    ## exposures of all extracted signatures which matches to
    ## this input signature.
    for (gtSigName in gtSigsNames) {
      matchedExtrSigIndex <- which(sigMatch$match1[,1] == gtSigName)

      if (length(matchedExtrSigIndex) > 0)
        ## 1 or more extracted signatures match to gtSigName in match1
        matchedExtrSigName <- rownames(sigMatch$match1)[matchedExtrSigIndex]
      else ## No extracted signatures match to gtSigName
        matchedExtrSigName <- NULL

      for (index in 1:ncol(inferredExposures)) { ## index refers to which tumor we are scrutinizing
        ## Each cycle traverses one tumor, and calculate the absolute difference
        ## between its inferred exposures and ground-truth exposures.
        gtExposureOneTumor <- gtExposures[gtSigName,index]
        inferredExposureOneTumor <- ifelse(length(matchedExtrSigIndex) > 0,
                                       yes = sum(inferredExposures[matchedExtrSigName,index]),
                                       no = 0)
        exposureDiff[gtSigName,1] <- exposureDiff[gtSigName,1] + gtExposureOneTumor
        exposureDiff[gtSigName,2] <- exposureDiff[gtSigName,2] + inferredExposureOneTumor
        exposureDiff[gtSigName,3] <- exposureDiff[gtSigName,3] +
          abs(gtExposureOneTumor - inferredExposureOneTumor)
      }
    }

    ## Only after the cycle, the exposureDiff[,c(1,3)] has been fixed.
    ## The Manhattan distance should normalize against the sum of exposures
    ## of one signature, not sum of exposures of all signatures.
    ## This can prevent the underestimation of discrepancy between inferred
    ## exposures and ground-truth exposures.
    exposureDiff[,4] <- exposureDiff[,3] / exposureDiff[,1]

    return(exposureDiff)
  }


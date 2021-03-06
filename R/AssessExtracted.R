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

    ## aggregated: a data.frame for aggregated exposure measures.
    {
      ## Initialize an empty data.frame for aggregated exposure difference
      aggregated <- data.frame(matrix(0,nrow = length(gtSigsNames),ncol = 4))
      rownames(aggregated) <- gtSigsNames
      colnames(aggregated) <- c("Ground.truth.exposure", ## Sum of all tumor's ground-truth exposure to gtSigName
                                "Inferred.exposure", ## Sum of all tumor's inferred exposure to gtSigName
                                "Aggregated.Manhattan.distance", ## Sum of absolute difference (L1-norm) of two exposures to gtSigName for each tumor
                                "Scaled.Aggregated.Manhattan.distance") ## Aggregated Manhattan Distance, scaled by the sum of
      ## ground-truth exposures to gtSigName in all tumors

      ## For each of the ground-truth signature, calculate the absolute difference
      ## between its input (ground-truth) exposure and its inferred exposure.
      ## Inferred exposure of a ground-truth signature equals to
      ## exposures of the extracted signature most similar to
      ## this ground-truth signature, and their cosine simlarity needs to
      ## be > 0.90.
      for (gtSigName in gtSigsNames) {
        matchedExtrSigIndex <- intersect(
          which(rownames(sigMatch$match2) == gtSigName),
          which(sigMatch$match2[,2] >= 0.9)
        )

        if (length(matchedExtrSigIndex) > 0)
          ## 1 extracted signature match most similar to gtSigName in match2
          matchedExtrSigName <- sigMatch$match2[matchedExtrSigIndex,1]
        else ## No extracted signatures match to gtSigName
          matchedExtrSigName <- NULL

        for (nTumor in 1:ncol(inferredExposures)) { ## nTumor refers to which tumor we are scrutinizing
          ## Each cycle traverses one tumor, and calculate the absolute difference
          ## between its inferred exposures and ground-truth exposures.
          gtExposureOneTumor <- gtExposures[gtSigName,nTumor]
          inferredExposureOneTumor <- ifelse(length(matchedExtrSigIndex) > 0,
                                             yes = sum(inferredExposures[matchedExtrSigName,nTumor]),
                                             no = 0)
          aggregated[gtSigName,1] <- aggregated[gtSigName,1] + gtExposureOneTumor
          aggregated[gtSigName,2] <- aggregated[gtSigName,2] + inferredExposureOneTumor
          aggregated[gtSigName,3] <- aggregated[gtSigName,3] +
            abs(gtExposureOneTumor - inferredExposureOneTumor)
        }
      }

      ## Only after the cycle, the aggregated[,c(1,3)] has been fixed.
      ## The Scaled aggregated Manhattan distance should normalize against
      ## the sum of exposures of one signature, not sum of exposures of all signatures.
      ## This can prevent the underestimation of discrepancy between inferred
      ## exposures and ground-truth exposures.
      aggregated[,4] <- aggregated[,3] / aggregated[,1]
    }

    ## separated: Manhattan distance and scaled Manhattan distance
    ## for each individual tumor:
    {
      separated = list()

      ## Initialize an empty data.frame for aggregated exposure difference
      for(gtSigName in gtSigsNames){
        separated[[gtSigName]] <- data.frame(
          "Ground.truth.exposure" = numeric(), ## Sum of all tumor's ground-truth exposure to gtSigsName
          "Inferred.exposure" = numeric(), ## Sum of all tumor's inferred exposure to gtSigsName
          "Manhattan.distance" = numeric(), ## Sum of absolute difference (L1-norm) of two exposure values for each tumor
          "Scaled.Manhattan.distance" = numeric()) ## Manhattan Distance, scaled by the ground-truth exposure of each tumor to gtSigName
      }


      ## exposures of ground-truth signature in all tumors

      for (gtSigName in gtSigsNames) {
        matchedExtrSigIndex <- intersect(
          which(rownames(sigMatch$match2) == gtSigName),
          which(sigMatch$match2[,2] >= 0.9)
        )

        if (length(matchedExtrSigIndex) > 0)
          ## 1 extracted signature match most similar to gtSigName in match2
          matchedExtrSigName <- sigMatch$match2[matchedExtrSigIndex,1]
        else ## No extracted signatures match to gtSigName
          matchedExtrSigName <- NULL

        for (nTumor in 1:ncol(inferredExposures)) { ## nTumor refers to which tumor we are scrutinizing
          ## Each cycle traverses one tumor, and calculate the absolute difference
          ## between its inferred exposures and ground-truth exposures.
          gtExposureOneTumor <- gtExposures[gtSigName,nTumor]
          inferredExposureOneTumor <- ifelse(length(matchedExtrSigIndex) > 0,
                                             yes = sum(inferredExposures[matchedExtrSigName,nTumor]),
                                             no = 0)
          values4OneTumor <- matrix(nrow = 1, ncol = 4)
          values4OneTumor[1,1] <- gtExposureOneTumor
          values4OneTumor[1,2] <- inferredExposureOneTumor
          values4OneTumor[1,3] <- abs(gtExposureOneTumor - inferredExposureOneTumor)
          values4OneTumor[1,4] <- values4OneTumor[1,3] / values4OneTumor[1,1]

          colnames(values4OneTumor) <- c(
            "Ground.truth.exposure",
            "Inferred.exposure",
            "Manhattan.distance",
            "Scaled.Manhattan.distance")

          separated[[gtSigName]] <- rbind(separated[[gtSigName]],values4OneTumor)
        }

        rownames(separated[[gtSigName]]) <- colnames(inferredExposures)
      }

    }

    exposureDiff <- list(
      aggregated = aggregated,
      separated = separated
    )

    return(exposureDiff)
  }


#' Run \code{MatchSigs2Directions}, then
#' plot its results and write them as .csv files.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general).
#
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure "Ground truth" exposures used generate the
#'   synthetic data from which \code{ex.sigs} were extracted.
#'
#' @return A list with the elements \code{averCosSim}, \code{match1},
#' \code{match2} as for \code{SigSetSimilarity}, with \code{match1}
#'  being matches for the the extracted signatures (\code{ex.sigs})
#'  and \code{match2} being the
#'  matches for the ground truth signatures (\code{gt.sigs}). The return list
#'  also echos the input arguments \code{ex.sigs} and \code{gt.sigs}.
#
#' @export
#' @family signature matching functions

MatchSigsAndRelabel <-
  function(ex.sigs, gt.sigs, exposure) {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }

    # IMPORTANT Remove signatures that are not present in
    # the exposure from which the synthetic data were
    # generated
    exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
    # Make sure we do not have an signatures in exposures that
    # are not in gt.sigs.
    stopifnot(
      setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
    gt.sigs <- gt.sigs[  , exposed.sig.names]

    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    ## Software-reported signatures with a best cosine similarity lower than 0.90
    ## is not considered an "extracted signature", it will rather be regarded as
    ## an artefact.
    true.match1 <- sim$match1
    true.match1 <- true.match1[true.match1$sim >= 0.9,]
    true.match2 <- sim$match2
    true.match2 <- true.match2[true.match2$sim >= 0.9,]

    sim$extracted.with.no.best.match <-
      setdiff(colnames(ex.sigs), true.match2$to)

    sim$ground.truth.with.no.best.match <-
      setdiff(colnames(gt.sigs), true.match1$to)
    # TODO(Steve) Review documentation / explanation. Note that
    # e.g. SBS29 might have a best match (BI_COMPOSITE_SBS18_P)
    # but no BI signatures has SBS29 as its best match
    #

    # TODO(Steve) Document the complexity below; mostly it deals
    # with setting up plotting that is easy(?) to interpret.
    labels <- character(ncol(ex.sigs))
    names(labels) <- colnames(ex.sigs)
    nums <- NumFromId(sim$match1$to)
    reordered.ex <- colnames(ex.sigs)[order(nums)]
    ex.sigs.x <- ex.sigs[ , order(nums),drop = FALSE]
    bestmatch.id <- sim$match1[reordered.ex, "to"]
    bestmatch.sim <- sim$match1[reordered.ex, "sim"]
    bestmatch.sim <- round(bestmatch.sim, digits=4)
    init.labels <-
      paste0(reordered.ex, " (", bestmatch.id, " ", bestmatch.sim, ")")
    names(init.labels) <- reordered.ex
    laggards <- setdiff(rownames(sim$match2), bestmatch.id)
    # Falling back to a loop here:
    for (lag in laggards) {
      my.ex.id  <- sim$match2[lag, "to"]
      my.ex.sim <- round(sim$match2[lag, "sim"], digits = 4)
      init.labels[my.ex.id] <-
        paste0(init.labels[my.ex.id],
               " (", lag, " ", my.ex.sim, ")")
    }
    colnames(ex.sigs.x) <- init.labels

    sim$ex.sigs <- ICAMS::as.catalog(
      ex.sigs.x,
      region = "genome",
      catalog.type = "counts.signature")
    sim$gt.sigs <- ICAMS::as.catalog(
      gt.sigs,
      region = "genome",
      catalog.type = "counts.signature")

    ## Calculate cosine similarity between all extracted signatures,
    ## and each of the ground-truth signatures.
    # E.g. First calculate the cosine similarity between ground-truth SBS5 and all
    # extracted signatures most similar to SBS5; then calculate the cosine similarity
    # between SBS1 and all extracted signatures most similar to SBS1.
    if(TRUE){ ## debug
      sim$cosSim <- list()

      gtSigNames <- rownames(sim$match2)
      exSigNames <- rownames(sim$match1)

      for(gtSigName in gtSigNames){
        ## sim$match1 denotes the ground-truth signature each extracted
        ## signature is most similar to, and their pairwise cosine similarity.
        tmp <- sim$match1
        ## In sim$match1, Find all extracted signatures similar to
        ## ground-truth signature "gtSigName".
        values <- tmp[which(tmp[,1] == gtSigName),2]

        if(is.nan(mean(values))) {
          ## None of the extracted signatures were most similar to "gtSigName"
          ## In this way, we go to sim$match2 instead, and find out
          ## the extracted signature "gtSigName" is most similar to.
          ##
          ## This is because there are some cases scenarios, two ground-truth
          ## signatures have been blended into 1 in extraction results.
          ##
          ## In this way, we can still study how an extracted signature is
          ## similar to multiple ground-truth signatures, and we can also
          ## compare the performance of different software packages in a
          ## more reasonable way.
          tmp <- sim$match2
          value <- tmp[gtSigName,2]
          sim$cosSim[[gtSigName]] <- value
        } else{
          ## There are some extracted signatures most similar to "gtSigName"
          ## Average cosine similarity of all extracted "gtSigName" signature
          sim$cosSim[[gtSigName]] <- mean(values)
        }
      }
    }

    invisible(sim)
  }

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
#' @return See \code{\link{MatchSigsAndRelabel}}
#'
#' @details Generates output files by calling
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           ground.truth.sigs,
           ground.truth.exposures) {
    ex.sigs <- ICAMS::ReadCatalog(extracted.sigs,
                                  region = "genome",
                                  catalog.type = "counts.signature")
    # read.extracted.sigs.fn(extracted.sigs)
    gt.sigs <- ICAMS::ReadCatalog(ground.truth.sigs, region = "genome",
                                  catalog.type = "counts.signature")
    # read.ground.truth.sigs.fn(ground.truth.sigs)
    exposure <- ReadExposure(ground.truth.exposures)
    # Rows are signatures, columns are samples.

    return(
      MatchSigsAndRelabel(ex.sigs, gt.sigs, exposure))
  }


#' @title Assess how well attributed exposures match input exposures
#'
#' We assume that in many cases attribution programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
#'
#' @param ground.truth.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param attributed.exp.path File containing mutation counts (exposures)
#' of synthetic tumors which are attributed to extracted or input signatures.
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
#' \code{Attributed.exposure}: sum of attributed exposures of
#' all tumors to all ground-truth signatures.
#' Here, attributed exposure of a tumor to a ground-truth
#' signature equals to the sum of the exposures of this tumor
#' to all extracted signatures which are most similar to
#' a ground-truth signature.
#' If there is no extracted signature resembling an ground-truth
#' signature, the attributed exposure of this ground-truth
#' signature will be \code{0}.
#'
#' \code{Absolute.difference}: sum of absolute difference between
#' ground-truth exposure and attributed exposure of all tumors
#' to all ground-truth signatures.
#'
#'
#' @details Generates output files by calling
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeExposures <-
  function(extracted.sigs,
           ground.truth.sigs,
           attributed.exp.path,
           ground.truth.exposures) {

    ## Bilaterally matching between ground-truth and extracted signatures
    sigMatch <- ReadAndAnalyzeSigs(extracted.sigs,
                                   ground.truth.sigs,
                                   ground.truth.exposures)


    ## Read in ground-truth and attributed exposures in ICAMS format
    gtExposures <- ReadExposure(ground.truth.exposures)
    attrExposures <- ReadExposure(attributed.exp.path)

    ## Names of ground-truth signatures
    gtSigsNames <- colnames(sigMatch$gt.sigs)

    ## Initialize an empty data.frame for exposure difference
    exposureDiff <- data.frame(matrix(0,nrow = length(gtSigsNames),ncol = 4))
    rownames(exposureDiff) <- gtSigsNames
    colnames(exposureDiff) <- c("Ground.truth.exposure", ## Sum of all tumor's ground-truth exposure to gtSigsName
                                "Attributed.exposure", ## Sum of all tumor's attributed exposure to gtSigsName
                                "Absolute.difference", ## Sum of absolute difference of two exposure values for each tumor
                                "Manhattan.distance") ## L1-difference betwen ground-truth exposure and attributed exposure.
    ## = Absolute.difference/Ground.truth.Exposure

    ## For each of the ground-truth signature, calculate the absolute difference
    ## between its input (ground-truth) exposure and its attributed exposure.
    ## Attributed exposure of a input signature equals to the sum of
    ## exposures of all extracted signatures which matches to
    ## this input signature.
    for (gtSigName in gtSigsNames) {
      matchedExtrSigIndex <- which(sigMatch$match1[,1] == gtSigName)

      if (length(matchedExtrSigIndex) > 0)
        ## 1 or more extracted signatures match to gtSigName in match1
        matchedExtrSigName <- rownames(sigMatch$match1)[matchedExtrSigIndex]
      else ## No extracted signatures match to gtSigName
        matchedExtrSigName <- NULL

      for (index in 1:ncol(attrExposures)) { ## index refers to which tumor we are scrutinizing
        ## Each cycle traverses one tumor, and calculate the absolute difference
        ## between its attributed exposures and ground-truth exposures.
        gtExposureOneTumor <- gtExposures[gtSigName,index]
        attrExposureOneTumor <- ifelse(length(matchedExtrSigIndex) > 0,
                                       yes = sum(attrExposures[matchedExtrSigName,index]),
                                       no = 0)
        exposureDiff[gtSigName,1] <- exposureDiff[gtSigName,1] + gtExposureOneTumor
        exposureDiff[gtSigName,2] <- exposureDiff[gtSigName,2] + attrExposureOneTumor
        exposureDiff[gtSigName,3] <- exposureDiff[gtSigName,3] +
          abs(gtExposureOneTumor - attrExposureOneTumor)
      }
    }

    ## Only after the cycle, the exposureDiff[,c(1,3)] has been fixed.
    ## The Manhattan distance should normalize against the sum of exposures
    ## of one signature, not sum of exposures of all signatures.
    ## This can prevent the underestimation of discrepancy between attributed
    ## exposures and ground-truth exposures.
    exposureDiff[,4] <- exposureDiff[,3] / exposureDiff[,1]

    return(exposureDiff)
  }


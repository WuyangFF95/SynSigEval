#' Run \code{\link[SynSigGen]{MatchSigs2Directions}}, then
#' plot its results and write them as .csv files.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general).
#
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure If \code{NULL}, then match
#'   \code{ex.sigs} against all signatures in \code{gt.sigs}.
#'   Otherwise this should be ground-truth exposures used generate the
#'   synthetic spectra from which \code{ex.sigs} were extracted.
#'   In this case we do not
#'   match to ground-truth signatures to that were not in the ground
#'   truth exposure.
#'
#' @return A list with the elements \code{averCosSim}, \code{match1},
#' \code{match2} as for \code{\link[SynSigGen]{MatchSigs2Directions}},
#' with \code{match1} being matches for the the extracted signatures
#' (\code{ex.sigs}) and \code{match2} being the matches for
#' the ground truth signatures (\code{gt.sigs}). The return list
#' also echos the input arguments \code{ex.sigs} and \code{gt.sigs}.
#'
#' @export
#'
#' @family signature matching functions

MatchSigsAndRelabel <-
  function(ex.sigs, gt.sigs, exposure = NULL) {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }

    if (!is.null(exposure)) {
      # Remove signatures that are not present in
      # the exposure from which the synthetic data were
      # generated
      exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
      # Make sure we do not have an signatures in exposures that
      # are not in gt.sigs.
      stopifnot(
        setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
      gt.sigs <- gt.sigs[  , exposed.sig.names]
    }

    sim <- SynSigGen::MatchSigs2Directions(ex.sigs, gt.sigs)

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
    nums <- SynSigGen::NumFromId(sim$match1$to)
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

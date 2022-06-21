

#' @title Assess how well inferred exposures match input exposures
#'
#' @description We assume that in many cases attribution programs will be run
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
#' \code{\link[ICAMSxtra]{TP_FP_FN_avg_sim}}
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
    gtExposures <- mSigAct::ReadExposure(
      ground.truth.exposures,check.names = F)
    inferredExposures <- mSigAct::ReadExposure(
      inferred.exp.path,check.names = F)

    ## Names of ground-truth signatures
    gtSigsNames <- colnames(sigMatch$gt.sigs)
    ## Names of spectra
    spectrumNames <- colnames(inferredExposures)

    ## I. Manhattan: Manhattan distance and scaled Manhattan distance
    ## for each individual spectrum:
    {
      Manhattan = list()
      ## Initialize an empty data.frame for Manhattan distance for each tumor
      for(spectrumName in spectrumNames){
        Manhattan[[spectrumName]] <- data.frame(
          "Ground.truth.exposure" = numeric(), ## ground-truth exposure to one signature in tumor "spectrumName"
          "Inferred.exposure" = numeric(), ## inferred exposure
          "Manhattan.distance" = numeric(), ## Absolute difference (L1-norm) of two exposure valuess
          "Scaled.Manhattan.distance" = numeric()) ## Manhattan Distance, divided by the ground-truth exposure of each tumor to all signatures
      }
      for (spectrumName in spectrumNames) {
        for (gtSigName in gtSigsNames) {
          matchedExtrSigIndex <- intersect(
            which(rownames(sigMatch$match2) == gtSigName),
            which(sigMatch$match2[,2] >= 0.9)
          )
          if (length(matchedExtrSigIndex) > 0)
            ## At least 1 extracted signature in match2 matches to gtSigName
            ## i.e. having highest cosine similarity.
            matchedExtrSigNames <- sigMatch$match2[matchedExtrSigIndex,1]
          else ## No extracted signatures match to gtSigName
            matchedExtrSigNames <- NULL
          gtExposureOneSig <- gtExposures[gtSigName,spectrumName]
          inferredExposureOneSig <- ifelse(length(matchedExtrSigIndex) > 0,
                                           yes = sum(inferredExposures[matchedExtrSigNames,spectrumName]),
                                           no = 0)
          values4OneSig <- matrix(nrow = 1, ncol = 4)
          values4OneSig[1,1] <- gtExposureOneSig
          values4OneSig[1,2] <- inferredExposureOneSig
          values4OneSig[1,3] <- abs(gtExposureOneSig - inferredExposureOneSig)
          ## Initializing 4-th column for Scaled Manhattan distance.
          values4OneSig[1,4] <- NA
          colnames(values4OneSig) <- c(
            "Ground.truth.exposure",
            "Inferred.exposure",
            "Manhattan.distance",
            "Manhattan.distance.Scaled.by.Spectrum")
          ## Scaled Manhattan distance at 4-th column.
          ## Manhattan distance of signature X / sum of ground-truth exposures to all signatures
          ## in tumor "spectrumName".
          Manhattan[[spectrumName]] <- rbind(Manhattan[[spectrumName]],values4OneSig)
        }
        Manhattan[[spectrumName]][,4] <- Manhattan[[spectrumName]][,3] / colSums(Manhattan[[spectrumName]])[1]
        rownames(Manhattan[[spectrumName]]) <- gtSigsNames
      }
    }

    ## II. SumOfManhattan: a data.frame for sum of Manhattan distance for each mutational spectrum.
    {
      ## Initialize an empty data.frame for sum of Manhattan distance
      SumOfManhattan <- data.frame(matrix(NA,nrow = length(spectrumNames),ncol = 4))
      rownames(SumOfManhattan) <- spectrumNames
      colnames(SumOfManhattan) <- c("Ground.truth.exposure", ## Sum of ground-truth exposure to all signatures' in one tumor
                                "Inferred.exposure", ##Sum of inferred exposure to all signatures' in one tumor
                                "Manhattan.distance", ## Sum of absolute difference (L1-norm) of ground-truth and inferred signatures to each of the signatures
                                "Manhattan.distance.Scaled.by.Spectrum") ## Manhattan distance, divided by total ground-truth mutations in one spectrum
      for (spectrumName in spectrumNames) {
        SumOfManhattan[spectrumName,] <- colSums(Manhattan[[spectrumName]])
      }
    }

    ## III. F1 score related measures:
	##
    ## PPV: Positive Predictive Value for one tumor.
    ## TPR: True Positive Rate for one tumor.
    ## F1 score: Harmonic mean of PPV and TPR.
    {
      PPV <- TPR <- F1 <- numeric(length = ncol(inferredExposures))
      names(PPV) <- names(TPR) <- names(F1) <- colnames(inferredExposures)
      for (spectrumName in spectrumNames) {
        ## Ground-truth positive signatures
        P <- names(which(gtExposures[,spectrumName] > 0))
        ## Discovered signatures
        D <- names(which(inferredExposures[,spectrumName] > 0))

        ## Calculate FP, FN and TP.
        FP <- setdiff(D,P)
        FN <- setdiff(P,D)
        TP <- setdiff(P,FN)

        ## PPV (precision) = #TP / (#TP + # FP)
        PPV[spectrumName] <- length(TP) / (length(TP) + length(FP))
        ## TPR (recall) = #TP / #P = #TP / (#TP + #FN)
        TPR[spectrumName] <- length(TP) / length(P)
        ## F1 score, harmonic mean of PPV (precision) and TPR (recall),
        ## also equals to #TP / [#TP + (#FP + #FN) / 2].
        F1[spectrumName] <-
          length(TP) / (length(TP) + (length(FP) + length(FN))/2 )
      }
    }

    ## Final return value - exposureDiff.
    exposureDiff <- list(
      Manhattan = Manhattan,
      SumOfManhattan = SumOfManhattan,
	  F1.measures = data.frame(
	   PPV = PPV,
	   TPR = TPR,
	   F1 = F1)
    )
    return(exposureDiff)
  }


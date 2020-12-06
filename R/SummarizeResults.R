CopyWithChecks <- function(from, to.dir, overwrite = FALSE) {
  if (!file.exists(from)) {
    warning("Cannot find", from, "\n\nSkipping\n\n")
  } else {
    copy.res <- file.copy(
      from = from, to = to.dir, overwrite = overwrite)
    if (!copy.res)
      cat("Copy from", from, "to directory", to.dir, "failed\n\n")
  }
}


#' Assess/evaluate one result sub-folder from a computational approach
#'
#' Note: For summarizing SigProExtractor or SignatureAnalyzer,
#' users should use SigProExtractor(SigProfiler-Python) v0.0.5.43+
#' and SignatureAnalyzer 2018-Apr-18.
#'
#' @param run.dir Lowest level path to result of a run. That is,
#' \code{top.dir}/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.1/. or
#' \code{top.dir}/sa.sa.96/Attr/YAPSA.results/seed.691/
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{run.dir}
#' which stores the software output.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' Usually, it refers to \code{sub.dir}, i.e. \code{run.dir}/../../../
#'
#' @param extracted.sigs.path Path to extracted sigs file, e.g.
#' \code{<run.dir>/SBS96/Selected_Solution/De_Novo_Solution/signatures.PCAWG.format.csv}.
#'
#' @param inferred.exp.path Path to inferred exposures file.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @param summary.folder.name The name of the folder containing summary results.
#' Usually, it equals to "summary".
#'
#' @importFrom utils write.csv capture.output sessionInfo
#'
#' @keywords internal

SummarizeSigOneSubdir <-
  function(run.dir,
           ground.truth.exposure.dir,
           extracted.sigs.path,
           inferred.exp.path = NULL,
           # TODO(Steve): copy this to the summary and do analysis on how much
           # extracted signature contributes to exposures.
           overwrite = FALSE,
           summary.folder.name = "summary") {

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(run.dir, "/", summary.folder.name)

    ## Analyze signature extraction similarity
    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
        ground.truth.exposures =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv")
      )

    ## Workout for not breaking the downstream summarizing code
    sigAnalysis$cosSim <- sigAnalysis$gt.mean.cos.sim
    sigAnalysis$gt.mean.cos.sim <- NULL

    if (dir.exists(outputPath)) {
      if (!overwrite) stop(outputPath, " already exists")
    }
    suppressWarnings(dir.create(outputPath))

    # Copies ground.truth exposures from second.level.dir
    # to outputPath == run.dir/<summary.folder.name>.
    CopyWithChecks(
      from = paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"),
      to.dir = outputPath,
      overwrite = TRUE)

    # Writes bi-directional matching and cos.sim calculation
    write.csv(sigAnalysis$match1,
              file = paste(outputPath,"match1.csv",sep = "/"))
    write.csv(sigAnalysis$match2,
              file = paste(outputPath,"match2.csv",sep = "/"))

    # Writes ground truth and extracted signatures
    # write.cat.fn(
    ICAMS::WriteCatalog(
      sigAnalysis$gt.sigs,
      paste(outputPath,"ground.truth.sigs.csv",sep = "/"),
    )
    # write.cat.fn
    ICAMS::WriteCatalog(
      sigAnalysis$ex.sigs,
      paste(outputPath,"extracted.sigs.csv",sep = "/"))

    # Dumps other outputs into "other.results.txt"
    capture.output(
      cat("Average cosine similarity\n"),
      sigAnalysis$averCosSim,
      cat("Average cosine similarity to each ground-truth signature\n"),
      sigAnalysis$cosSim,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis$ex.sigs),
      cat("\nsigAnalysis$extracted.with.no.best.match\n"),
      sigAnalysis$extracted.with.no.best.match,
      cat("\nsigAnalysis$ground.truth.with.no.best.match\n"),
      sigAnalysis$ground.truth.with.no.best.match,
      file = paste0(outputPath,"/other.results.txt"))

    ## Plot signatures as "counts.signatures" typed catalog
    ## Currently, ICAMS cannot plot COMPOSITE catalog.
    ## TODO(Wuyang): To add a ICAMS:::PlotCatalog.COMPOSITECatalog function

    # Output ground-truth sigs to a PDF file
    if("COMPOSITECatalog" %in% class(sigAnalysis$gt.sigs) == FALSE){
      ICAMS::PlotCatalogToPdf(sigAnalysis$gt.sigs,
                              paste0(outputPath,"/ground.truth.sigs.pdf"))
    }
    if("COMPOSITECatalog" %in% class(sigAnalysis$ex.sigs) == FALSE){
      # Output extracted sigs to a PDF file
      ICAMS::PlotCatalogToPdf(sigAnalysis$ex.sigs,
                              paste0(outputPath,"/extracted.sigs.pdf"))
    }

    ## Analyze signature attribution (a.k.a. exposure inference)
    # To be compatible with PCAWG project which only studies
    # signature extraction not signature attribution,
    # errors will not be thrown if !is.null(inferred.exp.path) == F.
    # Here we shouldn't use "exists("attritbuted.exp.path")" because
    # inferred.exp.path is defaulted to be NULL, but is always defined
    # therefore exists.
    if(!is.null(inferred.exp.path)) {

      if(file.exists(inferred.exp.path)) {
        exposureDiff <- ReadAndAnalyzeExposures(
          extracted.sigs = extracted.sigs.path,
          ground.truth.sigs =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
          inferred.exp.path = inferred.exp.path,
          ground.truth.exposures =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"))

        # Write results of exposure inference measures,
        # in aggregated format
        write.csv(exposureDiff$aggregated,
                  file = paste0(outputPath,"/aggregatedExposureDifference.csv"),
                  quote = T)
        # Write results of exposure inference measures,
        # in aggregated format for each tumor and each
        # ground-truth signature.
        for(gtSigName in names(exposureDiff$separated)){
          write.csv(exposureDiff$separated[[gtSigName]],
                    file = paste0(outputPath,"/separatedExposureDifference.",gtSigName,".csv"),
                    quote = T)
        }

        # Copy inferred exposures to summary folder.
        CopyWithChecks(inferred.exp.path,
                       paste0(outputPath,"/inferred.exposures.csv"),
                       overwrite = overwrite)
      }
      else {
        warning("Cannot find", inferred.exp.path, "\n\nSkipping\n\n")
      }
    }

    ## Log of system time and session info
    capture.output(Sys.time(), sessionInfo(),
                   file = paste0(outputPath,"/log.txt"))

    ## Save Signature extraction summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    save(sigAnalysis,
         file = paste0(outputPath,"/sigAnalysis.RDa"))
    ## Save exposure inference summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    if(file.exists(inferred.exp.path)){
      save(exposureDiff,
           file = paste0(outputPath,"/exposureDiff.RDa"))
    }

    invisible(sigAnalysis) # So we have something to check in tests
  }


#' Assess/evaluate multiple summarized runs for one dataset from one computational approach.
#'
#' Summarize results from each computational approach in \code{tool.dir}/\code{run.names}
#' (generated by running a computational approach),
#' combine them into \code{tool.dir}.
#'
#' @param datasetName Name of the dataset. (e.g. "S.0.1.Rsq.0.1").
#' Usually, it is has the same name as \code{basename(top.dir)}.
#'
#' @param toolName Name of computational approach. (e.g. "SigProExtractor")
#'
#' @param tool.dir Fourth level path from the \code{top.dir}. Expected to have
#' multiple runs with different names (e.g. "seed.1")
#' That is,
#' \code{top.dir}/sp.sp/ExtrAttr/sa.results/. or
#' \code{top.dir}/sa.sa.96/Attr/deconstructSigs.results/
#'
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{tool.dir}
#' which stores the software output.
#'
#' @param run.names A character vector records the list of \code{run.dir},
#' or fifth level directories from the dataset top-level folder.
#' E.g., c("seed.1","seed.691")
#'
#' @return A list contain values of measures measures in multiple runs: \itemize{
#' \item $averCosSim Cosine similarity
#' \item $truePos True Positives(TP): Ground-truth signatures which are active in
#' the spectra, and extracted.
#' \item $falseNeg False Negatives(FN): Ground-truth signatures not extracted.
#' \item $falsePos False Positives(FP): Signatures wrongly extracted, not resembling
#' any ground-truth signatures.
#' \item $TPR True positive rate (TPR, Sensitivity): TP / (TP + FN)
#' \item $PPV Positive predictive value (PPV): TP / (FP + TP)
#' \item $cosSim Average cosine similarity to each of the ground-truth signatures.
#' \item $AggManhattanDist Scaled Manhattan distance between ground-truth and inferred
#' exposures to each of the ground-truth signatures.
#' }
#' This list also contains \code{mean} and \code{sd}, and other
#' statistics of these measures in \itemize{
#' \item $fivenum
#' \item $fivenumMD
#' \item $meanSD
#' \item $meanSDMD
#' }
#'
#' @details Also writes multiple files into folder \code{tool.dir}.
#'
#' @importFrom utils write.csv capture.output sessionInfo
#'
#' @importFrom rlang .data
#'
#' @export
SummarizeMultiRuns <-
  function(datasetName,
           toolName,
           tool.dir,
           run.names){

    ## Indexes for signature extraction in multiple runs
    indexes <- c("averCosSim","falseNeg","falsePos",
                 "truePos","TPR","PPV")
    for(index in indexes) assign(index,numeric(0))
    cosSim <- list()


    for(runName in run.names){
      ## Load directories
      runDir <- paste0(tool.dir,"/",runName)
      summaryDir <- paste0(runDir,"/summary")
      sigAnalysisFile <- paste0(summaryDir,"/sigAnalysis.RDa")
      ## Add sigAnalysis <- NULL to please the Rcheck
      sigAnalysis <- NULL
      load(file = sigAnalysisFile)

      ## Names of ground-truth signatures
      ## Useful-in calculating the average of
      ## one-signature cosine similarity.
      gtSigNames <- rownames(sigAnalysis$match2)

      ## Concatenate average cosine similarity
      averCosSim <- c(averCosSim,sigAnalysis$averCosSim)

      ## Concatenate true positive, true negative and false positive signatures.
      falseNegNames <- sigAnalysis$ground.truth.with.no.best.match
      falsePosNames <- sigAnalysis$extracted.with.no.best.match
      truePosNames <- setdiff(gtSigNames,falseNegNames)
      falseNeg <- c(falseNeg,length(falseNegNames))
      falsePos <- c(falsePos,length(falsePosNames))
      truePos <- c(truePos, length(truePosNames))

      ## Concatenate TPR (True positive rate) and PPV (Positive predictive value)
      currentTPR <- length(truePosNames) / length(gtSigNames)
      currentPPV <- length(truePosNames) / (length(truePosNames) + length(falsePosNames))
      TPR <- c(TPR, currentTPR)
      PPV <- c(PPV, currentPPV)

      ## Concatenating one-signature cosine similarity
      for(gtSigName in gtSigNames) {
        if(is.null(cosSim[[gtSigName]]))
          cosSim[[gtSigName]] <- numeric(0)
        cosSim[[gtSigName]] <- c(cosSim[[gtSigName]],sigAnalysis$cosSim[[gtSigName]])
      }
    }

    ## Make every vector named by run names (e.g. "seed.1")
    names(averCosSim) <- run.names
    names(falseNeg) <- run.names
    names(falsePos) <- run.names
    names(truePos) <- run.names
    names(TPR) <- run.names
    names(PPV) <- run.names
    for(gtSigName in gtSigNames)
      names(cosSim[[gtSigName]]) <- run.names


    multiRun <- list()
    ## Save name of the computational approach and the dataset.
    multiRun$datasetName <- datasetName
    multiRun$toolName <- toolName
    ## Save extraction indexes on multiple runs
    multiRun$averCosSim <- averCosSim
    multiRun$falseNeg <- falseNeg
    multiRun$falsePos <- falsePos
    multiRun$truePos <- truePos
    multiRun$TPR <- TPR
    multiRun$PPV <- PPV
    ## Save one-signature cosine similarity on multiple runs
    multiRun$cosSim <- cosSim



    ## Calculate mean and SD for indexes of signature extraction
    multiRun$meanSD <- matrix(nrow = 6, ncol = 2)
    indexes <- c("averCosSim","falseNeg","falsePos",
                 "truePos","TPR","PPV")
    rownames(multiRun$meanSD) <- indexes
    colnames(multiRun$meanSD) <- c("mean","stdev")
    for(index in indexes){
      currentMean <- mean(multiRun[[index]])
      currentStdev <- stats::sd(multiRun[[index]])
      multiRun$meanSD[index,] <- c(currentMean, currentStdev)
    }

    ## Calculate fivenums for signature extraction
    multiRun$fivenum <- matrix(nrow = 6, ncol = 5)
    indexes <- c("averCosSim","falseNeg","falsePos",
                 "truePos","TPR","PPV")
    rownames(multiRun$fivenum) <- indexes
    colnames(multiRun$fivenum) <- c("min","lower-hinge","median","upperhinge","max")
    for(index in indexes){
      currentFiveNum <- stats::fivenum(multiRun[[index]])
      multiRun$fivenum[index,] <- currentFiveNum
    }

    ## Plot boxplot + beeswarm plot for signature extraction
    if(FALSE){

      titles <- c("averCosSim" = "Average cosine similarity",
                  "falseNeg" = "False negatives",
                  "falsePos" = "False positives",
                  "truePos" = "True positives",
                  "TPR" = "True positive rate (sensitivity)",
                  "PPV" = "Positive predictive value (PPV)")
      subtitles <- c("averCosSim" = "",
                     "falseNeg" = "Number of ground-truth signatures not extracted",
                     "falsePos" = "Number of signatures extracted, but different from ground-truth signatures",
                     "truePos" = "Number of ground-truth signatures extracted",
                     "TPR" = "#True Positives / (#True Positives + #False Negatives)",
                     "PPV" = "#True Positives / (#True Positives + #False Positives)")

      ## ggplot2 boxplot + beeswarm plot
      ggplotList <- list()
      for(index in indexes){
        indexNum <- which(index == indexes)
        ggplotList[[index]] <- ggplot2::ggplot(
          data.frame(value = multiRun[[index]],
                     indexName = index),
          ggplot2::aes(x = .data$indexName, y = .data$value))
        ggplotList[[index]] <- ggplotList[[index]] +
          ggplot2::ggtitle(titles[indexNum],subtitle = subtitles[indexNum])
        ggplotList[[index]] <- ggplotList[[index]] +
          ggplot2::geom_boxplot() +
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE, size = 0.3) +
          ## Restrict the decimal numbers of values of indexes to be 2
          ggplot2::scale_y_continuous(labels =function(x) sprintf("%.2f", x))
      }
      for(gtSigName in gtSigNames){
        ggplotList[[gtSigName]] <- ggplot2::ggplot(
          data.frame(value = multiRun$cosSim[[gtSigName]],
                     gtSigName = gtSigName),
          ggplot2::aes(x = .data$gtSigName, y = .data$value))
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
          ggplot2::ggtitle(label = paste0("Cosine similarity to signature ",gtSigName),
                           subtitle = paste0("Considers all extracted signatures resembling ", gtSigName))
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
          ggplot2::geom_boxplot() +
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE, size = 0.3) +
          ## Restrict the decimal numbers of values of indexes to be 2
          ggplot2::scale_y_continuous(labels =function(x) sprintf("%.2f", x))
      }

      ## Print high-resolution extraction indexes into one png file
      ## Only include extraction index plots
      ## in tempPlotList.
      if(FALSE){
        tempPlotList <- list()
        for(index in indexes){
          tempPlotList[[index]] <- ggplotList[[index]]
        }
        suppressMessages(
          ggplot2::ggsave(
            filename = paste0(tool.dir,"/boxplot.extraction.png"),
            plot = ggpubr::ggarrange(plotlist = tempPlotList),
            device = "png",
            dpi = 1000,
            limitsize = FALSE
          )
        )
      }

      ## Print extraction indexes into one pdf file
      grDevices::pdf(paste0(tool.dir,"/boxplot.extraction.measures.pdf"), pointsize = 1)
      for (index in indexes) print(ggplotList[[index]])
      grDevices::dev.off()


      ## Print high-resolution extraction indexes into one png file
      ## Only include one-signature cosine similarity plots
      ## in tempPlotList.
      if(FALSE){
        tempPlotList <- list()
        for(gtSigName in gtSigNames){
          tempPlotList[[gtSigName]] <- ggplotList[[gtSigName]]
        }
        suppressMessages(
          ggplot2::ggsave(
            filename = paste0(tool.dir,"/boxplot.onesig.cossim.png"),
            plot = ggpubr::ggarrange(plotlist = tempPlotList),
            device = "png",
            dpi = 1000,
            limitsize = FALSE
          )
        )
      }

      ## Print extraction indexes into one pdf file
      grDevices::pdf(paste0(tool.dir,"/boxplot.onesig.cossim.pdf"), pointsize = 1)
      for (gtSigName in gtSigNames) print(ggplotList[[gtSigName]])
      grDevices::dev.off()
    }


    ## Check whether runs of the computational approach
    ## involves exposure inferrence summary.
    exposureFlag <- TRUE
    for(runName in run.names){
      runDir <- paste0(tool.dir,"/",runName)
      summaryDir <- paste0(runDir,"/summary")
      exposureDiffFile <- paste0(summaryDir,"/exposureDiff.RDa")
      if(!file.exists(exposureDiffFile)){
        exposureFlag <- FALSE
        message("Skip summarizing scaled Manhattan distance...\n")
        break
      }
    }

    ## Summarize aggregated scaled Manhattan distance
    ## only if there are
    ## exposureDiff.Rda in all runDirs.
    if(exposureFlag){
      ## Read scaled aggregated Manhattan distances in multiple runs
      AggManhattanDist <- matrix(nrow = length(gtSigNames), ncol = length(run.names))
      rownames(AggManhattanDist) <- gtSigNames
      colnames(AggManhattanDist) <- run.names
      for(runName in run.names){
        runDir <- paste0(tool.dir,"/",runName)
        summaryDir <- paste0(runDir,"/summary")
        exposureDiffFile <- paste0(summaryDir,"/exposureDiff.RDa")
        ## Add exposureDiff <- NULL to please the R check
        exposureDiff <- NULL
        load(file = exposureDiffFile)
        AggManhattanDist[gtSigNames,runName] <- exposureDiff$aggregated[gtSigNames,"Scaled.Aggregated.Manhattan.distance"]
      }
      multiRun$AggManhattanDist <- AggManhattanDist

      ## Calculate mean and SD for aggregated  Manhattan distance
      meanSDAggMD <- matrix(nrow = length(gtSigNames), ncol = 2)
      rownames(meanSDAggMD) <- gtSigNames
      colnames(meanSDAggMD) <- c("mean","stdev")
      for(gtSigName in gtSigNames){
        meanSDAggMD[gtSigName,"mean"] <- mean(AggManhattanDist[gtSigName,])
        meanSDAggMD[gtSigName,"stdev"] <- stats::sd(AggManhattanDist[gtSigName,])
      }
      multiRun$meanSDAggMD <- meanSDAggMD

      ## Calculate fivenums for exposure inference Scaled Manhattan distance
      multiRun$fivenumAggMD <- matrix(nrow = length(gtSigNames), ncol = 5)
      rownames(multiRun$fivenumAggMD) <- gtSigNames
      colnames(multiRun$fivenumAggMD) <- c("min","lower-hinge","median","upperhinge","max")
      for(gtSigName in gtSigNames){
        multiRun$fivenumAggMD[gtSigName,] <- stats::fivenum(AggManhattanDist[gtSigName,])
      }
    }

    ## Summarize Manhattan distance scaled for
    ## individual tumors only if there are
    ## exposureDiff.Rda in all runDirs.
    if(exposureFlag){

      if(TRUE){
        ## Read scaled aggregated Manhattan distances in multiple runs
        meanSepMD <- matrix(nrow = length(gtSigNames), ncol = length(run.names))
        sdSepMD <- matrix(nrow = length(gtSigNames), ncol = length(run.names))
        rownames(meanSepMD) <- gtSigNames
        colnames(meanSepMD) <- run.names
        rownames(sdSepMD) <- gtSigNames
        colnames(sdSepMD) <- run.names

        for(runName in run.names){
          runDir <- paste0(tool.dir,"/",runName)
          summaryDir <- paste0(runDir,"/summary")
          exposureDiffFile <- paste0(summaryDir,"/exposureDiff.RDa")
          ## Add exposureDiff <- NULL to please the R check
          exposureDiff <- NULL
          load(file = exposureDiffFile)
          for(gtSigName in gtSigNames){
            meanSepMD[gtSigName,runName] <- mean(exposureDiff$separated[[gtSigName]][,"Scaled.Manhattan.distance"])
            sdSepMD[gtSigName,runName] <- sd(exposureDiff$separated[[gtSigName]][,"Scaled.Manhattan.distance"])
          }
        }
        multiRun$meanSepMD <- meanSepMD
        multiRun$sdSepMD <- sdSepMD
      }

    }


    ## Save data and results
    save(multiRun,file = paste0(tool.dir,"/multiRun.RDa"))
    write.csv(x = multiRun$meanSD,
              file = paste0(tool.dir,"/meanSD.csv"))
    write.csv(x = multiRun$fivenum,
              file = paste0(tool.dir,"/fivenum.csv"))
    if(exposureFlag){
      write.csv(x = multiRun$AggManhattanDist,
                file = paste0(tool.dir,"/Scaled.Aggregated.ManhattanDist.csv"))
      write.csv(x = multiRun$meanSDAggMD,
                file = paste0(tool.dir,"/meanSD.Scaled.Aggregated.Manhattan.dist.csv"))
      write.csv(x = multiRun$fivenumAggMD,
                file = paste0(tool.dir,"/fivenum.Scaled.Aggregated.Manhattan.dist.csv"))
      write.csv(x = multiRun$meanSepMD,
                file = paste0(tool.dir,"/mean.of.sep.Scaled.Manhattan.dist.csv"))
      write.csv(x = multiRun$sdSepMD,
                file = paste0(tool.dir,"/stdev.of.sep.Scaled.Manhattan.dist.csv"))
    }
    invisible(multiRun)
  }



#' Combine results for a single dataset, from different computational approaches.
#'
#' Summarize results from each computational approach in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#'
#' @param third.level.dir Third level path distinguishing de novo extraction
#' + attribution packages from attribution-only packages.
#' Examples:
#' \code{top.dir}/sp.sp/ExtrAttr/
#' \code{top.dir}/sa.sa/Attr/
#'
#' @param toolNames Names of computational approach. (e.g. "SigProExtractor")
#'
#' @param tool.dirnames Third level path from the \code{top.dir}. Expected to have
#' summarized results generated by \code{\link{SummarizeMultiRuns}}.
#' (multiRun.RDa, ManhattanDist.csv, meanSD.csv, meanSD.Manhattan.dist.csv)
#' Examples:
#' \code{"signeR.results"} (Under \code{third.level.dir} "ExtrAttr")
#' \code{"deconstructSigs.results"} (Under \code{third.level.dir} "Attr")
#'
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{tool.names}
#' which stores the software output.
#'
#' @param datasetGroup Numeric or character vector specifying the groups
#' each dataset belong to.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider slope as the group:
#' c("slope=0.1","slope=0.5","slope=1","slope=2","slope=5","slope=10")
#' Default: "Default"
#'
#' @param datasetGroupName Meaning or label of all datasetGroup.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider \code{"SBS1:SBS5 mutation count ratio"}
#' as the label of the \code{datasetGroup} slope.
#'
#' @param datasetSubGroup Optional. Numeric or character vector differentiating
#' datasets within each group.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider Pearson's R^2
#' as the subgroup:
#' c("Rsq=0.1","Rsq=0.2","Rsq=0.3","Rsq=0.6")
#' Default: Names of datasets, which are \code{basename(dataset.dirs)}
#'
#' @param datasetSubGroupName Optional. Meaning or label of all datasetSubGroup.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider \code{"Pearson's R squared"}
#' as the label of the \code{datasetSubGroup} Pearson's R^2.
#'
#' @return A list contain c(\code{mean},\code{sd}) of multiple runs:
#' Cosine similarity
#' True Positives(TP): Ground-truth signatures which are active in
#' the spectra, and extracted.
#' False Negatives(FN): Ground-truth signatures not extracted.
#' False Positives(FP): Signatures wrongly extracted, not resembling
#' any ground-truth signatures.
#' True positive rate (TPR, Sensitivity): TP / (TP + FN)
#' Positive predictive value (PPV, Precision): TP / (FP + TP)
#'
#' @details This function generates \code{multiTools.RDa} under
#' \code{third.level.dir}
#'
#' @importFrom utils write.csv capture.output sessionInfo
#'
#' @export
#'
SummarizeMultiToolsOneDataset <- function(
  third.level.dir,
  toolNames,
  tool.dirnames,
  datasetGroup,
  datasetGroupName,
  datasetSubGroup = NULL,
  datasetSubGroupName = NULL){

  multiTools <- list()
  combMeanSD <- NULL
  combMeanSDAggMD <- NULL

  for(toolNumber in 1:length(toolNames)){
    toolName <- toolNames[toolNumber]
    toolDirName <- tool.dirnames[toolNumber]
    toolPath <- paste0(third.level.dir,"/",toolDirName)
    ## Add multiRun <- NULL to please the R check
    multiRun <- NULL
    datasetName <- NULL
    load(paste0(toolPath,"/multiRun.RDa"))
    if(!is.null(datasetName)) {
      if(datasetName != multiRun$datasetName) {
        stop("Must provide results of different approaches on the SAME dataset.\n")
      }
    }
    datasetName <- multiRun$datasetName

    ## Combine multi-runs and multi-tools for each measure
    {
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","PPV")
      indexLabels <- c("averCosSim" = "Average cosine similarity of all signatures",
                       "falseNeg" = "Number of False negatives",
                       "falsePos" = "Number of False positives",
                       "truePos" = "Number of True positives",
                       "TPR" = "True positive rate",
                       "PPV" = "Positive predictive value")
      for(index in indexes){
        indexNum <- which(index == indexes)
        if(!exists("datasetSubGroup")) { # datasetSubGroup is not provided
          measure4OneTool <- data.frame(seed = names(multiRun[[index]]),
                                        value = multiRun[[index]],
                                        toolName = toolName,
                                        datasetName = multiRun$datasetName,
                                        datasetGroup = datasetGroup,
                                        stringsAsFactors = FALSE)
        } else {
          measure4OneTool <- data.frame(seed = names(multiRun[[index]]),
                                        value = multiRun[[index]],
                                        toolName = toolName,
                                        datasetName = multiRun$datasetName,
                                        datasetGroup = datasetGroup,
                                        datasetSubGroup = datasetSubGroup,
                                        stringsAsFactors = FALSE)
        }
        rownames(measure4OneTool) <- NULL
        ## Create a data.frame for each index,
        ## and summarize multi-Run, multiDataset values
        ## for each index.
        if(is.null(multiTools[[index]])){
          multiTools[[index]] <- data.frame()
        }
        multiTools[[index]] <- rbind(multiTools[[index]],measure4OneTool)
      }
    }

    ## meanSD contains mean and standard deviation
    ## for each extraction measure.
    {
      meanSD <- multiRun$meanSD
      colnames(meanSD) <- paste0(toolDirName,".", colnames(meanSD))
      if(is.null(meanSD)){
        combMeanSD <- meanSD
      } else{
        combMeanSD <- cbind(combMeanSD,meanSD)
      }
    }

    ## Combine multi-runs and multi-tools for
    ## one-signature cosine similarity.
    {
      gtSigNames <- names(multiRun$cosSim)
      multiTools$gtSigNames <- gtSigNames
      if(is.null(multiTools$cosSim)) multiTools$cosSim <- list()

      for(gtSigName in gtSigNames){
        if(!exists("datasetSubGroup")) {
          gtMeanCosSim4OneTool <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                                             value = multiRun$cosSim[[gtSigName]],
                                             toolName = toolName,
                                             datasetName = multiRun$datasetName,
                                             datasetGroup = datasetGroup,
                                             stringsAsFactors = FALSE)
        } else {
          gtMeanCosSim4OneTool <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                                             value = multiRun$cosSim[[gtSigName]],
                                             toolName = toolName,
                                             datasetName = multiRun$datasetName,
                                             datasetGroup = datasetGroup,
                                             datasetSubGroup = datasetSubGroup,
                                             stringsAsFactors = FALSE)
        }
        rownames(gtMeanCosSim4OneTool) <- NULL
        ## Create a data.frame for each ground-truth signature,
        ## and summarize multi-Run, multiDataset values
        ## for each ground-truth signature.
        if(is.null(multiTools$cosSim[[gtSigName]])){
          multiTools$cosSim[[gtSigName]] <- data.frame()
        }
        multiTools$cosSim[[gtSigName]] <- rbind(multiTools$cosSim[[gtSigName]],gtMeanCosSim4OneTool)
      }
    }

    ## Combine multi-runs and multi-tools for
    ## aggregated scaled Manhattan distance.
    if(!is.null(multiRun$AggManhattanDist)){
      ## Combine multi-runs and multi-tools for Manhattan
      ## distance of each ground-truth signature
      {
        if(is.null(multiTools$AggManhattanDist)) multiTools$AggManhattanDist <- list()
        for(gtSigName in gtSigNames){
          if(!exists("datasetSubGroup")) {
            gtAggManhattanDist4OneTool <- data.frame(seed = colnames(multiRun$AggManhattanDist),
                                                     value = multiRun$AggManhattanDist[gtSigName,],
                                                     toolName = toolName,
                                                     datasetName = multiRun$datasetName,
                                                     datasetGroup = datasetGroup,
                                                     stringsAsFactors = FALSE)
          } else{
            gtAggManhattanDist4OneTool <- data.frame(seed = colnames(multiRun$AggManhattanDist),
                                                     value = multiRun$AggManhattanDist[gtSigName,],
                                                     toolName = toolName,
                                                     datasetName = multiRun$datasetName,
                                                     datasetGroup = datasetGroup,
                                                     datasetSubGroup = datasetSubGroup,
                                                     stringsAsFactors = FALSE)
          }
          rownames(gtAggManhattanDist4OneTool) <- NULL
          ## Create a data.frame for each ground-truth signature,
          ## and summarize multi-Run, multiDataset values
          ## for each ground-truth signature.
          if(is.null(multiTools$AggManhattanDist[[gtSigName]])){
            multiTools$AggManhattanDist[[gtSigName]] <- data.frame()
          }
          multiTools$AggManhattanDist[[gtSigName]] <- rbind(multiTools$AggManhattanDist[[gtSigName]],gtAggManhattanDist4OneTool)
        }
      }

      ## meanSDAggMD contains mean and standard deviation
      ## for aggregated Scaled Manhattan distance between ground-truth exposures
      ## and inferred exposures for each ground-truth signature
      {
        meanSDAggMD <- multiRun$meanSDAggMD
        colnames(meanSDAggMD) <- paste0(toolDirName,".", colnames(meanSDAggMD))
        if(is.null(meanSDAggMD)){
          combMeanSDAggMD <- meanSDAggMD
        } else{
          combMeanSDAggMD <- cbind(combMeanSDAggMD,meanSDAggMD)
        }
      }
    }


    ## Combine multi-runs and multi-tools for
    ## mean of scaled Manhattan distance for each tumor.
    if(!is.null(multiRun$meanSepMD)){
      ## Combine multi-runs and multi-tools for Manhattan
      ## distance of each ground-truth signature
      if(is.null(multiTools$meanSepMD)) multiTools$meanSepMD <- list()
      for(gtSigName in gtSigNames){
        if(!exists("datasetSubGroup")) {
          gtmeanSepMD4OneTool <- data.frame(seed = colnames(multiRun$meanSepMD),
                                            value = multiRun$meanSepMD[gtSigName,],
                                            toolName = toolName,
                                            datasetName = multiRun$datasetName,
                                            datasetGroup = datasetGroup,
                                            stringsAsFactors = FALSE)
        } else{
          gtmeanSepMD4OneTool <- data.frame(seed = colnames(multiRun$meanSepMD),
                                            value = multiRun$meanSepMD[gtSigName,],
                                            toolName = toolName,
                                            datasetName = multiRun$datasetName,
                                            datasetGroup = datasetGroup,
                                            datasetSubGroup = datasetSubGroup,
                                            stringsAsFactors = FALSE)
        }
        rownames(gtmeanSepMD4OneTool) <- NULL
        ## Create a data.frame for each ground-truth signature,
        ## and summarize multi-Run, multiDataset values
        ## for each ground-truth signature.
        if(is.null(multiTools$meanSepMD[[gtSigName]])){
          multiTools$meanSepMD[[gtSigName]] <- data.frame()
        }
        multiTools$meanSepMD[[gtSigName]] <- rbind(multiTools$meanSepMD[[gtSigName]],gtmeanSepMD4OneTool)
      }
    }

    ## Combine multi-runs and multi-tools for
    ## standard deviation of scaled Manhattan distance for each tumor.
    if(!is.null(multiRun$sdSepMD)){
      ## Combine multi-runs and multi-tools for Manhattan
      ## distance of each ground-truth signature

      if(is.null(multiTools$sdSepMD)) multiTools$sdSepMD <- list()
      for(gtSigName in gtSigNames){
        if(!exists("datasetSubGroup")) {
          gtsdSepMD4OneTool <- data.frame(seed = colnames(multiRun$sdSepMD),
                                          value = multiRun$sdSepMD[gtSigName,],
                                          toolName = toolName,
                                          datasetName = multiRun$datasetName,
                                          datasetGroup = datasetGroup,
                                          stringsAsFactors = FALSE)
        } else{
          gtsdSepMD4OneTool <- data.frame(seed = colnames(multiRun$sdSepMD),
                                          value = multiRun$sdSepMD[gtSigName,],
                                          toolName = toolName,
                                          datasetName = multiRun$datasetName,
                                          datasetGroup = datasetGroup,
                                          datasetSubGroup = datasetSubGroup,
                                          stringsAsFactors = FALSE)
        }
        rownames(gtsdSepMD4OneTool) <- NULL
        ## Create a data.frame for each ground-truth signature,
        ## and summarize multi-Run, multiDataset values
        ## for each ground-truth signature.
        if(is.null(multiTools$sdSepMD[[gtSigName]])){
          multiTools$sdSepMD[[gtSigName]] <- data.frame()
        }
        multiTools$sdSepMD[[gtSigName]] <- rbind(multiTools$sdSepMD[[gtSigName]],gtmeanSepMD4OneTool)
      }
    }


  }

  multiTools$combMeanSD <- combMeanSD
  if(exists("combMeanSDAggMD")){
    multiTools$combMeanSDAggMD <- combMeanSDAggMD
  }
  multiTools$datasetName <- datasetName
  multiTools$datasetGroupName <- datasetGroupName
  multiTools$datasetSubGroupName <- datasetSubGroupName

  save(multiTools,file = paste0(third.level.dir,"/multiTools.RDa"))
  write.csv(x = multiTools$combMeanSD,
            file = paste0(third.level.dir,"/combined.meanSD.csv"))
  if(!is.null(multiTools$combMeanSDAggMD)){
    write.csv(x = multiTools$combMeanSDAggMD,
              file = paste0(third.level.dir,"/combined.meanSD.Aggregated.Manhattan.dist.csv"))
  }
  if(!is.null(multiTools$meanSepMD)){
    write.csv(x = multiTools$meanSepMD,
              file = paste0(third.level.dir,"/mean.of.sep.Scaled.Manhattan.dist.csv"))
  }
  if(!is.null(multiTools$sdSepMD)){
    write.csv(x = multiTools$sdSepMD,
              file = paste0(third.level.dir,"/stdev.of.sep.Scaled.Manhattan.dist.csv"))
  }
  invisible(multiTools)
}

#' Summarize results for multiple datasets, by different computational approaches.
#'
#' Summarize results of mutational signature extraction and exposure inferrence
#' by multiple computational approaches on multiple datasets. Before running this
#' function, make sure the summary file for each single data set
#'  \code{third.level.dir}/\code{multiTools.Rda} exists.
#'
#' \code{multiTools.Rda} is generated by \code{\link{SummarizeMultiToolsOneDataset}}).
#'
#' @param dataset.dirs Paths of top-level dataset directories trees you want
#' to investigate.
#' E.g. "./S.0.1.Rsq.0.1"
#'
#' @param second.third.level.dirname Name of the second.level.dir (e.g. "sp.sp")
#' and the third.level.dir (e.g. "ExtrAttr") to be investigated.
#'
#' Examples are: "sp.sp/ExtrAttr", "sa.sa.96/Attr"
#'
#' Note: \code{multiTools.RDa} are expected to be exist under
#' \code{dataset.dirs}/\code{second.third.level.dirname}.
#'
#' @param out.dir Path of the output directory.
#'
#' @param overwrite Whether to overwrite the contents in out.dir if
#' it already exists. (Default: FALSE)
#'
#' @importFrom rlang .data
#'
#' @importFrom utils write.csv
#'
#' @export
#'
SummarizeMultiToolsMultiDatasets <-
  function(dataset.dirs,
           second.third.level.dirname,
           out.dir,
           overwrite = FALSE){

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## For each measure,
    ## Create a data.frame integrating results of
    ## all runs and for all datasets
    {
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","PPV")
      indexLabels <- c("averCosSim" = "Average cosine similarity of all signatures",
                       "falseNeg" = "False negatives",
                       "falsePos" = "False positives",
                       "truePos" = "True positives",
                       "TPR" = "True positive rate (TPR, sensitivity)",
                       "PPV" = "Positive predictive value (PPV, precision)")
      indexNums <- length(indexes)
    }

    ## Summarizing different measures for extraction performance
    ## separately into list "FinalExtr".
    ## Showing individual values rather than
    ## only showing mean and standard deviation of multiple runs
    {
      FinalExtr <- list()
      toolNames <- character(0)
      for(index in indexes) {
        FinalExtr[[index]] <- data.frame()
      }
      FinalExtr$cosSim <- list()


      ## Combine extraction measures of different datasets:
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        datasetGroupName <- multiTools$datasetGroupName
        datasetSubGroupName <- multiTools$datasetSubGroupName


        ## Find tool names
        toolNames <- unique(multiTools[["averCosSim"]][,"toolName"])

        ## Bind values of measures in multiTools into FinalExtr.
        for(index in indexes){
          FinalExtr[[index]] <- rbind(FinalExtr[[index]],multiTools[[index]])
        }

        ## Bind values of cosine similarity in multiTools$cosSim into FinalExtr$cosSim
        gtSigNames <- multiTools$gtSigNames
        if(length(FinalExtr$cosSim) == 0){
          for(gtSigName in gtSigNames) {
            FinalExtr$cosSim[[gtSigName]] <- data.frame()
          }
        }
        for(gtSigName in gtSigNames){
          FinalExtr$cosSim[[gtSigName]] <- rbind(FinalExtr$cosSim[[gtSigName]],multiTools$cosSim[[gtSigName]])
        }
      }

      ## Calculate composite measure for each datasetDir.
      ## It equals to:
      ## True Positive Rate (TPR) + Positive Predictive Value (PPV)
      ## Cosine similarity to each of signature (SBS1 and SBS5 in SBS1-SBS5 paper)
      FinalExtr$compositeMeasure <- FinalExtr$TPR
      FinalExtr$compositeMeasure$value <- FinalExtr$TPR$value + FinalExtr$PPV$value
      for(gtSigName in gtSigNames){
        FinalExtr$compositeMeasure$value <- FinalExtr$compositeMeasure$value + FinalExtr$cosSim[[gtSigName]]$value
      }

    }

    ## Generating csv tables for extraction performance measure
    ## and cosine similarities.
    {

      ## Output combined extraction
      for(index in c(indexes,"compositeMeasure")){

        output <- FinalExtr[[index]]

        output <- output[,-4]
        colnames(output)[1] <- "Seed or run number"
        colnames(output)[2] <- paste0("Cosine similarity to ground-truth signature ",gtSigName)
        colnames(output)[3] <- "Name of computational approach"
        colnames(output)[4] <- datasetGroupName
        colnames(output)[5] <- datasetSubGroupName


        write.csv(output,
                  file = paste0(out.dir,"/",index,".csv"))
      }


      for(gtSigName in gtSigNames){

        output <- FinalExtr$cosSim[[gtSigName]]

        output <- output[,-4]
        colnames(output)[1] <- "Seed or run number"
        colnames(output)[2] <- paste0("Cosine similarity to ground-truth signature ",gtSigName)
        colnames(output)[3] <- "Name of computational approach"
        colnames(output)[4] <- datasetGroupName
        colnames(output)[5] <- datasetSubGroupName

        write.csv(output,
                  file = paste0(out.dir,"/cossim.to.",gtSigName,".csv"))
      }
    }


    ## For each extraction measures,
    ## merge values from multiple runs
    ## into one data.frame FinalExtr$<measure_name>
    ## and computational approaches for easier plotting.
    {
      ## Combine all extraction measurements, FinalExtr[[index]] into FinalExtr$Combined
      ## combined all one-signature cosine similarity, FinalExtr[[gtSigName]] into FinalExtr$Combined
      FinalExtr$combined <- data.frame()
      for(index in c("TPR","PPV")){
        plotDFOneMeasure <- data.frame(FinalExtr[[index]], indexLabel = indexLabels[index])
        FinalExtr$combined <- rbind(FinalExtr$combined,plotDFOneMeasure)
      }

      plotDFOneMeasure <- data.frame(FinalExtr$compositeMeasure, indexLabel = "Composite measure")
      FinalExtr$combined <- rbind(FinalExtr$combined,plotDFOneMeasure)

      for(gtSigName in gtSigNames){
        plotDFOneMeasure <- data.frame(FinalExtr$cosSim[[gtSigName]], indexLabel = paste0("Cosine similarity to ",gtSigName))
        FinalExtr$combined <- rbind(FinalExtr$combined,plotDFOneMeasure)
      }


      ## Convert FinalExtr$combined$datasetGroup and
      ## Let their levels follow gtools::mixedsort() fashion
      ## So that the order of the facet labels will be more reasonable for readers.
      FinalExtr$combined$datasetGroup <- factor(
        FinalExtr$combined$datasetGroup,
        levels = gtools::mixedsort(unique(FinalExtr$combined$datasetGroup)))

      if(!is.null(multiTools$datasetSubGroupName)) {
        FinalExtr$combined$datasetSubGroup <- factor(
          FinalExtr$combined$datasetSubGroup,
          levels = gtools::mixedsort(unique(FinalExtr$combined$datasetSubGroup)))
      }
    }



    ## Plot general pdf for all extraction measures
    ## Plot a general violin + beeswarm plot for multiple measures
    ## in all runs and in all datasets.
    {

      ## Plot a multi-facet ggplot for all measures and all runs.
      {
        ggplotList <- list()
        ## Generate a ggplot object based on FinalExtr$combined
        ggplotList$general <- ggplot2::ggplot(
          FinalExtr$combined,
          ggplot2::aes(x = .data$toolName, y = .data$value)) +
          ## Draw geom_violin and geom_quasirandom
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            #ggplot2::aes(fill = index),
            ## Maximize the violin plot width
            scale = "width",
            ## Make bandwidth larger
            #position = "dodge",
            #width = 1.2
            ## Hide outliers
            #outlier.shape = NA
          ) +
          #
          #ggbeeswarm::geom_quasirandom(
          #  groupOnX = TRUE, size = 0.3
          #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
          #) +
          ## Show median of the extraction measure distribution, as a solid dot.
          ggplot2::stat_summary(fun.y="median", geom="point", fill = "red", shape = 21) +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", fill="blue", shape = 23) +
          ## Change title for general violin + beeswarm plot
          ggplot2::ggtitle(label = "Extraction performance",
                           subtitle = "for all methods across all data sets and replicates") +
          ## Change axis titles
          ggplot2::labs(x = "Computational approach") +
          ## Remove axis.title.y (defaults to be "value", meaningless)
          ## Rotate axis.text.x 90 degrees,
          ## move axis.text.x right below the tick marks,
          ## and remove legends.
          ggplot2::theme(
            ## Remove axis.title.y
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(
              ## Rotate the axis.text.x
              angle = 90,
              ## move axis.text.x right below the tick marks
              hjust = 1,vjust = 0.5),
            ## Make font size of facet label smaller.
            strip.text = ggplot2::element_text(size = 10),
            ## remove legends.
            legend.position = "none") +
          ## Split the plot into multiple facets,
          ## according to different measures
          ggplot2::facet_wrap(
            ggplot2::vars(.data$indexLabel),
            ## Force facet_wrap to have 2 columns
            ncol = 2,
            scales = "free",
            ## Let facet label to be print on multiple lines
            labeller = ggplot2::label_wrap_gen(multi_line = T),
            ## Let facets be plotted horizontally
            dir = "h",
            ## Put facet label to the top
            strip.position = "top"
          ) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a multi-facet ggplot,
      ## facets are separated by measures and datasetGroup
      ## (in example, it refers to slope.)
      if(!is.null(multiTools$datasetSubGroupName)) {
        bys <- c("datasetGroup","datasetSubGroup")
      } else{
        bys <- c("datasetGroup")
      }

      for(by in bys)  {

        ## The value of "datasetGroupName" or "datasetSubGroupName"
        ## which is the caption of "datasetGroup"
        byCaption <- eval(parse(text = paste0("multiTools$",by,"Name")))

        ## Generate a ggplot object based on FinalExtr$combined
        ggplotList[[by]] <- ggplot2::ggplot(
          FinalExtr$combined,
          ggplot2::aes(x = .data$toolName, y = .data$value)) +
          ## Draw geom_violin and geom_quasirandom
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            #ggplot2::aes(fill = index),
            ## Maximize the violin plot width
            scale = "width",
            ## Make bandwidth larger
            #position = "dodge",
            #width = 1.2
            ## Hide outliers
            #outlier.shape = NA
          ) +
          #
          #ggbeeswarm::geom_quasirandom(
          #  groupOnX = TRUE, size = 0.3
          #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
          #) +
          ## Show median of the extraction measure distribution, as a solid dot.
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ## Change axis titles
          ggplot2::labs(x = "Computational approach") +
          ## Remove axis.title.y (defaults to be "value", meaningless)
          ## Rotate the axis.text.x (names of tools),
          ## move axis.text.x right below the tick marks
          ## and remove legends
          ggplot2::theme(
            ## Remove axis.title.y
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(
              ## Rotate the axis.text.x (names of tools)
              angle = 90,
              ## move axis.text.x right below the tick marks
              hjust = 1, vjust = 0.5),
            ## remove legends
            legend.position = "none") +
          ## Split the plot into multiple facets,
          ## according to different measures
          ggplot2::facet_grid(
            rows =  ggplot2::vars(.data$indexLabel),
            cols = eval(parse(text = paste0("ggplot2::vars(",by,")"))),
            scales = "free",
            ## Let facet label to be print on multiple lines
            labeller = ggplot2::label_wrap_gen(multi_line = T)
          ) +
          ## Make facet label font size smaller
          ggplot2::theme(strip.text.y = ggplot2::element_text(size = 4)) +
          ## Add title for general violin + beeswarm plot
          ggplot2::ggtitle(
            label = paste0("Measures of extraction performance as a function of"),
            subtitle = paste0(byCaption,".")) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(labels =function(x) sprintf("%.2f", x))
      }
      ## Plot violin + beeswarm plots in pdf format
      grDevices::pdf(paste0(out.dir,"/extraction.violins.pdf"),
                     #paper = "a4", ## A4 size
                     width = 7,
                     height = 10.5, ## Make height larger
                     pointsize = 1)
      for(by in names(ggplotList)){
        print(ggplotList[[by]])
      }
      grDevices::dev.off()
    }



    ## Plot general png and pdf for one-signature cosine similarity summary
    ## Plot a general violin + beeswarm plot for multiple signatures
    ## in all runs and in all datasets.
    {
      ## For ground-truth signature,
      ## Create a data.frame integrating results of
      ## all runs and for all datasets
      gtSigNames <- multiTools$gtSigNames
      sigNums <- length(gtSigNames)

      ## Combine all FinalExtr$cosSim[[gtSigName]] into FinalExtr$cosSimCombined
      FinalExtr$cosSimCombined <- data.frame()
      for(gtSigName in gtSigNames){
        plotDFOneMeasure <- data.frame(FinalExtr$cosSim[[gtSigName]], gtSigName = gtSigName)
        FinalExtr$cosSimCombined <- rbind(FinalExtr$cosSimCombined,plotDFOneMeasure)
      }
      ## Convert FinalExtr$combined$datasetGroup and
      ## Let their levels follow gtools::mixedsort() fashion
      ## So that the order of the facet labels will be more reasonable for readers.
      FinalExtr$cosSimCombined$datasetGroup <- factor(
        FinalExtr$cosSimCombined$datasetGroup,
        levels = gtools::mixedsort(unique(FinalExtr$cosSimCombined$datasetGroup)))

      if(!is.null(multiTools$datasetSubGroupName)) {
        FinalExtr$cosSimCombined$datasetSubGroup <- factor(
          FinalExtr$cosSimCombined$datasetSubGroup,
          levels = gtools::mixedsort(unique(FinalExtr$cosSimCombined$datasetSubGroup)))
      }


      ggplotList <- list()
      ## Plot a multi-facet ggplot for all gtSigNames and all runs.
      {
        ## Generate a ggplot object based on FinalExtr$combined
        ggplotList$general <- ggplot2::ggplot(
          FinalExtr$cosSimCombined,
          ggplot2::aes(x = .data$toolName, y = .data$value))
        ## Draw geom_violin and geom_quasirandom
        ggplotList$general <- ggplotList$general +
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            #ggplot2::aes(fill = gtSigName),
            ##
            scale = "width",
            ## Make bandwidth larger
            #position = "dodge",
            #width = 1.2
            ## Hide outliers
            #outlier.shape = NA
          ) +
          #
          #ggbeeswarm::geom_quasirandom(
          #  groupOnX = TRUE, size = 0.3
          #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
          #) +
          ## Show median of the cosine similarity distribution
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ## Add title for general violin + beeswarm plot
          ggplot2::ggtitle(label = "Average cosine similarity between ground-truth and extracted signatures",
                           subtitle = "for all computational approaches, ratios and correlation values.") +
          ## Change axis titles
          ggplot2::labs(x = "Computational approach",
                        y = "Cosine Similarity") +
          ## Rotate the axis.text.x (names of tools),
          ## move axis.text.x right below the tick marks
          ## and remove legends
          ggplot2::theme(axis.text.x = ggplot2::element_text(
            ## Rotate the axis.text.x (names of tools)
            angle = 90,
            ## move axis.text.x right below the tick marks
            hjust = 1, vjust = 0.5),
            ## remove legends.
            legend.position = "none") +
          ## Split the plot into multiple facets,
          ## according to different gtSigNames
          ggplot2::facet_wrap(
            ggplot2::vars(gtSigName),
            ## Force facet_wrap to have 2 columns
            ncol = 2,
            scales = "free",
            ## Let facets be plotted vertically
            dir = "v"
          ) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(
            ## For one-signature cosine similarity, set ylim from the minimum of cosine similarity value to 1
            limits = c(min(FinalExtr$cosSimCombined$value), 1),
            labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a multi-facet ggplot,
      ## facets are separated by gtSigNames and datasetGroup
      ## (in example, it refers to slope.)
      if(!is.null(multiTools$datasetSubGroupName)) {
        bys <- c("datasetGroup","datasetSubGroup")
      } else {
        bys <- c("datasetGroup")
      }

      for(by in bys)  {

        ## The value of "datasetGroupName" or "datasetSubGroupName"
        ## which is the caption of "datasetGroup"
        byCaption <- eval(parse(
          text = paste0("multiTools$",by,"Name")))

        ## Generate a ggplot object based on FinalExtr$combined
        ggplotList[[by]] <- ggplot2::ggplot(
          FinalExtr$cosSimCombined,
          ggplot2::aes(x = .data$toolName, y = .data$value)) +
          ## Draw geom_violin and geom_quasirandom
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            #ggplot2::aes(fill = gtSigName),
            ## Maximize the violin plot width
            scale = "width",
            ## Make bandwidth larger
            #position = "dodge",
            #width = 1.2
            ## Hide outliers
            #outlier.shape = NA
          ) +
          #ggbeeswarm::geom_quasirandom(
          #  groupOnX = TRUE, size = 0.3
          #  ## Need to add a single color (different from black)
          #  ## for all data points.
          #  , ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
          #) +
          ## Show median of the extraction measure distribution
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ## Add title for general violin + beeswarm plot
          ggplot2::ggtitle(
            label = paste0("Extraction cosine similarity as a function of"),
            subtitle = paste0(byCaption,".")) +
          ## Change axis titles
          ggplot2::labs(x = "Computational approach",
                        y = "Cosine Similarity") +
          ## Rotate the axis.text.x (names of tools),
          ## move axis.text.x right below the tick marks
          ## and remove legends
          ggplot2::theme(axis.text.x = ggplot2::element_text(
            ## Rotate the axis.text.x (names of tools)
            angle = 90,
            ## move axis.text.x right below the tick marks
            hjust = 1, vjust = 0.5),
            ## remove legends.
            legend.position = "none") +
          ## Split the plot into multiple facets,
          ## according to different gtSigNames
          ggplot2::facet_grid(rows = ggplot2::vars(gtSigName),
                              cols = eval(parse(text = paste0("ggplot2::vars(",by,")"))),
                              scales = "free") +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(
            ## For one-signature cosine similarity, set ylim from the minimum of cosine similarity value to 1
            limits = c(min(FinalExtr$cosSimCombined$value), 1),
            labels = function(x) sprintf("%.2f", x))
      }

      ## Plot violin + beeswarm plots in pdf format
      grDevices::pdf(paste0(out.dir,"/onesig.cossim.violins.pdf"), pointsize = 1)
      for(by in names(ggplotList)){
        print(ggplotList[[by]])
      }
      grDevices::dev.off()
    }

    ## Summarize aggregated scaled Manhattan distance only if
    ## multiTools$AggManhattanDist exists.
    {
      flagExposure <- TRUE
      ## Combine attribution assessment onto multiple sheets.
      ## Each sheet shows Scaled Manhattan distance for one mutational signature.
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))
        if(is.null(multiTools$AggManhattanDist)){
          flagExposure <- FALSE
          message("Skip summarizing scaled Manhattan distance...\n")
          break
        }
      }
    }

    ## Summarizing aggregated Scaled Manhattan distance results
    if(flagExposure){
      {
        FinalAttr <- list()
        FinalAttr$AggManhattanDist <- list()
        ## Combine attribution assessment onto multiple sheets.
        ## Each sheet shows Scaled Manhattan distance for one mutational signature.
        for(datasetDir in dataset.dirs){
          thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
          ## Add multiTools <- NULL to please R check
          multiTools <- NULL
          load(paste0(thirdLevelDir,"/multiTools.RDa"))

          gtSigNames <- multiTools$gtSigNames
          sigNums <- length(gtSigNames)

          if(length(FinalAttr$AggManhattanDist) == 0){
            for(gtSigName in gtSigNames) {
              FinalAttr$AggManhattanDist[[gtSigName]] <- data.frame()
            }
          }

          ## Combine Scaled Manhattan distance
          for(gtSigName in gtSigNames){
            FinalAttr$AggManhattanDist[[gtSigName]] <- rbind(
              FinalAttr$AggManhattanDist[[gtSigName]],
              multiTools$AggManhattanDist[[gtSigName]])
          }
        }

        ## For the purpose of SBS1-SBS5 paper,
        ## don't output summary tables for aggregated scaled Manhattan distance.
        if(FALSE){
          for(gtSigName in gtSigNames){
            output <- FinalAttr$AggManhattanDist[[gtSigName]]

            output <- output[,-4]
            colnames(output)[1] <- "Seed or run number"
            colnames(output)[2] <- paste0("Scaled distance of ",gtSigName)
            colnames(output)[3] <- "Name of computational approach"
            colnames(output)[4] <- datasetGroupName
            colnames(output)[5] <- datasetSubGroupName

            write.csv(output,
                      file = paste0(out.dir,"/Agg.ManhattanDist.",gtSigName,".csv"))
          }
        }

        ## Plot general png and pdf for attribution Scaled Manhattan distance summary
        ## Plot a general violin + beeswarm plot for multiple signatures
        ## in all runs and in all datasets.
        {

          ## Combine all FinalAttr$AggManhattanDist[[gtSigName]] into FinalAttr$AggManhattanDist$Combined
          FinalAttr$AggManhattanDist$combined <- data.frame()
          for(gtSigName in gtSigNames){
            plotDFOneMeasure <- data.frame(FinalAttr$AggManhattanDist[[gtSigName]], gtSigName = gtSigName)
            FinalAttr$AggManhattanDist$combined <- rbind(FinalAttr$AggManhattanDist$combined,plotDFOneMeasure)
          }

          ## Convert FinalAttr$AggManhattanDist$combined$datasetGroup and
          ## Let their levels follow gtools::mixedsort() fashion
          ## So that the order of the facet labels will be more reasonable for readers.
          FinalAttr$AggManhattanDist$combined$datasetGroup <- factor(
            FinalAttr$AggManhattanDist$combined$datasetGroup,
            levels = gtools::mixedsort(unique(FinalAttr$AggManhattanDist$combined$datasetGroup)))

          if(!is.null(multiTools$datasetSubGroupName)) {
            FinalAttr$AggManhattanDist$combined$datasetSubGroup <- factor(
              FinalAttr$AggManhattanDist$combined$datasetSubGroup,
              levels = gtools::mixedsort(unique(FinalAttr$AggManhattanDist$combined$datasetSubGroup)))
          }

          ggplotList <- list()
          ## Plot a multi-facet ggplot for all gtSigNames and all runs.
          {
            ## Generate a ggplot object based on FinalAttr$AggManhattanDist$combined
            ggplotList$general <- ggplot2::ggplot(
              FinalAttr$AggManhattanDist$combined,
              ggplot2::aes(x = .data$toolName, y = .data$value)) +
              ## Draw geom_violin and geom_quasirandom
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                #ggplot2::aes(fill = gtSigName),
                ## Maximize the violin plot width
                scale = "width",
                ## Make bandwidth larger
                #position = "dodge",
                #width = 1.2
                ## Hide outliers
                #outlier.shape = NA
              ) +
              #ggbeeswarm::geom_quasirandom(
              #  groupOnX = TRUE, size = 0.3
              #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
              #) +
              ## Show median of the Scaled Manhattan distance distribution
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Add title for general violin + beeswarm plot
              ggplot2::ggtitle(label = "Scaled aggregated Manhattan distance between inferred and ground-truth",
                               subtitle = "exposures for all computational approaches, ratios and correlation values.") +
              ## Change axis titles
              ggplot2::labs(x = "Computational approach",
                            y = "Scaled aggregated Manhattan distance") +
              ## Rotate the names of tools,
              ## move axis.text.x right below the tick marks
              ## and remove legends
              ggplot2::theme(axis.text.x = ggplot2::element_text(
                ## Rotate the axis.text.x (names of tools),
                angle = 90,
                ## move axis.text.x right below the tick marks
                hjust = 1, vjust = 0.5),
                ## remove legends.
                legend.position = "none") +
              ## Split the plot into multiple facets,
              ## according to different gtSigNames
              ggplot2::facet_wrap(
                ggplot2::vars(gtSigName),
                ## Force facet_wrap to have 2 columns
                ncol = 2,
                scales = "free",
                ## Let facets be plotted vertically
                dir = "v"
              ) +
              ## Restrict the decimal numbers of values of measures to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0, max(FinalAttr$AggManhattanDist$combined$value)),
                labels =function(x) sprintf("%.2f", x))
          }
          ## Plot a multi-facet ggplot,
          ## facets are separated by gtSigNames and datasetGroup
          ## (in example, it refers to slope.)
          if(!is.null(multiTools$datasetSubGroupName)) {
            bys <- c("datasetGroup","datasetSubGroup")
          } else {
            bys <- c("datasetGroup")
          }

          for(by in bys)  {

            ## The value of "datasetGroupName" or "datasetSubGroupName"
            ## which is the caption of "datasetGroup"
            byCaption <- eval(parse(
              text = paste0("multiTools$",by,"Name")))


            ## Generate a ggplot object based on FinalAttr$AggManhattanDist$combined
            ggplotList[[by]] <- ggplot2::ggplot(
              FinalAttr$AggManhattanDist$combined,
              ggplot2::aes(x = .data$toolName, y = .data$value))
            ## Draw geom_violin and geom_quasirandom
            ggplotList[[by]] <- ggplotList[[by]] +
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                #ggplot2::aes(fill = gtSigName),
                ## Maximize the violin plot width
                scale = "width"
                #,
                ## Make bandwidth larger
                #position = "dodge",
                #width = 1.2
                ## Hide outliers
                #outlier.shape = NA
              ) +
              #ggbeeswarm::geom_quasirandom(
              #  groupOnX = TRUE, size = 0.3
              #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
              #) +
              ## Show median of the Scaled Manhattan distance distribution
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Add title for general violin + beeswarm plot
              ggplot2::ggtitle(
                label = paste0("Scaled aggregated Manhattan distance summary plot as a function of "),
                subtitle = paste0("ground-truth signature names and ",byCaption,".")) +
              ## Change axis titles
              ggplot2::labs(x = "Computational approach",
                            y = "Scaled aggregated Manhattan distance") +
              ## Rotate the axis.text.x (names of tools),
              ## move axis.text.x right below the tick marks
              ## and remove legends
              ggplot2::theme(axis.text.x = ggplot2::element_text(
                ## Rotate the axis.text.x (names of tools),
                angle = 90,
                ## move axis.text.x right below the tick marks
                hjust = 1, vjust = 0.5),
                ## remove legends.
                legend.position = "none") +
              ## Split the plot into multiple facets,
              ## according to different gtSigNames
              ggplot2::facet_grid(rows =  ggplot2::vars(gtSigName),
                                  cols = eval(parse(text = paste0("ggplot2::vars(",by,")"))),
                                  scales = "free") +
              ## Restrict the decimal numbers of values of measures to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0, max(FinalAttr$AggManhattanDist$combined$value)),
                labels =function(x) sprintf("%.2f", x))
          }

          ## Plot violin + beeswarm plots in pdf format
          grDevices::pdf(paste0(out.dir,"/Aggregated.Scaled.Manhattan.Dist.violins.pdf"), pointsize = 1)
          for(by in names(ggplotList)){
            print(ggplotList[[by]])
          }
          grDevices::dev.off()
        }
      }


    }

    ## Summarizing results for mean and stdev of separated Manhattan distance
    if(flagExposure){

      fileNames = c(
        "meanSepMD" = "mean.of.sep.Scaled.Manhattan.dist",
        "sdSepMD" = "stdev.of.sep.Scaled.Manhattan.dist")
      titles = c(
        "meanSepMD" = "Mean of Manhattan distances of individual tumors",
        "sdSepMD" = "Standard deviation of Manhattan distances of individual tumors"
      )

      for(measure in c("meanSepMD", "sdSepMD"))
      {
        FinalAttr[[measure]] <- list()
        ## Combine attribution assessment onto multiple sheets.
        ## Each sheet shows Scaled Manhattan distance for one mutational signature.
        for(datasetDir in dataset.dirs){
          thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
          ## Add multiTools <- NULL to please R check
          multiTools <- NULL
          load(paste0(thirdLevelDir,"/multiTools.RDa"))

          gtSigNames <- multiTools$gtSigNames
          sigNums <- length(gtSigNames)

          if(length(FinalAttr[[measure]]) == 0){
            for(gtSigName in gtSigNames) {
              FinalAttr[[measure]][[gtSigName]] <- data.frame()
            }
          }

          ## Combine Scaled Manhattan distance
          for(gtSigName in gtSigNames){
            FinalAttr[[measure]][[gtSigName]] <- rbind(
              FinalAttr[[measure]][[gtSigName]],
              multiTools[[measure]][[gtSigName]])
          }
        }

        ## For the purpose of SBS1-SBS5 paper,
        ## don't output summary tables for scaled Manhattan distance.
        if(FALSE){
          for(gtSigName in gtSigNames){
            output <- FinalAttr[[measure]][[gtSigName]]

            output <- output[,-4]
            colnames(output)[1] <- "Seed or run number"
            colnames(output)[2] <- paste0("Scaled distance of ",gtSigName)
            colnames(output)[3] <- "Name of computational approach"
            colnames(output)[4] <- datasetGroupName
            colnames(output)[5] <- datasetSubGroupName

            write.csv(output,
                      file = paste0(out.dir,"/",fileNames[measure],".",gtSigName,".csv"))
          }
        }

        ## Plot general png and pdf for attribution Scaled Manhattan distance summary
        ## Plot a general violin + beeswarm plot for multiple signatures
        ## in all runs and in all datasets.
        {

          ## Combine all FinalAttr[[measure]][[gtSigName]] into FinalAttr[[measure]]$Combined
          FinalAttr[[measure]]$combined <- data.frame()
          for(gtSigName in gtSigNames){
            plotDFOneMeasure <- data.frame(FinalAttr[[measure]][[gtSigName]], gtSigName = gtSigName)
            FinalAttr[[measure]]$combined <- rbind(FinalAttr[[measure]]$combined,plotDFOneMeasure)
          }

          ## Convert FinalAttr[[measure]]$combined$datasetGroup and
          ## Let their levels follow gtools::mixedsort() fashion
          ## So that the order of the facet labels will be more reasonable for readers.
          FinalAttr[[measure]]$combined$datasetGroup <- factor(
            FinalAttr[[measure]]$combined$datasetGroup,
            levels = gtools::mixedsort(unique(FinalAttr[[measure]]$combined$datasetGroup)))

          if(!is.null(multiTools$datasetSubGroupName)) {
            FinalAttr[[measure]]$combined$datasetSubGroup <- factor(
              FinalAttr[[measure]]$combined$datasetSubGroup,
              levels = gtools::mixedsort(unique(FinalAttr[[measure]]$combined$datasetSubGroup)))
          }

          ggplotList <- list()
          ## Plot a multi-facet ggplot for all gtSigNames and all runs.
          {
            ## Generate a ggplot object based on FinalAttr[[measure]]$combined
            ggplotList$general <- ggplot2::ggplot(
              FinalAttr[[measure]]$combined,
              ggplot2::aes(x = .data$toolName, y = .data$value)) +
              ## Draw geom_violin and geom_quasirandom
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                #ggplot2::aes(fill = gtSigName),
                ## Maximize the violin plot width
                scale = "width",
                ## Make bandwidth larger
                #position = "dodge",
                #width = 1.2
                ## Hide outliers
                #outlier.shape = NA
              ) +
              #ggbeeswarm::geom_quasirandom(
              #  groupOnX = TRUE, size = 0.3
              #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
              #) +
              ## Show median of the Scaled Manhattan distance distribution
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Add title for general violin + beeswarm plot
              ggplot2::ggtitle(label = titles[measure],
                               subtitle = " between inferred and ground-truth exposures.") +
              ## Change axis titles
              ggplot2::labs(x = "Computational approach",
                            y = titles[measure]) +
              ## Rotate the names of tools,
              ## move axis.text.x right below the tick marks
              ## and remove legends
              ggplot2::theme(axis.text.x = ggplot2::element_text(
                ## Rotate the axis.text.x (names of tools),
                angle = 90,
                ## move axis.text.x right below the tick marks
                hjust = 1, vjust = 0.5),
                ## remove legends.
                legend.position = "none") +
              ## Split the plot into multiple facets,
              ## according to different gtSigNames
              ggplot2::facet_wrap(
                ggplot2::vars(gtSigName),
                ## Force facet_wrap to have 2 columns
                ncol = 2,
                scales = "free",
                ## Let facets be plotted vertically
                dir = "v"
              ) +
              ## Restrict the decimal numbers of values of measures to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0, max(FinalAttr[[measure]]$combined$value)),
                labels =function(x) sprintf("%.2f", x))
          }
          ## Plot a multi-facet ggplot,
          ## facets are separated by gtSigNames and datasetGroup
          ## (in example, it refers to slope.)
          if(!is.null(multiTools$datasetSubGroupName)) {
            bys <- c("datasetGroup","datasetSubGroup")
          } else {
            bys <- c("datasetGroup")
          }

          for(by in bys)  {

            ## The value of "datasetGroupName" or "datasetSubGroupName"
            ## which is the caption of "datasetGroup"
            byCaption <- eval(parse(
              text = paste0("multiTools$",by,"Name")))


            ## Generate a ggplot object based on FinalAttr[[measure]]$combined
            ggplotList[[by]] <- ggplot2::ggplot(
              FinalAttr[[measure]]$combined,
              ggplot2::aes(x = .data$toolName, y = .data$value))
            ## Draw geom_violin and geom_quasirandom
            ggplotList[[by]] <- ggplotList[[by]] +
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                #ggplot2::aes(fill = gtSigName),
                ## Maximize the violin plot width
                scale = "width"
                #,
                ## Make bandwidth larger
                #position = "dodge",
                #width = 1.2
                ## Hide outliers
                #outlier.shape = NA
              ) +
              #ggbeeswarm::geom_quasirandom(
              #  groupOnX = TRUE, size = 0.3
              #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
              #) +
              ## Show median of the Scaled Manhattan distance distribution
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Add title for general violin + beeswarm plot
              ggplot2::ggtitle(
                label = paste0(titles[measure]," as a function of "),
                subtitle = paste0("ground-truth signature names and ",byCaption,".")) +
              ## Change axis titles
              ggplot2::labs(x = "Computational approach",
                            y = titles[measure]) +
              ## Rotate the axis.text.x (names of tools),
              ## move axis.text.x right below the tick marks
              ## and remove legends
              ggplot2::theme(axis.text.x = ggplot2::element_text(
                ## Rotate the axis.text.x (names of tools),
                angle = 90,
                ## move axis.text.x right below the tick marks
                hjust = 1, vjust = 0.5),
                ## remove legends.
                legend.position = "none") +
              ## Split the plot into multiple facets,
              ## according to different gtSigNames
              ggplot2::facet_grid(rows =  ggplot2::vars(gtSigName),
                                  cols = eval(parse(text = paste0("ggplot2::vars(",by,")"))),
                                  scales = "free") +
              ## Restrict the decimal numbers of values of measures to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0, max(FinalAttr[[measure]]$combined$value)),
                labels =function(x) sprintf("%.2f", x))
          }

          ## Plot violin + beeswarm plots in pdf format
          grDevices::pdf(paste0(out.dir,"/",fileNames[measure],".",".pdf"), pointsize = 1)
          for(by in names(ggplotList)){
            print(ggplotList[[by]])
          }
          grDevices::dev.off()
        }
      }

    }




    FinalSummary <- list()
    FinalSummary$FinalExtr <- FinalExtr
    if(flagExposure) {
      FinalSummary$FinalAttr <- FinalAttr
    }

    save(FinalSummary,file = paste0(out.dir,"/FinalSummary.RDa"))

    invisible(FinalSummary)
  }


#' Combine results for multiple datasets, from one computational approaches.
#'
#' Summarize results from each computational approach in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#' @param dataset.dirs Paths of top-level dataset directories trees you want
#' to investigate.
#' E.g. "./S.0.1.Rsq.0.1"
#'
#' @param datasetGroup Numeric or character vector specifying the group
#' each dataset belong to.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider slope
#' (SBS1:SBS5 count ratio) as the group:
#' \code{c(0.1,0.5,1,2,5,10)}
#' Default: "Default"
#'
#' @param datasetGroupName Meaning or label of all datasetGroup.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider \code{"SBS1:SBS5 mutation count ratio"}
#' as the label of the \code{datasetGroup} slope.
#'
#' @param datasetSubGroup Numeric or character vector differentiating
#' datasets within each group.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider Pearson's R^2
#' as the subgroup:
#' c(0.1,0.2,0.3,0.6)
#' Default: Names of datasets, which are \code{basename(dataset.dirs)}
#'
#' @param datasetSubGroupName Meaning or label of all datasetSubGroup.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider \code{"Pearson's R squared"}
#' as the label of the \code{datasetSubGroup} Pearson's R^2.
#'
#' @param toolName Name of computational approach to be investigated
#' (e.g. "SigProExtractor")
#'
#' @param tool.dirname Name of the second.level.dir (e.g. "sp.sp"),
#' third.level.dir (e.g. "ExtrAttr") and tool.dir
#' (e.g. "SigProExtractor.results") to be investigated.
#'
#' One example: "sp.sp/ExtrAttr/SigProExtractor.results"
#'
#' Note: this function expects the summary generated by
#' \code{SummarizeSigOneSubdir} under \code{dataset.dirs}/\code{tool.dirname}
#'
#' @param out.dir Path of the output directory.
#'
#' @param overwrite Whether to overwrite the contents in out.dir if
#' it already exists. (Default: FALSE)
#'
#' @importFrom rlang .data
#'
#' @importFrom utils write.csv
#'
#' @export
#'
SummarizeOneToolMultiDatasets <-
  function(dataset.dirs,
           datasetGroup,
           datasetGroupName,
           datasetSubGroup = NULL,
           datasetSubGroupName = NULL,
           toolName,
           tool.dirname,
           out.dir,
           overwrite = FALSE){

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exists")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Wrap all datasets into one group, if datasetGroup is NULL.
    ## Re-order the dataset.group for better visualization of
    ## ggplot facets.
    {
      datasetNames <- basename(dataset.dirs)

      if(is.null(datasetGroup))
        datasetGroup <- rep("Default",length(dataset.dirs))
      datasetGroup <- factor(
        datasetGroup,
        levels = gtools::mixedsort(unique(datasetGroup)))
      names(datasetGroup) <- datasetNames

      if(!exists("datasetSubGroup"))
        datasetSubGroup <- datasetNames
      datasetSubGroup <- factor(
        datasetSubGroup,
        levels = gtools::mixedsort(unique(datasetSubGroup)))
      names(datasetSubGroup) <- datasetNames
    }

    ## Calculate summary tables for measures of extraction performance
    ## Need to calculate tables for All 6 measures
    {
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","PPV")
      indexLabels <- c("averCosSim" = "Average cosine similarity of all signatures",
                       "falseNeg" = "False negatives",
                       "falsePos" = "False positives",
                       "truePos" = "True positives",
                       "TPR" = "True positive rate (TPR, sensitivity)",
                       "PPV" = "Positive predictive value (PPV, precision)")
      subtitles <- c("averCosSim" = "",
                     "falseNeg" = "Number of missing ground-truth signatures",
                     "falsePos" = "Number of artefact signatures extracted, but different from ground-truth signatures",
                     "truePos" = "Number of extracted ground-truth signatures",
                     "TPR" = "True Positives / (True Positives + False Negatives)",
                     "PPV" = "True Positives / (True Positives + False Positives)")
      names(indexLabels) <- indexes
      names(subtitles) <- indexes
      indexNums <- length(indexes)

      ## Construct a summary list for storage
      OneToolSummary <- list()

      ## Combine each measurement for extraction performance for multiple datasets
      ## in multiple runs onto one summary table:
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
        toolName <- strsplit(basename(tool.dirname),".results")[[1]]
        ## Add multiRun <- NULL to please R check
        multiRun <- NULL
        load(paste0(thirdLevelDir,"/multiRun.RDa"))
        datasetName <- basename(datasetDir)
        for(index in indexes){
          if(FALSE){ ## debug
            measure4OneDataset <- data.frame(seed = names(multiRun[[index]]),
                                             index = index,
                                             value = multiRun[[index]],
                                             toolName = toolName,
                                             datasetName = datasetName,
                                             datasetGroup = datasetGroup[datasetName],
                                             datasetGroupName = datasetGroupName,
                                             datasetSubGroup = datasetSubGroup[datasetName],
                                             datasetSubGroupName = datasetSubGroupName,
                                             stringsAsFactors = FALSE)
          } else {
            measure4OneDataset <- data.frame(seed = names(multiRun[[index]]),
                                             value = multiRun[[index]],
                                             toolName = toolName,
                                             datasetGroup = datasetGroup[datasetName],
                                             datasetSubGroup = datasetSubGroup[datasetName],
                                             stringsAsFactors = FALSE)
          }
          rownames(measure4OneDataset) <- NULL

          ## Create a data.frame for each measure,
          ## and summarize multi-Run, multiDataset values
          ## for each measure.
          if(is.null(OneToolSummary[[index]])){
            OneToolSummary[[index]] <- data.frame()
          }
          OneToolSummary[[index]] <- rbind(OneToolSummary[[index]],measure4OneDataset)
        }
      }

      ## Calculate the stats (returned by summary()) of each extraction performance measure.
      OneToolSummary$stats <- list()
      for(index in indexes){
        currentStats <- summary(OneToolSummary[[index]][,"value"])
        OneToolSummary$stats[[index]] <- currentStats
      }

      ## For TPR (sensitivity), PPV and Number of False Negatives,
      ## calculate the proportion of 1.
      OneToolSummary$prop1 <- list()
      for(index in c("TPR","PPV","falseNeg")){
        currentProp <- length(which(OneToolSummary[[index]][,"value"] == 1)) / length(OneToolSummary[[index]][,"value"])
        OneToolSummary$prop1[[index]] <- currentProp
      }


    }

    ## Draw boxplot + beeswarm plot for extraction measures
    {
      ## Create a list to store ggplot2 boxplot + beeswarm plot objects
      ggplotList <- list()
      ## Plot a value~datasetSubGroup beeswarm for each measure.
      for(index in indexes){
        indexNum <- which(indexes == index)
        ## ggplot2::ggplot() sets coordinates
        ggplotList[[index]] <- ggplot2::ggplot(
          OneToolSummary[[index]],
          ## Make sure that only one x-label is shown in one small facet.
          #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
          ggplot2::aes(x = .data$toolName, y = .data$value)
        )
        ## Add facets
        ggplotList[[index]] <- ggplotList[[index]] +
          ggplot2::facet_grid(
            rows = ggplot2::vars(datasetSubGroup),
            cols = ggplot2::vars(datasetGroup),
            ## Move x facet labels to the right,
            ## This is to let the facet labels correspond to axis.title.
            switch = "x") +
          ## Draw boxplots and beeswarm plots on multi-facets.
          ## Draw geom_violin
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            ## Maximize the violin plot width
            scale = "width"
          ) +
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ## Draw beeswarm plot
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                       ## Make dot size smaller
                                       size = 0.3
                                       #,
                                       ## Remove differentiated colors for beeswarm dots
                                       ## Set groups for the filling functionalities to differentiate
                                       #ggplot2::aes(color = .data$datasetGroup)
          ) +
          ## Change filling color
          ggplot2::scale_fill_brewer(palette = "Greys") +
          ## Change titles
          ## and change axis titles.
          ## ggplot2::labs() has stronger function than ggplo2::ggtitle.
          ggplot2::labs(
            ## Add title for value~datasetSubGroup beeswarm plot,
            title = paste0(toolName,": ",indexLabels[index]),
            subtitle = subtitles[index],
            ## Change title of y axis (axis.title.y) into measure info (same as title)
            y = indexLabels[index],
            ## Change title of x axis to "Pearson's R squared"
            x = "Pearson's R squared") +
          ## Change title of legend to datasetGroupName
          ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
          ## Change axis.text and tickmarks
          ggplot2::theme(
            ## Remove axis.text.x
            axis.text.x = ggplot2::element_blank(),
            ## Remove tick marks on x axis (axis.ticks.x)
            axis.ticks.x = ggplot2::element_blank(),
            ## Remove entire legend
            legend.position = "none"
          ) +
          ## Restrict the decimal numbers of values of measures (y) to be 2
          ggplot2::scale_y_continuous(
            labels =function(x) sprintf("%.2f", x),
            ## Add a secondary axis title on the top of the plot
            ## Showing axis label indicating facets
            sec.axis = ggplot2::dup_axis(
              breaks = NULL, ## Don't show axis breaks
              labels = NULL, ## Don't show axis tickmarks
              name = "SBS1:SBS5 mutation count ratio")
          )
      }


      ## Output multiple extraction measures in a pdf file
      grDevices::pdf(paste0(out.dir,"/boxplot.onetool.extraction.measures.pdf"), pointsize = 1)
      for(index in indexes)
        suppressMessages(suppressWarnings(print(ggplotList[[index]])))
      grDevices::dev.off()
    }


    ## Summarize one-signature cosine similarity for one tool.
    {
      OneToolSummary$cosSim <- list()

      ## Combine one-signature cosine similarity from different spectra datasets
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
        toolName <- strsplit(basename(tool.dirname),".results")[[1]]
        ## Add multiRun <- NULL to please R check
        multiRun <- NULL
        load(paste0(thirdLevelDir,"/multiRun.RDa"))
        gtSigNames <- names(multiRun$cosSim)
        sigNums <- length(gtSigNames)
        datasetName <- basename(datasetDir)

        for(gtSigName in gtSigNames){

          if(FALSE){ # debug
            gtMeanCosSim4OneDataset <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                                                  gtSigName = gtSigName,
                                                  value = multiRun$cosSim[[gtSigName]],
                                                  toolName = toolName,
                                                  datasetName = datasetName,
                                                  datasetGroup = datasetGroup[datasetName],
                                                  datasetGroupName = datasetGroupName,
                                                  datasetSubGroup = datasetSubGroup[datasetName],
                                                  datasetSubGroupName = datasetSubGroupName,
                                                  stringsAsFactors = FALSE)
          } else {
            gtMeanCosSim4OneDataset <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                                                  value = multiRun$cosSim[[gtSigName]],
                                                  toolName = toolName,
                                                  datasetGroup = datasetGroup[datasetName],
                                                  datasetSubGroup = datasetSubGroup[datasetName],
                                                  stringsAsFactors = FALSE)
          }
          rownames(gtMeanCosSim4OneDataset) <- NULL

          ## Create a data.frame for each measure,
          ## and summarize multi-Run, multiDataset values
          ## for each measure.
          if(is.null(OneToolSummary$cosSim[[gtSigName]])){
            OneToolSummary$cosSim[[gtSigName]] <- data.frame()
          }
          OneToolSummary$cosSim[[gtSigName]] <- rbind(OneToolSummary$cosSim[[gtSigName]],gtMeanCosSim4OneDataset)
        }
      }

      ## Calculate the stats (returned by summary()) of one-signature cosine similarity.
      OneToolSummary$stats$cosSim <- list()
      for(gtSigName in gtSigNames){
        currentStats <- summary(OneToolSummary$cosSim[[gtSigName]][,"value"])
        OneToolSummary$stats$cosSim[[gtSigName]] <- currentStats
      }



      ## Combine multiple one-signature cosine similarity data.frame
      ## into OneToolSummary$cosSim$combined
      OneToolSummary$cosSim$combined <- data.frame()
      for(gtSigName in gtSigNames){
        gtMeanCosSim4AllDatasets <- data.frame(OneToolSummary$cosSim[[gtSigName]],
                                               stringsAsFactors = FALSE)
        rownames(gtMeanCosSim4AllDatasets) <- NULL

        if(nrow(OneToolSummary$cosSim$combined) == 0 |
           ncol(OneToolSummary$cosSim$combined) == 0 |
           is.null(dim(OneToolSummary$cosSim$combined)) ) {
          OneToolSummary$cosSim$combined <- gtMeanCosSim4AllDatasets
        } else {
          OneToolSummary$cosSim$combined <-
            rbind(OneToolSummary$cosSim$combined,gtMeanCosSim4AllDatasets)
        }
      }

    }
    ## Plot one-signature cosine similarity boxplot + beeswarm plot for one tool
    { ## debug
      ## Create a list to store ggplot2 boxplot + beeswarm plot objects
      ggplotList$cosSim <- list()
      ## Plot a value~datasetSubGroup beeswarm plot for each signature.
      for(gtSigName in gtSigNames){
        sigNum <- which(gtSigNames == gtSigName)
        ggplotList$cosSim[[gtSigName]] <- ggplot2::ggplot(
          OneToolSummary$cosSim[[gtSigName]],
          ## Make sure that only one x-label is shown in one small facet.
          #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
          ggplot2::aes(x = .data$toolName, y = .data$value)
        )
        ## Add facets
        ggplotList$cosSim[[gtSigName]] <- ggplotList$cosSim[[gtSigName]] +
          ggplot2::facet_grid(
            rows = ggplot2::vars(datasetSubGroup),
            cols = ggplot2::vars(datasetGroup),
            ## Move x facet labels to the right,
            ## This is to let the facet labels correspond to axis.title.
            switch = "x") +
          ## Draw beeswarm plots on multiple facets
          ## Draw geom_violin
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            ## Maximize the violin plot width
            scale = "width",
            ## Hide outliers
            #outlier.shape = NA
          ) +
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ## Draw beeswarm plot
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                       size = 0.3 ## Make dot size smaller
                                       #,
                                       ## Remove differentiated colors for beeswarm dots
                                       ## Set groups for the filling functionalities to differentiate
                                       #ggplot2::aes(color = .data$datasetGroup)
          ) +
          ## Change filling color
          ggplot2::scale_fill_brewer(palette = "Greys") +
          ## Change axis.text and tickmarks
          ggplot2::theme(
            ## Remove axis.text.x
            axis.text.x = ggplot2::element_blank(),
            ## Remove tick marks on x axis (axis.ticks.x)
            axis.ticks.x = ggplot2::element_blank(),
            ## Remove entire legend
            legend.position = "none"
          ) +
          ## Add titles
          ggplot2::labs(
            ## Add title for value~datasetSubGroup beeswarm plot
            title = paste0(toolName,": Average cosine similarity between signature ",gtSigName),
            subtitle = paste0("and all extracted signatures resembling ",gtSigName),
            ## Change title of y axis (axis.title.y) into gtSigName info (same as title)
            y = paste0("Cosine similarity to signature ",gtSigName),
            ## Change title of x axis to "Pearson's R squared"
            x = "Pearson's R squared") +
          ## Change title of legend to datasetGroupName
          ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
          ## Restrict the decimal numbers of values of measures (y) to be 2
          ggplot2::scale_y_continuous(
            ## For one-signature cosine similarity, set ylim from the minimum of Manhattan distance value to 1.
            limits = c(min(OneToolSummary$cosSim$combined$value),1),
            labels =function(x) sprintf("%.2f", x),
            ## Add a secondary axis title on the top of the plot
            ## Showing axis label indicating facets
            sec.axis = ggplot2::dup_axis(
              breaks = NULL, ## Don't show axis breaks
              labels = NULL, ## Don't show axis tickmarks
              name = "SBS1:SBS5 mutation count ratio"))
      }


      ## Output multiple extraction measures in a pdf file
      grDevices::pdf(paste0(out.dir,"/boxplot.onetool.onesig.cossim.pdf"), pointsize = 1)
      for(gtSigName in gtSigNames)
        suppressMessages(suppressWarnings(print(ggplotList$cosSim[[gtSigName]])))
      grDevices::dev.off()
    }


    ## Summarize scaled Manhattan distance only if
    ## scaled Manhattan distance data exists in object multRun
    exposureFlag <- TRUE
    {
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
        toolName <- strsplit(basename(tool.dirname),".results")[[1]]
        ## Add multiRun <- NULL to please R check
        multiRun <- NULL
        load(paste0(thirdLevelDir,"/multiRun.RDa"))
        if(is.null(multiRun$AggManhattanDist)){
          exposureFlag <- FALSE
          message("Skip summarizing scaled Manhattan distance...\n")
          break
        }
      }
    }

    ## Summarize aggregated scaled Manhattan distance for one tool.
    if(exposureFlag){
      ## Summarize aggregated scaled Manhattan distance for one tool.
      {
        OneToolSummary$AggManhattanDist <- list()

        for(datasetDir in dataset.dirs){
          thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
          toolName <- strsplit(basename(tool.dirname),".results")[[1]]
          ## Add multiRun <- NULL to please R check
          multiRun <- NULL
          load(paste0(thirdLevelDir,"/multiRun.RDa"))
          gtSigNames <- rownames(multiRun$AggManhattanDist)
          sigNums <- length(gtSigNames)
          datasetName <- basename(datasetDir)

          for(gtSigName in gtSigNames){
            if(FALSE) { # debug
              gtAggManhattanDist4OneDataset <- data.frame(seed = colnames(multiRun$AggManhattanDist),
                                                          gtSigName = gtSigName,
                                                          value = multiRun$AggManhattanDist[gtSigName,],
                                                          toolName = toolName,
                                                          datasetName = datasetName,
                                                          datasetGroup = datasetGroup[datasetName],
                                                          datasetGroupName = datasetGroupName,
                                                          datasetSubGroup = datasetSubGroup[datasetName],
                                                          datasetSubGroupName = datasetSubGroupName,
                                                          stringsAsFactors = FALSE)
            } else {
              gtAggManhattanDist4OneDataset <- data.frame(seed = colnames(multiRun$AggManhattanDist),
                                                          value = multiRun$AggManhattanDist[gtSigName,],
                                                          toolName = toolName,
                                                          datasetGroup = datasetGroup[datasetName],
                                                          datasetSubGroup = datasetSubGroup[datasetName],
                                                          stringsAsFactors = FALSE)
            }
            rownames(gtAggManhattanDist4OneDataset) <- NULL

            ## Create a data.frame for each measure,
            ## and summarize multi-Run, multiDataset values
            ## for each measure.
            if(is.null(OneToolSummary$AggManhattanDist[[gtSigName]])){
              OneToolSummary$AggManhattanDist[[gtSigName]] <- data.frame()
            }
            OneToolSummary$AggManhattanDist[[gtSigName]] <- rbind(OneToolSummary$AggManhattanDist[[gtSigName]],gtAggManhattanDist4OneDataset)
          }
        }

        ## Combine multiple ground-truth signature Manhattan-distance data.frame
        ## into OneToolSummary$AggManhattanDist$combined.
        OneToolSummary$AggManhattanDist$combined <- data.frame()
        for(gtSigName in gtSigNames){
          gtAggManhattanDist4AllDatasets <- data.frame(OneToolSummary$AggManhattanDist[[gtSigName]],
                                                       stringsAsFactors = FALSE)
          rownames(gtAggManhattanDist4AllDatasets) <- NULL

          if(nrow(OneToolSummary$AggManhattanDist$combined) == 0 |
             ncol(OneToolSummary$AggManhattanDist$combined) == 0 |
             is.null(dim(OneToolSummary$AggManhattanDist$combined)) ) {
            OneToolSummary$AggManhattanDist$combined <- gtAggManhattanDist4AllDatasets
          } else {
            OneToolSummary$AggManhattanDist$combined <-
              rbind(OneToolSummary$AggManhattanDist$combined,gtAggManhattanDist4AllDatasets)
          }
        }

      }
      ## Plot aggregated scaled Manhattan distance violin plot
      ## + beeswarm plot for one tool
      { ## debug
        ## Create a list to store ggplot2 boxplot + beeswarm plot objects
        ggplotList$AggManhattanDist <- list()
        ## Plot a value~datasetSubGroup beeswarm plot for each signature.
        for(gtSigName in gtSigNames){
          sigNum <- which(gtSigNames == gtSigName)
          ggplotList$AggManhattanDist[[gtSigName]] <- ggplot2::ggplot(
            OneToolSummary$AggManhattanDist[[gtSigName]],
            ## Make sure that only one x-label is shown in one small facet.
            #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
            ggplot2::aes(x = .data$toolName, y = .data$value)
          )
          ## Add facets
          ggplotList$AggManhattanDist[[gtSigName]] <- ggplotList$AggManhattanDist[[gtSigName]] +
            ggplot2::facet_grid(
              rows = ggplot2::vars(datasetSubGroup),
              cols = ggplot2::vars(datasetGroup),
              ## Move x facet labels to the right,
              ## This is to let the facet labels correspond to axis.title.
              switch = "x") +
            ## Draw beeswarm plots on multiple facets
            ## Draw geom_violin
            ggplot2::geom_violin(
              ## Change filling color to white
              fill = "#FFFFFF",
              ## Maximize the violin plot width
              scale = "width",
              ## Hide outliers
              #outlier.shape = NA
            ) +
            ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
            ## Show mean of the extraction meaasure distribution, as a blue diamond.
            ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
            ## Draw beeswarm plot
            ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                         size = 0.3 ## Make dot size smaller
                                         ,
                                         ## Remove differentiated colors for beeswarm dots
                                         ## Set groups for the filling functionalities to differentiate
                                         #ggplot2::aes(color = .data$datasetGroup)
            ) +
            ## Change filling color
            ggplot2::scale_fill_brewer(palette = "Greys") +
            ## Change axis.text and tickmarks
            ggplot2::theme(
              ## Remove axis.text.x
              axis.text.x = ggplot2::element_blank(),
              ## Remove tick marks on x axis (axis.ticks.x)
              axis.ticks.x = ggplot2::element_blank(),
              ## Remove entire legend
              legend.position = "none"
            ) +
            ## Change titles
            ggplot2::labs(
              ## Add title for value~datasetSubGroup beeswarm plot
              title = paste0(toolName,": Scaled Manhattan distance of ",gtSigName," exposure"),
              subtitle = "Between ground-truth exposure and inferred exposure",
              ## Change title of y axis (axis.title.y) same as gtSigName info (same as title)
              y = paste0("Scaled aggregated Manhattan distance of ",gtSigName," exposure"),
              ## Change title of x axis to "Pearson's R squared"
              x = "Pearson's R squared") +
            ## Change title of legend to datasetGroupName
            ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
            ## Restrict the decimal numbers of values of measures (y) to be 2
            ggplot2::scale_y_continuous(
              ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
              limits = c(0,max(OneToolSummary$AggManhattanDist$combined$value)),
              ## Restrict the decimal numbers of values of measures (y) to be 2
              labels =function(x) sprintf("%.2f", x),
              ## Add a secondary axis title on the top of the plot
              ## Showing axis label indicating facets
              sec.axis = ggplot2::dup_axis(
                breaks = NULL, ## Don't show axis breaks
                labels = NULL, ## Don't show axis tickmarks
                name = "SBS1:SBS5 mutation count ratio"))
        }


        ## Output multiple extraction measures in a pdf file
        grDevices::pdf(paste0(out.dir,"/boxplot.onetool.aggregated.Manhattan.dist.pdf"), pointsize = 1)
        for(gtSigName in gtSigNames)
          suppressMessages(suppressWarnings(print(ggplotList$AggManhattanDist[[gtSigName]])))
        grDevices::dev.off()
      }
    }

    if(TRUE){
      ## Summarize mean of scaled Manhattan distance
      ## separated for individual tumors for each tool.
      if(exposureFlag){
        {
          OneToolSummary$meanSepMD <- list()

          for(datasetDir in dataset.dirs){
            thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
            toolName <- strsplit(basename(tool.dirname),".results")[[1]]
            ## Add multiRun <- NULL to please R check
            multiRun <- NULL
            load(paste0(thirdLevelDir,"/multiRun.RDa"))
            gtSigNames <- rownames(multiRun$meanSepMD)
            sigNums <- length(gtSigNames)
            datasetName <- basename(datasetDir)

            for(gtSigName in gtSigNames){
              if(FALSE) { # debug
                gtmeanSepMD4OneDataset <- data.frame(seed = colnames(multiRun$meanSepMD),
                                                     gtSigName = gtSigName,
                                                     value = multiRun$meanSepMD[gtSigName,],
                                                     toolName = toolName,
                                                     datasetName = datasetName,
                                                     datasetGroup = datasetGroup[datasetName],
                                                     datasetGroupName = datasetGroupName,
                                                     datasetSubGroup = datasetSubGroup[datasetName],
                                                     datasetSubGroupName = datasetSubGroupName,
                                                     stringsAsFactors = FALSE)
              } else {
                gtmeanSepMD4OneDataset <- data.frame(seed = colnames(multiRun$meanSepMD),
                                                     value = multiRun$meanSepMD[gtSigName,],
                                                     toolName = toolName,
                                                     datasetGroup = datasetGroup[datasetName],
                                                     datasetSubGroup = datasetSubGroup[datasetName],
                                                     stringsAsFactors = FALSE)
              }
              rownames(gtmeanSepMD4OneDataset) <- NULL

              ## Create a data.frame for each measure,
              ## and summarize multi-Run, multiDataset values
              ## for each measure.
              if(is.null(OneToolSummary$meanSepMD[[gtSigName]])){
                OneToolSummary$meanSepMD[[gtSigName]] <- data.frame()
              }
              OneToolSummary$meanSepMD[[gtSigName]] <- rbind(OneToolSummary$meanSepMD[[gtSigName]],gtmeanSepMD4OneDataset)
            }
          }

          OneToolSummary$meanSepMD$combined <- data.frame()
          for(gtSigName in gtSigNames){
            gtmeanSepMD4AllDatasets <- data.frame(OneToolSummary$meanSepMD[[gtSigName]],
                                                  stringsAsFactors = FALSE)
            rownames(gtmeanSepMD4AllDatasets) <- NULL

            if(nrow(OneToolSummary$meanSepMD$combined) == 0 |
               ncol(OneToolSummary$meanSepMD$combined) == 0 |
               is.null(dim(OneToolSummary$meanSepMD$combined)) ) {
              OneToolSummary$meanSepMD$combined <- gtmeanSepMD4AllDatasets
            } else {
              OneToolSummary$meanSepMD$combined <-
                rbind(OneToolSummary$meanSepMD$combined,gtmeanSepMD4AllDatasets)
            }
          }

        }
        ## Plot mean of separated scaled Manhattan distance violin plot
        ## + beeswarm plot for one tool
        { ## debug
          ## Create a list to store ggplot2 boxplot + beeswarm plot objects
          ggplotList$meanSepMD <- list()
          ## Plot a value~datasetSubGroup beeswarm plot for each signature.
          for(gtSigName in gtSigNames){
            sigNum <- which(gtSigNames == gtSigName)
            ggplotList$meanSepMD[[gtSigName]] <- ggplot2::ggplot(
              OneToolSummary$meanSepMD[[gtSigName]],
              ## Make sure that only one x-label is shown in one small facet.
              #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
              ggplot2::aes(x = .data$toolName, y = .data$value)
            )
            ## Add facets
            ggplotList$meanSepMD[[gtSigName]] <- ggplotList$meanSepMD[[gtSigName]] +
              ggplot2::facet_grid(
                rows = ggplot2::vars(datasetSubGroup),
                cols = ggplot2::vars(datasetGroup),
                ## Move x facet labels to the right,
                ## This is to let the facet labels correspond to axis.title.
                switch = "x") +
              ## Draw beeswarm plots on multiple facets
              ## Draw geom_violin
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                ## Maximize the violin plot width
                scale = "width",
                ## Hide outliers
                #outlier.shape = NA
              ) +
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Draw beeswarm plot
              ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                           size = 0.3 ## Make dot size smaller
                                           ,
                                           ## Remove differentiated colors for beeswarm dots
                                           ## Set groups for the filling functionalities to differentiate
                                           #ggplot2::aes(color = .data$datasetGroup)
              ) +
              ## Change filling color
              ggplot2::scale_fill_brewer(palette = "Greys") +
              ## Change axis.text and tickmarks
              ggplot2::theme(
                ## Remove axis.text.x
                axis.text.x = ggplot2::element_blank(),
                ## Remove tick marks on x axis (axis.ticks.x)
                axis.ticks.x = ggplot2::element_blank(),
                ## Remove entire legend
                legend.position = "none"
              ) +
              ## Change titles
              ggplot2::labs(
                title = paste0(toolName,": Mean of scaled Manhattan distance of "),
                subtitle = paste0(gtSigName," exposure in individual tumors"),
                ## Change title of y axis (axis.title.y) same as gtSigName info (same as title)
                y = paste0("mean(separated Manhattan distance of ",gtSigName," exposure)"),
                ## Change title of x axis to "Pearson's R squared"
                x = "Pearson's R squared") +
              ## Change title of legend to datasetGroupName
              ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
              ## Restrict the decimal numbers of values of measures (y) to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0,max(OneToolSummary$meanSepMD$combined$value)),
                ## Restrict the decimal numbers of values of measures (y) to be 2
                labels =function(x) sprintf("%.2f", x),
                ## Add a secondary axis title on the top of the plot
                ## Showing axis label indicating facets
                sec.axis = ggplot2::dup_axis(
                  breaks = NULL, ## Don't show axis breaks
                  labels = NULL, ## Don't show axis tickmarks
                  name = "SBS1:SBS5 mutation count ratio"))
          }


          ## Output multiple extraction measures in a pdf file
          grDevices::pdf(paste0(out.dir,"/boxplot.onetool.mean.of.sep.Scaled.Manhattan.dist.pdf"), pointsize = 1)
          for(gtSigName in gtSigNames)
            suppressMessages(suppressWarnings(print(ggplotList$meanSepMD[[gtSigName]])))
          grDevices::dev.off()
        }
      }


      ## Summarize standard deviation of scaled Manhattan distance
      ## separated for individual tumors for each tool.
      if(exposureFlag){
        {
          OneToolSummary$sdSepMD <- list()

          for(datasetDir in dataset.dirs){
            thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
            toolName <- strsplit(basename(tool.dirname),".results")[[1]]
            ## Add multiRun <- NULL to please R check
            multiRun <- NULL
            load(paste0(thirdLevelDir,"/multiRun.RDa"))
            gtSigNames <- rownames(multiRun$sdSepMD)
            sigNums <- length(gtSigNames)
            datasetName <- basename(datasetDir)

            for(gtSigName in gtSigNames){
              if(FALSE) { # debug
                gtsdSepMD4OneDataset <- data.frame(seed = colnames(multiRun$sdSepMD),
                                                   gtSigName = gtSigName,
                                                   value = multiRun$sdSepMD[gtSigName,],
                                                   toolName = toolName,
                                                   datasetName = datasetName,
                                                   datasetGroup = datasetGroup[datasetName],
                                                   datasetGroupName = datasetGroupName,
                                                   datasetSubGroup = datasetSubGroup[datasetName],
                                                   datasetSubGroupName = datasetSubGroupName,
                                                   stringsAsFactors = FALSE)
              } else {
                gtsdSepMD4OneDataset <- data.frame(seed = colnames(multiRun$sdSepMD),
                                                   value = multiRun$sdSepMD[gtSigName,],
                                                   toolName = toolName,
                                                   datasetGroup = datasetGroup[datasetName],
                                                   datasetSubGroup = datasetSubGroup[datasetName],
                                                   stringsAsFactors = FALSE)
              }
              rownames(gtsdSepMD4OneDataset) <- NULL

              ## Create a data.frame for each measure,
              ## and summarize multi-Run, multiDataset values
              ## for each measure.
              if(is.null(OneToolSummary$sdSepMD[[gtSigName]])){
                OneToolSummary$sdSepMD[[gtSigName]] <- data.frame()
              }
              OneToolSummary$sdSepMD[[gtSigName]] <- rbind(OneToolSummary$sdSepMD[[gtSigName]],gtsdSepMD4OneDataset)
            }
          }

          OneToolSummary$sdSepMD$combined <- data.frame()
          for(gtSigName in gtSigNames){
            gtsdSepMD4AllDatasets <- data.frame(OneToolSummary$sdSepMD[[gtSigName]],
                                                stringsAsFactors = FALSE)
            rownames(gtsdSepMD4AllDatasets) <- NULL

            if(nrow(OneToolSummary$sdSepMD$combined) == 0 |
               ncol(OneToolSummary$sdSepMD$combined) == 0 |
               is.null(dim(OneToolSummary$sdSepMD$combined)) ) {
              OneToolSummary$sdSepMD$combined <- gtsdSepMD4AllDatasets
            } else {
              OneToolSummary$sdSepMD$combined <-
                rbind(OneToolSummary$sdSepMD$combined,gtsdSepMD4AllDatasets)
            }
          }

        }
        ## Plot standard deviation of separated scaled Manhattan distance
        ## violin plot + beeswarm plot for one tool
        { ## debug
          ## Create a list to store ggplot2 boxplot + beeswarm plot objects
          ggplotList$sdSepMD <- list()
          ## Plot a value~datasetSubGroup beeswarm plot for each signature.
          for(gtSigName in gtSigNames){
            sigNum <- which(gtSigNames == gtSigName)
            ggplotList$sdSepMD[[gtSigName]] <- ggplot2::ggplot(
              OneToolSummary$sdSepMD[[gtSigName]],
              ## Make sure that only one x-label is shown in one small facet.
              #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
              ggplot2::aes(x = .data$toolName, y = .data$value)
            )
            ## Add facets
            ggplotList$sdSepMD[[gtSigName]] <- ggplotList$sdSepMD[[gtSigName]] +
              ggplot2::facet_grid(
                rows = ggplot2::vars(datasetSubGroup),
                cols = ggplot2::vars(datasetGroup),
                ## Move x facet labels to the right,
                ## This is to let the facet labels correspond to axis.title.
                switch = "x") +
              ## Draw beeswarm plots on multiple facets
              ## Draw geom_violin
              ggplot2::geom_violin(
                ## Change filling color to white
                fill = "#FFFFFF",
                ## Maximize the violin plot width
                scale = "width",
                ## Hide outliers
                #outlier.shape = NA
              ) +
              ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
              ## Show mean of the extraction meaasure distribution, as a blue diamond.
              ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
              ## Draw beeswarm plot
              ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                           size = 0.3 ## Make dot size smaller
                                           ,
                                           ## Remove differentiated colors for beeswarm dots
                                           ## Set groups for the filling functionalities to differentiate
                                           #ggplot2::aes(color = .data$datasetGroup)
              ) +
              ## Change filling color
              ggplot2::scale_fill_brewer(palette = "Greys") +
              ## Change axis.text and tickmarks
              ggplot2::theme(
                ## Remove axis.text.x
                axis.text.x = ggplot2::element_blank(),
                ## Remove tick marks on x axis (axis.ticks.x)
                axis.ticks.x = ggplot2::element_blank(),
                ## Remove entire legend
                legend.position = "none"
              ) +
              ## Change titles
              ggplot2::labs(
                ## Add title for value~datasetSubGroup beeswarm plot
                title = paste0(toolName,": Standard deviation of scaled Manhattan distance of "),
                subtitle = paste0(gtSigName," exposure in individual tumors"),
                ## Change title of y axis (axis.title.y) same as gtSigName info (same as title)
                y = paste0("sd(separated Manhattan distance of ",gtSigName," exposure)"),
                ## Change title of x axis to "Pearson's R squared"
                x = "Pearson's R squared") +
              ## Change title of legend to datasetGroupName
              ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
              ## Restrict the decimal numbers of values of measures (y) to be 2
              ggplot2::scale_y_continuous(
                ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
                limits = c(0,max(OneToolSummary$sdSepMD$combined$value)),
                ## Restrict the decimal numbers of values of measures (y) to be 2
                labels =function(x) sprintf("%.2f", x),
                ## Add a secondary axis title on the top of the plot
                ## Showing axis label indicating facets
                sec.axis = ggplot2::dup_axis(
                  breaks = NULL, ## Don't show axis breaks
                  labels = NULL, ## Don't show axis tickmarks
                  name = "SBS1:SBS5 mutation count ratio"))
          }


          ## Output multiple extraction measures in a pdf file
          grDevices::pdf(paste0(out.dir,"/boxplot.onetool.stdev.of.sep.Scaled.Manhattan.dist.pdf"), pointsize = 1)
          for(gtSigName in gtSigNames)
            suppressMessages(suppressWarnings(print(ggplotList$sdSepMD[[gtSigName]])))
          grDevices::dev.off()
        }
      }


    }


    ## Write Summary tables for extraction measures
    for(index in indexes){
      output <- OneToolSummary[[index]]

      ## Change "value” to label of measure.
      colnames(output)[1] <- "Seed or run number"
      colnames(output)[2] <- indexLabels[index]
      colnames(output)[3] <- "Name of computational approach"
      colnames(output)[4] <- datasetGroupName
      colnames(output)[5] <- datasetSubGroupName

      write.csv(output,
                file = paste0(out.dir,"/",index,".csv"),
                quote = F, row.names = F)
    }

    ## Write Summary tables for signature cosine similarity.
    for(gtSigName in gtSigNames){
      output <- OneToolSummary[[index]]

      ## Change "value” to label of measure.
      colnames(output)[1] <- "Seed or run number"
      colnames(output)[2] <- paste0("Cosine similarity to ground-truth signature ",gtSigName)
      colnames(output)[3] <- "Name of computational approach"
      colnames(output)[4] <- datasetGroupName
      colnames(output)[5] <- datasetSubGroupName

      write.csv(output,
                file = paste0(out.dir,"/cossim.to.",gtSigName,".csv"),
                quote = F, row.names = F)
    }

    ## Write stat summary information into a text file.
    utils::capture.output(OneToolSummary$stats,file = paste0(out.dir,"/stats.txt"))
    utils::capture.output(OneToolSummary$prop1,file = paste0(out.dir,"/prop1.txt"))


    ## Add datasetGroupName and datasetSubGroupName into OneToolSummary
    OneToolSummary$datasetGroupName <- datasetGroupName
    OneToolSummary$datasetSubGroupName <- datasetSubGroupName

    save(OneToolSummary, file = paste0(out.dir,"/OneToolSummary.RDa"))
    invisible(OneToolSummary)
  }

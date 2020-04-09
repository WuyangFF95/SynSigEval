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


#' Assess/evaluate one result subfolder from a software package
#'
#' Note: For summarizing sigproextractor or SignatureAnalyzer,
#' users should use sigproextractor(SigProfiler-Python) v0.0.5.43+
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
#' @keywords internal
#'
#' @importFrom utils capture.output sessionInfo

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

        # Write results of exposure inference analysis
        write.csv(exposureDiff,
                  file = paste0(outputPath,"/exposureDifference.csv"),
                  quote = T)

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


#' Assess/evaluate multiple summarized runs for one dataset from one software package.
#'
#' Summarize results from each software package in \code{tool.dir}/\code{run.names}
#' (generated by running a software package),
#' combine them into \code{tool.dir}.
#'
#' @param datasetName Name of the dataset. (e.g. "S.0.1.Rsq.0.1").
#' Usually, it is has the same name as \code{basename(top.dir)}.
#'
#' @param toolName Name of software package. (e.g. "sigproextractor")
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
#' @return A list contain c(\code{mean},\code{sd}) of multiple runs:
#' Cosine similarity
#' True Positives(TP): Ground-truth signatures which are active in
#' the spectra, and extracted.
#' False Negatives(FN): Ground-truth signatures not extracted.
#' False Positives(FP): Signatures wrongly extracted, not resembling
#' any ground-truth signatures.
#' True positive rate (TPR, Sensitivity): TP / (TP + FN)
#' False discovery rate (FDR): FP / (FP + TP)
#'
#' @details Also writes multiple files into folder \code{tool.dir}:
#'
#'
#' @importFrom utils capture.output sessionInfo
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
                 "truePos","TPR","FDR")
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

      ## Concatenate TPR (True positive rate) and FDR (False discovery rate)
      currentTPR <- length(truePosNames) / length(gtSigNames)
      currentFDR <- length(falsePosNames) / (length(truePosNames) + length(falsePosNames))
      TPR <- c(TPR, currentTPR)
      FDR <- c(FDR, currentFDR)

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
    names(FDR) <- run.names
    for(gtSigName in gtSigNames)
      names(cosSim[[gtSigName]]) <- run.names


    multiRun <- list()
    ## Save name of the software package and the dataset.
    multiRun$datasetName <- datasetName
    multiRun$toolName <- toolName
    ## Save extraction indexes on multiple runs
    multiRun$averCosSim <- averCosSim
    multiRun$falseNeg <- falseNeg
    multiRun$falsePos <- falsePos
    multiRun$truePos <- truePos
    multiRun$TPR <- TPR
    multiRun$FDR <- FDR
    ## Save one-signature cosine similarity on multiple runs
    multiRun$cosSim <- cosSim



    ## Calculate mean and SD for indexes of signature extraction
    multiRun$meanSD <- matrix(nrow = 6, ncol = 2)
    indexes <- c("averCosSim","falseNeg","falsePos",
                 "truePos","TPR","FDR")
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
                 "truePos","TPR","FDR")
    rownames(multiRun$fivenum) <- indexes
    colnames(multiRun$fivenum) <- c("min","lower-hinge","median","upperhinge","max")
    for(index in indexes){
      currentFiveNum <- stats::fivenum(multiRun[[index]])
      multiRun$fivenum[index,] <- currentFiveNum
    }

    ## Plot boxplot + beeswarm plot for signature extraction
    if(FALSE){
      titles <- c("Average cosine similarity",
                  "False negatives",
                  "False positives",
                  "True positives",
                  "True positive rate (sensitivity)",
                  "False discovery rate (FDR)")
      subtitles <- c("","Number of ground-truth signatures not extracted",
                     "Number of signatures extracted, but different from ground-truth signatures",
                     "Number of ground-truth signatures extracted",
                     "True Positives / (True Positives + False Negatives)",
                     "False Positives / (True Positives + False Positives)")

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

    ## Indexes for exposure inference in multiple runs
    ManhattanDist <- matrix(nrow = length(gtSigNames), ncol = length(run.names))
    rownames(ManhattanDist) <- gtSigNames
    colnames(ManhattanDist) <- run.names
    for(runName in run.names){
      runDir <- paste0(tool.dir,"/",runName)
      summaryDir <- paste0(runDir,"/summary")
      exposureDiffFile <- paste0(summaryDir,"/exposureDiff.RDa")
      ## Add exposureDiff <- NULL to please the R check
      exposureDiff <- NULL
      load(file = exposureDiffFile)
      ManhattanDist[gtSigNames,runName] <- exposureDiff[gtSigNames,"Manhattan.distance"]
    }
    multiRun$ManhattanDist <- ManhattanDist

    ## Calculate mean and SD for indexes of exposure inference
    meanSDMD <- matrix(nrow = length(gtSigNames), ncol = 2)
    rownames(meanSDMD) <- gtSigNames
    colnames(meanSDMD) <- c("mean","stdev")
    for(gtSigName in gtSigNames){
      meanSDMD[gtSigName,"mean"] <- mean(ManhattanDist[gtSigName,])
      meanSDMD[gtSigName,"stdev"] <- stats::sd(ManhattanDist[gtSigName,])
    }
    multiRun$meanSDMD <- meanSDMD

    ## Calculate fivenums for exposure inference Scaled Manhattan distance
    multiRun$fivenumMD <- matrix(nrow = length(gtSigNames), ncol = 5)
    rownames(multiRun$fivenumMD) <- gtSigNames
    colnames(multiRun$fivenumMD) <- c("min","lower-hinge","median","upperhinge","max")
    for(gtSigName in gtSigNames){
      multiRun$fivenumMD[gtSigName,] <- stats::fivenum(ManhattanDist[gtSigName,])
    }

    ## Plot boxplot + beeswarm plot for exposure inference
    if(FALSE){
      ggplotList <- list()
      for(gtSigName in gtSigNames){
        ggplotList[[gtSigName]] <- ggplot2::ggplot(
          data.frame(value = ManhattanDist[gtSigName,],
                     gtSigName = gtSigName),
          ggplot2::aes(x = .data$gtSigName, y = .data$value))
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
          ggplot2::ggtitle(paste0("L1-difference of exposure of signature ",gtSigName))
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
          ggplot2::geom_boxplot() +
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE, size = 0.3) +
          ## Restrict the decimal numbers of values of indexes to be 2
          ggplot2::scale_y_continuous(labels =function(x) sprintf("%.2f", x))
      }


      ## Print high-resolution extraction indexes into one single png file
      if(FALSE){
        tempPlotList <- list()
        for(gtSigName in gtSigNames) {
          tempPlotList[[gtSigName]] <- ggplotList[[gtSigName]]
        }

        suppressMessages(
          ggplot2::ggsave(
            filename = paste0(tool.dir,"/boxplot.Manhattan.Dist.png"),
            plot = ggpubr::ggarrange(plotlist = tempPlotList),
            device = "png",
            dpi = 1000,
            limitsize = FALSE
          )
        )
      }

      ## Print extraction indexes into one pdf file
      grDevices::pdf(paste0(tool.dir,"/boxplot.attribution.measures.pdf"), pointsize = 1)
      for(gtSigName in gtSigNames) print(ggplotList[[gtSigName]])
      grDevices::dev.off()
    }


    ## Save data and results
    save(multiRun,file = paste0(tool.dir,"/multiRun.RDa"))
    write.csv(x = multiRun$ManhattanDist,
              file = paste0(tool.dir,"/ManhattanDist.csv"))
    write.csv(x = multiRun$meanSD,
              file = paste0(tool.dir,"/meanSD.csv"))
    write.csv(x = multiRun$meanSDMD,
              file = paste0(tool.dir,"/meanSD.Manhattan.dist.csv"))
    write.csv(x = multiRun$fivenum,
              file = paste0(tool.dir,"/fivenum.csv"))
    write.csv(x = multiRun$fivenumMD,
              file = paste0(tool.dir,"/fivenum.Manhattan.dist.csv"))
    invisible(multiRun)
  }



#' Combine results for a single dataset, from different software packages.
#'
#' Summarize results from each software package in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#'
#' @param third.level.dir Third level path distinguishing de-novo extraction
#' + attribution packages from attribution-only packages.
#' Examples:
#' \code{top.dir}/sp.sp/ExtrAttr/
#' \code{top.dir}/sa.sa/Attr/
#'
#' @param toolNames Names of software package. (e.g. "sigproextractor")
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
#' @param datasetGroup Numeric or character vector specifying the group
#' each dataset belong to.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider slope as the group:
#' c("slope=0.1","slope=0.5","slope=1","slope=2","slope=5","slope=10")
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
#' c("Rsq=0.1","Rsq=0.2","Rsq=0.3","Rsq=0.6")
#' Default: Names of datasets, which are \code{basename(dataset.dirs)}
#'
#' @param datasetSubGroupName Meaning or label of all datasetSubGroup.
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
#' False discovery rate (FDR): FP / (FP + TP)
#'
#' @export
#'
#' @importFrom utils capture.output sessionInfo
SummarizeMultiToolsOneDataset <- function(
  third.level.dir,
  toolNames,
  tool.dirnames,
  datasetGroup,
  datasetGroupName,
  datasetSubGroup,
  datasetSubGroupName){

  multiTools <- list()
  combMeanSD <- NULL
  combMeanSDMD <- NULL

  for(toolNumber in 1:length(toolNames)){
    toolName <- toolNames[toolNumber]
    toolDirName <- tool.dirnames[toolNumber]
    toolPath <- paste0(third.level.dir,"/",toolDirName)
    ## Add multiRun <- NULL to please the R check
    multiRun <- NULL
    datasetName <- NULL
    load(paste0(toolPath,"/multiRun.RDa"))

    ## Combine multi-runs and multi-tools for each index
    {
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "Number of False negatives",
                       "Number of False positives",
                       "Number of True positives",
                       "True positive rate",
                       "False discovery rate")
      for(index in indexes){
        indexNum <- which(index == indexes)
        tmp <- data.frame(seed = names(multiRun[[index]]),
                          index = index,
                          indexLabel = indexLabels[indexNum],
                          value = multiRun[[index]],
                          toolName = toolName,
                          datasetName = multiRun$datasetName,
                          datasetGroup = datasetGroup,
                          datasetGroupName = datasetGroupName,
                          datasetSubGroup = datasetSubGroup,
                          datasetSubGroupName = datasetSubGroupName,
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL
        ## Create a data.frame for each index,
        ## and summarize multi-Run, multiDataset values
        ## for each index.
        if(is.null(multiTools[[index]])){
          multiTools[[index]] <- data.frame()
        }
        multiTools[[index]] <- rbind(multiTools[[index]],tmp)
      }
    }

    ## Combine multi-runs and multi-tools for
    ## one-signature cosine similarity.
    {
      gtSigNames <- rownames(multiRun$ManhattanDist)
      multiTools$gtSigNames <- gtSigNames
      if(is.null(multiTools$cosSim)) multiTools$cosSim <- list()

      for(gtSigName in gtSigNames){
        tmp <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                          gtSigName = gtSigName,
                          label = paste0("Cosine similarity to signature ",gtSigName),
                          value = multiRun$cosSim[[gtSigName]],
                          toolName = toolName,
                          datasetName = multiRun$datasetName,
                          datasetGroup = datasetGroup,
                          datasetGroupName = datasetGroupName,
                          datasetSubGroup = datasetSubGroup,
                          datasetSubGroupName = datasetSubGroupName,
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL
        ## Create a data.frame for each ground-truth signature,
        ## and summarize multi-Run, multiDataset values
        ## for each ground-truth signature.
        if(is.null(multiTools$cosSim[[gtSigName]])){
          multiTools$cosSim[[gtSigName]] <- data.frame()
        }
        multiTools$cosSim[[gtSigName]] <- rbind(multiTools$cosSim[[gtSigName]],tmp)
      }
    }


    ## Combine multi-runs and multi-tools for Manhattan
    ## distance of each ground-truth signature
    {
      if(is.null(multiTools$ManhattanDist)) multiTools$ManhattanDist <- list()
      for(gtSigName in gtSigNames){
        tmp <- data.frame(seed = colnames(multiRun$ManhattanDist),
                          gtSigName = gtSigName,
                          value = multiRun$ManhattanDist[gtSigName,],
                          toolName = toolName,
                          datasetName = multiRun$datasetName,
                          datasetGroup = datasetGroup,
                          datasetGroupName = datasetGroupName,
                          datasetSubGroup = datasetSubGroup,
                          datasetSubGroupName = datasetSubGroupName,
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL
        ## Create a data.frame for each ground-truth signature,
        ## and summarize multi-Run, multiDataset values
        ## for each ground-truth signature.
        if(is.null(multiTools$ManhattanDist[[gtSigName]])){
          multiTools$ManhattanDist[[gtSigName]] <- data.frame()
        }
        multiTools$ManhattanDist[[gtSigName]] <- rbind(multiTools$ManhattanDist[[gtSigName]],tmp)
      }
    }

    ## meanSDMD contains mean and standard deviation
    ## for each extraction index.
    {
      meanSD <- multiRun$meanSD
      colnames(meanSD) <- paste0(toolDirName,".", colnames(meanSD))
      if(is.null(meanSD)){
        combMeanSD <- meanSD
      } else{
        combMeanSD <- cbind(combMeanSD,meanSD)
      }
    }

    ## meanSDMD contains mean and standard deviation
    ## for Scaled Manhattan distance between ground-truth exposures
    ## and inferred exposures for each ground-truth signature
    {
      meanSDMD <- multiRun$meanSDMD
      colnames(meanSDMD) <- paste0(toolDirName,".", colnames(meanSDMD))
      if(is.null(meanSDMD)){
        combMeanSDMD <- meanSDMD
      } else{
        combMeanSDMD <- cbind(combMeanSDMD,meanSDMD)
      }
    }
  }

  multiTools$combMeanSD <- combMeanSD
  multiTools$combMeanSDMD <- combMeanSDMD
  multiTools$datasetGroupName <- datasetGroupName
  multiTools$datasetSubGroupName <- datasetSubGroupName

  save(multiTools,file = paste0(third.level.dir,"/multiTools.RDa"))
  write.csv(x = multiTools$combMeanSD,
            file = paste0(third.level.dir,"/combined.meanSD.csv"))
  write.csv(x = multiTools$combMeanSDMD,
            file = paste0(third.level.dir,"/combined.meanSD.Manhattan.dist.csv"))
  invisible(multiTools)
}

#' Combine results for multiple datasets, from different software packages.
#'
#' Summarize results from each software package in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
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

    ## For each index,
    ## Create a data.frame integrating results of
    ## all runs and for all datasets
    {
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "False negatives",
                       "False positives",
                       "True positives",
                       "True positive rate (TPR, sensitivity)",
                       "False discovery rate (FDR)")
      indexNums <- length(indexes)
    }

    ## Summarizing different measures for extraction performance
    ## separately into list "FinalExtr".
    ## Showing individual values rather than
    ## only showing mean and standard deviation of multiple runs
    {
      FinalExtr <- list()
      toolNames <- character(0)

      ## Combine extraction assessment onto 7 sheets:
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        ## Find tool names
        toolNames <- unique(multiTools[["averCosSim"]][,"toolName"])

        if(length(FinalExtr) == 0){
          for(index in indexes) {
            FinalExtr[[index]] <- data.frame()
          }
        }
        current <- list()
        for(index in indexes){
          current[[index]] <- multiTools[[index]]
          #rownames(current[[index]]) <- datasetDir
          FinalExtr[[index]] <- rbind(FinalExtr[[index]],current[[index]])
        }
      }
      for(index in indexes){
        write.csv(FinalExtr[[index]],
                  file = paste0(out.dir,"/",index,".csv"))
      }
    }


    ## For each extraction measures,
    ## merge values from multiple runs
    ## into one data.frame plotDFList$<measure_name>
    ## and software approaches for easier plotting.
    {
      plotDFList <- list()

      ## The first index to plot are truePos and FDR.
      indexes <- c("truePos","FDR")
      indexLabels<- c("Number of true positives",
                      "False discovery rate")
      for(index in indexes){
        plotDFList[[index]] <- data.frame()
      }

      ## For each dataset, combine the index values into plotDFList[[index]]
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        for(index in indexes){
          plotDFList[[index]] <- rbind(plotDFList[[index]],multiTools[[index]])
        }

        ## The next to plot are one-signature cosine similarity into plotDFList.
        ## In SBS1-SBS5-correlation study, refer to average SBS1 and average SBS5 similarity.
        gtSigNames <- names(multiTools$cosSim)
        for(gtSigName in gtSigNames){
          if(is.null(plotDFList[[gtSigName]])){
            plotDFList[[gtSigName]] <- data.frame()
          }
          plotDFList[[gtSigName]] <- rbind(plotDFList[[gtSigName]],multiTools$cosSim[[gtSigName]])
        }

        ## The third to plot is the overall one-dimension measure
        ## equals to:
        ## Average of (1-FDR) in 20 runs +
        ## Average of one-signature cosine similarity in 20 runs +
        ## Proportion of True Positives = Total number of ground-truth signatures.
        if(is.null(plotDFList[["compositeMeasure"]])){
          plotDFList[["compositeMeasure"]] <- data.frame()
        }

        ## Calculate values for each tool/software package.
        toolNames <- unique(multiTools$averCosSim$toolName)
        for(toolName in toolNames){

          ## Calculate the value for mean(1 - FDR) for all runs of each tool on this dataset.
          rowNum <- which(multiTools$FDR[,"toolName"] == toolName)
          currentAver1MinusFDR <- 1 - mean(multiTools$FDR[rowNum,"value"])

          ## Calculate the value for mean(one signature cosine similarity) for all runs of each tool on this dataset.
          currentAverOneSigCosSim <- list()
          for(gtSigName in gtSigNames){
            rowNum <- which(multiTools$cosSim[[gtSigName]][,"toolName"] == toolName)
            currentAverOneSigCosSim[[gtSigName]] <- mean(multiTools$cosSim[[gtSigName]][rowNum,"value"])
          }

          if(FALSE){ ## Considering whether to calculate all-extracted-proportion, or the mean of .
            ## Calculate the value for mean(one signature cosine similarity) for all runs of each tool on this dataset.
            rowNum <- which(multiTools$truePos[,"toolName"] == toolName)
            currentPropAllExtracted <-
              sum(multiTools$truePos[rowNum,"value"] == length(multiTools$gtSigNames)) / length(rowNum)
          } else {
            ## Calculate the value for mean(truePos) / total number of ground-truth signatures for all runs of each tool on this dataset.
            rowNum <- which(multiTools$truePos[,"toolName"] == toolName)
            currentAverTruePosProp <-   mean(multiTools$truePos[rowNum,"value"]) / length(multiTools$gtSigNames)

          }

          ## Sum up all the measurements, as the overall measurement
          if(FALSE){
            currentMeasureValue <- sum(c(currentAver1MinusFDR , unlist(currentAverOneSigCosSim) , currentPropAllExtracted))
          } else{
            currentMeasureValue <- sum(c(currentAver1MinusFDR , unlist(currentAverOneSigCosSim) , currentAverTruePosProp))
          }

          ## Add current composite measure to the plotDFList$compositeMeasure
          currentMeasureDF <- data.frame(
            "seed" = "all",
            "index" = "compositeMeasure",
            "indexLabel" = "Composite measure",
            "value" = currentMeasureValue,
            "toolName" = toolName,
            "datasetName" = multiTools$averCosSim$datasetName[1],
            "datasetGroup" = multiTools$averCosSim$datasetGroup[1],
            "datasetGroupName" = multiTools$datasetGroupName,
            "datasetSubGroup" = multiTools$averCosSim$datasetSubGroup[1],
            "datasetSubGroupName" = multiTools$datasetSubGroupName)

          plotDFList$compositeMeasure <- rbind(plotDFList$compositeMeasure,currentMeasureDF)
        }
      }

      ## Change column names of plotDFList[[gtSigName]] (one-signature cosine similarity DFs)
      for(gtSigName in gtSigNames){
        colnames(plotDFList[[gtSigName]])[2] <- "index"
        colnames(plotDFList[[gtSigName]])[3] <- "indexLabel"
      }

      ## Combine all extraction measurements, plotDFList[[index]] into plotDFList$Combined
      ## combined all one-signature cosine similarity, plotDFList[[gtSigName]] into plotDFList$Combined
      plotDFList$combined <- data.frame()
      for(index in c(indexes,gtSigNames,"compositeMeasure")){
        plotDFList$combined <- rbind(plotDFList$combined,plotDFList[[index]])
      }

      ## Convert plotDFList$combined$datasetGroup and
      ## Let their levels follow gtools::mixedsort() fashion
      ## So that the order of the facet labels will be more reasonable for readers.
      plotDFList$combined$datasetGroup <- factor(
        plotDFList$combined$datasetGroup,
        levels = gtools::mixedsort(unique(plotDFList$combined$datasetGroup)))

      plotDFList$combined$datasetSubGroup <- factor(
        plotDFList$combined$datasetSubGroup,
        levels = gtools::mixedsort(unique(plotDFList$combined$datasetSubGroup)))
    }



    ## Plot general pdf for all extraction measures
    ## Plot a general violin + beeswarm plot for multiple measures
    ## in all runs and in all datasets.
    {

      ## Plot a multi-facet ggplot for all measures and all runs.
      {
        ggplotList <- list()
        ## Generate a ggplot object based on plotDFList$combined
        ggplotList$general <- ggplot2::ggplot(
          plotDFList$combined,
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
          ggplot2::labs(x = "Software package") +
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
      for(by in c("datasetGroup","datasetSubGroup"))  {

        ## The value of "datasetGroupName" or "datasetSubGroupName"
        ## which is the caption of "datasetGroup"
        byCaption <- eval(parse(text = paste0("multiTools$",by,"Name")))

        ## Generate a ggplot object based on plotDFList$combined
        ggplotList[[by]] <- ggplot2::ggplot(
          plotDFList$combined,
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
          ggplot2::labs(x = "Software package") +
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



    ## Write a table for extraction measures.
    if(FALSE){ ## Remove redundant measures
      indexes <- c("averCosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "False negatives",
                       "False positives",
                       "True positives",
                       "True positive Rate (TPR, sensitivity)",
                       "False discovery Rate (FDR)")
    } else{
      indexes <- c("averCosSim","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "False discovery rate (FDR)")
    }


    ## Summarizing one-signature cosine similarity results.
    {
      FinalExtr$cosSim <- list()
      ## Combine assessment onto multiple sheets.
      ## Each sheet shows cosine similarity to one mutational signature.
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        gtSigNames <- rownames(multiTools$combMeanSDMD)
        sigNums <- length(gtSigNames)

        if(length(FinalExtr$cosSim) == 0){
          for(gtSigName in gtSigNames) {
            FinalExtr$cosSim[[gtSigName]] <- data.frame()
          }
        }

        ## Combine one-signature cosine similarity
        current <- list()
        for(gtSigName in gtSigNames){
          current[[gtSigName]] <- multiTools$cosSim[[gtSigName]]
          FinalExtr$cosSim[[gtSigName]] <- rbind(FinalExtr$cosSim[[gtSigName]],current[[gtSigName]])
        }
      }
      for(gtSigName in gtSigNames){
        write.csv(FinalExtr$cosSim[[gtSigName]],
                  file = paste0(out.dir,"/onesig.cossim.",gtSigName,".csv"))
      }
    }


    ## Plot general png and pdf for one-signature cosine similarity summary
    ## Plot a general violin + beeswarm plot for multiple signatures
    ## in all runs and in all datasets.
    {
      plotDFList <- list()

      ## For ground-truth signature,
      ## Create a data.frame integrating results of
      ## all runs and for all datasets
      gtSigNames <- multiTools$gtSigNames
      sigNums <- length(gtSigNames)

      for(gtSigName in gtSigNames){
        plotDFList[[gtSigName]] <- data.frame()
      }

      ## For each dataset, combine the gtSigName values into plotDFList[[gtSigName]]
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        for(gtSigName in gtSigNames){
          plotDFList[[gtSigName]] <- rbind(plotDFList[[gtSigName]],multiTools$cosSim[[gtSigName]])
        }
      }

      ## Combine all plotDFList[[gtSigName]] into plotDFList$Combined
      plotDFList$combined <- data.frame()
      for(gtSigName in gtSigNames){
        plotDFList$combined <- rbind(plotDFList$combined,plotDFList[[gtSigName]])
      }

      ggplotList <- list()
      ## Plot a multi-facet ggplot for all gtSigNames and all runs.
      {
        ## Generate a ggplot object based on plotDFList$combined
        ggplotList$general <- ggplot2::ggplot(
          plotDFList$combined,
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
                           subtitle = "for all software packages, ratios and correlation values.") +
          ## Change axis titles
          ggplot2::labs(x = "Software package",
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
            limits = c(min(plotDFList$combined$value), 1),
            labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a multi-facet ggplot,
      ## facets are separated by gtSigNames and datasetGroup
      ## (in example, it refers to slope.)
      for(by in c("datasetGroup","datasetSubGroup"))  {

        ## The value of "datasetGroupName" or "datasetSubGroupName"
        ## which is the caption of "datasetGroup"
        byCaption <- eval(parse(
          text = paste0("multiTools$",by,"Name")))

        ## Generate a ggplot object based on plotDFList$combined
        ggplotList[[by]] <- ggplot2::ggplot(
          plotDFList$combined,
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
          ggplot2::labs(x = "Software package",
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
            limits = c(min(plotDFList$combined$value), 1),
            labels = function(x) sprintf("%.2f", x))
      }

      ## Plot violin + beeswarm plots in pdf format
      grDevices::pdf(paste0(out.dir,"/onesig.cossim.violins.pdf"), pointsize = 1)
      for(by in names(ggplotList)){
        print(ggplotList[[by]])
      }
      grDevices::dev.off()
    }


    ## Summarizing attribution Scaled Manhattan distance results.
    {
      FinalAttr <- list()
      ## Combine attribution assessment onto multiple sheets.
      ## Each sheet shows Scaled Manhattan distance for one mutational signature.
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        gtSigNames <- rownames(multiTools$combMeanSDMD)
        sigNums <- length(gtSigNames)

        if(length(FinalAttr) == 0){
          for(gtSigName in gtSigNames) {
            FinalAttr[[gtSigName]] <- data.frame()
          }
        }

        ## Combine Scaled Manhattan distance
        current <- list()
        for(gtSigName in gtSigNames){
          current[[gtSigName]] <- multiTools$combMeanSDMD[gtSigName,,drop = F]
          FinalAttr[[gtSigName]] <- rbind(FinalAttr[[gtSigName]],current[[gtSigName]])
        }
      }
      for(gtSigName in gtSigNames){
        write.csv(FinalAttr[[gtSigName]],
                  file = paste0(out.dir,"/ManhattanDist.",gtSigName,".csv"))
      }
    }
    ## Plot general png and pdf for attribution Scaled Manhattan distance summary
    ## Plot a general violin + beeswarm plot for multiple signatures
    ## in all runs and in all datasets.
    {
      plotDFList <- list()

      ## For ground-truth signature,
      ## Create a data.frame integrating results of
      ## all runs and for all datasets
      gtSigNames <- multiTools$gtSigNames
      sigNums <- length(gtSigNames)

      for(gtSigName in gtSigNames){
        plotDFList[[gtSigName]] <- data.frame()
      }

      ## For each dataset, combine the gtSigName values into plotDFList[[gtSigName]]
      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
        ## Add multiTools <- NULL to please R check
        multiTools <- NULL
        load(paste0(thirdLevelDir,"/multiTools.RDa"))

        for(gtSigName in gtSigNames){
          plotDFList[[gtSigName]] <- rbind(plotDFList[[gtSigName]],multiTools$ManhattanDist[[gtSigName]])
        }
      }

      ## Combine all plotDFList[[gtSigName]] into plotDFList$Combined
      plotDFList$combined <- data.frame()
      for(gtSigName in gtSigNames){
        plotDFList$combined <- rbind(plotDFList$combined,plotDFList[[gtSigName]])
      }

      ggplotList <- list()
      ## Plot a multi-facet ggplot for all gtSigNames and all runs.
      {
        ## Generate a ggplot object based on plotDFList$combined
        ggplotList$general <- ggplot2::ggplot(
          plotDFList$combined,
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
          ggplot2::ggtitle(label = "Scaled Manhattan distance between inferred and grond-truth exposures",
                           subtitle = "for all software packages, ratios and correlation values.") +
          ## Change axis titles
          ggplot2::labs(x = "Software package",
                        y = "Scaled Manhattan distance") +
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
            limits = c(0, max(plotDFList$combined$value)),
            labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a multi-facet ggplot,
      ## facets are separated by gtSigNames and datasetGroup
      ## (in example, it refers to slope.)
      for(by in c("datasetGroup","datasetSubGroup"))  {

        ## The value of "datasetGroupName" or "datasetSubGroupName"
        ## which is the caption of "datasetGroup"
        byCaption <- eval(parse(
          text = paste0("multiTools$",by,"Name")))


        ## Generate a ggplot object based on plotDFList$combined
        ggplotList[[by]] <- ggplot2::ggplot(
          plotDFList$combined,
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
            label = paste0("Scaled Manhattan distance summary plot as a function of "),
            subtitle = paste0("ground-truth signature names and ",byCaption,".")) +
          ## Change axis titles
          ggplot2::labs(x = "Software package",
                        y = "Scaled Manhattan distance") +
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
            limits = c(0, max(plotDFList$combined$value)),
            labels =function(x) sprintf("%.2f", x))
      }

      ## Plot violin + beeswarm plots in pdf format
      grDevices::pdf(paste0(out.dir,"/Manhattan.Dist.violins.pdf"), pointsize = 1)
      for(by in names(ggplotList)){
        print(ggplotList[[by]])
      }
      grDevices::dev.off()
    }

    FinalSummary <- list()
    FinalSummary$FinalExtr <- FinalExtr
    FinalSummary$FinalAttr <- FinalAttr

    save(FinalSummary,file = paste0(out.dir,"/FinalSummary.RDa"))

    invisible(FinalSummary)
  }


#' Combine results for multiple datasets, from one software packages.
#'
#' Summarize results from each software package in \code{third.level.dir}/\code{tool.dirnames}
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
#' @param toolName Name of software package to be investigated
#' (e.g. "sigproextractor")
#'
#' @param tool.dirname Name of the second.level.dir (e.g. "sp.sp"),
#' third.level.dir (e.g. "ExtrAttr") and tool.dir
#' (e.g. "sigproextractor.results") to be investigated.
#'
#' One example: "sp.sp/ExtrAttr/sigproextractor.results"
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
#' @export
#'
SummarizeOneToolMultiDatasets <-
  function(dataset.dirs,
           datasetGroup = NULL,
           datasetGroupName,
           datasetSubGroup = NULL,
           datasetSubGroupName,
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

      if(is.null(datasetSubGroup))
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
                   "truePos","TPR","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "False negatives",
                       "False positives",
                       "True positives",
                       "True positive rate (TPR, sensitivity)",
                       "False discovery rate (FDR)")
      subtitles <- c("",
                     "Number of missing ground-truth signatures",
                     "Number of artefact signatures extracted, but different from ground-truth signatures",
                     "Number of extracted ground-truth signatures",
                     "True Positives / (True Positives + False Negatives)",
                     "False Positives / (True Positives + False Positives)")
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
        ## Omit redundant measures
        if(FALSE){
          indexes <- rownames(multiRun$meanSD)
          indexNums <- length(indexes)
        }
        datasetName <- basename(datasetDir)
        for(index in indexes){
          tmp <- data.frame(seed = names(multiRun[[index]]),
                            index = index,
                            value = multiRun[[index]],
                            toolName = toolName,
                            datasetName = datasetName,
                            datasetGroup = datasetGroup[datasetName],
                            datasetGroupName = datasetGroupName,
                            datasetSubGroup = datasetSubGroup[datasetName],
                            datasetSubGroupName = datasetSubGroupName,
                            stringsAsFactors = FALSE)

          rownames(tmp) <- NULL

          ## Create a data.frame for each measure,
          ## and summarize multi-Run, multiDataset values
          ## for each measure.
          if(is.null(OneToolSummary[[index]])){
            OneToolSummary[[index]] <- data.frame()
          }
          OneToolSummary[[index]] <- rbind(OneToolSummary[[index]],tmp)
        }
      }

      ## Combine tables of different measures for extraction performance
      ## into OneToolSummary$extraction.
      OneToolSummary[["extraction"]] <- data.frame()
      for(index in indexes){
        tmp <- data.frame(OneToolSummary[[index]],
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL

        if(nrow(OneToolSummary[["extraction"]]) == 0 |
           ncol(OneToolSummary[["extraction"]]) == 0 |
           is.null(dim(OneToolSummary[["extraction"]])) ) {
          OneToolSummary[["extraction"]] <- tmp
        } else {
          OneToolSummary[["extraction"]] <-
            rbind(OneToolSummary[["extraction"]],tmp)
        }
      }

      ## Calculate the stats of each extraction performance measure.
      OneToolSummary$stats <- list()
      for(index in indexes){
        currentStats <- summary(OneToolSummary[[index]][,"value"])
        OneToolSummary$stats[[index]] <- currentStats
      }

      ## For TPR (sensitivity) and Number of False Negatives,
      ## calculate the proportion of 1.
      OneToolSummary$prop1 <- list()
      for(index in c("TPR","falseNeg")){
        currentProp <- length(which(OneToolSummary[[index]][,"value"] == 1)) / length(OneToolSummary[[index]][,"value"])
        OneToolSummary$prop1[[index]] <- currentProp
      }


    }

    ## Draw boxplot + beeswarm plot for extraction measures
    {
      ## Only average cosine similarity, one-signature cosine similarity and FDR are plotted.
      indexes <- c("averCosSim","FDR")
      indexLabels <- c("Average cosine similarity of all signatures",
                       "False discovery rate (FDR)")
      subtitles <- c("",
                     "False Positives / (True Positives + False Positives)")

      ## Designate titles and subtitles for each page
      titles <- indexLabels
      names(titles) <- indexes
      names(subtitles) <- indexes
      ylabels <- titles


      ## Create a list to store ggplot2 boxplot + beeswarm plot objects
      ggplotList <- list()
      ## Plot a general boxplot + beeswarm plot for multiple measures
      if(FALSE){
        ggplotList[["general"]] <- ggplot2::ggplot(
          OneToolSummary[["extraction"]],
          ggplot2::aes(x = .data$toolName, y = .data$value))
        ggplotList[["general"]] <- ggplotList[["general"]] +
          ## Draw boxplot
          ggplot2::geom_boxplot(
            ## Change filling color to white
            fill = "#FFFFFF",
            ## Maximize the violin plot width
            scale = "width"
            ## Hide outliers
            #outlier.shape = NA
          ) +
          ## Draw beeswarm plot
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE, ## Make repetitive points with the same Y to dodge on X axis
                                       size = 0.3, ## Make dot size smaller
          ) +
          ## If fill is set to single color, disable scale_fill_manual
          #ggplot2::scale_fill_manual(
          #  values = grDevices::topo.colors(length(indexes)))
          ## Add title for general boxplot + beeswarm plot
          ggplot2::labs(
            title = paste0(toolName,": Summary plot for extraction measures")
          ) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(
            labels =function(x) sprintf("%.2f", x)
          )
      }
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
            title = paste0(toolName,": ",titles[index]),
            subtitle = subtitles[index],
            ## Change title of y axis (axis.title.y) into measure info (same as title)
            y = titles[index],
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


      ## Output measures in a png file
      for(index in indexes){
        ## Need to suppress warning,
        ## because for each facet, there was only one possible x
        ## and the layout$ in ggplot2:::build.ggplot() will throw a warning:
        ##  Warning in min(x) : no non-missing arguments to min; returning Inf
        ##  Warning in max(x) : no non-missing arguments to max; returning -Inf
        ##
        ## This warning is not an error, but it is annoying.
        suppressMessages(
          suppressWarnings(
            ggplot2::ggsave(paste0(out.dir,"/boxplot.onetool.extraction.",index,".png"),
                            plot = ggplotList[[index]], device = "png", dpi = 1000,limitsize = FALSE)
          )
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
          tmp <- data.frame(seed = names(multiRun$cosSim[[gtSigName]]),
                            gtSigName = gtSigName,
                            value = multiRun$cosSim[[gtSigName]],
                            toolName = toolName,
                            datasetName = datasetName,
                            datasetGroup = datasetGroup[datasetName],
                            datasetGroupName = datasetGroupName,
                            datasetSubGroup = datasetSubGroup[datasetName],
                            datasetSubGroupName = datasetSubGroupName,
                            stringsAsFactors = FALSE)
          rownames(tmp) <- NULL

          ## Create a data.frame for each measure,
          ## and summarize multi-Run, multiDataset values
          ## for each measure.
          if(is.null(OneToolSummary$cosSim[[gtSigName]])){
            OneToolSummary$cosSim[[gtSigName]] <- data.frame()
          }
          OneToolSummary$cosSim[[gtSigName]] <- rbind(OneToolSummary$cosSim[[gtSigName]],tmp)
        }
      }

      ## Calculate the stats of one-signature cosine similarity.
      OneToolSummary$stats$cosSim <- list()
      for(gtSigName in gtSigNames){
        currentStats <- summary(OneToolSummary$cosSim[[gtSigName]][,"value"])
        OneToolSummary$stats$cosSim[[gtSigName]] <- currentStats
      }



      ## Combine multiple one-signature cosine similarity data.frame
      ## into OneToolSummary$cosSim$combined
      OneToolSummary$cosSim$combined <- data.frame()
      for(gtSigName in gtSigNames){
        tmp <- data.frame(OneToolSummary$cosSim[[gtSigName]],
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL

        if(nrow(OneToolSummary$cosSim$combined) == 0 |
           ncol(OneToolSummary$cosSim$combined) == 0 |
           is.null(dim(OneToolSummary$cosSim$combined)) ) {
          OneToolSummary$cosSim$combined <- tmp
        } else {
          OneToolSummary$cosSim$combined <-
            rbind(OneToolSummary$cosSim$combined,tmp)
        }
      }

    }
    ## Plot one-signature cosine similarity boxplot + beeswarm plot for one tool
    { ## debug
      ## Create a list to store ggplot2 boxplot + beeswarm plot objects
      ggplotList <- list()
      ## Plot a general boxplot + beeswarm plot for multiple measures
      if(FALSE){
        ggplotList[["general"]] <- ggplot2::ggplot(
          OneToolSummary$cosSim$combined,
          ggplot2::aes(x = .data$toolName, y = .data$value))
        ## Draw boxplot + beeswarm plot
        ggplotList[["general"]] <- ggplotList[["general"]] +
          ## Draw geom_violin
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            ## Maximize the violin plot width
            scale = "width"
            ## Hide outliers
            #outlier.shape = NA
          ) +
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                       size = 0.3, ## Make dot size smaller
                                       #position = ggplot2::position_dodge(0.9),
          ) +
          ## Change filling colors
          ## If fill is set to single color, disable scale_fill_manual
          #ggplot2::scale_fill_manual(
          #  values = grDevices::topo.colors(length(indexes)))
          ## Add title for general boxplot + beeswarm plot
          ggplot2::labs(title = paste0(toolName,": Summary plot for one-signature cosine similarity")) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(
            ## For one-signature cosine similarity, set ylim from the minimum of Manhattan distance value to 1.
            limits = c(min(OneToolSummary$cosSim$combined$value),1),
            labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a value~datasetSubGroup beeswarm plot for each signature.
      for(gtSigName in gtSigNames){
        sigNum <- which(gtSigNames == gtSigName)
        ggplotList[[gtSigName]] <- ggplot2::ggplot(
          OneToolSummary$cosSim[[gtSigName]],
          ## Make sure that only one x-label is shown in one small facet.
          #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
          ggplot2::aes(x = .data$toolName, y = .data$value)
        )
        ## Add facets
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
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


      ## Output average cosine similarity in high resolution png file
      for(gtSigName in gtSigNames){
        suppressMessages(
          suppressWarnings(
            ggplot2::ggsave(filename = paste0(out.dir,"/boxplot.onetool.",gtSigName,".onesig.cossim.png"),
                            plot = ggplotList[[gtSigName]], device = "png", dpi = 1000,limitsize = FALSE)
          )
        )
      }

      ## Output multiple extraction measures in a pdf file
      grDevices::pdf(paste0(out.dir,"/boxplot.onetool.onesig.cossim.pdf"), pointsize = 1)
      for(gtSigName in gtSigNames)
        suppressMessages(suppressWarnings(print(ggplotList[[gtSigName]])))
      grDevices::dev.off()
    }



    ## Summarize attribution scaled Manhattan distance for one tool.
    {
      OneToolSummary$ManhattanDist <- list()

      for(datasetDir in dataset.dirs){
        thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
        toolName <- strsplit(basename(tool.dirname),".results")[[1]]
        ## Add multiRun <- NULL to please R check
        multiRun <- NULL
        load(paste0(thirdLevelDir,"/multiRun.RDa"))
        gtSigNames <- rownames(multiRun$ManhattanDist)
        sigNums <- length(gtSigNames)
        datasetName <- basename(datasetDir)

        for(gtSigName in gtSigNames){
          tmp <- data.frame(seed = colnames(multiRun$ManhattanDist),
                            gtSigName = gtSigName,
                            value = multiRun$ManhattanDist[gtSigName,],
                            toolName = toolName,
                            datasetName = datasetName,
                            datasetGroup = datasetGroup[datasetName],
                            datasetGroupName = datasetGroupName,
                            datasetSubGroup = datasetSubGroup[datasetName],
                            datasetSubGroupName = datasetSubGroupName,
                            stringsAsFactors = FALSE)
          rownames(tmp) <- NULL

          ## Create a data.frame for each measure,
          ## and summarize multi-Run, multiDataset values
          ## for each measure.
          if(is.null(OneToolSummary$ManhattanDist[[gtSigName]])){
            OneToolSummary$ManhattanDist[[gtSigName]] <- data.frame()
          }
          OneToolSummary$ManhattanDist[[gtSigName]] <- rbind(OneToolSummary$ManhattanDist[[gtSigName]],tmp)
        }
      }

      ## Combine multiple ground-truth signature Manhattan-distance data.frame
      ## into OneToolSummary$ManhattanDist$combined.
      OneToolSummary$ManhattanDist$combined <- data.frame()
      for(gtSigName in gtSigNames){
        tmp <- data.frame(OneToolSummary$ManhattanDist[[gtSigName]],
                          stringsAsFactors = FALSE)
        rownames(tmp) <- NULL

        if(nrow(OneToolSummary$ManhattanDist$combined) == 0 |
           ncol(OneToolSummary$ManhattanDist$combined) == 0 |
           is.null(dim(OneToolSummary$ManhattanDist$combined)) ) {
          OneToolSummary$ManhattanDist$combined <- tmp
        } else {
          OneToolSummary$ManhattanDist$combined <-
            rbind(OneToolSummary$ManhattanDist$combined,tmp)
        }
      }

    }
    ## Plot normalized Manhattan distance violin plot + beeswarm plot for one tool
    { ## debug
      ## Create a list to store ggplot2 boxplot + beeswarm plot objects
      ggplotList <- list()
      ## Plot a general boxplot + beeswarm plot for Manhttan distance of multiple signatures
      if(FALSE){
        ggplotList[["general"]] <- ggplot2::ggplot(
          OneToolSummary$ManhattanDist$combined,
          ggplot2::aes(x = .data$toolName, y = .data$value))
        ## Draw boxplot + beeswarm plot
        ggplotList[["general"]] <- ggplotList[["general"]] +
          ## Draw geom_violin
          ggplot2::geom_violin(
            ## Change filling color to white
            fill = "#FFFFFF",
            ## Maximize the violin plot width
            scale = "width"
            ## Hide outliers
            #outlier.shape = NA
          ) +
          ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
          ## Show mean of the extraction meaasure distribution, as a blue diamond.
          ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
          ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                       size = 0.3, ## Make dot size smaller
                                       #position = ggplot2::position_dodge(0.9),
          ) +
          ## Change filling colors
          ## If fill is set to single color, disable scale_fill_manual
          #ggplot2::scale_fill_manual(
          #  values = grDevices::topo.colors(length(indexes)))
          ## Add title for general boxplot + beeswarm plot
          ggplot2::labs(title = paste0(toolName,": Summary plot for Scaled Manhattan distance")) +
          ## Restrict the decimal numbers of values of measures to be 2
          ggplot2::scale_y_continuous(
            ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
            limits = c(0,max(OneToolSummary$ManhattanDist$combined$value)),
            labels =function(x) sprintf("%.2f", x))
      }
      ## Plot a value~datasetSubGroup beeswarm plot for each signature.
      for(gtSigName in gtSigNames){
        sigNum <- which(gtSigNames == gtSigName)
        ggplotList[[gtSigName]] <- ggplot2::ggplot(
          OneToolSummary$ManhattanDist[[gtSigName]],
          ## Make sure that only one x-label is shown in one small facet.
          #ggplot2::aes(x = .data$datasetGroup, y = .data$value)
          ggplot2::aes(x = .data$toolName, y = .data$value)
        )
        ## Add facets
        ggplotList[[gtSigName]] <- ggplotList[[gtSigName]] +
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
            y = paste0("Scaled Manhattan distance of ",gtSigName," exposure"),
            ## Change title of x axis to "Pearson's R squared"
            x = "Pearson's R squared") +
          ## Change title of legend to datasetGroupName
          ggplot2::guides(color = ggplot2::guide_legend(title = datasetGroupName)) +
          ## Restrict the decimal numbers of values of measures (y) to be 2
          ggplot2::scale_y_continuous(
            ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
            limits = c(0,max(OneToolSummary$ManhattanDist$combined$value)),
            ## Restrict the decimal numbers of values of measures (y) to be 2
            labels =function(x) sprintf("%.2f", x),
            ## Add a secondary axis title on the top of the plot
            ## Showing axis label indicating facets
            sec.axis = ggplot2::dup_axis(
              breaks = NULL, ## Don't show axis breaks
              labels = NULL, ## Don't show axis tickmarks
              name = "SBS1:SBS5 mutation count ratio"))
      }


      ## Output average cosine similarity in high resolution png file
      for(gtSigName in gtSigNames){
        suppressMessages(
          suppressWarnings(
            ggplot2::ggsave(filename = paste0(out.dir,"/boxplot.onetool.",gtSigName,".Manhattan.dist.png"),
                            plot = ggplotList[[gtSigName]], device = "png", dpi = 1000,limitsize = FALSE)
          )
        )
      }

      ## Output multiple extraction measures in a pdf file
      grDevices::pdf(paste0(out.dir,"/boxplot.onetool.Manhattan.dist.pdf"), pointsize = 1)
      for(gtSigName in gtSigNames)
        suppressMessages(suppressWarnings(print(ggplotList[[gtSigName]])))
      grDevices::dev.off()
    }


    ## Write Summary tables for measures
    for(index in indexes){
      write.csv(OneToolSummary[[index]],
                file = paste0(out.dir,"/",index,".csv"),
                quote = F, row.names = F)
    }

    ## Write stat summary information into a text file.
    utils::capture.output(OneToolSummary$stats,file = paste0(out.dir,"/stats.txt"))
    utils::capture.output(OneToolSummary$prop1,file = paste0(out.dir,"/prop1.txt"))

    OneToolSummary$datasetGroupName <- datasetGroupName
    OneToolSummary$datasetSubGroupName <- datasetSubGroupName

    save(OneToolSummary, file = paste0(out.dir,"/OneToolSummary.RDa"))
    invisible(OneToolSummary)
  }

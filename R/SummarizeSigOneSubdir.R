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
#' @param run.dir A directory which contains output of computational approach in one
#' run on a specific dataset, possibly with a specified seed. E.g.
#' \code{2b.Full_output_K_as_2/hdp.results/S.0.1.Rsq.0.1/seed.1/}.
#'
#' This code depends on a conventional directory structure documented
#' in \code{NEWS.md}.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures. Should
#' contain a file named \code{ground.truth.syn.exposures.csv}.
#' In PCAWG paper: \code{run.dir/../../}
#' In SBS1-SBS5 paper: \code{0.Input_datasets/S.0.1.Rsq.0.1/}
#'
#' @param extracted.sigs.path Path to extracted sigs file, For most computational approaches,
#' it's in:
#' \code{<run.dir>/extracted.signatures.csv}.
#'
#' For SigProExtractor, it's in:
#' \code{<run.dir>/<catalog.type>/Selected_Solution/De_Novo_Solution/De_Novo_Solution_Signatures_<catalog.type>.txt}
#' , and is converted to
#' \code{<run.dir>/extracted.signatures.csv} by
#' wrapper function \code{\link{SummarizeSigProExtractor}}.
#'
#' For SignatureAnalyzer, it's in:
#'  \code{<run.dir>/sa.output.sigs.csv}
#'
#' For EMu, helmsman.NMF, MultiModalMuSig.LDA and MultiModalMuSig.CTM, their extracted signature
#' files are not standard \code{ICAMS} catalog format. These non-standard extracted signatures
#' are converted to \code{extracted.signatures.csv}, a standard \code{ICAMS} catalog by their
#' respective wrapper functions:
#' (\code{\link{ReadEMuCatalog}}, \code{\link{helmsmanCatalog2ICAMS}}, \code{\link{MMCatalog2ICAMS}}).
#' They are dumped into \code{extracted.sigs.csv} same as most approaches.
#'
#' @param inferred.exp.path Path to inferred exposures file. For most computational approaches,
#' it's in:
#' \code{<run.dir>/inferred.exposures.csv}.
#'
#' For SigProExtractor, it's in
#' \code{<run.dir>/<catalog.type>/Selected_Solution/De_Novo_Solution/De_Novo_Solution_Activities_<catalog.type>.txt}
#' , and is converted to
#' \code{<run.dir>/inferred.exposures.csv} by
#' wrapper function \code{\link{SummarizeSigProExtractor}}.
#'
#' For SignatureAnalyzer, it's in:
#'  \code{<run.dir>/sa.output.exp.csv}
#'
#' For EMu, helmsman.NMF, MultiModalMuSig.LDA and MultiModalMuSig.CTM, their extracted signature
#' files are not standard \code{ICAMSxtra} catalog format. These non-standard extracted signatures
#' are converted to \code{extracted.signatures.csv}, a standard \code{ICAMSxtra} exposures by their
#' respective wrapper functions:
#' (\code{\link{ReadEMuExposureFile}}, \code{\link{ReadhelmsmanExposure}}, \code{\link{ReadExposureMM}}).
#' They are dumped into \code{inferred.exposures.csv} same as most approaches.
#'
#'
#' @param summarize.exp Whether to summarize exposures when the file specified
#' by \code{inferred.exp.path} exists.
#'
#' @param overwrite If TRUE overwrite and files in existing \code{run.dir/summary} folder.
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
           summarize.exp = TRUE,
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
      cat("Cosine similarity to each ground-truth signature\n"),
      sigAnalysis$cosSim,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis$ex.sigs),
      cat("\nFalse positive signatures\n"),
      sigAnalysis$extracted.with.no.best.match,
      cat("\nFalse negative signatures\n"),
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
    if(!is.null(inferred.exp.path) & summarize.exp) {

      if(file.exists(inferred.exp.path)) {
        exposureDiff <- ReadAndAnalyzeExposures(
          extracted.sigs = extracted.sigs.path,
          ground.truth.sigs =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
          inferred.exp.path = inferred.exp.path,
          ground.truth.exposures =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"))

        # Write results of sum of Manhattan distance
        write.csv(exposureDiff$SumOfManhattan,
                  file = paste0(outputPath,"/SumOfManhattan.csv"),
                  quote = T)

        ## Write TPR, PPV, and F1 sscore
        ## for exposure inference measures.
        write.csv(exposureDiff$F1.measures,
                  file = paste0(outputPath,"/F1.measures.csv"),
                  quote = T)

        # Write results of exposure inference measures,
        # in aggregated format for each tumor and each
        # ground-truth signature.
        for(spectrumName in names(exposureDiff$Manhattan)){
          ## Replace all characters unsuitable for filenames.
          cleanedName <- fs::path_sanitize(spectrumName, replacement = ".")

          write.csv(exposureDiff$Manhattan[[spectrumName]],
                    file = paste0(outputPath,"/Manhattan.",cleanedName,".csv"),
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
                   file = paste0(outputPath,"/assessment.sessionInfo.txt"))

    ## Save Signature extraction summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    save(sigAnalysis,
         file = paste0(outputPath,"/sigAnalysis.RDa"))
    ## Save exposure inference summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    if(file.exists(inferred.exp.path) & summarize.exp){
      save(exposureDiff,
           file = paste0(outputPath,"/exposureDiff.RDa"))
    }

    invisible(sigAnalysis) # So we have something to check in tests
  }


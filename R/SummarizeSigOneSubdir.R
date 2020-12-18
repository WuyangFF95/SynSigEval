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


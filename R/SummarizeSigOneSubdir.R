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
#' @param run.dir A directory which contains output of computational approach
#' in one run on a specific dataset, possibly with a specified seed. E.g.
#' \code{2b.Full_output_K_as_2/hdp.results/S.0.1.Rsq.0.1/seed.1/}.
#'
#' This code depends on a conventional directory structure documented
#' in \code{NEWS.md}.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' Should contain a file named \code{ground.truth.syn.exposures.csv}.
#' In PCAWG paper: \code{run.dir/../../}
#' In SBS1-SBS5 paper: \code{0.Input_datasets/S.0.1.Rsq.0.1/}
#'
#' @param extracted.sigs.path Path to extracted sigs file.
#' For most computational approaches, it's in:
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
#' For EMu, helmsman.NMF, MultiModalMuSig.LDA and MultiModalMuSig.CTM,
#' their extracted signature files are not standard \code{ICAMS} catalog format.
#' These non-standard extracted signatures are converted to
#' \code{extracted.signatures.csv}, a standard \code{ICAMS} catalog
#' by their respective wrapper functions: \itemize{
#' \item \code{\link{ReadEMuCatalog}},
#' \item \code{\link{helmsmanCatalog2ICAMS}}, and
#' \item \code{\link{MMCatalog2ICAMS}}
#' }
#' They are dumped into \code{extracted.sigs.csv} same as most approaches.
#'
#' @param inferred.exp.path Path to inferred exposures file.
#' For most computational approaches, it's in:
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
#' For EMu, helmsman.NMF, MultiModalMuSig.LDA and MultiModalMuSig.CTM,
#' their extracted signature files are not standard \code{ICAMSxtra} catalog format.
#' These non-standard extracted signatures are converted to
#' \code{extracted.signatures.csv}, a standard \code{ICAMSxtra} exposures
#' by their respective wrapper functions: \itemize{
#' \item \code{\link{ReadEMuExposureFile}}
#' \item \code{\link{ReadhelmsmanExposure}}, and
#' \item \code{\link{ReadExposureMM}}
#' }
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
#' @param export.Manhattan.each.spectrum Whether to export csv files for Manhattan
#' distance of each mutational spectrum.
#'
#' @importFrom utils write.csv capture.output sessionInfo
#' @importFrom gtools mixedsort
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
           summary.folder.name = "summary",
           export.Manhattan.each.spectrum = FALSE) {

    # Specify and create output directory -------------------------------------

    outputPath <- paste0(run.dir, "/", summary.folder.name)
    if (dir.exists(outputPath)) {
      if (!overwrite) stop(outputPath, " already exists")
    }
    suppressWarnings(dir.create(outputPath))



    # Calculate signature extraction measures ---------------------------------

    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
        ground.truth.exposures =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv")
      )

     # Temporary workaround.
     # table does not have a correct column name.
     # replace "gt,sig" with "gt.sig"
     colnames(sigAnalysis$table)[2] <- "gt.sig"

     sigAnalysis$TP.ex.names <- sigAnalysis$table$ex.sig
     sigAnalysis$TP.gt.names <- sigAnalysis$table$gt.sig
     sigAnalysis$FP.names <-
       setdiff(colnames(sigAnalysis$ex.sigs), sigAnalysis$table$ex.sig)
     sigAnalysis$FN.names <-
       setdiff(colnames(sigAnalysis$gt.sigs), sigAnalysis$table$gt.sig)

     sigAnalysis$averCosSim <- sigAnalysis$avg.cos.sim
     sigAnalysis$avg.cos.sim <- NULL
     sigAnalysis$PPV <- sigAnalysis$TP / (sigAnalysis$TP + sigAnalysis$FP)
     sigAnalysis$TPR <- sigAnalysis$TP / (sigAnalysis$TP + sigAnalysis$FN)



     # Calculate best cosine similarity for each gt sig -----------------------

     # gt sigs in sigAnalysis$table are gt sigs with best match.
     # For these sigs, we directly record the cosSim values in $table.
     cosSim_TP_gt <- numeric(0)
     if(nrow(sigAnalysis$table) > 0) {
       cosSim_TP_gt <- sigAnalysis$table$sim
       names(cosSim_TP_gt) <- sigAnalysis$table$gt.sig
     }

     # For gt sigs which do not have a best match (and thus FN sigs),
     # their best cosine similarity values are calculated in cosSim_FN
     cosSim_FN <- numeric(0)
     if(sigAnalysis$FN > 0) {
       cosSim_FN <- numeric(sigAnalysis$FN)
       names(cosSim_FN) <- sigAnalysis$FN.names
       for(gt_sig_name in sigAnalysis$FN.names) {
         cosSim_FN[gt_sig_name] <- max(sigAnalysis$sim.matrix[, gt_sig_name])
       }
     }
     sigAnalysis$cosSim <- c(cosSim_TP_gt, cosSim_FN)





    # Export results summary from sigAnalysis ---------------------------------

    # Copy ground.truth exposures from second.level.dir
    # to outputPath == run.dir/<summary.folder.name>.
    CopyWithChecks(
      from = paste0(ground.truth.exposure.dir, "/ground.truth.syn.exposures.csv"),
      to.dir = outputPath,
      overwrite = TRUE)

    # Export bi-lateral matching table between extracted and ground-truth sigs
    write.csv(sigAnalysis$table,
              file = paste(outputPath,"match.ex.to.gt.csv",sep = "/"),
              row.names = F)
    # Export full cosine similarity table
    write.csv(sigAnalysis$sim.matrix,
              file = paste(outputPath,"full.cossims.ex.to.gt.csv",sep = "/"),
              row.names = T)

    # Export ground-truth sigs with non-zero ground-truth exposures
    ICAMS::WriteCatalog(
      sigAnalysis$gt.sigs,
      paste(outputPath,"ground.truth.sigs.csv",sep = "/"),
    )
    # Export extracted signatures with new names
    ex.sigs.renamed <- RelabelExSigs(sigAnalysis)
    sigAnalysis$ex.sigs <- ex.sigs.renamed
    ICAMS::WriteCatalog(
      sigAnalysis$ex.sigs,
      paste(outputPath,"extracted.sigs.csv",sep = "/"))

    # Exports other outputs into "other.results.txt"
    capture.output(
      cat("Average cosine similarity, only best matches are considered\n"),
      sigAnalysis$averCosSim,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis$ex.sigs),
      cat("\nTrue positive ground-truth signatures\n"),
      mixedsort(sigAnalysis$TP.gt.names),
      cat("\nTrue positive extracted signatures\n"),
      mixedsort(sigAnalysis$TP.ex.names),
      cat("\nFalse positive signatures\n"),
      mixedsort(sigAnalysis$FP.names),
      cat("\nFalse negative signatures\n"),
      mixedsort(sigAnalysis$FN.names),
      file = paste0(outputPath,"/other.results.txt"))



    # Plot signatures to PDF files --------------------------------------------

    # Currently, ICAMS cannot plot COMPOSITE catalog.
    # TODO(Wuyang): To add a ICAMS:::PlotCatalog.COMPOSITECatalog function

    # Plot ground-truth sigs with non-zero ground-truth exposures
    if("COMPOSITECatalog" %in% class(sigAnalysis$gt.sigs) == FALSE){
      ICAMS::PlotCatalogToPdf(sigAnalysis$gt.sigs,
                              paste0(outputPath,"/ground.truth.sigs.pdf"))
    }
    # Plot renamed and sorted extracted sigs
    if("COMPOSITECatalog" %in% class(sigAnalysis$ex.sigs) == FALSE){
      ICAMS::PlotCatalogToPdf(sigAnalysis$ex.sigs,
                              paste0(outputPath,"/extracted.sigs.pdf"))
    }



    # Summarize signature attribution (a.k.a. exposure inference) -------------

    # Only when inferred.exp.path is not empty
    # AND flag "summarize.exp" is set as TRUE
    if(!is.null(inferred.exp.path) && summarize.exp == TRUE) {

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

        # Write TPR, PPV, and F1 score
        # for exposure inference measures.
        write.csv(exposureDiff$F1.measures,
                  file = paste0(outputPath,"/F1.measures.csv"),
                  quote = T)

        # Write results of exposure inference measures,
        # in aggregated format for each tumor and each
        # ground-truth signature.
        if(export.Manhattan.each.spectrum){
          for(spectrumName in names(exposureDiff$Manhattan)){
            ## Replace all characters unsuitable for filenames.
            cleanedName <- fs::path_sanitize(spectrumName, replacement = ".")

            write.csv(exposureDiff$Manhattan[[spectrumName]],
                      file = paste0(outputPath,"/Manhattan.",cleanedName,".csv"),
                      quote = T)
          }
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



    # Save and return result summary lists, and log file ----------------------

    # Log of system time and session info
    capture.output(Sys.time(), sessionInfo(),
                   file = paste0(outputPath,"/assessment.sessionInfo.txt"))

    # Save Signature extraction summary in an RDa file,
    # for reuse in SummarizeMultiRuns().
    save(sigAnalysis,
         file = paste0(outputPath,"/sigAnalysis.RDa"))
    # Save exposure inference summary in an RDa file,
    # for reuse in SummarizeMultiRuns().
    if(file.exists(inferred.exp.path) & summarize.exp){
      save(exposureDiff,
           file = paste0(outputPath,"/exposureDiff.RDa"))
    }
    # Returned for checking in tests
    invisible(sigAnalysis)
  }


#' Assess/evaluate results from packages which can do
#' BOTH extraction and attribution,
#' excluding SigProfiler-Python and SignatureAnalyzer.
#'
#' Packages including but not limited to:
#' hdp, MutationalPatterns, sigfit,
#' signeR, SomaticSignatures.
#'
#' @inheritParams SummarizeSigOneSubdir
#'
#' @export
#'
#' @importFrom ICAMS WriteCatalog ReadCatalog
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
SummarizeSigOneExtrAttrSubdir <-
  function(run.dir,
           ground.truth.exposure.dir = paste0(run.dir,"/../../../"),
           summarize.exp = TRUE,
           overwrite = FALSE,
           summary.folder.name = "summary",
           export.Manhattan.each.spectrum = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if SigProExtractor updates!
    inputPath <- run.dir
    stopifnot(dir.exists(inputPath))

    # Specify the path of extracted signatures in ICAMS csv format.
    extracted.sigs.path <- paste0(inputPath,"/extracted.signatures.csv")

    # Specify the path of inferred exposures in SynSig csv format.
    inferred.exp.path <- paste0(inputPath,"/inferred.exposures.csv")

    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under run.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        run.dir = inputPath,
        ground.truth.exposure.dir = ground.truth.exposure.dir,
        extracted.sigs.path = extracted.sigs.path,
        inferred.exp.path = inferred.exp.path,
        summarize.exp = summarize.exp,
        overwrite = overwrite,
        summary.folder.name = summary.folder.name,
        export.Manhattan.each.spectrum = export.Manhattan.each.spectrum)

    invisible(retval) # So we can test without looking at a file.
  }






#' Assess/evaluate results from packages which can
#' ONLY do exposure attribution.
#'
#' Packages including but not limited to:
#' deconstructSigs, YAPSA.
#'
#' Here, we excluded SignatureEstimation. Although it is also
#' a package with only attribution, but it has two attribution
#' algorithms. Therefore the naming of the results are slightly
#' different from the other two packages.
#'
#' @inheritParams SummarizeSigOneSubdir
#'
#' @param run.dir Lowest level path to results, e.g.
#' \code{<top.dir>}\code{/sa.sa.96/Attr/YAPSA.results/seed.1/}
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. For packages which can do both extraction and attribution,
#' we expect two files, \code{ground.truth.signatures.csv}
#' and \code{inferred.exposures.csv} are in the folder.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' It defaults to be \code{sub.dir}, i.e. \code{run.dir}/../../
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS WriteCatalog ReadCatalog
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
SummarizeSigOneAttrSubdir <-
  function(run.dir,
           ground.truth.exposure.dir = paste0(run.dir,"/../../../"),
           overwrite = FALSE,
           summary.folder.name = "summary",
           export.Manhattan.each.spectrum = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if SigProExtractor updates!
    inputPath <- run.dir
    stopifnot(dir.exists(inputPath))

    # Specify the path of extracted signatures in ICAMS csv format.
    ground.truth.sigs.path <- paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv")

    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under run.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        run.dir = inputPath,
        ground.truth.exposure.dir = ground.truth.exposure.dir,
        extracted.sigs.path = ground.truth.sigs.path,
        inferred.exp.path = paste0(inputPath,"/inferred.exposures.csv"),
        summarize.exp = TRUE,
        overwrite = overwrite,
        summary.folder.name = summary.folder.name,
        export.Manhattan.each.spectrum = export.Manhattan.each.spectrum)

    invisible(retval) # So we can test without looking at a file.
  }

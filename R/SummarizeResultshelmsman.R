#' Assess/evaluate results from SigProfiler-python (a.k.a. SigProExtractor)
#' Assessment is restricted to v0.0.5.43,
#' because different version has different folder structure.
#'
#' @param run.dir Lowest level path to results, e.g.
#' \code{<top.dir>}\code{/sa.sa.96/ExtrAttr/SigProExtractor.results/seed.1/}
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory
#' \code{<run.dir>}\code{/SBS96} which
#' stores SigProfiler results.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' Usually, it refers to \code{sub.dir}, i.e. \code{run.dir}/../../../
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @param hierarchy Whether the user have enabled hierarchy = True
#' when running SigProExtractor.
#' specifying True or False into SigProExtractor will cause the program
#' to generate different folder structure. (Default: \code{FALSE})
#'
#' @export
#'
#' @importFrom ICAMS WriteCatalog ReadCatalog
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
SummarizeSigOnehelmsmanSubdir <-
  function(run.dir,
           ground.truth.exposure.dir = paste0(run.dir,"/../../../"),
           overwrite = FALSE,
           hierarchy = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if SigProExtractor updates!
    inputPath <- paste0(run.dir)
    stopifnot(dir.exists(inputPath))

    # Read in extracted signatures in SigProExtractor txt format,
    # and convert it to ICAMS csv format.
    # Need special function to read in extracted signatures
    # Converted signatures will be included in the /summary folder.
    extractedSigs <- helmsmanCatalog2ICAMS(
      cat = paste0(inputPath,"/H_loadings.txt"),
      region = "unknown",
      catalog.type = "counts.signature")
    ## extracted signatures need to be normalized.
    for(sigName in colnames(extractedSigs)){
      extractedSigs[,sigName] <- extractedSigs[,sigName] / sum(extractedSigs[,sigName])
    }

    extracted.sigs.path <- paste0(run.dir,"/extracted.signatures.csv")
    ICAMS::WriteCatalog(extractedSigs, extracted.sigs.path)

    # Read in inferred exposures in SP format,
    # and convert it into our internal format
    inferred.exp.path.helmsman.format <-
      paste0(inputPath,"/W_components.txt")
    rawExposure <- ReadhelmsmanExposure(
      inferred.exp.path.helmsman.format,
      check.names = FALSE)

    spectra <- ICAMS::ReadCatalog(
      file = paste0(ground.truth.exposure.dir,"/ground.truth.syn.catalog.csv"),
      catalog.type = "counts",
      strict = FALSE)
    exposureCounts <- rawExposure
    ## The sum of exposure of each spectrum needs to
    ## be normalized to the total number of mutations
    ## in each spectrum.
    for(sample in colnames(exposureCounts)){
      exposureCounts[,sample] <- rawExposure[,sample] / sum(rawExposure[,sample]) * sum(spectra[,sample])
    }

    inferred.exp.path <- paste0(run.dir,"/inferred.exposures.csv")
    ICAMSxtra::WriteExposure(exposureCounts,inferred.exp.path)


    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under run.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        run.dir = run.dir,
        ground.truth.exposure.dir = ground.truth.exposure.dir,
        extracted.sigs.path = extracted.sigs.path,
        inferred.exp.path = inferred.exp.path,
        overwrite = overwrite)

    invisible(retval) # So we can test without looking at a file.
  }

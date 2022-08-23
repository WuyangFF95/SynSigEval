#' Assess/evaluate results from helmsman.NMF
#'
#' @inheritParams SummarizeSigOneSubdir
#'
#' @param run.dir A directory which contains output of helmsman.NMF in one
#' run on a specific dataset, possibly with a specified seed. E.g.
#' \code{2b.Full_output_K_as_2/helmsman.NMF.results/S.0.1.Rsq.0.1/seed.1/}.
#'
#' This code depends on a conventional directory structure documented
#' in \code{NEWS.md}.
#'
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
           summarize.exp = TRUE,
           overwrite = FALSE) {

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
    ## The sum of exposure of each spectrum needs to
    ## be normalized to the total number of mutations
    ## in each spectrum.
    rawExposure <- ReadhelmsmanExposure(
      inferred.exp.path.helmsman.format,
      check.names = FALSE)
    spectra <- ICAMS::ReadCatalog(
      file = paste0(ground.truth.exposure.dir,"/ground.truth.syn.catalog.csv"),
      catalog.type = "counts",
      strict = FALSE)
    exposureCounts <- rawExposure
    for(sample in colnames(exposureCounts)){
      exposureCounts[,sample] <- rawExposure[,sample] / sum(rawExposure[,sample]) * sum(spectra[,sample])
    }

    inferred.exp.path <- paste0(run.dir,"/inferred.exposures.csv")
    mSigHdp::WriteExposure(exposureCounts,inferred.exp.path)


    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under run.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        run.dir = run.dir,
        ground.truth.exposure.dir = ground.truth.exposure.dir,
        extracted.sigs.path = extracted.sigs.path,
        inferred.exp.path = inferred.exp.path,
        summarize.exp = summarize.exp,
        overwrite = overwrite)

    invisible(retval) # So we can test without looking at a file.
  }

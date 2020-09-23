#' @title SynSigEval
#'
#' @description Assess the performance of mutational-signature analysis
#' programs Using catalogs of synthetic mutational spectra created
#' by package \code{SynSigGen}.
#'
#' @section Overview:
#'
#' The main focus is generating synthetic catalogs of mutational
#' spectra (mutations in tumors) based on known mutational signature
#' profiles and attributions (assignment of exposures to tumors) in
#' the PCAWG7 data. We call this kind of synthetic
#' data broadly "reality-based" synthetic
#' data.  The package also has a set of functions that
#' generate random mutational signature profiles and then create
#' synthetic catalogs based on these random signature profiles. We
#' call this kind of synthetic data "random" synthetic data, while
#' pointing out that much depends on the distributions from which
#' the random signature profiles and attributions are generated.
#'
#' Typical workflow for generating catalogs of "reality-based" synthetic
#' mutational spectra is as follows.
#'
#' \preformatted{
#'
#' In \code{SynSigGen}:
#'
#' Input (based on SignatureAnalyzer or SigProfiler analysis of PCAWG tumors)
#'   A, matrix of attributions (signatures x samples)
#'   S, mutational signature profiles (mutation type x signature)
#'
#' P <- GetSynSigParamsFromExposures(A, ...)
#'
#' synthetic.exposures <- GenerateSyntheticExposures(P, ...)
#'
#' synthetic.spectra <- CreateAndWRiteCatalog(S, synthetic.exposures, ...)
#'
#' In \code{SynSigEval}:
#'
#' T <- Signatures extracted by SignatureAnalzer, SigProfiler,
#' or othercomputational approaches on synthetic.spectra
#'
#' SummarizeResults(T, S, synthetic.exposures, ...)
#'
#' }
#'
#' @section (In SynSigGen) Creating Synthetic Mutational Catalogs:
#'
#' These functions create synthetic mutational catalogs based
#' on parameters derived from signature profiles
#' and attributions (exposures).
#'
#' @section Summarize results (of signature extraction):
#'
#' Relevant functions are: \enumerate{
#'
#' \item \code{\link{SummarizeSigProExtractor}}
#' \item \code{\link{SignatureAnalyzerSummarizeTopLevel}}
#' \item \code{\link{SignatureAnalyzerSummarizeSBS1SBS5}}
#' \item \code{\link{SummarizeSigOneExtrAttrSubdir}}
#' \item \code{\link{SummarizeSigOneAttrSubdir}}
#' \item \code{\link{SummarizeSigOnehelmsmanSubdir}}
#' \item \code{\link{SummarizeSigOneSigProSSSubdir}}
#'
#' }
#'
#' @section Comparing two sets of mutational signatures:
#'
#' Functions for comparing mutational signatures and
#' sets of mutational signatures. Often we will be interested
#' in comparing signature profiles extracted from synthetic data to the
#' ground-truth signature profiles.
#'
#' \code{\link{Match1Sig}},
#' \code{\link{MatchSigs1Direction}},
#' \code{\link{MatchSigs2Directions}},
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @section Folder structure for SynSigEval v0.2:
#'
#' Summary function will fit to the new 5-level folder structure:
#'
#' First Level - \code{top.level.dir}: dataset folder (e.g. "S.0.1.Rsq.0.1", "syn.pancreas").
#' All spectra datasets under any top.level.dir have the same exposure.
#'
#' Second Level - \code{ground.truth.exposure.dir}: spectra folder: (e.g. "sp.sp", "sa.sa.96").
#' All spectra datasets under any second.level.dir have the same signature and
#' the same exposure counts.
#'
#' Third Level - \code{third.level.dir}: It can be ("Attr") for storing results of packages
#' which can only do exposure attribution of known signatures ("Attr");
#' it can also be ("ExtrAttr"), folder to store results of software packages which
#' can do de-novo extraction and following attribution.
#'
#' Fourth Level - \code{tool.dir}: The results of a software package
#' (e.g. "SigProExtractor.results","SignatureEstimation.QP.results").
#' Under this level, \code{tool.dir} may contain multiple \code{run.dir},
#' each is a run of the software package using a specific number of seed.
#'
#' Fifth level - \code{run.dir}: contains results from a run of the software package
#' using a specific number of seed. (e.g. "seed.1")
#'
#'
#' @docType package
#' @name SynSigEval

NULL

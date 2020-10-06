#' @title SynSigEval
#'
#' @description Assess the performance of two steps in mutational signature
#' analysis: \itemize{
#' \item signature extraction
#' \item exposure inferrence (a.k.a. signature attribution)
#' } by computational approaches, using catalogs of synthetic mutational
#' spectra created by package \code{SynSigGen}.
#'
#'
#' @section Input:
#'
#' \code{SynSigEval} requires the input data listed below: \enumerate{
#' \item E, matrix of synthetic exposures (signatures x samples)
#'
#' \item S, mutational signature profiles (mutation type x signature)
#'
#' \item synthetic.spectra, synthetic mutational spectra with known
#' ground-truth mutational signature profiles (S) and exposures
#' (synthetic.exposures). It can be created from \code{SynSigGen}.
#'
#' \item T, signatures extracted by SignatureAnalzer, SigProfiler,
#' or other computational approaches on \code{synthetic.spectra}.
#' For attribution-only approaches, T=S.
#'
#' \item F, exposures inferred by computational approaches on
#' \code{synthetic.spectra}.
#' }
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
#' @section Summarize results:
#'
#' \enumerate{
#'
#' \item Summarize results in fifth-level \code{run.dir}:
#'
#' Relevant functions are: \itemize{
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
#' \item Summarize results of multiple runs by a computational approach on one spectra data set:
#'
#' \code{SummarizeMultiRuns}
#'
#' \item Summarize results of multiple computational approaches on one spectra data set:
#'
#' \code{SummarizeMultiToolsOneDataset}
#'
#' \item Summarize results of multiple computational approaches on multiple spectra data sets:
#'
#' \code{SummarizeMultiToolsMultiDatasets}
#'
#' }
#'
#' @docType package
#' @name SynSigEval

NULL

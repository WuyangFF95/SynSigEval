#' Summarize results (SBS96, DBS, ID or COMPOSITE) from SignatureAnalyzer
#'
#' @param run.dir Lowest level path to results, that is
#' \code{<top.dir>}\code{/sa.sa.96/sa.results/},
#' \code{<top.dir>}\code{/sp.sp/sa.results/},
#' \code{<top.dir>}\code{/sa.sa.COMPOSITE/sa.results/}, or
#' \code{<top.dir>}\code{/sp.sa.COMPOSITE/sa.results/}.
#'
#' We expect \code{run.dir} contain the best-run subdirectory (e.g. "best.run").
#' The name of the subdirectory needs to be given to \code{which.run} parameter.
#'
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' It defaults to be \code{sub.dir}, i.e. \code{run.dir}/../../
#'
#' @param which.run Name of subdirectory under \code{run.dir}
#' containing the run to summarize.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @keywords internal
#'
#' @importFrom ICAMS WriteCatalog ReadCatalog PlotCatalogToPdf
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSASubdir <-
  function(run.dir,
           ground.truth.exposure.dir = paste0(run.dir,"/../../"),
           which.run = "/best.run/",
           overwrite = FALSE) {
    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(run.dir, "/", which.run)

    if (!dir.exists(inputPath)) stop(inputPath, "does not exist")

    retval <-
      SummarizeSigOneSubdir(
        run.dir = run.dir,
        ground.truth.exposure.dir = ground.truth.exposure.dir,
        extracted.sigs.path = paste0(inputPath,"/sa.output.sigs.csv"),
        inferred.exp.path = paste0(inputPath,"/sa.output.exp.csv"),
        overwrite = overwrite)

    invisible(retval)
  }

#' Summarize all subdirectories of SignatureAnalyzer results
#' on a major dataset.
#'
#' This function depends on a particular directory structure: see
#' argument \code{top.level.dir}. This function finds the
#' best of multiple SignatureAnalyzer extraction runs and summarizes
#' the comparison of the best run with the ground truth.
#'
#' @param top.level.dir Path to top level directory, which
#' must contain the following subdirectories:
#' \itemize{
#' \item \code{sa.sa.96/sa.results/}
#' \item \code{sp.sp/sa.results/}
#' \item \code{sa.sa.COMPOSITE/sa.results/}
#' \item \code{sp.sa.COMPOSITE/sa.results/}
#' }
#' Each of the directories must contain
#' additional subdirectories, one for each SignatureAnalyzer
#' run, names \code{sa.run.<n>}, where <n> is an integer
#' (string of digits).
#'
#' @param overwrite If TRUE overwrite existing summary files.
#'
#' @export

SignatureAnalyzerSummarizeTopLevel <-
  function(top.level.dir, overwrite = FALSE) {
    stopifnot(dir.exists(top.level.dir))

    assign("last.warning", NULL, envir = baseenv())

    options(warn = 2) # Warnings treated as errors

    sa.sa.96.dir <- paste0(top.level.dir, "/sa.sa.96/sa.results")
    stopifnot(dir.exists(sa.sa.96.dir))
    sp.sp.dir <- paste0(top.level.dir, "/sp.sp/sa.results")
    stopifnot(dir.exists(sp.sp.dir))
    sa.sa.COMPOSITE.dir <-
      paste0(top.level.dir, "/sa.sa.COMPOSITE/sa.results")
    stopifnot(dir.exists(sa.sa.COMPOSITE.dir))
    sp.sa.COMPOSITE.dir <-
      paste0(top.level.dir, "/sp.sa.COMPOSITE/sa.results")
    stopifnot(dir.exists(sp.sa.COMPOSITE.dir))

    CopyBestSignatureAnalyzerResult(sa.sa.96.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sp.sp.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sa.sa.COMPOSITE.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sp.sa.COMPOSITE.dir, overwrite = overwrite)

    retval <-
      list(sa.sa.96 =
             SummarizeSigOneSASubdir(
               sa.sa.96.dir, overwrite = overwrite),
           sp.sp =
             SummarizeSigOneSASubdir(
               sp.sp.dir, overwrite = overwrite),
           sa.sa.COMPOSITE =
             SummarizeSigOneSASubdir(
               sa.sa.COMPOSITE.dir, overwrite = overwrite),
           sp.sa.COMPOSITE =
             SummarizeSigOneSASubdir(
               sp.sa.COMPOSITE.dir, overwrite = overwrite))

    capture.output(print(retval), file = paste0(top.level.dir, "/retval.txt"))
    invisible(retval)
  }

#' Summarize all sub-directories of SignatureAnalyzer results
#' on the correlated SBS1 / SBS5.
#'
#' This is special-purpose function to summarize results
#' from one in-silico experiment that examines how well
#' signatures can be extracted from synthetic tumors with
#' correlated SBS1 and SBS5.
#'
#' @param top.level.dir Path to top level directory.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
SignatureAnalyzerSummarizeSBS1SBS5 <-
  function(top.level.dir, overwrite = FALSE) {
    stopifnot(dir.exists(top.level.dir))

    assign("last.warning", NULL, envir = baseenv())

    options(warn = 2) # Warnings treated as errors

    subdirs <-
      c("S.0.1.Rsq.0.1", "S.0.1.Rsq.0.2", "S.0.1.Rsq.0.3", "S.0.1.Rsq.0.6",
        "S.0.5.Rsq.0.1", "S.0.5.Rsq.0.2", "S.0.5.Rsq.0.3", "S.0.5.Rsq.0.6",
        "S.1.Rsq.0.1",   "S.1.Rsq.0.2",   "S.1.Rsq.0.3",   "S.1.Rsq.0.6")

    Summarize1 <- function(pseudo.top.level) { # pseudo.third.level is e.g. S.0.1.Rsq.0.1
      sp.sp.path <- paste0(top.level.dir, "/", pseudo.top.level, "/sp.sp/")
      if (!dir.exists(sp.sp.path)) stop(sp.sp.path, "does not exist")
      sa.results.path <- paste0(sp.sp.path, "sa.results/")
      if (!dir.exists(sa.results.path)) stop(sa.results.path, "does not exist")
      CopyBestSignatureAnalyzerResult(sa.results.path, overwrite = overwrite)
      SummarizeSigOneSASubdir(sa.results.path, overwrite = overwrite)
    }

    retval <- lapply(subdirs, Summarize1)

    capture.output(print(retval), file = paste0(top.level.dir, "/retval.txt"))
    invisible(retval)
}

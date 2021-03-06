% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummarizeResultsSA.R
\name{SummarizeSigOneSASubdir}
\alias{SummarizeSigOneSASubdir}
\title{Summarize results (SBS96, DBS, ID or COMPOSITE) from SignatureAnalyzer}
\usage{
SummarizeSigOneSASubdir(
  run.dir,
  ground.truth.exposure.dir = paste0(run.dir, "/../../"),
  which.run = "/best.run/",
  summarize.exp = TRUE,
  overwrite = FALSE
)
}
\arguments{
\item{run.dir}{Lowest level path to results, for example in PCAWG paper:
\code{<top.dir>/sa.sa.96/sa.results/},
\code{<top.dir>/sp.sp/sa.results/},
\code{<top.dir>/sa.sa.COMPOSITE/sa.results/}, or
\code{<top.dir>/sp.sa.COMPOSITE/sa.results/}.

In SBS1-SBS5 paper:
\code{2b.Full_output_K_as_2/SignatureAnalyzer.results/S.0.1.Rsq.0.1/seed.1/}


We expect \code{run.dir} contain the best-run subdirectory (e.g. "best.run").
The name of the subdirectory needs to be given to \code{which.run} parameter.

Here, \code{<top.dir>} refers to a top-level directory which contains the
full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
This code depends on a conventional directory structure documented
elsewhere.}

\item{ground.truth.exposure.dir}{Folder which stores ground-truth exposures. Should
contain a file named \code{ground.truth.syn.exposures.csv}.
In PCAWG paper: \code{run.dir/../../}
In SBS1-SBS5 paper: \code{0.Input_datasets/S.0.1.Rsq.0.1/}}

\item{which.run}{Name of subdirectory under \code{run.dir}
containing the run to summarize.}

\item{summarize.exp}{Whether to summarize exposures when the
\code{<run.dir>/<which.run>/sa.output.exp.csv} exists.}

\item{overwrite}{If TRUE overwrite and files in existing \code{run.dir/summary} folder.}
}
\description{
Summarize results (SBS96, DBS, ID or COMPOSITE) from SignatureAnalyzer
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummarizeResultsSP.R
\name{SummarizeSigProExtractor}
\alias{SummarizeSigProExtractor}
\title{Summarize SigProfiler results in the sa.sa.96 and/or sp.sp subdirectories.}
\usage{
SummarizeSigProExtractor(
  top.dir,
  sub.dir = c("sa.sa.96", "sp.sp"),
  overwrite = FALSE
)
}
\arguments{
\item{top.dir}{The top directory of a conventional data structure containing
at least one of the subdirectories: sa.sa.96/sp.results and sp.sp/sp.results;
see further documentation elsewhere.}

\item{sub.dir}{The subdirectory under \code{top.dir}, and containing a folder
named sp.results. By default, it contains both \code{c("sa.sa","sp.sp")}.
But you should specify \code{sub.dir = "sp.sp"} for \code{top.dir} with only
the \code{sp.sp} subdirectory
(as is the case for the correlated SBS1-and-SBS5-containing data sets).}

\item{overwrite}{If TRUE overwrite and files in existing \code{run.dir/summary} folder.}
}
\description{
Summarize SigProfiler results in the sa.sa.96 and/or sp.sp subdirectories.
}
\details{
Results are put in standardized subdirectories of \code{top.dir}.
}

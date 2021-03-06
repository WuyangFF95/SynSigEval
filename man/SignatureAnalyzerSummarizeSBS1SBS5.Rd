% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummarizeResultsSA.R
\name{SignatureAnalyzerSummarizeSBS1SBS5}
\alias{SignatureAnalyzerSummarizeSBS1SBS5}
\title{Summarize all sub-directories of SignatureAnalyzer results
on the correlated SBS1 / SBS5.}
\usage{
SignatureAnalyzerSummarizeSBS1SBS5(
  top.level.dir,
  summarize.exp = TRUE,
  overwrite = FALSE
)
}
\arguments{
\item{top.level.dir}{Path to top level directory.}

\item{summarize.exp}{Whether to summarize exposures when the file specified
by \code{inferred.exp.path} exists.}

\item{overwrite}{If TRUE overwrite existing directories and files.}
}
\description{
This is special-purpose function to summarize results
from one in-silico experiment that examines how well
signatures can be extracted from synthetic tumors with
correlated SBS1 and SBS5.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AssessExtracted.R
\name{ReadAndAnalyzeSigs}
\alias{ReadAndAnalyzeSigs}
\title{Assess how well extracted signatures match input signatures

We assume that in many cases extraction programs will be run
outside of R on file inputs and will generate fill outputs.}
\usage{
ReadAndAnalyzeSigs(extracted.sigs, ground.truth.sigs, ground.truth.exposures)
}
\arguments{
\item{extracted.sigs}{Path to file containing the extracted signature profiles.}

\item{ground.truth.sigs}{File containing signature profiles from which the
synthetic data were generated.}

\item{ground.truth.exposures}{File containing the exposures from which
the synthetic catalogs were generated.  This file is used to restrict
assessment to only those signatures in \code{ground.truth.sigs}
that were actually represented in the exposures.}
}
\value{
See \code{\link[ICAMSxtra]{MatchSigsAndRelabel}}
}
\description{
Assess how well extracted signatures match input signatures

We assume that in many cases extraction programs will be run
outside of R on file inputs and will generate fill outputs.
}
\details{
Generates output files by calling
\code{\link[ICAMSxtra]{MatchSigsAndRelabel}}
}

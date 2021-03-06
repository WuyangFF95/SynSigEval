% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMuinteraction.R
\name{CreateEMuOutput}
\alias{CreateEMuOutput}
\title{Prepare input file for EMu from a
EMu formatted catalog file.}
\usage{
CreateEMuOutput(
  catalog,
  out.dir = paste0(dirname(catalog), "/ExtrAttr/EMu.results"),
  overwrite = FALSE
)
}
\arguments{
\item{catalog}{a catalog in ICAMS format. It can be
a .csv file, or a matrix or data.frame.
Usually, it refers to \code{"ground.truth.syn.catalog.csv"}.}

\item{out.dir}{Directory that will be created for the output;
abort if it already exists. Usually, the \code{out.dir} will
be a \code{EMu.results} folder directly under the folder
storing \code{catalog}.}

\item{overwrite}{If TRUE, overwrite existing output}
}
\value{
\code{invisible(catalog)},
original catalog in EMu format
}
\description{
Prepare input file for EMu from a
EMu formatted catalog file.
}
\details{
Creates folder named \code{EMu.results} containing catalogs
in EMu-formatted catalogs: Rows are signatures;
the first column is the name of the mutation type, while the remaining
columns are samples (tumors).
These EMu-formatted catalogs will the input when running EMu program
later on compiled binary.
}

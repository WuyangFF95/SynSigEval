% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummarizeMultiRuns.R
\name{SummarizeMultiRuns}
\alias{SummarizeMultiRuns}
\title{Assess/evaluate multiple summarized runs on one dataset by one computational approach.}
\usage{
SummarizeMultiRuns(datasetName, toolName, resultPath, run.names)
}
\arguments{
\item{datasetName}{Name of the dataset. (e.g. "S.0.1.Rsq.0.1").
Usually, it is has the same name as \code{basename(top.dir)}.}

\item{toolName}{Name of computational approach. (e.g. "SigProExtractor")}

\item{resultPath}{Path expected to have multiple result folders each
named as \code{run.names} (e.g. "seed.1").
The example \code{resultPath} is
S.0.1.Rsq.0.1/sp.sp/ExtrAttr/hdp.results/ in old folder structure, or
3a.Original_output_K_unspecified/hdp/S.0.1.Rsq.0.1 in new folder structure.}

\item{run.names}{A character vector records the list of directories which are under
\code{resultPath} and contain results of computational approach, and a \code{summary} folder generated
by \code{\link{SummarizeSigOneExtrAttrSubdir}}.}
}
\value{
A list contain values of measures measures in multiple runs: \itemize{
\item $averCosSim Cosine similarity
\item $truePos True Positives(TP): Ground-truth signatures which are active in
the spectra, and extracted.
\item $falseNeg False Negatives(FN): Ground-truth signatures not extracted.
\item $falsePos False Positives(FP): Signatures wrongly extracted, not resembling
any ground-truth signatures.
\item $TPR True positive rate (TPR, Sensitivity): TP / (TP + FN)
\item $PPV Positive predictive value (PPV): TP / (FP + TP)
\item $cosSim Average cosine similarity to each of the ground-truth signatures.
\item $NumSigsSimilar Number of signatures resembling each of ground-truth signatures,
with pairwise cosine similarity > 0.9.
\item $AggManhattanDist (if exposures of signatures were inferred) Scaled
Manhattan distance between ground-truth and inferred
exposures to each of the ground-truth signatures.
}
This list also contains \code{mean} and \code{sd}, and other
statistics of these measures in \itemize{
\item $fivenum - summary generated by \code{\link[stats]{fivenum}} - columns of this table
refer to Tukey's five number summary for each extraction measure across all runs: \itemize{

\item \code{min} - minimum
\item \code{lower-hinge} - first quartile. Serve as the lower-hinge of the box-whisker plot.
\item \code{median} - median of measure across all runs.
\item \code{upper-hinge} - third quartile. Serve as the upper-hinge of the box-whisker plot.
\item \code{max} - maximum

}

\item $fivenumMD - Tukey's five number summary for aggregately-scaled Manhattan distance.
\item $meanSD - mean and standard deviation for extraction measures.
\item $meanSDMD - mean and standard deviation for aggregately-scaled Manhattan distance.
}
}
\description{
Summarize results from each computational approach in \code{resultPath}/\code{run.names}
(generated by running a computational approach),
combine them into \code{resultPath}.
}
\details{
Also writes multiple files into folder \code{resultPath}.
}

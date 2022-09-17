
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigEval

<!-- badges: start -->

[![R-CMD-check](https://github.com/WuyangFF95/SynSigEval/workflows/R-CMD-check/badge.svg?branch=main)](https://github.com/WuyangFF95/SynSigEval/actions?query=workflow%3AR-CMD-check+branch%3Amain)
<!-- badges: end -->

Assess the performance of mutational-signature analysis programs using
catalogs of synthetic mutational spectra created by package
[**`SynSigGen`**](https://github.com/steverozen/SynSigGen).

This version (0.4.0) and previous version (0.3.1) are suitable for
assessing the extraction accuracy in paper *mSigHdp: hierarchical
Dirichlet processes for mutational signature discovery* - Liu et
al. (2022).

For assessing extraction accuracy on data sets presented by paper
*Accuracy of Mutational Signature Software on Correlated Signatures* -
Wu et al. (2022), please proceed to version
[**`0.2.2`**](https://github.com/WuyangFF95/SynSigEval/tree/0.2.2).

Check NEWS.md for differences between version 0.4.0 and previous
versions.

## Installation

Before installation, prerequisites in
[**`Bioconductor`**](https://www.bioconductor.org/) needs to be
installed:

``` r
install.packages("BiocManager")
BiocManager::install(
  c("Biostrings", "BSgenome", "GenomeInfoDb", "GenomicRanges")
)
```

Install the development version of **`SynSigEval`** from GitHub with the
R command line:

``` r
install.packages("remotes")
# Stable version
remotes::install_github("WuyangFF95/SynSigEval", ref = "v0.4.0-branch")
# Develop version
remotes::install_github("WuyangFF95/SynSigEval", ref = "main")
```

## Reference manual

<https://github.com/WuyangFF95/SynSigEval/blob/v0.4.0-branch/data-raw/SynSigEval_0.4.0.pdf>

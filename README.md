
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigEval

<!-- badges: start -->

[![R-CMD-check](https://github.com/WuyangFF95/SynSigEval/workflows/R-CMD-check/badge.svg)](https://github.com/WuyangFF95/SynSigEval/actions)
<!-- badges: end -->

Assess the performance of mutational-signature analysis programs using
catalogs of synthetic mutational spectra created by package
[SynSigGen](https://github.com/steverozen/SynSigGen).

## Installation

Before installation, prerequisites in
[Bioconductor](https://www.bioconductor.org/) needs to be installed:

``` r
install.packages("BiocManager")
BiocManager::install(c("Biostrings","BSgenome","GenomeInfoDb","GenomicRanges"))
```

Install the development version of SynsigEval from
[GitHub](https://github.com/) with the R command line:

``` r
install.packages("devtools")
devtools::install_github("WuyangFF95/SynSigEval")
```

## Reference manual

<https://github.com/WuyangFF95/SynSigEval/blob/master/data-raw/SynSigEval_0.2.3.pdf>

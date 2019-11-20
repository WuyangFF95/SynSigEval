#' Plot the SBS96 part of a SignatureAnalyzer COMPOSITE signature or catalog
#'
#' @param catalog Catalog or signature matrix
#'
#' @param name Name of file to print to.
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @keywords internal

Plot96PartOfComposite <- function(catalog, name, type = "density") {
  catalog.type <- attr(catalog,"catalog.type")
  region <- attr(catalog,"region")
  cat1536 <- ICAMS::as.catalog(
    catalog[1:1536, ],
    catalog.type = catalog.type,
    region = region)
  cat96 <- ICAMS::Collapse1536CatalogTo96(cat1536)
  all.0 <- which(colSums(cat96) == 0)
  if (length(all.0) > 0 ) {
    cat96[ , all.0] <- 1
    cn <- colnames(cat96)
    cn[all.0] <- paste(cn[all.0], "WARNING all 0")
    colnames(cat96) <- cn
  }
  ICAMS::PlotCatalogToPdf(
    catalog = cat96/sum(cat96),
    file = name)
}

#' Plot the a SignatureAnalyzer COMPOSITE signature or catalog into separate pdfs
#'
#' @param catalog Catalog or signature matrix
#'
#' @param filename.header Contain path and the beginning part of the file name.
#' The name of the pdf files will be:
#' \code{filename.header}.SNS.96.pdf
#' \code{filename.header}.SNS.1536.pdf
#' \code{filename.header}.DNS.78.pdf
#' \code{filename.header}.ID.83.pdf
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @param id A vector containing the identifiers of the samples
#' or signatures in \code{catalog}.
#'
#' @importFrom ICAMS Collapse1536CatalogTo96 PlotCatalogToPdf
#'
#' @export
PlotCatCOMPOSITE <- function(catalog, filename.header, type, id = colnames(catalog)) {

  ## Read in COMPOSITE catalogue
  test.COMPOSITE.sigs <-
    ICAMS::ReadCatalog(catalog)

  ## Check
  stopifnot(nrow(test.COMPOSITE.sigs) == 1697)
  # TODO WUYang: check whether the base context is in correct order

  ## Subsetting COMPOSITE catalogue
  catalog.list <- SplitCatCOMPOSITE(test.COMPOSITE.sigs)

  ## Plot using ICAMS embedded plotting function
  ICAMS::PlotCatalogToPdf(
    catalog.list$SBS1536,
    filename = paste0(filename.header,".SNS.1536.pdf"))
  ICAMS::PlotCatalogToPdf(
    catalog.list$DBS78,
    filename = paste0(filename.header,".DNS.78.pdf"))
  ICAMS::PlotCatalogToPdf(
    catalog.list$ID83,
    filename = paste0(filename.header,".ID.83.pdf"))
}

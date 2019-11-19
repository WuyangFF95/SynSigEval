# ReadWriteCOMPOSITECatalogs.R


#' Split COMPOSITE (SNS1536+DBS78+ID83) catalogues
#' in ICAMS format into 3 individual catalogs.
#' @param catalog Input catalog, can be a .csv file or matrix
#' in ICAMS COMPOSITE format.
#' @return a list, containing 3 catalog matrices in MultiModalMuSig format.
#' Each matrix contains SNS1536, DBS78 and ID83 information, respectively.
SplitCatCOMPOSITE <- function(catalog) {

  # Read COMPOSITE catalog. Either from file or matrix-like
  stopifnot(is.character(catalog) | is.data.frame(catalog) | is.matrix(catalog))
  if(class(catalog) == "character")
    catMatrix <- ReadCatCOMPOSITE(catalog)
  else
    catMatrix <- catalog

  # Split COMPOSITE catalog to 3 catalogues
  # (1 SNS1536, 1 DNS78, 1 Indel83)
  # use function in ReadWriteCatalogs.R
  catList <- list("SNS1536" = catMatrix[1:1536,],
                  "DNS78" = catMatrix[1537:1614,],
                  "ID83" = catMatrix[1615:1697,])

  return(catList)
}

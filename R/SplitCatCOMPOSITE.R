#' Split COMPOSITE (SNS1536+DBS78+ID83) catalogs
#' in ICAMS format into 3 individual catalogs.
#' @param catalog Input catalog, can be a .csv file or matrix
#' in ICAMS COMPOSITE format.
#' @return a list, containing 3 catalog matrices in MultiModalMuSig format.
#' Each matrix contains SNS1536, DBS78 and ID83 information, respectively.
SplitCatCOMPOSITE <- function(catalog) {

  # Read COMPOSITE catalog. Either from file or matrix-like
  stopifnot(is.character(catalog) | is.data.frame(catalog) | is.matrix(catalog))
  if (methods::is(catalog, "character"))
    catMatrix <- ICAMS::ReadCatalog(catalog)
  else
    catMatrix <- catalog

  ref.genome <- attr(catMatrix,"ref.genome")
  catalog.type <- attr(catMatrix,"catalog.type")
  region <- attr(catMatrix,"region")


  # Split COMPOSITE catalog to 3 catalogues
  # (1 SNS1536, 1 DNS78, 1 Indel83)
  # use function in ReadWriteCatalogs.R
  catList <- list("SBS1536" = catMatrix[1:1536,],
                  "DBS78" = catMatrix[1537:1614,],
                  "ID83" = catMatrix[1615:1697,])

  # Subsetting will generate a matrix without "Catalog" type.
  # Use ICAMS::as.catalog() to fix it.
  lapply(catList,ICAMS::as.catalog,
         ref.genome = ref.genome,
         catalog.type = catalog.type,
         region = region)

  return(catList)
}

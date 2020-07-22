# MultiModalMuSigInteraction.R
# Interacting functions for running MultiModalMuSig julia package

#' Convert Catalogs from ICAMS format to MM format
#'
#' @param catalog A catalog matrix in ICAMS format. (SNS/DNS/ID)
#'
#' @return a catalog matrix in MultiModalMuSig format.
#'
#' @export
ICAMSCatalog2MM <- function(catalog) {
  # Read catalog. From matrix-like
  stopifnot(is.data.frame(catalog) | is.matrix(catalog))

  catalog <- data.frame(
    "term" = rownames(catalog),
    catalog,
    check.names = FALSE)
  return(catalog)
}



#' Convert Catalogs (File or Matrix) from MM format to ICAMS format
#'
#' @param cat Input catalog, can be a tab-delimited file
#' or matrix in MultiModalMuSig format.
#'
#' @param region Catalog region. Can be a specific genomic
#' or exonic region, or "unknown".
#' Default: "unknown"
#'
#' @param catalog.type Is the catalog a signature catalog,
#' or a spectrum catalog?
#' Default: "counts.signature"
#'
#' @return a catalog matrix in ICAMS format.
#'
#' @export
MMCatalog2ICAMS <- function(
  cat,
  region = "unknown",
  catalog.type = "counts.signature") {

  # Read MM-formatted catalog. Either from file or matrix-like
  stopifnot(is.character(cat) | is.data.frame(cat) | is.matrix(cat))
  if(class(cat) == "character") {
    catMatrix <- read.table(
      file = cat, header = T,
      sep = "\t", row.names = 1,
      as.is = T, check.names = FALSE)
  } else {
    catMatrix <- cat
    rownames(catMatrix) <- cat[,1]
    catMatrix <- catMatrix[,-1]
  }
  # For some SNS96 catalogs,
  # rownames of catMatrix may be in MM format
  # (A[C->A]G) rather than in ICAMS format (ACGA).
  # In this case, the format will be changed.
  if(nrow(catMatrix)==96 & grepl("->",rownames(catMatrix)[1])){

    oldrowNames <- rownames(catMatrix)
    newrowNames <- character(length(oldrowNames))
    before <- character(length(oldrowNames))
    ref <- character(length(oldrowNames))
    after <- character(length(oldrowNames))
    var <- character(length(oldrowNames))

    for(ii in 1:length(oldrowNames)){
      before[ii] <- strsplit(oldrowNames[ii],"\\[")[[1]][1]

      ref[ii] <- strsplit(oldrowNames[ii],"->")[[1]][1]
      ref[ii] <- strsplit(ref[ii],"\\[")[[1]][2]

      after[ii] <- strsplit(oldrowNames[ii],"\\]")[[1]][2]

      var[ii] <- strsplit(oldrowNames[ii],"->")[[1]][2]
      var[ii] <- strsplit(ref[ii],"\\]")[[1]][1]

      newrowNames[ii] <- paste0(before[ii],ref[ii],after[ii],var[ii])
    }

    rownames(catMatrix) <- newrowNames
  }

  catalog <- ICAMS::as.catalog(
    catMatrix,
    region = region,
    catalog.type = catalog.type)

  return(catalog)
}


#' Read Catalog files in MM format
#' @param exposureFile Input exposure file, can be a tab-delimited
#' text file in MultiModalMuSig format.
#' @return a exposure matrix in ICAMS format.
#'
#' @export
ReadExposureMM <- function(exposureFile){
  exposure <- read.table(
    exposureFile, header = T,
    sep = "\t", row.names = 1,
    as.is = T, check.names = FALSE)
  return(exposure)
}


#' Prepare input file for MultiModalMuSig from a
#' MultiModalMuSig formatted catalog file.
#'
#' @param catalog a catalog in ICAMS format. It can be
#' a .csv file, or a matrix or data.frame.
#' Usually, it refers to \code{"ground.truth.syn.catalog.csv"}.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exists. Usually, the \code{out.dir} will
#' be a \code{MultiModalMuSig.results} folder directly under the
#' folder storing \code{catalog}.
#'
#' @param overwrite If TRUE, overwrite existing output
#'
#' @return \code{invisible(catMatrix)},
#' original catalog in MultiModalMuSig format
#'
#' @details Creates folder named \code{MultiModalMuSig.results} containing catalogs
#' in MultiModalMuSig-formatted catalogs: Rows are signatures;
#' the first column is the name of the mutation type, while the remaining
#' columns are samples (tumors).
#' These MM-formated catalogs will the input when running MultiModalMuSig program
#' later on Julia platform.
#'
#' @export
#'
#' @importFrom utils capture.output
CreateMultiModalMuSigOutput <-
  function(catalog,
           out.dir = paste0(dirname(catalog),"/ExtrAttr/MultiModalMuSig.results"),
           overwrite = FALSE) {

  ## If catalog is a string of file path
  if(is.character(catalog)){
    ## Read in catalog matrix using ICAMS::ReadCatalog.
    catMatrix <- ICAMS::ReadCatalog(catalog, strict = FALSE)
    ## Convert catalog to MM format
    catMatrix <- ICAMSCatalog2MM(catMatrix)
    ## Fetch the name of catalog file without extension
    oldFileName <- tools::file_path_sans_ext(basename(catalog))
  } else if(is.data.frame(catalog) | is.matrix(catalog)){
    ## Convert catalog to MM format
    catMatrix <- ICAMSCatalog2MM(catalog)
    ## Fetch the name of catalog file
    oldFileName <- "ground.truth.syn.catalog"
  }

  ## Create out.dir
  dir.create(out.dir,recursive = T)

	## Dump catMatrix into out.dir
	newFileName <- paste0(out.dir,"/",oldFileName,".tsv")
  utils::write.table(catMatrix, file = newFileName,
              sep = "\t", quote = F, row.names = F)

  invisible(catMatrix)
}



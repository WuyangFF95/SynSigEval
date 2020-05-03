## Load required packages
library(ICAMS)
library(SynSigEval)



## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)

## Create helmsman-formatted input-catalog under
## <dataset.top.level.dir>/sp.sp/ExtrAttr/helmsman.results
for(datasetName in datasetNames){
  for(index in 1:20){
    CreateEMuOutput(catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
                    out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/EMu.results",
                                    "/run.",index),
                    overwrite = T)
  }
}

## Run EMu in C++
## Using Run.EMu.py
reticulate::py_run_file("Run.EMu.py")


## After running, convert EMu-formatted tsv files
## to SynSigEval/ICAMS-formatted csv files.
for(datasetName in datasetNames){
  for(nrun in 1:20){
    ## Grep the names of EMu-formatted files.
    ## When the K is different, the output file would be different.
    ## Extracted signature file name:
    ## _{K}_ml_spectra.txt
    ## Attributed exposure file name:
    ## _{K}_assigned.txt
    resultDir <-
      paste0(datasetName,
            "/sp.sp/ExtrAttr/EMu.results/",
            "run.",nrun,"/")

    files <- list.files(resultDir)

    signatureFile <- files[grep(pattern = "_ml_spectra.txt",x = files)]
    exposureFile <- files[grep(pattern = "_assigned.txt",x = files)]


    ## Convert signatures
    signatures <- SynSigEval::ReadEMuCatalog(
      paste0(resultDir,"/",signatureFile),
      mutTypes = ICAMS::catalog.row.order$SBS96,
      sampleOrSigNames = NULL,
      region = "unknown",
      catalog.type = "counts.signature")
    ICAMS::WriteCatalog(
      signatures,
      paste0(resultDir,"/extracted.signatures.csv"))

    ## Convert exposures
    exposure <- SynSigEval::ReadEMuExposureFile(
      exposureFile = paste0(resultDir,"/",exposureFile),
      sigNames = NULL,
      sampleNames = paste0("TwoCorreSigsGen::",1:500))
    SynSigGen::WriteExposure(
      exposure,
      paste0(resultDir,"/attributed.exposures.csv"))

  }
}



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

## Create MultimodalMuSig-formatted input-catalog under
## <dataset.top.level.dir>/sp.sp/ExtrAttr/MultimodalMuSig.LDA.results and
## <dataset.top.level.dir>/sp.sp/ExtrAttr/MultimodalMuSig.CTM.results
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    CreateMultiModalMuSigOutput(
      catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results",
                       "/seed.",seedInUse),
      overwrite = T)
    CreateMultiModalMuSigOutput(
      catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results",
                       "/seed.",seedInUse),
      overwrite = T)
  }
}

## Run MultiModalMuSig in Julia Language
## Using Run.MultiModalMuSig.jl

## After running, convert MultiModalMuSig-formatted tsv files
## to SynSigEval/ICAMS-formatted csv files.
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Convert signatures
    signatures <- SynSigEval::MMCatalog2ICAMS(
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,
             "/extracted.signatures.tsv"),
      region = "unknown",
      catalog.type = "counts.signature")
    ICAMS::WriteCatalog(
      signatures,
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,
             "/extracted.signatures.csv"))

    ## Convert exposures
    exposure <- SynSigEval::ReadExposureMM(
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,
             "/inferred.exposures.tsv"))
    SynSigGen::WriteExposure(
      exposure,
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,
             "/inferred.exposures.csv"))

  }
}

for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Convert signatures
    signatures <- SynSigEval::MMCatalog2ICAMS(
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,
             "/extracted.signatures.tsv"),
      region = "unknown",
      catalog.type = "counts.signature")
    ICAMS::WriteCatalog(
      signatures,
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,
             "/extracted.signatures.csv"))

    ## Convert exposures
    exposure <- SynSigEval::ReadExposureMM(
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,
             "/inferred.exposures.tsv"))
    SynSigGen::WriteExposure(
      exposure,
      paste0(datasetName,
             "/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,
             "/inferred.exposures.csv"))
  }
}

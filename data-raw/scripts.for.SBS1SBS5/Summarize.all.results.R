require(ICAMS)
require(SynSigEval)

## Specify default options
options(stringsAsFactors = F)

## Specify dataset names
slopes <- c("0.1","0.5","1","2","10")
Rsqs <- c("0.1","0.2","0.3","0.6")
datasetNames <- character(0)

for(slope in slopes){
  for(Rsq in Rsqs){
    datasetNames <- c(datasetNames,
                      paste0("S.",slope,".Rsq.",Rsq))
  }
}
## Specify tool Names
# R-based tools which can do both extraction and attribution,
# excluding SignatureAnalyzer (due to special folder structure)
# and maftools (its seed is hard-coded)
RBasedExtrAttrToolNames <- c("hdp","MutationalPatterns","sigfit.EMu","sigfit.NMF","signeR")
# Python or other language based tools.
# excluding maftools (seed is fixed) and EMu (cannot designate seed)
otherExtrAttrToolNames <- c("MultiModalMuSig.CTM","MultiModalMuSig.LDA")
# Tools can only do attribution
attrToolNames <-
  c("decompTumor2Sig","deconstructSigs","mSigAct",
    "MutationalPatterns","mutSignatures",
    "SignatureEstimation.QP","SignatureEstimation.SA",
    "YAPSA","sigfit.NMF","sigfit.EMu")


## Specify seeds used in analysis.
## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)

## Before summarizing, convert MultiModalMuSig-formatted catalog
## to ICAMS-formatted catalog.
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Do for MultiModalMuSig.LDA
    cat <- MMCatalog2ICAMS(
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,"/extracted.signatures.tsv"),
      region = "genome",
      catalog.type = "counts.signature")
    ICAMS::WriteCatalog(
      cat,
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,"/extracted.signatures.csv"))
    expo <- ReadExposureMM(paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
                           "seed.",seedInUse,"/attributed.exposures.tsv"))
    WriteExposure(
      expo,
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/",
             "seed.",seedInUse,"/attributed.exposures.csv"))

    ## Do for MultiModalMuSig.CTM
    cat <- MMCatalog2ICAMS(
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,"/extracted.signatures.tsv"),
      region = "genome",
      catalog.type = "counts.signature")
    ICAMS::WriteCatalog(
      cat,
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,"/extracted.signatures.csv"))
    expo <- ReadExposureMM(paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
                                  "seed.",seedInUse,"/attributed.exposures.tsv"))
    WriteExposure(
      expo,
      paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/",
             "seed.",seedInUse,"/attributed.exposures.csv"))
  }
}

## Part I: Run Summarize functions in SynSigEval
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Summarize R-based Extraction and attribution tools.
    for(extrAttrToolName in RBasedExtrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrToolName,
                            ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize R-based attribution-only tools.
    for(attrToolName in attrToolNames){
      SynSigEval::SummarizeSigOneAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/Attr/",attrToolName,
                            ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize sigproSS (a Python-based attribution-only tool)
    SynSigEval::SummarizeSigOneSigProSSSubdir(
      run.dir = paste0(datasetName,"/sp.sp/Attr/SigProSS.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize non-R Extraction and attribution tools.
    for(extrAttrToolName in otherExtrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrToolName,
                         ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize SigProExtractor
    SynSigEval::SummarizeSigOneSigProExtractorSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttr/SigProExtractor.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize helmsman
    SynSigEval::SummarizeSigOnehelmsmanSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttr/helmsman.NMF.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize SignatureAnalyzer
    {
      SynSigEval::CopyBestSignatureAnalyzerResult(
        sa.results.dir = paste0(datasetName,
                                "/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.",
                                seedInUse,"/"),
        overwrite = T)

      SynSigEval:::SummarizeSigOneSASubdir(
        run.dir = paste0(datasetName,
                         "/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.",
                         seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        which.run = "/best.run/",
        overwrite = T)
    }
  }
  ## Summarize maftools
  SynSigEval::SummarizeSigOneExtrAttrSubdir(
    run.dir = paste0(datasetName,
                     "/sp.sp/ExtrAttr/maftools.results/seed.123456","/"),
    ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
    overwrite = T)
  ## Summarize EMu
  for(nrun in 1:20){
    SynSigEval::SummarizeSigOneExtrAttrSubdir(
      run.dir = paste0(datasetName,"/sp.sp/ExtrAttr/","EMu",
                       ".results/run.",nrun,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
  }
}

## Part II: Write summary table for 20 seeds for each tool with each dataset
otherExtrAttrToolNames <- c("helmsman.NMF","MultiModalMuSig.CTM","MultiModalMuSig.LDA","SigProExtractor","SignatureAnalyzer")
for(datasetName in datasetNames){
  ## For each dataset, summarize 20 runs
  ## using different seeds by EMu
  for(extrAttrToolName in c(
    RBasedExtrAttrToolNames,otherExtrAttrToolNames)){
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = extrAttrToolName,
      tool.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrToolName,
                          ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
  ## For each dataset, summarize 20 runs
  ## (without seeds) by EMu
  SynSigEval::SummarizeMultiRuns(
    datasetName = datasetName,
    toolName = "EMu",
    tool.dir = paste0(datasetName,"/sp.sp/ExtrAttr/","EMu",
                      ".results/"),
    run.names = paste0("run.",1:20))
  ## For each dataset, summarize 1 run by maftools
  SynSigEval::SummarizeMultiRuns(
    datasetName = datasetName,
    toolName = "maftools",
    tool.dir = paste0(datasetName,"/sp.sp/ExtrAttr/","maftools",
                      ".results/"),
    run.names = paste0("seed.","123456"))
  for(attrToolName in c(attrToolNames,"SigProSS")){
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = attrToolName,
      tool.dir = paste0(datasetName,"/sp.sp/Attr/",attrToolName,
                        ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
}

## Part III: Write Summary table of multiple tools for each dataset
for(slope in slopes){
  for(Rsq in Rsqs){
    datasetName <- paste0("S.",slope,".Rsq.",Rsq)
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/ExtrAttr/"),
      toolName = c(RBasedExtrAttrToolNames,otherExtrAttrToolNames,
                   "EMu","maftools"),
      tool.dirnames = paste0(c(RBasedExtrAttrToolNames,otherExtrAttrToolNames,
                               "EMu","maftools"),".results"),
      datasetGroup = Rsq,
      datasetGroupName = "Pearson's R^2",
      datasetSubGroup = slope,
      datasetSubGroupName = "SBS1:SBS5 mutation count ratio")
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/Attr/"),
      toolName = c(attrToolNames,"SigProSS"),
      tool.dirnames = paste0(c(attrToolNames,"SigProSS"),".results"),
      datasetGroup =  Rsq,
      datasetGroupName = "Pearson's R^2",
      datasetSubGroup = slope,
      datasetSubGroupName = "SBS1:SBS5 mutation count ratio"
      )
  }
}
## Part IV: Generate a summary table and boxplot for results
## of multiple datasets, from each separate tool.
datasetGroup <- rep(c("0.1","0.2","0.3","0.6"),5)
names(datasetGroup) <- datasetNames
datasetSubGroup <- rep(c("0.1","0.5","1","2","10"),each = 4)
names(datasetSubGroup) <- datasetNames


for(toolName in c(
  RBasedExtrAttrToolNames,otherExtrAttrToolNames,
  "maftools","EMu")){
  SummarizeOneToolMultiDatasets(
    dataset.dirs = datasetNames,
    datasetGroup = datasetGroup,
    datasetGroupName = "Pearson's R^2",
    datasetSubGroup = datasetSubGroup,
    datasetSubGroupName = "SBS1:SBS5 mutation count ratio",
    toolName = toolName,
    tool.dirname = paste0("sp.sp/ExtrAttr/",toolName,".results/"),
    out.dir = paste0("FinalToolWiseSummary/ExtrAttr/",toolName,"/"),
    overwrite = T)
}
for(toolName in c(attrToolNames,"SigProSS")){
  SummarizeOneToolMultiDatasets(
    dataset.dirs = datasetNames,
    datasetGroup = datasetGroup,
    datasetGroupName = "Pearson's R^2",
    datasetSubGroup = datasetSubGroup,
    datasetSubGroupName = "SBS1:SBS5 mutation count ratio",
    toolName = toolName,
    tool.dirname = paste0("sp.sp/Attr/",toolName,".results/"),
    out.dir = paste0("FinalToolWiseSummary/Attr/",toolName,"/"),
    overwrite = T)
}



## Part V: Generate a combined summary table for results
## of multiple datasets, from multiple tools
FinalExtrAttr <- SummarizeMultiToolsMultiDatasets(dataset.dirs = datasetNames,
                       second.third.level.dirname = "sp.sp/ExtrAttr",
                       out.dir = "./FinalExtrAttrSummary", overwrite = T)
FinalAttr <- SummarizeMultiToolsMultiDatasets(dataset.dirs = datasetNames,
                       second.third.level.dirname = "sp.sp/Attr",
                       out.dir = "./FinalAttrSummary", overwrite = T)



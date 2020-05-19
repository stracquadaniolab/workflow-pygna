####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")
requiredPackages = c("BiocManager","SummarizedExperiment", "TCGAbiolinks", "org.Hs.eg.db", "recount", "TCGAutils", "limma", "biomaRt")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) {
    res = tryCatch({
      install.packages(p)  
    }, error = function(e) {
      BiocManager::install(p) 
    })
  }
  library(p,character.only = TRUE)
}

elaborateTcga <- function(query) {
  # Start the elaboration on TCGA (there are only tumor samples)
  print("Downloading data from TCGA")
  GDCdownload(query)
  experiment <- GDCprepare(query = query)
  log["PT"] = count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[2])
  log["NT"] = count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[3])
  eset.tcga.cancer = assay(experiment)
  print("Normal tissue found in TCGA")
  print("Performing data normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = eset.tcga.cancer, geneInfo =  geneInfoHT)
  print("Performing data filtering")
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("NT"))
  samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("TP"))
  DEG.ucs <- TCGAanalyze_DEA( mat1 = dataFilt[,samplesNT],
                              mat2 = dataFilt[,samplesTP],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 1,
                              logFC.cut = 0,
                              method = "glmLRT")
  return (DEG.ucs)
}


getTissue <-function(project) {
  data = NULL
  switch (project,
    "TCGA-DLBC" ={data = list("blood","Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Normal")},    #Lymphoma
    "TCGA-LUSC" ={data = list("lung", "Primary Tumor", "Solid Tissue Normal")},                                             #Lung
    "TCGA-LAML" ={data = list("bone_marrow","Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Normal")},    #Leukemia
    "TCGA-LCML" ={data = list("blood","Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Normal")},    #Myelogenous Leukemia
    "TCGA-PRAD" ={data = list("prostate", "Primary Tumor", "Solid Tissue Normal")},                                         #Prostate
    "TCGA-BRCA" ={data = list("breast", "Primary Tumor", "Solid Tissue Normal")}                                            #Breast
    )
  return(data)
}

DATAFOLDER= snakemake@params[["folder"]]
PROJECT = snakemake@params[["name"]]
OUTPUTFILE= snakemake@output[[1]]
LOGFILE = snakemake@output[[2]]
###################################

# Let's query TCGA to see if it has both TP and NP

downloadFromTcga = tryCatch({query =  GDCquery(project = PROJECT,
                                      data.category = "Transcriptome Profiling",
                                      data.type = "Gene Expression Quantification",
                                      workflow.type = "HTSeq - Counts")
          results = query[[1]][[1]]
          sampleTypes = results[c("sample_type")]
          sampleTypes = unique(sampleTypes)
          if ("Primary Tumor" %in% sampleTypes$sample_type) 
            if ("Solid Tissue Normal" %in% sampleTypes$sample_type)
              downloadFromTcga = query
            else 
              downloadFromTcga = "Part"
            else 
              downloadFromTcga = FALSE
          
}, error = function(e) {
  return(FALSE)
})

print(PROJECT)
log = data.frame(matrix(ncol = 6, nrow = 1))
colnames(log) = c("Project", "Tissue", "PT", "PTSource", "NT", "NTSource")
log["Project"] = PROJECT

if (typeof(downloadFromTcga) == "list") {
  # TP from TCGA - NT from TCGA
  log["PTSource"] = "TCGA-GDC"
  log["NTSource"] = "TCGA-GDC"
  DEG.ucs = elaborateTcga(query = query)
} else {
  # NT from GTEX
  data = getTissue(PROJECT)
  log["Tissue"] = data
  log["NTSource"] = "GTEX"
  tissue = data[[1]]
  GTEX = paste("GTEX_",tissue,sep="")
  print("No normal-tissue from TCGA: downloading NT data from GTEX")
  ucs.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
  log["NT"] = ncol(ucs.recount.gtex[[GTEX]])
  eset.gtex<-assays(scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE))$counts
  rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))

  if (downloadFromTcga=="Part") {
    # TP from TCGA with GDC
    print("Downloading data from TCGA using GDC Query")
    log["PTSource"] = "TCGA-GDC"
    GDCdownload(query)
    experiment <- GDCprepare(query = query)
    log["PT"] = count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[2])
    eset.tcga.cancer = assay(experiment)
  } else {
    # TP from TCGA with recount
    print("Downloading data from TCGA using Recount")
    log["PTSource"] = "Recount"
    TCGA = paste("TCGA_",tissue,sep="")
    ucs.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue=tissue)
    log["PT"] = ncol(ucs.recount.tcga[[TCGA]])
    eset.tcga<-assays(scale_counts(ucs.recount.tcga[[TCGA]], round = TRUE))$counts
    colnames(eset.tcga)<-colData(ucs.recount.tcga[[TCGA]])$gdc_cases.samples.portions.analytes.aliquots.submitter_id
    rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))
    eset.tcga.cancer<-eset.tcga[,which(colData(ucs.recount.tcga[[TCGA]])$gdc_cases.samples.sample_type==data[[2]])]
  }
  ##merging data by row names
  print("Performing data preparaton")
  dataPrep.ucs<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0)
  rownames(dataPrep.ucs)<-dataPrep.ucs$Row.names
  dataPrep.ucs$Row.names<-NULL
  print("Performing data normalization")
  dataNorm.ucs <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs, geneInfo = geneInfoHT)
  print("Performing data filtering")
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm.ucs, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  DEG.ucs <- TCGAanalyze_DEA( mat1 = dataFilt[,colnames(eset.gtex)],
                              mat2 = dataFilt[,colnames(eset.tcga.cancer)],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 1 ,
                              logFC.cut = 0,
                              method = "glmLRT")
}
print("-----Dataset creation: Start-----")
dataset = DEG.ucs
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']= tryCatch({
  mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')},
  error = function(e) {
    mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ENSEMBL')}
)
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3 & dataset$FDR<0.01)
OUTPUTFILE = lapply(OUTPUTFILE, tolower)
OUTPUTFILE = str_replace(OUTPUTFILE,"-","_")
write.csv(dataset, OUTPUTFILE)
print("-----Dataset creation: Done-----")

write.csv(log,LOGFILE)

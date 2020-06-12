####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")
requiredPackages = c("BiocManager","SummarizedExperiment", "TCGAbiolinks", "org.Hs.eg.db", "recount", "TCGAutils", "limma", "biomaRt", "stringr")
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
  
  GDCdownload(query, DATAFOLDER)
  experiment <- GDCprepare(query = query)
  
  dataPrep <- TCGAanalyze_Preprocessing(object = experiment, datatype = "HTSeq - Counts")
  
  log["PT"] <<- count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[2])
  log["NT"] <<- count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[3])
  
  tcgaDf.cancer = assay(experiment)
  print("Normal tissue found in TCGA")
  print("Performing data normalization")
  
  dataNorm <- TCGAanalyze_Normalization(tabDF = tcgaDf.cancer, geneInfo =  geneInfoHT, method= "gcContent")
  print("Performing data filtering")
  
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  
  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("NT"))
  samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("TP"))
  DEG <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                         mat2 = dataFilt[,samplesTP],
                         Cond1type = "Normal",
                         Cond2type = "Tumor",
                         fdr.cut = 1,
                         logFC.cut = 0,
                         method = "glmLRT")
  return (DEG)
}


getTissue <-function(project) {
  data = NULL
  switch (project,
    "TCGA-DLBC" ={data = list("blood","Primary Tumor", "Blood Derived Normal")},                                            #Lymphoma
    "TCGA-LCML" ={data = list("bone_marrow","Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Normal")},    #Myelogenous Leukemia
    "TCGA-LAML" ={data = list("bone_marrow","Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Normal")},    #Leukemia
    "TCGA-PRAD" ={data = list("prostate", "Primary Tumor", "Solid Tissue Normal")},                                         #Prostate
    "TCGA-BRCA" ={data = list("breast", "Primary Tumor", "Solid Tissue Normal")},                                           #Breast
    "TCGA-LUSC" ={data = list("lung", "Primary Tumor", "Solid Tissue Normal")},                                             #Lung
    "TCGA-BLCA" ={data = list("bladder","Primary Tumor", "Solid Tissue Normal")}                                            #Bladder
    )
  return(data)
}

PROJECT = snakemake@params[["name"]]
DATAFOLDER = snakemake@params[["raw_data"]]
OUTPUTFILE= snakemake@output[[1]]
LOGFILE = snakemake@output[[2]]
###################################
#PROJECT ="TCGA-BLCA"
#OUTPUTFILE = "blca.csv"
#LOGFILE = "log.csv"

PROJECT = toupper(PROJECT)
PROJECT = str_replace(PROJECT,"_","-")
data = getTissue(PROJECT)


# Let's query TCGA to see if it has both TP and NP
downloadFromTcga = tryCatch({query =  GDCquery(project = PROJECT,
                                      data.category = "Transcriptome Profiling",
                                      data.type = "Gene Expression Quantification",
                                      workflow.type = "HTSeq - Counts")
          results = query[[1]][[1]]
          sampleTypes = results[c("sample_type")]
          sampleTypes = unique(sampleTypes)
          if (data[[2]] %in% sampleTypes$sample_type) 
            if (data[[3]] %in% sampleTypes$sample_type)
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
tissue = data[[1]]
log["Tissue"] = tissue

if (typeof(downloadFromTcga) == "list") {
  # TP from TCGA - NT from TCGA
  log["PTSource"] = "TCGA-GDC"
  log["NTSource"] = "TCGA-GDC"
  DEG = elaborateTcga(query = query)
} else {
  
  # NT from GTEX
  log["NTSource"] = "Recount-GTEX"
  print("No normal-tissue from TCGA: downloading NT data from GTEX")
  if (PROJECT=="TCGA-LAML") {
    GTEX = paste("GTEX_","blood",sep="")
    log["GTEX-Tissue"] = "blood"
    recountGtex<-TCGAquery_recount2(project="GTEX", tissue="blood")
  } else {
    GTEX = paste("GTEX_",tissue,sep="")  
    log["GTEX-tissue"] = tissue
    recountGtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
  }
  log["NT"] = ncol(recountGtex[[GTEX]])
  gtexNt<-assays(scale_counts(recountGtex[[GTEX]], round = TRUE))$counts
  rownames(gtexNt) <- gsub("\\..*", "", rownames(gtexNt))
 
  # TP from TCGA with recount
  print("Downloading data from TCGA using Recount")
  log["PTSource"] = "Recount-TCGA"
  if (PROJECT == "TCGA-DLBC") {
    log["Recount-Tissue"] = "lymph_nodes"
    TCGA = paste("TCGA_","lymph_nodes",sep="")
    recountTcga<-TCGAquery_recount2(project="TCGA", tissue="lymph_nodes")
  } else {
    log["Recount-Tissue"] = tissue
    TCGA = paste("TCGA_",tissue,sep="")
    recountTcga<-TCGAquery_recount2(project="TCGA", tissue=tissue)
  }
  
  if (downloadFromTcga =="Part") {
    if (PROJECT!= "TCGA-DLBC"){
      # TP from TCGA with GDC
      print("Downloading data from TCGA using GDC Query")
      GDCdownload(query)
      experiment <- GDCprepare(query = query)
      log["PT"] = count(experiment@colData@listData[["sample_type"]] == getTissue(PROJECT)[2])
      barcodes = experiment@colData@listData[["submitter_id"]]
      log["PTbarcodefilter"] = "True"
      
      filter = colData(recountTcga[[TCGA]])$gdc_cases.submitter_id %in% barcodes
      recountTcga<-recountTcga[[TCGA]][,filter]
      log["PT"] = ncol(recountTcga)
      tcgaDf<-assays(scale_counts(recountTcga, round = TRUE))$counts
      colnames(tcgaDf)<-colData(recountTcga)$gdc_cases.samples.portions.analytes.aliquots.submitter_id
      rownames(tcgaDf) <- gsub("\\..*", "", rownames(tcgaDf))
      tcgaDf.cancer<-tcgaDf[,which(colData(recountTcga)$gdc_cases.samples.sample_type==data[[2]])]
    } else {
      log["PT"] = ncol(recountTcga[[TCGA]])
      tcgaDf<-assays(scale_counts(recountTcga[[TCGA]], round = TRUE))$counts
      colnames(tcgaDf)<-colData(recountTcga[[TCGA]])$gdc_cases.samples.portions.analytes.aliquots.submitter_id
      rownames(tcgaDf) <- gsub("\\..*", "", rownames(tcgaDf))
      tcgaDf.cancer<-tcgaDf[,which(colData(recountTcga[[TCGA]])$gdc_cases.samples.sample_type==data[[2]])]
    }
  } else {
    log["PT"] = ncol(recountTcga[[TCGA]])
    tcgaDf<-assays(scale_counts(recountTcga[[TCGA]], round = TRUE))$counts
    colnames(tcgaDf)<-colData(recountTcga[[TCGA]])$gdc_cases.samples.portions.analytes.aliquots.submitter_id
    rownames(tcgaDf) <- gsub("\\..*", "", rownames(tcgaDf))
    tcgaDf.cancer<-tcgaDf[,which(colData(recountTcga[[TCGA]])$gdc_cases.samples.sample_type==data[[2]])]
  }
  
  ##merging data by row names
  print("Performing data preparaton")
  dataPrep.ucs<-merge(as.data.frame(gtexNt), as.data.frame(tcgaDf.cancer), by=0)
  rownames(dataPrep.ucs)<-dataPrep.ucs$Row.names
  dataPrep.ucs$Row.names<-NULL
  print("Performing data normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs, geneInfo = geneInfoHT)
  print("Performing data filtering")
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  DEG <- TCGAanalyze_DEA( mat1 = dataFilt[,colnames(gtexNt)],
                              mat2 = dataFilt[,colnames(tcgaDf.cancer)],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 1 ,
                              logFC.cut = 0,
                              method = "glmLRT")
}
print("-----Dataset creation: Start-----")
dataset = DEG
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']= tryCatch({
  mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')},
  error = function(e) {
    mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ENSEMBL')}
)
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3 & dataset$FDR<0.01)
log["significant"] = count(c(dataset["significant"]==1))
log["not_significant"] = count(c(dataset["significant"]==0))
log["total"] = nrow(dataset)
write.csv(dataset, OUTPUTFILE)
print("-----Dataset creation: Done-----")
write.csv(log,LOGFILE)

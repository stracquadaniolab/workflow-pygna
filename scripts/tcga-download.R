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

getTissue <-function(project) {
  tissue = NULL
  switch (project,
    "TCGA-DLBC" ={tissue ="blood"},    #Lymphoma
    "TCGA-LUSC" ={tissue ="lung"},     #Lung
    "TCGA-LAML" ={tissue ="blood"},    #Leukemia
    "TCGA-LCML" ={tissue ="blood"},    #Myelogenous Leukemia
    "TCGA-PRAD" ={tissue ="prostate"}, #Prostate
    "TCGA-BRCA" ={tissue ="breast"}    #Breast
    )
  return(tissue)
}

DATAFOLDER= snakemake@params[["folder"]]
PROJECT = snakemake@params[["name"]]
OUTPUTFILE= snakemake@output[[1]]

###################################
alternate = c("TCGA-DLBC","TCGA-LAML","TCGA-LCML", "TCGA-LUSC")
tissue=getTissue(PROJECT)
GTEX = paste("GTEX_",tissue,sep="")
TCGA = paste("TCGA_",tissue,sep="")

print(PROJECT)
# Start the elaboration on TCGA (there are only tumor samples)
print("Downloading data from TCGA")
query<- GDCquery(project = PROJECT,
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "HTSeq - Counts")
GDCdownload(query)
experiment <- GDCprepare(query = query)
eset.tcga.cancer = assay(experiment)

# Start the elaboration for the GTEX (which is common to all cases, as it contains normal cells)
if (PROJECT %in% alternate) {
  print("No normal-tissue from TCGA: downloading data from GTEX")
  ucs.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
  SE.ucs.recount.gtex <- ucs.recount.gtex[[GTEX]]
  eset.gtex<-assays(scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE))$counts
  rse_scaled <- scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE)
  rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
  ##merging data by row names
  print("Performing data preparaton")
  dataPrep.ucs<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0)
  rownames(dataPrep.ucs)<-dataPrep.ucs$Row.names
  dataPrep.ucs$Row.names<-NULL
  print("Performing data normalization")
  dataNorm.ucs <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs, geneInfo = geneInfoHT)
  print("Performing data filtering")
  dataFilt.ucs <- TCGAanalyze_Filtering(tabDF = dataNorm.ucs, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  DEG.ucs <- TCGAanalyze_DEA( mat1 = dataFilt.ucs[,colnames(eset.gtex)],
                              mat2 = dataFilt.ucs[,colnames(eset.tcga.cancer)],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 1 ,
                              logFC.cut = 0,
                              method = "glmLRT")
} else {
  print("Normal tissue found in TCGA")
  print("Performing data normalization")
  dataNorm <- TCGAanalyze_Normalization(tabDF = eset.tcga.cancer, geneInfo =  geneInfoHT)
  print("Performing data filtering")
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
  print("Performing DEA")
  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("NT"))
  samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("TP"))
  DEG.ucs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                              mat2 = dataFilt[,samplesTP],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 1,
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
write.csv(dataset, OUTPUTFILE)
print("-----Dataset creation: Done-----")

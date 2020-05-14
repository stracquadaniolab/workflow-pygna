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
convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}

getTissue <-function(project) {
  tissue = NULL
  switch (project,
    "TCGA-DLBC" ={tissue ="blood"},
    "TCGA-LUSC" ={tissue ="blood"},
    "TCGA-LAML" ={tissue ="blood"},
    "TCGA-LCML" ={tissue ="blood"},
    "TCGA-PRAD" ={tissue ="prostate"},
    "TCGA-BRCA" ={tissue ="breast"}
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
# Start the elaboration for the GTEX (which is common to all cases, as it contains normal cells)
print("Downloading data from GTEX")
ucs.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
SE.ucs.recount.gtex <- ucs.recount.gtex[[GTEX]]
eset.gtex<-assays(scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE))$counts
rse_scaled <- scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE)
rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))

# Start the elaboration on TCGA (there are only tumor samples)
if (PROJECT %in% alternate) {
  print("Downloading data from TCGA [alternate version]")
  query<- GDCquery(project = PROJECT,
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts")
  GDCdownload(query)
  experiment <- GDCprepare(query = query)
  eset.tcga.cancer = assay(experiment)
} else {
  print("Downloading data from TCGA")
  ucs.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue=tissue)  
  SE.ucs.recount.tcga <- ucs.recount.tcga[[TCGA]]
  eset.tcga<-assays(scale_counts(ucs.recount.tcga[[TCGA]], round = TRUE))$counts
  colnames(eset.tcga)<-colData(ucs.recount.tcga[[TCGA]])$gdc_cases.samples.portions.analytes.aliquots.submitter_id
  rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))
  eset.tcga.cancer<-eset.tcga[,which(colData(ucs.recount.tcga[[TCGA]])$gdc_cases.samples.sample_type=="Primary Tumor")]
}


##merging data by row names
print("Performing data preparaton")
dataPrep.ucs<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0)

rownames(dataPrep.ucs)<-dataPrep.ucs$Row.names
dataPrep.ucs$Row.names<-NULL
print("Performing data normalization")
dataNorm.ucs <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs,
                                          geneInfo = geneInfoHT,
                                          method = "gcContent")
print("Performing data filtering")
dataFilt.ucs <- TCGAanalyze_Filtering(tabDF = dataNorm.ucs,
                                      method = "quantile", 
                                      qnt.cut =  0.25)
print("Performing DEA")
DEG.ucs <- TCGAanalyze_DEA( mat1 = dataFilt.ucs[,colnames(eset.gtex)],
                            mat2 = dataFilt.ucs[,colnames(eset.tcga.cancer)],
                            metadata =FALSE,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 1 ,
                            logFC.cut = 0,
                            method = "glmLRT")


print("-----Dataset creation: Start-----")
dataset = DEG.ucs
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']=mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3 & dataset$FDR<0.01)
write.csv(dataset, OUTPUTFILE)
print("-----Dataset creation: Done-----")

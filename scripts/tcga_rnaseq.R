#options(repos="https://cran.rstudio.com")
requiredPackages = c("BiocManager","SummarizedExperiment", "TCGAbiolinks", "org.Hs.eg.db", "recount", "TCGAutils")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) 
    install.packages(p)
  library(p,character.only = TRUE)
}

convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}

#DATAFOLDER= snakemake@params[["folder"]]
#PROJECT = snakemake@params[["name"]]
#OUTPUTFILE= snakemake@output[[1]]

PROJECT = "TCGA-BRCA"
OUTPUTFILE="results.csv"

sampleType= c("Primary Tumor","Solid Tissue Normal")

print("Downloading data from TCGA...")
print(PROJECT)

query <- GDCquery(project = PROJECT,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
experiment <- GDCprepare(query = query)

alternate = c("TCGA-DLBC","TCGA-LAML","TCGA-LCML", "TCGA-BRCA")

if (PROJECT %in% alternate) {
  tissue="breast"
  GTEX = paste("GTEX_",tissue,sep="")
  
  ucs.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
  SE.ucs.recount.gtex <- ucs.recount.gtex[GTEX]
  
  eset.gtex<-assays(scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE))$counts
  
  rse_scaled <- scale_counts(ucs.recount.gtex[[GTEX]], round = TRUE)
  summary(colSums(assays(rse_scaled)$counts)) / 1e6
  rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
  
  eset.tcga.cancer<-assay(experiment)
  eset.tcga.normal<-eset.gtex
  dataPrep.ucs<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)
  
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs,
                                            geneInfo = geneInfoHT,
                                            method = "gcContent")
} else {
  print("normalization of genes")
  dataNorm <- TCGAanalyze_Normalization(tabDF = assay(experiment), geneInfo =  geneInfo)
  
  print("selection of normal samples 'NT'")
  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                     typesample = c("NT"))
}  

print("quantile filter of genes")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
print("selection of tumor samples 'TP'")
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))
print("Diff.expr.analysis (DEA)")
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 1,
                            logFC.cut = 0,
                            method = "glmLRT")


print("-----Dataset creation: Start-----")
dataset = dataDEGs
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']=mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3 & dataset$FDR<0.01)
write.csv(dataset, OUTPUTFILE)
print("-----Dataset creation: Done-----")
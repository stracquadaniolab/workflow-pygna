#options(repos="https://cran.rstudio.com")

requiredPackages = c("BiocManager","SummarizedExperiment", "TCGAbiolinks", "org.Hs.eg.db")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

DATAFOLDER= snakemake@params[["folder"]]
PROJECT = snakemake@params[["name"]]
OUTPUTFILE= snakemake@output[[1]]

#filename = paste(DATAFOLDER,"blcaExp.rda",sep = "")
if (PROJECT %in% c("TCGA-LAML","TCGA-LCML")) {
  sampleType= c("Primary Blood Derived Cancer - Peripheral Blood","Blood Derived Normal")
} else
  sampleType= c("Primary Tumor","Solid Tissue Normal")

print("Downloading data from TCGA...")
print(PROJECT)
query <- GDCquery(project = PROJECT,
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = sampleType)
print(query)
GDCdownload(query)
experiment <- GDCprepare(query = query)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = assay(experiment), geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
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

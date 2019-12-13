#options(repos="https://cran.rstudio.com")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

DATAFOLDER= snakemake@params[[1]]
OUTPUTFILE= snakemake@output[[1]]

#filename = paste(DATAFOLDER,"blcaExp.rda",sep = "")

print("Downloading data from TCGA...")
query <- GDCquery(project = "TCGA-BLCA",
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
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
print("Dataset creation: Start")
dataset = dataDEGs
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']=mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3 & dataset$FDR<0.01)
write.csv(dataset, OUTPUTFILE)

print("Dataset creation: Done")

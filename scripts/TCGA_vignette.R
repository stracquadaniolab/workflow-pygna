options(repos="https://cran.rstudio.com")
install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

DATAFOLDER= snakemake@input[[1]]
OUTPUTFILE= snakemake@output[[1]]

filename = paste(DATAFOLDER,"brcaExp.rda",sep = "")

print("Downloading data from TCGA...")
query <- GDCquery(project = "TCGA-BRCA",
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query)
brca <- GDCprepare(query = query, save = TRUE, save.filename = filename)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = assay(brca), geneInfo =  geneInfo)

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
                            fdr.cut = 0.01 ,
                            logFC.cut = 0,
                            method = "glmLRT")
print("Dataset creation: Start")
dataset = dataDEGs
symbols = unlist(c(row.names(dataset)))
dataset['genes.Entrezid']=mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
dataset = dataset[order(dataset$logFC),]
dataset["significant"] = as.double(abs(dataset$logFC)>=3)
write.csv(dataset, OUTPUTFILE)

print("Dataset creation: Done")

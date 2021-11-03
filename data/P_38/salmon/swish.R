library(SummarizedExperiment)
library(BiocFileCache)
library(org.Mm.eg.db)
library(fishpond)
library(tximeta)
library(DESeq2)



dir <- "/home/ch219793/pipelines/bulk_rna_seq/selwyn_38/salmon/"
setwd(dir)

coldata <- read.table("salmon_sample_annot.txt", sep="\t",
                      header = T)
names(coldata) <- c('names', 'condition', 'type', 'files')

file_list <- list.files("salmon_GRCm38.p6")
coldata$files <- file.path(dir, "salmon_GRCm38.p6/", file_list, 
                                  "quant.sf")

all(file.exists(coldata$files))

#
se <- tximeta(coldata, type = 'salmon')
assayNames(se)
head(rownames(se))
gse <- summarizeToGene(se)
dds <- DESeqDataSet(gse, ~condition)


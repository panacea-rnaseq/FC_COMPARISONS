wd <- '/Users/sbunga/gitHub/panacea-rnaseq/FC_COMPARISONS/'
setwd(wd)
source(paste0(wd, '/functions.R'))


## Primary Mouse samples
dir <- paste0(wd, '/data/P_38/salmon/')
coldata <- read.table(paste0(dir, "salmon_sample_annot.txt"), sep="\t",
                      header = T)
names(coldata) <- c('names', 'condition', 'type', 'files')
file_str <- '/home/ch219793/pipelines/bulk_rna_seq/selwyn_38/salmon//salmon_vM25-report/quants/'
coldata$files <- str_replace(coldata$files, file_str, '')
coldata$files <- paste0('./data/P_38/salmon/salmon_GRCm38.p6/', coldata$files)
all(file.exists(coldata$files))
se <- tximeta(coldata, type = 'salmon', skipSeqinfo = T)
gse <- summarizeToGene(se)
nrow(gse@assays@data$abundance)
nrow(gse@assays@data$counts)

head(rownames(se))
coldata <- read.csv(paste0(dir, "sample_annot.txt"),
                    sep="\t",row.names="sample")

colnames(gse@assays@data$counts) <- row.names(coldata)
colnames(gse@assays@data$abundance) <- row.names(coldata)
cts_pMm <- gse@assays@data$abundance
len_pMm <- round(gse@assays@data$counts)

cts_symbol_pMm <- as.data.frame(ensembl2symbol(cts_pMm, 'Mouse'))
len_symbol_pMm <- as.data.frame(ensembl2symbol(len_pMm, 'Mouse'))
cts_symbol_pMm$Mm_genes <- row.names(cts_symbol_pMm) 
len_symbol_pMm$Mm_genes <- row.names(len_symbol_pMm) 


write.table(row.names(cts_symbol_pMm), quote = F, 
            file = paste0(dir, "../mgi_symbols.txt"), 
            row.names = F, col.names = F)


## Human iPSC samples
dir <- paste0(wd, '/data/P_39/salmon/')

coldata <- read.table(paste0(dir, "salmon_sample_annot.txt"), sep="\t",
                      header = T)
names(coldata) <- c('names', 'condition', 'type')
files_str <- './data/P_39/salmon/salmon_human_v35-report/quants/'
coldata$files <- paste0(files_str, list.files(files_str), '/quant.sf')
all(file.exists(coldata$files))

se <- tximeta(coldata, type = 'salmon', skipSeqinfo = T)
gse <- summarizeToGene(se)
nrow(gse@assays@data$counts)
nrow(gse@assays@data$counts)
head(rownames(se))

pasAnno <- paste0(dir,"/salmon_sample_annot.txt",sep="")
coldata <- read.csv(pasAnno,sep="\t",row.names="Sample")
colnames(gse@assays@data$counts) <- row.names(coldata)
colnames(gse@assays@data$abundance) <- row.names(coldata)
cts_pHs <- gse@assays@data$abundance
len_pHs <- round(gse@assays@data$counts)
cts_symbol_pHs <- as.data.frame(ensembl2symbol(cts_pHs, 'Human'))
len_symbol_pHs <- as.data.frame(ensembl2symbol(len_pHs, 'Human'))


## Human DRG samples (RAY LAB)
File_1 <- read.table('./data/P_55/File1_Sheet1.txt', sep = "\t", header = T)
conditions_table <- read.csv("./data/P_55/sample_annot.txt",
                             sep="\t",row.names="Sample")

remove_col <- c(13:length(colnames(File_1)))
human_drg_raw <- File_1[, -remove_col]
human_drg_raw <- human_drg_raw[,c(1,4,6:ncol(human_drg_raw))]

cts <- human_drg_raw[-as.integer(row.names(human_drg_raw[duplicated(human_drg_raw$h_gene_short_name), ])), ]
rownames(cts) <- cts$h_ensembl_id_ver

# Clean the table
cts <- cts[, c(2, 3:9)]
conditions_input <- conditions_table
rnames <- c('DRG', 'Spinal_cord', 'Brain_1_Nucleus_accumbens',
            'Brain_1_Caudate_nucleus', 'Brain_2_Hippocampus',
            'Brain_1_Frontal_cortex', 'Heart')
rownames(conditions_input) <- rnames
colnames(cts)[2:8] <- rnames
cts_symbol_dHs <- cts
row.names(cts_symbol_dHs) <- cts$h_gene_short_name
cts_symbol_dHs$h_gene_short_name <- NULL


## Human and Mouse DRG samples
human_mouse_p10 <- readxl::read_xls('./data/P_10/human_mouse_gene_foldchange.xls')

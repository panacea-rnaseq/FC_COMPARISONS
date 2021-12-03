wd <- './gitHub/panacea-rnaseq/FC_COMPARISONS/'
setwd(wd)

### Load Libraries
library(DESeq2)
library(stringr)
library(tximeta)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

# create output directory
dir.create('./final_outputs', showWarnings = F)

# Pre-process the files
source('./pre_process_files.R')
source('./functions.R')

# Read the human mouse orthologue table
hs_mm <- read.table('~/.biomart/human_mouse_ortho.tsv', sep='\t', header = T)
colnames(hs_mm) <- c('Hs_genes', 'Mm_genes')

# Preparing human dataframe
cts_symbol_pHs$Hs_symbol <- row.names(cts_symbol_pHs)
len_symbol_pHs$Hs_symbol <- row.names(len_symbol_pHs)
colnames(len_symbol_pHs)[1:15] <- paste(colnames(len_symbol_pHs)[1:15], 
                                        'read', sep='_')

col_len <- length(colnames(cts_symbol_pHs))
Hs_full_df <- merge(cts_symbol_pHs, len_symbol_pHs, by='Hs_symbol')
Hs_full_df <- cts_symbol_pHs[, c(col_len, 1:(col_len-1))]
row.names(Hs_full_df) <- Hs_full_df$Hs_symbol
col_length <- length(colnames(len_symbol_pHs))

# Creating a  list of human iPSC samples (P_39)
human_iPSC <- list(
  'Nociceptor_4_weeks' = prepare_df(Hs_full_df, len_symbol_pHs, 'Nocicep.*4Week'),
  'Nociceptor_8weeks' = prepare_df(Hs_full_df, len_symbol_pHs, 'Nocicep.*8Week'),
  'Motor' = prepare_df(Hs_full_df, len_symbol_pHs, 'Motor'),
  'Cardio' = prepare_df(Hs_full_df, len_symbol_pHs, 'Cardiomy'),
  'Cortical' = prepare_df(Hs_full_df, len_symbol_pHs, 'Cortical')
)

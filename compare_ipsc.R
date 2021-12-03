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


# Human iPSC samples vs Human iPSC

dir.create(paste0('./final_outputs/P_39_Human_iPSC/'),
           showWarnings = F,
           recursive = T)
WD <- paste0('./final_outputs/P_39_Human_iPSC/')
nociceptor_cols <- c("Nociceptor_4_weeks", "Nociceptor_8weeks")
other_cols <- c("Motor", "Cardio", "Cortical" )

for(hs in nociceptor_cols){
  hs_sample <- hs
  hs_sample_df <- as.data.frame(human_iPSC[hs])
  cnames <- colnames(hs_sample_df)
  colnames(hs_sample_df)[which(cnames ==na.omit(str_extract(cnames, '.*Hs_symbol')))] <- 'Hs_symbol'
  dir.create(paste0(WD, hs_sample), showWarnings = F)
  
  for(other_samples in other_cols){
    other_samples_df <- as.data.frame(human_iPSC[other_samples])
    cnames <- colnames(other_samples_df)
    colnames(other_samples_df)[which(cnames == na.omit(str_extract(cnames,'.*Hs_symbol')))] <- 'Hs_symbol'
    df <- merge(hs_sample_df, other_samples_df, by='Hs_symbol')
    df$FoldChange <- (df[,paste0(hs, ".Hs_Avg_TPM")]+1)/(df[,paste0(other_samples, ".Hs_Avg_TPM")]+1)
    
    write.csv(df, paste0(WD, hs_sample,'/', hs_sample, '_vs_', other_samples,'.csv'), 
              col.names = T, row.names = F)
    
    
    n <- nrow(df)
    # up regulated top 10
    up_reg <- df[order(df$FoldChange, decreasing = T), ][1:round(n*10/100), ]
    
    # down regulated top 10
    down_reg <- df[order(df$FoldChange, decreasing = F), ][1:round(n*10/100), ]
    
    write.csv(up_reg, paste0(WD, hs_sample,'/', hs_sample, '_vs_', other_samples,'_up_reg.csv'), 
              col.names = T, row.names = F)
    write.csv(down_reg, paste0(WD, hs_sample,'/', hs_sample, '_vs_', other_samples,'_down_reg.csv'), 
              col.names = T, row.names = F)
    
    }
}



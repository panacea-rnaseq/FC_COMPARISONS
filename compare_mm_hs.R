library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(stringr)
library(tximeta)

# Run pre-processing
source('./pre_process_files.R')

## Comparing Primary mouse and human iPSC
# Preparing mouse dataframe
hs_mm <- read.table('~/.biomart/human_mouse_ortho.tsv', sep='\t', header = T)
colnames(hs_mm) <- c('Hs_genes', 'Mm_genes')
pMm_full_df <- merge(x=cts_symbol_pMm, y=hs_mm, by.x='Mm_genes')
col_len <- length(colnames(pMm_full_df))

# rearranging the order
pMm_full_df <- pMm_full_df[, c(1,col_len, 2:(col_len-1))]


# Creating a list of all the samples from primary mouse
Mm_DRGs <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df),
                                                                 "DRG.*")))]
len_Mm_DRGs <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "DRG.*")))]

colnames(len_Mm_DRGs)[2:ncol(len_Mm_DRGs)] <- paste(colnames(len_Mm_DRGs)[2:ncol(len_Mm_DRGs)], 
                                                    'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "DRG.*")
Mm_Avg_TPM <- rowMeans(Mm_DRGs[, na.omit(cnames)])
Mm_DRGs <- merge(Mm_DRGs, len_Mm_DRGs, by = 'Mm_genes')
Mm_DRGs$Mm_Avg_TPM <- Mm_Avg_TPM



Mm_SC <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df), "Spinal_.*")))]

len_Mm_SC <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "Spinal_.*")))]

colnames(len_Mm_SC)[2:ncol(len_Mm_SC)] <- paste(colnames(len_Mm_SC)[2:ncol(len_Mm_SC)], 
                                                    'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "Spinal_.*")
Mm_Avg_TPM <- rowMeans(Mm_SC[, na.omit(cnames)])
Mm_SC <- merge(Mm_SC, len_Mm_SC, by = 'Mm_genes')
Mm_SC$Mm_Avg_TPM <- Mm_Avg_TPM


Mm_Brain1 <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df), "Brain_1.*")))]

len_Mm_Brain1 <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "Brain_1.*")))]

colnames(len_Mm_Brain1)[2:ncol(len_Mm_Brain1)] <- paste(colnames(len_Mm_Brain1)[2:ncol(len_Mm_Brain1)], 
                                                        'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "Brain_1.*")
Mm_Avg_TPM <- rowMeans(Mm_Brain1[, na.omit(cnames)])
Mm_Brain1 <- merge(Mm_Brain1, len_Mm_Brain1, by = 'Mm_genes')
Mm_Brain1$Mm_Avg_TPM <- Mm_Avg_TPM


Mm_Brain2 <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df), "Brain_2.*")))]

len_Mm_Brain2 <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "Brain_2.*")))]

colnames(len_Mm_Brain2)[2:ncol(len_Mm_Brain2)] <- paste(colnames(len_Mm_Brain2)[2:ncol(len_Mm_Brain2)], 
                                                        'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "Brain_2.*")   
Mm_Avg_TPM <- rowMeans(Mm_Brain2[, na.omit(cnames)])
Mm_Brain2 <- merge(Mm_Brain2, len_Mm_Brain2, by = 'Mm_genes')
Mm_Brain2$Mm_Avg_TPM <- Mm_Avg_TPM


Mm_Brain3 <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df), "Brain_3.*")))]

len_Mm_Brain3 <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "Brain_3.*")))]

colnames(len_Mm_Brain3)[2:ncol(len_Mm_Brain3)] <- paste(colnames(len_Mm_Brain3)[2:ncol(len_Mm_Brain3)], 
                                                        'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "Brain_3.*")  
Mm_Avg_TPM <- rowMeans(Mm_Brain3[, na.omit(cnames)])
Mm_Brain3 <- merge(Mm_Brain3, len_Mm_Brain3, by = 'Mm_genes')
Mm_Brain3$Mm_Avg_TPM <- Mm_Avg_TPM 


Mm_Heart <- 
  pMm_full_df[, c(1,2,which(colnames(pMm_full_df) == str_extract(colnames(pMm_full_df), "Heart.*")))]

len_Mm_Heart <- 
  len_symbol_pMm[, c(19, which(colnames(len_symbol_pMm) == str_extract(colnames(len_symbol_pMm), "Heart.*")))]

colnames(len_Mm_Heart)[2:ncol(len_Mm_Heart)] <- paste(colnames(len_Mm_Heart)[2:ncol(len_Mm_Heart)], 
                                                        'read', sep="_")
cnames <- str_extract(colnames(pMm_full_df), "Heart.*")
Mm_Avg_TPM <- rowMeans(Mm_Heart[, na.omit(cnames)])
Mm_Heart <- merge(Mm_Heart, len_Mm_Heart, by = 'Mm_genes')
Mm_Heart$Mm_Avg_TPM <- Mm_Avg_TPM



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

# Creating a list of all the samples from human iPSC
Hs_cortical <- 
  Hs_full_df[, c(1,which(colnames(Hs_full_df) == str_extract(colnames(Hs_full_df), "Cortical.*")))]

len_Hs_cortical <- 
  len_symbol_pHs[, c(col_length, which(colnames(len_symbol_pHs) == 
                                         str_extract(colnames(len_symbol_pHs), 
                                                     "Cortical.*")))]

cnames <- str_extract(colnames(Hs_full_df), "Cortical.*")
Hs_Avg_TPM <- rowMeans(Hs_cortical[, na.omit(cnames)])
Hs_cortical$Hs_Avg_TPM <- Hs_Avg_TPM
Hs_cortical <- merge(Hs_cortical, len_Hs_cortical, by = 'Hs_symbol')


Hs_motor <- 
  Hs_full_df[, c(1,which(colnames(Hs_full_df) == str_extract(colnames(Hs_full_df), 
                                                             "Motor.*")))]
len_Hs_motor <- 
  len_symbol_pHs[, c(col_length, which(colnames(len_symbol_pHs) == 
                                         str_extract(colnames(len_symbol_pHs), 
                                                     "Motor.*")))]
cnames <- str_extract(colnames(Hs_full_df), "Motor.*")
Hs_Avg_TPM <- rowMeans(Hs_motor[, na.omit(cnames)])
Hs_motor$Hs_Avg_TPM <- Hs_Avg_TPM
Hs_motor <- merge(Hs_motor, len_Hs_motor, by = 'Hs_symbol')



Hs_nociceptor_8week <- 
  Hs_full_df[, c(1,which(colnames(Hs_full_df) == str_extract(colnames(Hs_full_df), 
                                                             "Nocicep.*8Week.*")))]

len_Hs_nociceptor_8week <- 
  len_symbol_pHs[, c(col_length, which(colnames(len_symbol_pHs) == 
                                         str_extract(colnames(len_symbol_pHs), 
                                                     "Nocicep.*8Week.*")))]
cnames <- str_extract(colnames(Hs_full_df), "Nocicep.*8Week.*")

Hs_Avg_TPM <- rowMeans(Hs_nociceptor_8week[, na.omit(cnames)])
Hs_nociceptor_8week$Hs_Avg_TPM <- Hs_Avg_TPM
Hs_nociceptor_8week <- merge(Hs_nociceptor_8week, len_Hs_nociceptor_8week,
                             by = 'Hs_symbol')




Hs_nociceptor_4week <- 
  Hs_full_df[, c(1,which(colnames(Hs_full_df) == str_extract(colnames(Hs_full_df), 
                                                             "Nocicep.*4Week.*")))]

len_Hs_nociceptor_4week <- 
  len_symbol_pHs[, c(col_length, which(colnames(len_symbol_pHs) == 
                                         str_extract(colnames(len_symbol_pHs), 
                                                     "Nocicep.*4Week.*")))]

cnames <- str_extract(colnames(Hs_full_df), "Nocicep.*4Week.*")
Hs_Avg_TPM <- rowMeans(Hs_nociceptor_4week[, na.omit(cnames)])
Hs_nociceptor_4week$Hs_Avg_TPM <- Hs_Avg_TPM
Hs_nociceptor_4week <- merge(Hs_nociceptor_4week, len_Hs_nociceptor_4week,
                             by = 'Hs_symbol')


Hs_cardio <- 
  Hs_full_df[, c(1,which(colnames(Hs_full_df) == str_extract(colnames(Hs_full_df), 
                                                             "Cardiomy.*")))]

len_Hs_cardio <- 
  len_symbol_pHs[, c(col_length, which(colnames(len_symbol_pHs) == 
                                         str_extract(colnames(len_symbol_pHs), 
                                                     "Cardiomy.*")))]

cnames <- str_extract(colnames(Hs_full_df), "Cardiomy.*")
Hs_Avg_TPM <- rowMeans(Hs_cardio[, na.omit(cnames)])
Hs_cardio$Hs_Avg_TPM <- Hs_Avg_TPM
Hs_cardio <- merge(Hs_cardio, len_Hs_cardio, by = 'Hs_symbol')



## Creating a list of all samples from Human and Mouse DRGs
# Mouse list
#mDRG_full_df <- merge(x=cts_symbol_mDRG, y=mDRG_hgnc_genes, by.x='Mm_genes')
#col_len <- length(colnames(mDRG_full_df))
#mDRG_full_df <- mDRG_full_df[, c(col_len, 1:(col_len-1))]

#cnames <- str_extract(colnames(mDRG_full_df), "mDRG.*")
#colnames(len_symbol_mDRG)[1:3] <- paste(colnames(len_symbol_mDRG)[1:3], 'read', sep="_")
#mDRG_full_df <- merge(mDRG_full_df, len_symbol_mDRG, by = 'Mm_genes')
#mDRG_full_df$mDRG_Avg_TPM <- rowMeans(mDRG_full_df[, na.omit(cnames)])
#colnames(mDRG_full_df)[2] <- 'Hs_symbol'


# Human list
#col_len <- length(colnames(cts_symbol_hDRG))
#hDRG_full_df <- cts_symbol_hDRG[, c(col_len, 1:(col_len-1))]
#cnames <- str_extract(colnames(hDRG_full_df), "hDRG.*")
#colnames(len_symbol_hDRG)[1:3] <- paste(colnames(len_symbol_hDRG)[1:3], 'read', sep="_")
#hDRG_full_df <- merge(hDRG_full_df, len_symbol_hDRG, by = 'Hs_symbol')
#hDRG_full_df$hDRG_Avg_TPM <- rowMeans(hDRG_full_df[, na.omit(cnames)])



## Creating a list of all the samples from Human DRG's
# Preparing human drgs
cts_symbol_dHs$Hs_symbol <- row.names(cts_symbol_dHs)
col_len <- length(colnames(cts_symbol_dHs))
dHs_full_df <- cts_symbol_dHs[, c(col_len, 1:(col_len-1))]


dHs_drg <- 
  dHs_full_df[, c(1, which(colnames(dHs_full_df) == 
                             str_extract(colnames(dHs_full_df), "DRG.*")))]
dHs_drg$dHs_Avg_TPM <- dHs_drg$DRG

dHs_spinal_cord <- 
  dHs_full_df[, c(1, which(colnames(dHs_full_df) == 
                             str_extract(colnames(dHs_full_df), "Spinal_.*")))]
dHs_spinal_cord$dHs_Avg_TPM <- dHs_spinal_cord$Spinal_cord


dHs_brain_1 <- 
  dHs_full_df[, c(1, which(colnames(dHs_full_df) == 
                             str_extract(colnames(dHs_full_df), "Brain_1.*")))]

dHs_brain_1$dHs_Avg_TPM <- rowMeans(dHs_brain_1[, c(2:length(dHs_brain_1))])



dHs_brain_2 <- 
  dHs_full_df[, c(1, which(colnames(dHs_full_df) == 
                             str_extract(colnames(dHs_full_df), "Brain_2.*")))]
dHs_brain_2$dHs_Avg_TPM <- dHs_brain_2$Brain_2_Hippocampus


dHs_Heart <- 
  dHs_full_df[, c(1, which(colnames(dHs_full_df) == 
                             str_extract(colnames(dHs_full_df), "Heart.*")))]
dHs_Heart$dHs_Avg_TPM <- dHs_Heart$Heart




# Creating a list of primary mouse samples (P_38)
primary_mouse_list <- list(
  'Spinal_cord' = Mm_SC,
  'Heart' = Mm_Heart,
  'DRG' = Mm_DRGs,
  'Brain1' = Mm_Brain1,
  'Brain2' = Mm_Brain2,
  'Brain3' = Mm_Brain3
)

# Creating a  list of human iPSC samples (P_39)
human_iPSC <- list(
  'Nociceptor_4_weeks' = Hs_nociceptor_4week,
  'Nociceptor_8weeks' = Hs_nociceptor_8week,
  'Motor' = Hs_motor,
  'Cardio' = Hs_cardio,
  'Cortical' = Hs_cortical
)


# Creating a list of human DRG samples (P_52)
human_drg <- list(
  'DRG' = dHs_drg,
  'Spinal_cord' = dHs_spinal_cord,
  'Brain_1' = dHs_brain_1,
  'Brain_2' = dHs_brain_2,
  'Heart' = dHs_Heart
)

# Creating a list for human and mouse DRG (P 10)
mouse_drg_p10 <- list(
  'mDRG' = human_mouse_p10[,c(2,5,14,15,16,17,18,19,20)]
)
human_drg_p10 <- list(
  'hDRG' = human_mouse_p10[,c(2,7,8,9,10,11,12,13)]
)
colnames(mouse_drg_p10$mDRG)[2] <- 'Mm_genes'
colnames(mouse_drg_p10$mDRG)[9] <- 'mDRG_Avg_TPM'
colnames(human_drg_p10$hDRG)[1] <- 'Hs_symbol'
colnames(human_drg_p10$hDRG)[8] <- 'Avg_TPM'

## Primary Human VS ALL
dir.create(paste0(wd,'/outputs/P_39_Primary_Human'),
           showWarnings = F,
           recursive = T)


# Human iPSC samples vs Primary mouse
hs_count = 1
mm_count = 1

dir.create(paste0(wd,'/outputs/P_39_Human_iPSC/Human_iPSC_vs_Primary_Mouse'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd,'/outputs/P_39_Human_iPSC/Human_iPSC_vs_Primary_Mouse/')
for(hs in human_iPSC){
  hs_sample <- names(human_iPSC)[hs_count]
  dir.create(paste0(WD, hs_sample), showWarnings = F)
  for(mm in primary_mouse_list){
    mm_sample <- names(primary_mouse_list)[mm_count]
    colnames(mm)[2] <- 'Hs_symbol'
    df <- merge(hs, mm, by='Hs_symbol')
    df$FoldChange <- (df$Hs_Avg_TPM+1)/(df$Mm_Avg_TPM+1)
    write.csv(df, paste0(WD, hs_sample,'/', hs_sample, '_vs_', mm_sample,'.csv'), 
              col.names = T, row.names = F)
    mm_count = mm_count+1
  }
  mm_count=1
  hs_count = hs_count+1
}




## Primary mouse vs ALL
dir.create(paste0(wd,'/outputs/P_38_Primary_Mouse'),
           showWarnings = F,
           recursive = T)

# Primary mouse vs human iPSC samples
mm_count = 1
hs_count = 1
dir.create(paste0(wd,'/outputs/P_38_Primary_Mouse/Primary_mouse_vs_Human_iPSC'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd,'/outputs/P_38_Primary_Mouse/Primary_mouse_vs_Human_iPSC/')
for(mm in primary_mouse_list){
  mm_sample <- names(primary_mouse_list)[mm_count]
  colnames(mm)[2] <- 'Hs_symbol'
  dir.create(paste0(WD, mm_sample), showWarnings = F)
  for(hs in human_iPSC){
    hs_sample <- names(human_iPSC)[hs_count]
    df <- merge(mm, hs, by='Hs_symbol')
    df$FoldChange <- (df$Mm_Avg_TPM+1)/(df$Hs_Avg_TPM+1)
    write.csv(df, paste0(WD, mm_sample,'/', mm_sample, '_vs_', hs_sample,'.csv'), 
              col.names = T, row.names = F)
    hs_count = hs_count+1
    }
  hs_count=1
  mm_count = mm_count+1
}


# Primary mouse vs human DRG samples
mm_count = 1
hs_count = 1
dir.create(paste0(wd,'/outputs/P_38_Primary_Mouse/Primary_mouse_vs_Human_DRG'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd, '/outputs/P_38_Primary_Mouse/Primary_mouse_vs_Human_DRG/')
for(mm in primary_mouse_list){
  mm_sample <- names(primary_mouse_list)[mm_count]
  colnames(mm)[2] <- 'Hs_symbol'
  dir.create(paste0(WD, mm_sample), showWarnings = F)
  for(hs in human_drg){
    hs_sample <- names(human_drg)[hs_count]
    df <- merge(mm, hs, by='Hs_symbol')
    df$FoldChange <- (df$Mm_Avg_TPM+1)/(df$dHs_Avg_TPM+1)
    write.csv(df, paste0(WD, mm_sample,'/', mm_sample, '_vs_', hs_sample,'.csv'), 
              col.names = T, row.names = F)
    hs_count = hs_count+1
  }
  hs_count=1
  mm_count = mm_count+1
}




# Human DRG vs Human iPSC samples
dHs_count = 1
hs_count = 1
dir.create(paste0(wd,'/outputs/P_52_Human_DRG_Ray_Lab/Human_DRG_vs_Human_iPSC'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd,'/outputs/P_52_Human_DRG_Ray_Lab/Human_DRG_vs_Human_iPSC/')
for(dh in human_drg){
  dHs_sample <- names(human_drg)[dHs_count]
  print(dHs_sample)
  dir.create(paste0(WD, dHs_sample), 
             showWarnings = F)
  for(hs in human_iPSC){
    hs_sample <- names(human_iPSC)[hs_count]
    df <- merge(dh, hs, by='Hs_symbol')
    df$FoldChange <- (df$dHs_Avg_TPM+1)/(df$Hs_Avg_TPM+1)
    write.csv(df, paste0(WD, dHs_sample,'/', dHs_sample, '_vs_', hs_sample,'.csv'), 
              col.names = T, row.names = F)
    hs_count = hs_count+1
  }
  hs_count=1
  dHs_count = dHs_count+1
}





## Human DRG in-house vs All
# Human DRG in-house vs Primary mouse
hs_count = 1
mm_count = 1


dir.create(paste0(wd,'/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Primary_Mouse'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd, '/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Primary_Mouse')
for(hs in human_drg_p10){
  hs_sample <- names(human_drg_p10)[hs_count]
  dir.create(paste0(WD, '/', hs_sample), showWarnings = F)
  for(mm in primary_mouse_list){
    mm_sample <- names(primary_mouse_list)[mm_count]
    colnames(mm)[2] <- 'Hs_symbol'
    df <- merge(hs, mm, by='Hs_symbol')
    df$FoldChange <- (df$Avg_TPM+1)/(df$Mm_Avg_TPM+1)
    write.csv(df, paste0(WD, '/', hs_sample,'/', hs_sample, '_vs_', mm_sample,'.csv'), 
              col.names = T, row.names = F)
    mm_count = mm_count+1
  }
  mm_count=1
  hs_count = hs_count+1
}

# Human DRG in-house vs Human iPSC
hs_count = 1
ipsc_count = 1

dir.create(paste0(wd,'/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Human_iPSC'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd, '/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Human_iPSC')
for(hs in human_drg_p10){
  hs_sample <- names(human_drg_p10)[hs_count]
  dir.create(paste0(WD, '/', hs_sample), showWarnings = F)
  for(mm in human_iPSC){
    mm_sample <- names(human_iPSC)[ipsc_count]
    df <- merge(hs, mm, by='Hs_symbol')
    df$FoldChange <- (df$Avg_TPM+1)/(df$Hs_Avg_TPM+1)
    write.csv(df, paste0(WD, '/', hs_sample,'/', hs_sample, '_vs_', mm_sample,'.csv'), 
              col.names = T, row.names = F)
    ipsc_count = ipsc_count+1
  }
  ipsc_count=1
  hs_count = hs_count+1
}



# Human DRG in-house vs Human DRG (RAY LAB)
hs_count = 1
human_count = 1

dir.create(paste0(wd,'/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Human_DRG_Ray_Lab'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd, '/outputs/P_10_Human_DRG_IN_HOUSE/Human_DRG_vs_Human_DRG_Ray_Lab')
for(hs in human_drg_p10){
  hs_sample <- names(human_drg_p10)[hs_count]
  colnames(hs)[c(1,8)] <- c('Hs_symbol', 'hDRG_Avg_TPM')
  dir.create(paste0(WD, '/', hs_sample), showWarnings = F)
  for(mm in human_drg){
    mm_sample <- names(human_drg)[human_count]
    df <- merge(hs, mm, by='Hs_symbol')
    df$FoldChange <- (df$hDRG_Avg_TPM+1)/(df$dHs_Avg_TPM+1)
    write.csv(df, paste0(WD, '/', hs_sample,'/', hs_sample, '_vs_', mm_sample,'.csv'), 
              col.names = T, row.names = F)
    human_count = human_count+1
  }
  human_count=1
  hs_count = hs_count+1
}



# DRG Mouse in-house vs Primary Mouse
mDRG_count = 1
mm_count = 1

dir.create(paste0(wd,'/outputs/P_10_Mouse_DRG_IN_HOUSE/Mouse_DRG_in_house_vs_Primary_Mouse'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd, '/outputs/P_10_Mouse_DRG_IN_HOUSE/Mouse_DRG_in_house_vs_Primary_Mouse')
for(mDRG in mouse_drg_p10){
  mDRG_sample <- names(mouse_drg_p10)[mDRG_count]
  dir.create(paste0(WD, '/', mDRG_sample), showWarnings = F)
  for(mm in primary_mouse_list){
    mm_sample <- names(primary_mouse_list)[mm_count]
    df <- merge(mDRG, mm, by='Mm_genes')
    df$FoldChange <- (df$mDRG_Avg_TPM+1)/(df$Mm_Avg_TPM+1)
    write.csv(df, paste0(WD, '/', mDRG_sample,'/', mDRG_sample, '_vs_', mm_sample,'.csv'), 
              col.names = T, row.names = F)
    mm_count = mm_count+1
  }
  mm_count=1
  mDRG_count = mDRG_count+1
}


## Primary Human VS ALL
dir.create(paste0(wd,'/outputs/P_39_Primary_Human'),
           showWarnings = F,
           recursive = T)

# Human iPSC samples vs Human iPSC
hs_count = 1
mm_count = 1

dir.create(paste0(wd,'/outputs/P_39_Human_iPSC/Human_iPSC_vs_Human_iPSC'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd,'/outputs/P_39_Human_iPSC/Human_iPSC_vs_Human_iPSC/')
for(hs in human_iPSC){
  hs_sample <- names(human_iPSC)[hs_count]
  dir.create(paste0(WD, hs_sample), showWarnings = F)
  for(mm in human_iPSC){
    if(all(names(mm) %in% names(hs)) == T){
      mm_count = mm_count+1
      next
    }else{
      mm_sample <- names(human_iPSC)[mm_count]
      #colnames(mm)[2] <- 'Hs_symbol'
      colnames(mm)[which(colnames(mm) == 'Hs_Avg_TPM')] <- 'Hs_Avg__TPM'
      df <- merge(hs, mm, by='Hs_symbol')
      df$FoldChange <- (df$Hs_Avg_TPM+1)/(df$Hs_Avg__TPM+1)
      write.csv(df, paste0(WD, hs_sample,'/', hs_sample, '_vs_', mm_sample,'.csv'), 
                col.names = T, row.names = F)
      mm_count = mm_count+1
    }
  }
  mm_count=1
  hs_count = hs_count+1
}




# Human DRG vs Human DRG
dHs_count = 1
hs_count = 1
dir.create(paste0(wd,'/outputs/P_52_Human_DRG_Ray_Lab/Human_DRG_vs_Human_DRG'),
           showWarnings = F,
           recursive = T)
WD <- paste0(wd,'/outputs/P_52_Human_DRG_Ray_Lab/Human_DRG_vs_Human_DRG/')
for(dh in human_drg){
  dHs_sample <- names(human_drg)[dHs_count]
  print(dHs_sample)
  dir.create(paste0(WD, dHs_sample), 
             showWarnings = F)
  for(hs in human_drg){
    if(all(names(dh) %in% names(hs)) == T){
      hs_count = hs_count+1
      next
    }else{
      hs_sample <- names(human_drg)[hs_count]
      print(hs_sample)
      colnames(hs)[which(colnames(hs) == 'dHs_Avg_TPM')] <- 'dHs_Avg__TPM'
      df <- merge(dh, hs, by='Hs_symbol')
      df$FoldChange <- (df$dHs_Avg_TPM+1)/(df$dHs_Avg__TPM+1)
      write.csv(df, paste0(WD, dHs_sample,'/', dHs_sample, '_vs_', hs_sample,'.csv'), 
                col.names = T, row.names = F)
      hs_count = hs_count+1
    }
  }
  hs_count=1
  dHs_count = dHs_count+1
}

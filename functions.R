prepare_df <- function(full_df, len_full_df, kword){
  full_df_subset <- 
    full_df[, c(1,which(colnames(full_df) == str_extract(colnames(full_df), paste0(kword, ".*"))))]
  
  len_Hs <- 
    len_full_df[, c(col_length, which(colnames(len_full_df) == 
                                           str_extract(colnames(len_full_df), 
                                                       paste0(kword, ".*"))))]
  
  cnames <- str_extract(colnames(full_df), paste0(kword, ".*"))
  Hs_Avg_TPM <- rowMeans(full_df_subset[, na.omit(cnames)])
  full_df_subset$Hs_Avg_TPM <- Hs_Avg_TPM
  full_df_subset <- merge(full_df_subset, len_Hs, by = 'Hs_symbol')
}


ensembl2symbol <-function(cts, sp){
  ens <- rownames(cts)
  ens %<>% as.character %>% gsub("\\.\\d+", "", .)
  
  if(sp == 'Human'){
    symbols <- mapIds(org.Hs.eg.db, keys = ens,
                      column = c('SYMBOL'), keytype = 'ENSEMBL')
  }
  if(sp == 'Mouse'){
    symbols <- mapIds(org.Mm.eg.db, keys = ens,
                      column = c('SYMBOL'), keytype = 'ENSEMBL')
  }
  #symbols
  
  rownames(cts) <-symbols
  keep <- !is.na(rownames(cts))
  cts <- cts[keep,]
  return(cts)
}



# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
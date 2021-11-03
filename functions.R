

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
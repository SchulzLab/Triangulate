args <- commandArgs(trailingOnly= T)
#file.format <- args[1] ## Either: 1) _Roadmap or 2) and empty string
#data.name <- args[2]
#data.path <- args[3]

feat.input <- args[1]
resp.input <- args[2]
feat.output <- args[3]
resp.output <- args[4]
############ Functions ############


###################################
###################################

library(biomaRt)

###################################
###################################


#feat <- read.csv(paste(data.path, "/scMTL_", data.name, "_feature_reduced_var", file.format, ".txt", sep= ""))
#response <- read.csv(paste(data.path, "/scMTL_", data.name, "_response_reduced_var", file.format, ".txt", sep= ""))

feat <- read.csv(feat.input)
response <- read.csv(resp.input)

rownames(feat) <- feat$X
rownames(response) <- response$X

feat <- feat[, -1]
response <- response[, -1]
##
print(dim(feat))
print(dim(response))

tfs <- colnames(feat)


###############################

#mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=84)) ## It used to work like this, but at some point the functions started having issues with the version argument
mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", host= "http://Mar2016.archive.ensembl.org")) ## The host provieded here corresponds to the version 84 of Ensembl database

G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"), values= tfs, mart= mart)


complex.tfs <- which(sapply(seq(length(tfs)), function(i)grepl("\\.\\.", tfs[i])) == T)
complex.tf.cnt <- length(complex.tfs)

complex.tfs.tokenized <- sapply(seq(complex.tf.cnt), function(i)strsplit(tfs[i], "\\.\\."))

tfs.tokens <- unique(unlist(sapply(seq(complex.tf.cnt), function(i)unlist(strsplit(tfs[i], "\\.\\.")))))
G_list_tokens <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"), values= tfs.tokens, mart= mart)

mapping <- unique(rbind(G_list, G_list_tokens))

gene.ids <- rownames(response)
gene.ids.df <- data.frame(ensembl_gene_id= gene.ids)
overlapping.genes <- merge(mapping, gene.ids.df, by= "ensembl_gene_id")

threshold <- floor(ncol(response) * .9)
remove.tfs <- NULL
for(i in seq(nrow(overlapping.genes))){
  lowly.expr <- which(response[overlapping.genes$ensembl_gene_id[i], ] <= 1)
  if(length(lowly.expr) >= threshold)
    remove.tfs <- c(remove.tfs, i)
}

gene.hits <- which(rownames(response) %in% overlapping.genes$ensembl_gene_id[remove.tfs])
if(length(gene.hits) > 0)
  response <- response[-gene.hits, ]

## TFs to be removed from the features
removableTFs <- which(tfs %in% overlapping.genes$hgnc_symbol[remove.tfs]);
for(i in seq(length(remove.tfs))){
  for(j in seq(length(complex.tfs.tokenized))){
    hit <- which(complex.tfs.tokenized[[j]] == overlapping.genes$hgnc_symbol[remove.tfs[i]]);
    if(length(hit) > 0)
      removableTFs <- c(removableTFs, complex.tfs[j])
  }
}

removableTFs <- unique(removableTFs)
print(removableTFs)
if(length(removableTFs) > 0)
  feat <- feat[-gene.hits, -removableTFs]

print(dim(feat))
print(dim(response))

#write.csv(feat, file= paste(data.path, "/scMTL_", data.name, "_feature_reduced_var", file.format, ".txt", sep= ""))
#write.csv(response, file= paste(data.path, "/scMTL_", data.name, "_response_reduced_var", file.format, ".txt", sep= ""))
write.csv(feat, file= feat.output)
write.csv(response, file= resp.output)

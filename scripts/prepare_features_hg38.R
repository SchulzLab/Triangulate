###############################################
## This script prepares the feature data corresponding to the genes provided in the gene expression data. The arguments are:
## 1) path to the input file containing the gene expression measurements in single cell (either RData from the monocle run or a matrix of genes (rows) and cells (columns)),
## 2) parameter *monocle_param* that indicates whether the monocle object was created from **scMTL\_pipeline\_part1.sm**,
## 3) path to the feature file,
## 4) path to the response file,
## 5) path to where the MST resulted from the monocle run should be saved.
###############################################
library(monocle)
args <- commandArgs(trailingOnly= T)
rdata.file <- args[1]
is.monocle <- args[2]
feature.outfile <- args[3]
response.outfile <- args[4]
mst.outfile <- args[5]
is.monocle <- ifelse(tolower(is.monocle) == "true" || tolower(is.monocle) == "t" , T, F)

feature_type <- strsplit(strsplit(feature.outfile, "scMTL_")[[1]][2], "_")[[1]][3]
sample_type <- strsplit(strsplit(feature.outfile, "scMTL_")[[1]][2], "_")[[1]][1]

ATAC_file_prefix <- "JAMM_singleCell_ATAC_CEL_PHH/" 

if(feature_type == "static"){
  tepic <- read.table("Static_hg38_PromFeatures_2kb_Gene_View.txt", header= T)
  if(sample_type == "HSMM"){
    tepic <- read.table("Single_Cell_TF_Annotation_2kb_Static_TEPIC_03_13_18_10_04_02_177323167_Affinties_2kb_Gene_windows.txt", header= T)
  }
}else if(feature_type == "epigenetic"){
  if(sample_type == "EndoMT" || sample_type == "EndoMT_avg_mat"){
    tepic <- read.table("HUVEC_hg38_3kb_DNase_TEPIC_04_08_19_17_43_21_129943828_Peak_Features_Affinity_Gene_View_Filtered.txt", header= T)
  }else{
    tepic <- read.table("01_HepG2_DNase_TEPIC_05_03_17_15_04_25_Three_Peak_Based_Features_Decay_Affinity_Gene_View_Filtered_Integrated.txt", header= T)
    tepic <- tepic[, -ncol(tepic)]
  }
}else if(feature_type == "dynamic"){
  tepic <- read.table("HepG2_ChIP_Extended_hg38_3kb_Gene_View.txt", header= T, stringsAsFactors= F)
  ## remove duplicates
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacHLC"){
  tepic <- read.table("JAMM_ATAC_HLC_3kb_TEPIC_07_05_19_12_53_25_421001953_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacAll"){
  tepic <- read.table("JAMM_ATAC_All_3kb_TEPIC_07_05_19_12_53_25_421001953_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacHLCR2"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_CEL_R2_HLC_909_3kb_TEPIC_05_03_19_09_33_28_248716596_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacHLCR3"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_CEL_R3_HLC_912_3kb_TEPIC_05_03_19_09_27_48_499437728_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacHLCR4"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_CEL_R4_HLC_915_3kb_TEPIC_05_03_19_09_21_45_061489181_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacPHHD2"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_PHH_AFJ_D2_1288_3kb_TEPIC_05_03_19_10_36_18_937398358_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacPHHD1"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_PHH_DJJ_D1_1287_3kb_TEPIC_05_06_19_13_06_20_046536308_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}else if(feature_type == "atacPHHD3"){
  tepic <- read.table(paste(ATAC_file_prefix, "/JAMM_ATAC_PHH_IAN_D3_1289_3kb_TEPIC_05_03_19_09_19_30_576286132_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt", sep= ""), header= T, stringsAsFactors= F)
  dupl.genes <- duplicated(tepic$geneID)
  tepic <- tepic[!dupl.genes, ]
}



print("done reading tepic data")
if(is.monocle){
  load(rdata.file)
  print("done loading the monocle file")
  ensgs <- sapply(seq(nrow(fData(HSMM))), function(i)strsplit(as.character(fData(HSMM)$gene_short_name[i]), split="\\.")[[1]][1])

  hit.idx <- sapply(seq(length(ensgs)), function(i)which(tepic$geneID == ensgs[i]))

  tepic_reordered <- tepic[unlist(hit.idx), ]
  missing.genes <- which(sapply(seq(length(hit.idx)), function(i)length(hit.idx[[i]])) == 0)
  if(length(missing.genes) != 0){
    exprs_filtered <- exprs(HSMM)[-missing.genes, ]
    ensgs <- ensgs[-missing.genes]
  }else{
    exprs_filtered <- exprs(HSMM)
  }

  monocle_MST <- HSMM_PT@auxOrderingData[[HSMM_PT@dim_reduce_type]]$pr_graph_cell_proj_tree
  monocle_MST_adj <- igraph::as_adjacency_matrix(monocle_MST)

  print("writing the results to file...")
  write.table(as.matrix(monocle_MST_adj), file= mst.outfile, row.names= pData(HSMM)$timepoint, col.names= pData(HSMM)$timepoint)
  print("done writing MST")

  write.table(tepic_reordered, feature.outfile, col.names= T, row.names= F)
  print("done writing features")

  write.table(exprs_filtered, response.outfile, col.names= pData(HSMM)$timepoint, row.names= ensgs)
  print("done writing response")
}else{

  file_type <- substr(rdata.file, nchar(rdata.file) - 2, nchar(rdata.file))
  if(file_type == "csv"){
    expr <- read.csv(rdata.file)
    rownames(expr) <- expr$X
    expr <- expr[, -1]
  }else{
    expr <- read.table(rdata.file)
  }
  ensgs <- sapply(seq(nrow(expr)), function(i)strsplit(rownames(expr)[i], split= "\\.")[[1]][1])
  hit.idx <- sapply(seq(length(ensgs)), function(i)which(tepic$geneID == ensgs[i]))

  tepic_reordered <- tepic[unlist(hit.idx), ]
  missing.genes <- which(sapply(seq(length(hit.idx)), function(i)length(hit.idx[[i]])) == 0)
  if(length(missing.genes) != 0){
    exprs_filtered <- expr[-missing.genes, ]
    ensgs <- ensgs[-missing.genes]
  }else{
    exprs_filtered <- expr
  }
  print("writing the results to file...")
  ## Since monocle is not run, the MST info is not available, but for the snakemake pipeline there should be some output.
  write.table("No Monocle Run -> Empty MST", file= mst.outfile)

  write.table(tepic_reordered, feature.outfile, col.names= T, row.names= F)
  print("done writing features")

  write.table(exprs_filtered, response.outfile, row.names= ensgs)
  print("done writing response")
}

library(monocle)
library(stringr)

args <- commandArgs(trailingOnly= T)
input <- args[1]
output1 <- args[2]
output2 <- args[3]

is.count <- which(tolower(unlist(strsplit(output1, "_"))) == "count") > 0
if(length(is.count) == 0)
  is.count <- F
is.csv <- F
if(str_count(input, "\\.csv$") > 0)
  is.csv <- T

## create a cellDataSet object
#HSMM_expr_matrix <- read.table("GSE52529_fpkm_matrix.txt") ## The original HSMM data downloaded from the repo


#HSMM_expr_matrix <- read.table("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/G_MTL/monocle/scMTL_HSMM_response.txt") ## The filtered data excluding the X & Y chr genes to be consistent with TEPIC provides
if(!is.csv){
  HSMM_expr_matrix <- read.table(input) ## The filtered data excluding the X & Y chr genes to be consistent with TEPIC provides
}else{
  HSMM_expr_matrix <- read.csv(input)
  rownames(HSMM_expr_matrix) <- HSMM_expr_matrix$X
  HSMM_expr_matrix <- HSMM_expr_matrix[, -1]
}

## This section aims to address the Inf values found in the iPCS dataset delivered by Kathrin. She said that some genes with read count of 1 turned out to end up having Inf TPM values. So, they should be set to 0.
inf.idx <- which(is.infinite(as.matrix(HSMM_expr_matrix)) == T)
if(length(inf.idx) > 0){
  HSMM_expr_matrix_corrected <- do.call(data.frame, lapply(HSMM_expr_matrix, function(x) replace(x, is.infinite(x), 0)))
  rownames(HSMM_expr_matrix_corrected) <- rownames(HSMM_expr_matrix)
  HSMM_expr_matrix <- HSMM_expr_matrix_corrected
}
####


HSMM_expr_matrix_original <- HSMM_expr_matrix

cols_split <- sapply(seq(ncol(HSMM_expr_matrix)),function(i) strsplit(colnames(HSMM_expr_matrix)[i],"_"))

HSMM_sample_sheet <- data.frame(timepoint= colnames(HSMM_expr_matrix))
HSMM_sample_sheet <- data.frame(timepoint= colnames(HSMM_expr_matrix), Hours= sapply(seq(length(cols_split)), function(i) cols_split[[i]][1]))

#HSMM_sample_sheet <- data.frame(timepoint= seq(colnames(HSMM_expr_matrix)))
HSMM_gene_annotation <- data.frame(gene_short_name= rownames(HSMM_expr_matrix))

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
colnames(HSMM_expr_matrix) <- 1:ncol(HSMM_expr_matrix)
rownames(HSMM_expr_matrix) <- 1:nrow(HSMM_expr_matrix)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),phenoData = pd, featureData = fd)

## Choosing a distribution for your data
#HSMM <- newCellDataSet(count_matrix,phenoData = pd,featureData = fd,expressionFamily=tobit())


## Estimate size factors and dispersions
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

#HSMM <- detectGenes(HSMM, min_expr = 0.1)
#print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))

## distribution of mRNA totals across the cells
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

if(!is.count) # A reasonable threshold for the TPM or FPKM normalized gene expression
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                               2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                               2*sd(log10(pData(HSMM)$Total_mRNAs)))

#pdf(paste(output, ".pdf", sep= ""))
pdf(output1)
qplot(Total_mRNAs, data = pData(HSMM), color = Hours, geom = "density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)

## This section classifies the cells according to the given marker genes, e.g, MYF5 and ANPEP
if(F){
## classifying and counting cells
MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "ENSG00000111049")) ## The gene symbol is "MYF5", ENSG00000111049.3
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ENSG00000166825")) ## The gene symbol is "ANPEP", ENSG00000166825.9


cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func =
                       function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
                    { x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

HSMM <- classifyCells(HSMM, cth, 0.1)

table(pData(HSMM)$CellType)
writeLines(text= as.character(pData(HSMM)$CellType), "HSMM_cellType_annotation.txt")
}

###############
## 
if(lower_bound < upper_bound){
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                     pData(HSMM)$Total_mRNAs < upper_bound]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
}
## distribution
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
library(reshape)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
ylab.name <- ifelse(is.count, "counts", "relative counts")
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab(paste("Standardized log(", ylab.name,")", sep= "")) +
ylab("Density")

HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
pheatmap::pheatmap(cellPairwiseDistances(HSMM),main="Cell PairWise Distance")
HSMM_PT <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "Hours")

dev.off()

#save(HSMM, HSMM_PT, file= paste(output, ".RData", sep= ""))
save(HSMM, HSMM_PT, file= output2)
monocle_MST <- HSMM_PT@auxOrderingData[[HSMM_PT@dim_reduce_type]]$pr_graph_cell_proj_tree
monocle_MST_adj <- igraph::as_adjacency_matrix(monocle_MST)
output <- strsplit(output1, ".pdf")[[1]][1]
write.table(as.matrix(monocle_MST_adj), file= paste(output, "MST_adj.txt", sep= ""), row.names= pData(HSMM)$timepoint, col.names= pData(HSMM)$timepoint)

if(F){
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, "Myoblast", classify_func =
                         function(x) { x[MYF5_id,] >= 1 })
  cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
                      { x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

  HSMM <- classifyCells(HSMM, cth, 0.1)
  ###
  table(pData(HSMM)$CellType)
  ###
  pie <- ggplot(pData(HSMM),
                aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
  pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  ## Clustering cells using marker genes
  marker_diff <- markerDiffTable(HSMM[expressed_genes,],cth,residualModelFormulaStr = "~Media + num_genes_expressed",cores = 1)

  candidate_clustering_genes <-
        row.names(subset(marker_diff, qval < 0.01))
  marker_spec <-
      calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
  head(selectTopMarkers(marker_spec, 3))

  semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
  HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
  plot_ordering_genes(HSMM)

  plot_pc_variance_explained(HSMM, return_all = F)

  HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 3,
                            norm_method = 'log',
                              reduction_method = 'tSNE',
                                residualModelFormulaStr = "~Media + num_genes_expressed",
                                  verbose = T)
  HSMM <- clusterCells(HSMM, num_clusters = 2)
  plot_cell_clusters(HSMM, 1, 2, color = "CellType")

  ## Constructing Single Cell Trajectories

  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,method = 'DDRTree')
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = "Hours")
  plot_cell_trajectory(HSMM_myo, color_by = "State")


  GM_state <- function(cds){
      if (length(unique(pData(cds)$State)) > 1){
            T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
      return(as.numeric(names(T0_counts)[which
                                  (T0_counts == max(T0_counts))]))
    } else {
          return (1)
      }
  }
  HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
  plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")

  plot_cell_trajectory(HSMM_myo, color_by = "State") +
      facet_wrap(~State, nrow = 1)
}

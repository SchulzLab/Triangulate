# Triangulate
TRee guIded estimAtioN of siNgle cell reGULATion

Triangulate provides a platform for linking regulatory elements to gene expression in single cells. Given a feature matrix consisting of estimated TF activities for each gene, and a response matrix consisting of gene expression measurements at single cells, Triangulate is able to infer the TF-to-cell activities through building a mult-task learning model. The figure below illustrates the input matrices for the multi-task learning framework. The estimated coefficients obtained from training the model on the TF data can be used to interpret the activity of each TF in individual cells.
![Triangulate](https://github.com/SchulzLab/Triangulate/blob/master/images/triangulate.001.png)

We divided the procedure into three parts, which are implremented using the snakemake workflow.
Throughout this workflow three wildcards are defined: *datasets*, *imputation\_status*, and *feature\_type*. The values that we used for *datasets* in our study are StemNet for the HLC/PHH cells and HSMM for the skeletal muscle myoblast cells.
For the *imputation\_status* we used imputed and notImputed. And finally, for the *feature\_type* wildcard {static epigenetic dynamic} values were used.
## Part 1: data preparation
This step invokes an R scripts in order to generate the desired response matrix of gene expression data. The file run\_monocle\_tuorial\.R takes three arguments, 1) path to the csv file containing TPM converted expression values, 2) path to the file where the output of _Monocle_ plots are created, 3) path to the RData object where the filtered expression data obtained from applying the filtereing steps suggested by the Monocle's tutorial should be saved.
scMTL\_pipeline\_part1.sm contains the snakemake workflow for the part 1 of analysis.
## Part 2: 
In this part, first the appropriate feature matrix is created using the wildcards {static epigenetic dynamic} indicating the type of regulatory feature required for the abalysis. Afterwards, a series of filtering steps are applied on both feature and response matrices. 
removeRedundantGenes\_varianceBased.R is invoked on the reduced gene expression data from the RData object obtained in part 1. Based on the 3rd quartile of variance computed for each gene overall its expression in single cells, a cutoff is used to discard the genes that exhibit low variance in their expression profile.
removeZeroExprTFs.R removes the TFs that their corresponding gene expression is smaller than 1 for more than 90\% of cells.
scMTL\_pipeline\_part2.sm contains the snakemake workflow for the part 1 of analysis.
## Part 3:
In this last step of the workflow, we train our multi-task learning model on the feature and response matrices prepared from the previous parts. The results are stored in the path provided to the Rscript file: *run\_TGGLasso.R*. From the RData object saved by this script, one can load the data partitioned into training and test sets, via the partition variable (partitiona$test$x for feature and partitiona$test$y for response of the test partition). The coefficients of model can be accessed via TGL.model$B and TGL.model$intercept. For instance, in order to obtain the prediction on the training and test data, one can use the following command:

pred.train <- cbind(1, x.train) %\*% rbind(TGL.model$intercept, TGL.model$B)</br>

pred.test <- cbind(1, x.test) %\*% rbind(TGL.model$intercept, TGL.model$B)<br/>


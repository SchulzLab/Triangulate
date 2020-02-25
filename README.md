# Table of contents
1. [Triangulate](#Triangulate)
2. [Case study](#case-study)
# Triangulate
TRee guIded estimAtioN of siNgle cell reGULATion

Triangulate provides a platform for linking regulatory elements to gene expression in single cells. Given a feature matrix consisting of estimated TF activities for each gene, and a response matrix consisting of gene expression measurements at single cell level, Triangulate is able to infer the TF-to-cell activities through building a multi-task learning model. The figure below illustrates the input matrices for the multi-task learning framework. The estimated coefficients obtained from training the model on the TF data can be used to interpret the activity of each TF in individual cells.
![Triangulate](https://github.com/SchulzLab/Triangulate/blob/master/images/triangulate.001.png)

We divided the procedure into three parts, which are implemented using the snakemake workflow.
Throughout this workflow three wildcards are defined: *datasets*, *imputation\_status*, and *feature\_type*. The values that we used for *datasets* in our study are StemNet for the HLC/PHH cells and HSMM for the skeletal muscle myoblast cells.
For the *imputation\_status* we used imputed and notImputed. And finally, for the *feature\_type* wildcard {static epigenetic dynamic} values were used.
## Part 1: data preparation
This step invokes an R scripts in order to generate the desired response matrix of gene expression data. The file *run\_monocle\_tuorial\.R* takes three arguments, 1) path to the csv file containing TPM converted expression values, 2) path to the file where the output of _Monocle_ plots are created, 3) path to the RData object where the filtered expression data obtained from applying the filtering steps suggested by the Monocle's tutorial should be saved.

**scMTL\_pipeline\_part1.sm** in the scripts folder contains the snakemake workflow for the part 1 of analysis.
## Part 2: data filtering
In this part, first the appropriate feature matrix is created using the wildcards {static epigenetic dynamic} indicating the type of regulatory feature required for the analysis. Afterwards, a series of filtering steps are applied on both feature and response matrices. 
*removeRedundantGenes\_varianceBased.R* is invoked on the reduced gene expression data from the RData object obtained in part 1. Based on the 3rd quartile of variance computed for each gene overall its expression in single cells, a cutoff is used to discard the genes that exhibit low variance in their expression profile.
*removeZeroExprTFs.R* removes the TFs that their corresponding gene expression is smaller than 1 for more than 90\% of cells.

**scMTL\_pipeline\_part2.sm** in the scripts folder contains the snakemake workflow for the part 2 of analysis.
## Part 3: model training
In this last step of the workflow, we train our multi-task learning model on the feature and response matrices prepared from the previous parts. The results are stored in the path provided to the Rscript file: *run\_TGGLasso.R*. From the RData object saved by this script, one can load the data partitioned into training and test sets, via the partition variable (partition$test$x for feature and partition$test$y for response of the test partition). The coefficients of model can be accessed via TGL.model$B and TGL.model$intercept. For instance, in order to obtain the prediction on the training and test data, one can use the following command:
```{r}
pred.train <- cbind(1, x.train) %*% rbind(TGL.model$intercept, TGL.model$B)

pred.test <- cbind(1, x.test) %*% rbind(TGL.model$intercept, TGL.model$B)
```
**scMTL\_pipeline\_part3.sm** in the scripts folder contains the snakemake workflow for the part 3 of analysis.

# Case study <a name="case-study"></a>
To demonstrate a usage of Triangulate, we provided a snakemake file that given the processed and filtered feature and response matrices, it runs the tree-guided MTL model on the static features, for notImputed expression of the HLC/PHH cells (StemNet).
As the first step, the user needs to clone the git repository onto their desired repository:
```console
git clone https://github.com/SchulzLab/Triangulate.git
```
It is required to have the following R packages installed for the compilation of the snakemake file:
* [LinearMTL]: <github.com/tohein/LinearMTL>
* parallel
* doParallel
* monocle

To build the Triangulate model for this case study, the following command should be used in bash. But before that make sure to decompress the scMTL\_StemNet\_notImputed\_static\_feature\_doubleReduced.zip file in the data directory.
```console
snakemake -s scripts/scMTL_case_study.sm
```
After a successful execution of the command above, there should exist an Rdata file named **scMTL\_StemNet\_notImputed\_static_TGGLasso.RData** in the directory the command was executed. By loading this file into R, one can generate the coefficient heat map. But since this matrix contains 683 many TFs, the illustration of all those TFs can result in an unappealing and cluttered plot. Therefore, we prefer to show only those TFs that have higher coefficient values compared to others.
```{r}
load("scMTL_StemNet_notImputed_static_TGGLasso.RData")
library(pheatmap)
pheatmap(TGL.model$B) ## unfavorable
## Select the top TFs
top_TFs <- which(rowSums(abs(TGL.model$B)) > .5)
pheatmap(TGL.model$B[top_TFs, ])
```
![top\_TF\_coefficients](https://github.com/SchulzLab/Triangulate/blob/master/images/topTFs_coef.png)

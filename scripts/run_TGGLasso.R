library(LinearMTL)
library(parallel)
library(doParallel)
library(ggplot2)
library(monocle)
registerDoParallel(8)
set.seed(1)

max_iter <- 2000
max_iter <- 200
##############################
##############################
### Functions
generate.random.data <- function(x, y){
  rand.idx <- sample(nrow(x))
  x.random <- x[rand.idx, ]
  y.random <- y
  return(list(x= x.random, y= y.random))
}
get.randomTree <- function(y){
  cell.names <- colnames(y)
  cell.cnt <- length(cell.names)
  adj.mat <- matrix(0, ncol= cell.cnt, nrow= cell.cnt)
  random.root <- sample(cell.cnt)[2] ## for the purpose of trying another random case, I'm changing sample(cell.cnt)[1] to sample(cell.cnt)[2]
  adj.mat[random.root, ] <- 1
  adj.mat[, random.root] <- 1

  return(adj.mat)
}

##############################
get.randomTree <- function(mst){
  random.mst <- mst
  group.idx <- which(rowSums(mst) > 1 & rowSums(mst) < ncol(mst))
  for(i in seq(length(group.idx))){
    random.mst[group.idx[i], ] <- 0
    rand.cells <- sample(ncol(mst))[seq(sum(mst[group.idx[i], ]))]
    random.mst[group.idx[i], rand.cells] <- 1
  }
  return(random.mst)

}

##############################
get.baselineTree <- function(y){
  starTree <- diag(ncol(y)) # creates an identity matrix of the size of the cells
  starTree <- rbind(starTree, 1)
  return(starTree)

}

##############################
get.rmse <- function(y, yp){
  rmse <- sqrt(1/(length(y)) * sum((y - yp)^2))
  return(rmse)
}

##############################
data.partition <- function(x,y,percent){
  #partitions the whole dataset into percent% training set and (1-percent)% test set
  N <- dim(x)[1]
  train.cnt <- floor(N*percent)
  train <- list(x=x[1:train.cnt,],y=y[1:train.cnt, ])
  test <- list(x=x[(train.cnt+1):N,],y=y[(train.cnt+1):N, ])
  return(list(train=train,test=test))
}

##############################
##############################
library(glmnet)
indiv.elasticnet <- function(x, y){
  alphas <- seq(0, 1, 0.1)
  cv.elastic <- list()
  el <- list()
  pred <- list()

  for(i in seq(length(alphas))){
      cv.elastic[[i]] <- cv.glmnet(x= x, y= y, alpha= alphas[i])
      el[[i]] <- glmnet(x= x, y= y, alpha= alphas[i], lambda= cv.elastic[[i]]$lambda.min)
      pred[[i]] <- predict(el[[i]], newx= test.x.norm)
  }

  all.cvms <- sapply(seq(length(cv.elastic)), function(i) min(cv.elastic[[i]]$cvm))
  best.el.cv <- cv.elastic[[which.min(all.cvms)]]
}

##############################
train.indiv.models <- function(x,y,alpha_seq){
  K <- ncol(y)
  all.pred.res <- list()
  all.glmnet.res <- list()
  for(i in seq(K)){
    cv.glmnet_res <- mclapply(alpha_seq,function(alpha){cv.glmnet(x=as.matrix(x),y=as.numeric(y[,i]),family="gaussian",alpha=alpha,parallel=T)},mc.cores=10)
    best.error <- Inf
    for(alpha in seq(length(cv.glmnet_res))){
      cvm.best <- min(cv.glmnet_res[[alpha]]$cvm)
      if(best.error > cvm.best){
        best.model <- cv.glmnet_res[[alpha]]
        best_alpha <- alpha_seq[alpha]
        best.error <- cvm.best
      }
    }
    glmnet_res <- glmnet(x=as.matrix(x),y=as.numeric(y[,i]),family="gaussian",lambda=best.model$lambda.min,alpha=best_alpha)
    pred.res <- predict(glmnet_res,as.matrix(x),s=best.model$lambda.min,alpha=best_alpha)
    all.glmnet.res[[i]] <- glmnet_res
    all.pred.res[[i]] <- pred.res
  }
  return(list(all.pred.res=all.pred.res,all.glmnet.res=all.glmnet.res))
}

##############################
cv.folds <- function(x, y, k){
  fold.size <- ceiling(nrow(x) / k)
  fold.idx <- list()
  for(i in seq(k)){
    fold.idx[[i]] <- seq((i - 1) * fold.size + 1, i * fold.size)
    if(fold.idx[[i]][fold.size] > nrow(x))
      fold.idx[[i]] <- seq((i - 1) * fold.size + 1, nrow(x))
  }
  return(fold.idx)
}

##############################
cv.TGL.core <- function(x, y, mst, weights, lambdas, fold.idx, lambda.idx){
  k <- length(fold.idx)
  pred.acc <- vector(mode= "numeric", length= k)
  for(i in seq(k)){
    train <- list(x= x[-fold.idx[[i]], ], y= y[-fold.idx[[i]], ])
    validation <- list(x= x[fold.idx[[i]], ], y= y[fold.idx[[i]], ])

    train.mean <- list(x= colMeans(train$x), y= colMeans(train$y))
    train.sd <- list(x= apply(train$x, 2, FUN= sd), y= apply(train$y, 2, FUN= sd))

    ### If the variance is zero, the scale function results in NaN due to division by zero. Therefore, I'll just leave that variable unscaled, by setting the corresponding mean and sd to 0 and 1, respectively.
    train.sd.zero <- list(x= which(train.sd$x == 0), y= which(train.sd$y == 0))
    
    if(length(train.sd.zero$x != 0)){
      train.sd$x[train.sd.zero$x] <- 1
      train.mean$x[train.sd.zero$x] <- 0
    }
    if(length(train.sd.zero$y != 0)){
      train.sd$y[train.sd.zero$y] <- 1
      train.mean$y[train.sd.zero$y] <- 0
    }

    train.scaled <- list(x= scale(train$x, center= train.mean$x, scale= train.sd$x), y= scale(train$y, center= train.mean$y, scale= train.sd$y))
    validation.scaled <- list(x= scale(validation$x, center= train.mean$x, scale= train.sd$x), y= scale(validation$y, center= train.mean$y, scale= train.sd$y))


    model <- TreeGuidedGroupLasso(X= train.scaled$x, Y= train.scaled$y, groups= mst, weights= weights, lambda= lambdas[lambda.idx], mu= 1e-6, standardize= F, max.iter= max_iter)
    coefs <- rbind(model$intercept, model$B);
    preds <- cbind(1, validation.scaled$x) %*% coefs
    models <- model

    pred.acc[i] <- get.rmse(preds, validation.scaled$y)
  }
  return(pred.acc)
}

##############################
cv.TGL.parallel <- function(x, y, mst, weights, lambdas, fold.idx){
  library(parallel)
  library(doParallel)
  registerDoParallel(8)

  models <- list()
  preds <- list()
  k <- length(fold.idx)
  pred.acc.mat <- matrix(NA, ncol= length(lambdas), nrow= k)

  pred.acc <- mclapply(seq(length(lambdas)), function(lambda.idx){cv.TGL.core(x, y, mst, weights, lambdas, fold.idx, lambda.idx)}, mc.cores= 10)

  for(i in seq(length(pred.acc))){
    pred.acc.mat[, i] <- pred.acc[[i]]
  }

  return(pred.acc.mat)
}

##############################
cv.TGL <- function(x, y, mst, weights, lambdas, fold.idx){
  k <- length(fold.idx)
  models <- list()
  preds <- list()
  lambda.idx <- 1
  pred.acc <- matrix(NA, ncol= length(lambdas), nrow= k)
  for(lambda in lambdas){
    for(i in seq(k)){
      train <- list(x= x[-fold.idx[[i]], ], y= y[-fold.idx[[i]], ])
      validation <- list(x= x[fold.idx[[i]], ], y= y[fold.idx[[i]], ])

      train.mean <- list(x= colMeans(train$x), y= colMeans(train$y))
      train.sd <- list(x= apply(train$x, 2, FUN= sd), y= apply(train$y, 2, FUN= sd))

      ### If the variance is zero, the scale function results in NaN due to division by zero. Therefore, I'll just leave that variable unscaled, by setting the corresponding mean and sd to 0 and 1, respectively.
      train.sd.zero <- list(x= which(train.sd$x == 0), y= which(train.sd$y == 0))
      
      if(length(train.sd.zero$x != 0)){
        train.sd$x[train.sd.zero$x] <- 1
        train.mean$x[train.sd.zero$x] <- 0
      }
      if(length(train.sd.zero$y != 0)){
        train.sd$y[train.sd.zero$y] <- 1
        train.mean$y[train.sd.zero$y] <- 0
      }

      train.scaled <- list(x= scale(train$x, center= train.mean$x, scale= train.sd$x), y= scale(train$y, center= train.mean$y, scale= train.sd$y))
      validation.scaled <- list(x= scale(validation$x, center= train.mean$x, scale= train.sd$x), y= scale(validation$y, center= train.mean$y, scale= train.sd$y))

      model <- TreeGuidedGroupLasso(X= train.scaled$x, Y= train.scaled$y, groups= mst, weights= weights, lambda= lambda, mu= 1e-6, standardize= F)
      coefs <- rbind(model$intercept, model$B);
      preds <- cbind(1, validation.scaled$x) %*% coefs
      models <- model

      pred.acc[i, lambda.idx] <- get.rmse(preds, validation.scaled$y)

    }
    lambda.idx <- lambda.idx + 1
  }
  return(pred.acc)
}
##############################
##############################

args <- commandArgs(trailingOnly= T) #"random" ## monocle
tree.type <- "hc"
#experiment.type <- feature_type
input.x <- args[1]
input.y <- args[2]
output <- args[3]
#x <- read.table("scMTL_HSMM_feature.txt", header= T, row.names= 1)
#y <- read.table("scMTL_HSMM_response.txt", header= T, row.names= 1)

print("reading data...")
x <- read.csv(input.x)
y <- read.csv(input.y)

rownames(x) <- x$X
x <- x[, -1]

rownames(y) <- y$X
y <- y[, -1]

x <- log2(1 + x)
y <- log2(1 + y)

print("done reading data!")


##############################
print(c("tree.type", tree.type))

if(tree.type == "hc"){
  tree.res <- BuildTreeHC(y)
  mst <- tree.res$groups
}

##############################
##### Shuffle the data
shuffle.idx <- sample(nrow(x))

x <- x[shuffle.idx, ]
y <- y[shuffle.idx, ]

#####################
## Partition into training and test sets

percent <- .6
partition <- data.partition(x, y, percent)

########################


if(tree.type == "hc" || tree.type == "random_hc" || tree.type == "monocle_hc" || tree.type == "random_monocle_hc"){
  weights <- tree.res$weights
}else if(tree.type == "starTree"){
  weights <- vector(mode= "numeric", length= nrow(mst))
  weights <- weights + 1
}
### CV
weights.matrix <- matrix(1, ncol= ncol(mst), nrow= length(weights))

print("Starting the RunGroupCrossvalidation...")
## This function is faulty. it takes days to run and in the end it crashes!
#mtl.res.cv <- RunGroupCrossvalidation(X= log2(1 + x), Y= log2(1 + y), groups= mst, weights.matrix= weights.matrix, lambda.vec= seq(0, .3, .01), mu= 1e-6, num.threads = 20)
print("Done running RunGroupCrossvalidation!")

#lambda <- 0.3

print("Starting TreeGuidedGroupLasso for a lambda search grid")

lambdas <- seq(0, 1, .05)
#lambdas <- 1 #StemNet run's last stage was aborted but I could see that the best lambda selected was 1, so now I just wanna resume that and get the final results
fold.idx <- cv.folds(partition$train$x, partition$train$y, 5)
pred.acc <- cv.TGL.parallel(partition$train$x, partition$train$y, mst, weights, lambdas, fold.idx)
print("Done with the cv.TGL.parallel! Starting to train the final model...")

best.lambda <- lambdas[which.min(colMeans(pred.acc))]
print(c("best.lambda", best.lambda))

############
train.mean <- list(x= colMeans(partition$train$x), y= colMeans(partition$train$y))
train.sd <- list(x= apply(partition$train$x, 2, FUN= sd), y= apply(partition$train$y, 2, FUN= sd))

### If the variance is zero, the scale function results in NaN due to division by zero. Therefore, I'll just leave that variable unscaled, by setting the corresponding mean and sd to 0 and 1, respectively.
train.sd.zero <- list(x= which(train.sd$x == 0), y= which(train.sd$y == 0))

if(length(train.sd.zero$x != 0)){
  train.sd$x[train.sd.zero$x] <- 1
  train.mean$x[train.sd.zero$x] <- 0
}
if(length(train.sd.zero$y != 0)){
  train.sd$y[train.sd.zero$y] <- 1
  train.mean$y[train.sd.zero$y] <- 0
}

TGL.model <- TreeGuidedGroupLasso(X= scale(partition$train$x, center= train.mean$x, scale= train.sd$x), Y= scale(partition$train$y, center= train.mean$y, scale= train.sd$y), groups= mst, weights= weights, lambda= best.lambda, mu= 1e-6, standardize= F)

print("Done running TreeGuidedGroupLasso!")

#save(partition, TGL.model, mst, pred.acc, lambdas, best.lambda, file= paste("scMTL_", data.name, "_TGGLasso_", tree.type, "_",experiment.type ,"_weights_corrected_max_iter_200.RData", sep= ""))
save(partition, TGL.model, mst, pred.acc, lambdas, best.lambda, file= output)

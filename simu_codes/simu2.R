rm(list=ls())
dir.file <- "D:\\Working\\Research paper\\cooperateProjects\\KernelPCA\\simu"
setwd(dir.file)
source("helpfunc.R")
library(SKPCA)

#install.packages("kpcaIG")
# case 1 ---------------------------

data("iris")
str(iris)


dim.PCs <- 2
nsample <- 1000
set.seed(1)
idx_resample <- sample( nrow(iris), nsample,replace = TRUE)
iris.upsample <- iris[idx_resample, ]
n <- nrow(iris.upsample)


methodNames <- c("SKPCA-FS-L", "KPCA-Permute-L", "SPCA", "RSPCA",
                 "SKPCA-FS-G","KPCA-Permute-G", 'KPCA-IG-L', "KPCA-IG-G")
n_methods <- length(methodNames)
N <- 10
metricList <- list(Tr=matrix(NA,N, n_methods), Dist=matrix(NA, N, n_methods),
                   ARI=matrix(NA,N, n_methods), NMI=matrix(NA, N, n_methods),
                   F1score=matrix(NA,N, n_methods), ACC = matrix(NA,N, n_methods),
                   R2=matrix(NA,N, n_methods),timeMat = matrix(NA, N, n_methods))
for(i in 1:length(metricList)) colnames(metricList[[i]]) <- methodNames


for(i in 1:N){
  # i <- 1
  message('i = ', i)
  Y <- factor(iris.upsample$Species); K <- length(unique(Y))
  set.seed(i)
  X.u <- matrix(rnorm(n*10), n, 10)
  X.f <- as.matrix(iris.upsample[,1:4]) + matrix(rnorm(n*4,  sd=0.1), n, 4)
  X <- cbind(X.f, X.u)
  p <- ncol(X)
  X.s <- scale(X)
  w0 <- rep(0, p); w0[1:4] <- 1
  relist.l <- SKPCA.FS(X.s, dim.PCs=2, num.fea=4, dim.sketch=5)
  str(relist.l)
  plot(relist.l$A, col=Y)

  ## Metrics: TPR, TNR, F1-score
  set.seed(1)
  fit.g <- kmeans(relist.l$A, K)
  metricList$ARI[i,1] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,1] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  # metricList$Tr[i,1] <- GFM::measurefun(relist.l$A, H)
  # metricList$Dist[i,1] <- dis_fun(relist.l$A, H)

  metricList$timeMat[i,1] <- relist.l$time.use
  metricList$F1score[i,1] <-classific_metric(relist.l$w, w0, type='F1score')
  metricList$ACC[i,1] <- classific_metric(relist.l$w, w0, type='ACC')
  # classific_metric(relist.l$w, w0, type='Precision')
  # classific_metric(relist.l$w, w0, type='Recall')
  dat <- NULL
  metricList$R2[i,1] <-get_r2_mcfadden(relist.l$A, factor(Y))

  relist.g <- SKPCA.FS(X.s, dim.PCs=2, num.fea=4, dim.sketch=4,  kern.type = 'gaussian',
                       kernel.use = rbfdot(sigma = 0.001))
  str(relist.g)
  ## Metrics: TPR, TNR, F1-score
  set.seed(1)
  fit.g <- kmeans(relist.g$A, K)
  metricList$ARI[i,5] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,5] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  # metricList$Tr[i,5] <- GFM::measurefun(relist.g$A, H)
  # metricList$Dist[i,5] <- dis_fun(relist.g$A, H)
  metricList$timeMat[i,5] <- relist.g$time.use
  metricList$F1score[i,5] <- classific_metric(relist.g$w, w0, type='F1score')
  metricList$ACC[i,5] <- classific_metric(relist.g$w, w0, type='ACC')
  # classific_metric(relist.g$w, w0, type='Precision')
  # classific_metric(relist.g$w, w0, type='Recall')
  #
  metricList$R2[i,5] <-get_r2_mcfadden(relist.g$A, factor(Y))
  sapply(metricList, colMeans, na.rm=T)



  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4)
  hw <- rep(0, p); hw[fit.permute$fea.select] <- 1
  # GFM::measurefun(fit.permute$PCs, H)
  metricList$R2[i,2] <- get_r2_mcfadden(fit.permute$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.permute$PCs, K)
  # metricList$Tr[i,2] <- GFM::measurefun(fit.permute$PCs, H)
  # metricList$Dist[i,2] <- dis_fun(fit.permute$PCs, H)
  metricList$timeMat[i,2] <- fit.permute$time.use
  metricList$ARI[i,2] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,2] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,2] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,2] <-classific_metric(hw, w0, type='ACC')


  tic <- proc.time()
  fit.spca <- SPCA.fit(X, dim.PCs=2, num.fea=4, method=c("SPCA"),deltaK=0,
                       logalpha.start=-2, logalpha.end = -0.01, log.step=0.01) #  "RSPCA"
  toc <- proc.time()

  hw <- rep(0, p); hw[rowSums(abs(fit.spca$loadings))>0] <- 1
  metricList$R2[i,3] <-get_r2_mcfadden(fit.spca$scores, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.spca$scores, K)
  # metricList$Tr[i,3] <- GFM::measurefun(fit.spca$scores, H)
  # metricList$Dist[i,3] <- dis_fun(fit.spca$scores, H)
  metricList$timeMat[i,3] <- toc[3] - tic[3]
  metricList$ARI[i,3] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,3] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,3] <-classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,3] <-classific_metric(hw, w0, type='ACC')


  tic <- proc.time()
  fit.rspca <- SPCA.fit(X, dim.PCs=2, num.fea=4, method=c("RSPCA"),,deltaK=0,
                        logalpha.start=-2, logalpha.end = -0.01, log.step=0.01) #  "RSPCA"
  toc <- proc.time()
  metricList$R2[i,4] <-get_r2_mcfadden(fit.rspca$scores, factor(Y))
  hw <- rep(0, p); hw[rowSums(abs(fit.rspca$loadings))>0] <- 1
  set.seed(1)
  fit.g <- kmeans(fit.rspca$scores, K)
  # metricList$Tr[i,4] <- GFM::measurefun(fit.rspca$scores, H)
  # metricList$Dist[i,4] <- dis_fun(fit.rspca$scores, H)
  metricList$timeMat[i,4] <- toc[3] - tic[3]
  metricList$ARI[i,4] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,4] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,4] <-classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,4] <-classific_metric(hw, w0, type='ACC')


  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4, kern.type ="gaussian.radial.basis")
  hw <- rep(0, p); hw[fit.permute$fea.select] <- 1
  # GFM::measurefun(fit.permute$PCs, H)
  metricList$R2[i,6] <- get_r2_mcfadden(fit.permute$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.permute$PCs, K)
  metricList$timeMat[i,6] <- fit.permute$time.use
  metricList$ARI[i,6] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,6] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,6] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,6] <-classific_metric(hw, w0, type='ACC')

  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = 'polydot')
  hw <- rep(0, p); hw[res.IG$fea.select] <- 1
  # GFM::measurefun(fit.permute$PCs, H)
  metricList$R2[i,7] <- get_r2_mcfadden(res.IG$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(res.IG$PCs,  K)
  metricList$timeMat[i,7] <- res.IG$time.use
  metricList$ARI[i,7] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,7] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,7] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,7] <-classific_metric(hw, w0, type='ACC')

  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = "rbfdot")
  hw <- rep(0, p); hw[res.IG$fea.select] <- 1
  # GFM::measurefun(fit.permute$PCs, H)
  metricList$R2[i,8] <- get_r2_mcfadden(res.IG$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(res.IG$PCs,  K)
  metricList$timeMat[i,8] <- res.IG$time.use
  metricList$ARI[i,8] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,8] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,8] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,8] <-classific_metric(hw, w0, type='ACC')

  sapply(metricList, colMeans, na.rm=T)

}
sapply(metricList, colMeans, na.rm=T)
save(metricList, file='simu2_dim.PCs2_nfea4_2025_1_23.rds')




# Rerun KPCA-Permute and KPCA-IG using post-selected features -------------

pos_kpca.permute.fit <- function(X_select, dim.PCs,kern.type= c("linear", "gaussian.radial.basis")){
  library(mixKernel)
  tic <- proc.time()
  if(kern.type == 'gaussian.radial.basis'){
    phychem.kernel.use <- compute.kernel(X_select, kernel.func = kern.type, sigma=1e-4)
  }else{
    phychem.kernel.use <- compute.kernel(X_select, kernel.func = kern.type)
  }
  
  # perform a KPCA
  kernel.use.pca.result <- kernel.pca(phychem.kernel.use, ncomp = dim.PCs)
  toc <- proc.time()
  return(list(PCs=kernel.use.pca.result$x, time.use=toc[3]-tic[3]))
  
}

post_kpcaIG.fit <- function(X_select, dim.PCs, kern.type = c('polydot', "rbfdot")){
  library(kpcaIG)
  kern.type <- match.arg(kern.type)
  tic <- proc.time()
  # compute one kernel.use for the psychem dataset
  if(kern.type == 'rbfdot'){
    kpar <- list(sigma=1e-4)
  }else{
    kpar <- list()
  }
  kpca_tan <- kernelpca(X_select,  kernel = kern.type, features = dim.PCs, kpar = kpar)
  toc <- proc.time()
  return(list(PCs=kpca_tan@pcv, time.use=toc[3]-tic[3]))
} 

for(i in 1:N){
  # i <- 2
  message('i = ', i)
  Y <- factor(iris.upsample$Species); K <- length(unique(Y))
  set.seed(i)
  X.u <- matrix(rnorm(n*10), n, 10)
  X.f <- as.matrix(iris.upsample[,1:4]) + matrix(rnorm(n*4,  sd=0.1), n, 4)
  X <- cbind(X.f, X.u)
  p <- ncol(X)
  X.s <- scale(X)
  w0 <- rep(0, p); w0[1:4] <- 1
  
  m.use <- round(sqrt(nrow(X)))
  
  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4)
  hw <- rep(0, p); hw[fit.permute$fea.select] <- 1
  fea.index <- fit.permute$fea.select
  # fit.permute.post.l <- pos_kpca.permute.fit((X[,fea.index]), dim.PCs, kern.type = c('linear'))
  fit.permute.post.l <- SKPCA(X[,fea.index], d=dim.PCs, kern= polydot(), m=m.use,  Stype = 'subGaussian', use.R=TRUE)
  fit.permute.post.l$PCs <- fit.permute.post.l$A
  
  dat <- NULL
  # GFM::measurefun(fit.permute$PCs, H)
  metricList$R2[i,2] <- get_r2_mcfadden(fit.permute.post.l$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.permute.post.l$PCs, K)
  # fit.g <- kmeans(X[,fea.index], K)
  metricList$timeMat[i,2] <- fit.permute$time.use
  metricList$ARI[i,2] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,2] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,2] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,2] <-classific_metric(hw, w0, type='ACC')
  
  
  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4, kern.type ="gaussian.radial.basis")
  hw <- rep(0, p); hw[fit.permute$fea.select] <- 1
  # fit.permute.post.g <- pos_kpca.permute.fit(X[,fit.permute$fea.select], dim.PCs, kern.type = c('gaussian.radial.basis'))
  fit.permute.post.g <- SKPCA(X[,fit.permute$fea.select], d=dim.PCs, kern= rbfdot(sigma = 0.001), m=m.use,  Stype = 'subGaussian', use.R=TRUE)
  fit.permute.post.g$PCs <- fit.permute.post.g$A
  metricList$R2[i,6] <- get_r2_mcfadden(fit.permute.post.g$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.permute.post.g$PCs, K)
  metricList$timeMat[i,6] <- fit.permute$time.use
  metricList$ARI[i,6] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,6] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,6] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,6] <-classific_metric(hw, w0, type='ACC')
  
  
  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = 'polydot')
  hw <- rep(0, p); hw[res.IG$fea.select] <- 1
  fea.index <- res.IG$fea.select
  #fit.ig.post.l <- post_kpcaIG.fit(X[,fea.index], dim.PCs, kern.type = c('polydot'))
  fit.ig.post.l <- SKPCA(X[,fea.index], d=dim.PCs, kern= polydot(), m=m.use,  Stype = 'subGaussian', use.R=T)
  
  fit_prin <- princomp(X[,fea.index])
  fit.ig.post.l$PCs <- fit.ig.post.l$A
  
  
  metricList$R2[i,7] <- get_r2_mcfadden(fit.ig.post.l$PCs, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.ig.post.l$PCs,  K)
  metricList$timeMat[i,7] <- res.IG$time.use
  metricList$ARI[i,7] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,7] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,7] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,7] <-classific_metric(hw, w0, type='ACC')
  
  
  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = "rbfdot")
  hw <- rep(0, p); hw[res.IG$fea.select] <- 1
  fea.index <- res.IG$fea.select
  # fit.ig.post.g <- post_kpcaIG.fit(X[,fea.index], dim.PCs, kern.type = c('rbfdot'))
  fit.ig.post.g <- SKPCA(X[,fea.index], d=dim.PCs, kern= rbfdot(sigma = 0.001), m=m.use,  Stype = 'subGaussian', use.R=TRUE)
  fit.ig.post.g$PCS <- fit.ig.post.g$A
  metricList$R2[i,8] <- get_r2_mcfadden(fit.ig.post.g$PCS, factor(Y))
  set.seed(1)
  fit.g <- kmeans(fit.ig.post.g$PCS,  K)
  metricList$timeMat[i,8] <- res.IG$time.use
  metricList$ARI[i,8] <- cluster_metric(fit.g$cluster, Y, type='ARI')
  metricList$NMI[i,8] <- cluster_metric(fit.g$cluster, Y, type='NMI')
  metricList$F1score[i,8] <- classific_metric(hw, w0, type='F1score')
  metricList$ACC[i,8] <-classific_metric(hw, w0, type='ACC')
  
  sapply(metricList, colMeans, na.rm=T)
  
}
sapply(metricList, colMeans, na.rm=T)
save(metricList, file='simu2_dim.PCs2_nfea4_2025_02_16.rds')


# Compare the computation time by varying nsample -------------------------
methodNames <- c("SKPCA-FS-L", "SKPCA-FS-G", "KPCA-Permute-L", "SPCA", "RSPCA",
                 "KPCA-Permute-G", 'KPCA-IG-L', "KPCA-IG-G")
nsample.vec <- c(200, 400, 800, 1600, 3200, 6400)
n_methods <- length(methodNames)
N <- 10 # 10
timeList <- list(time.array=array(dim = c(length(nsample.vec), n_methods, N)))
row.names(timeList[[1]]) <- nsample.vec; colnames(timeList[[1]]) <- methodNames

for(ni in 5:length(nsample.vec)){
  # ni <- 1
  message("ni = ", ni, "/", length(nsample.vec))
  nsample <- nsample.vec[ni]
  idx_resample <- sample( nrow(iris), nsample,replace = TRUE)
  iris.upsample <- iris[idx_resample, ]
  n <- nrow(iris.upsample)
  
  Y <- factor(iris.upsample$Species); K <- length(unique(Y))
  set.seed(ni)
  X.u <- matrix(rnorm(n*10), n, 10)
  X.f <- as.matrix(iris.upsample[,1:4]) + matrix(rnorm(n*4,  sd=0.1), n, 4)
  X <- cbind(X.f, X.u)
  
  for(i in 1:N){
    # i <- 1
    message('i = ', i)
    try({
      relist.l <- SKPCA.FS(X, dim.PCs=2, num.fea=4, dim.sketch = 40 , use.R=TRUE)
      timeList$time.array[ni,1, i] <- relist.l$time.use
    },silent=T)
    
    try({
      relist.g <- SKPCA.FS(X, dim.PCs=2, num.fea=4,  kern.type = 'gaussian',
                           kernel.use = rbfdot(sigma = 0.001), use.R=TRUE)
      timeList$time.array[ni,2, i] <- relist.g$time.use
    },silent=T)
    
    
    
    fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4)
    timeList$time.array[ni,3, i] <- fit.permute$time.use



    fit.spca <- SPCA.fit(X, dim.PCs=2, num.fea=4, method=c("SPCA"),deltaK=0,
                         logalpha.start=-2, logalpha.end = -0.01, log.step=0.01) #  "RSPCA"
    timeList$time.array[ni,4, i] <- fit.spca$time.use


    fit.rspca <- SPCA.fit(X, dim.PCs=2, num.fea=4, method=c("RSPCA"),,deltaK=0,
                          logalpha.start=-2, logalpha.end = -0.01, log.step=0.01) #  "RSPCA"

    timeList$time.array[ni,5, i] <- fit.rspca$time.use


    fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4, kern.type ="gaussian.radial.basis")
    timeList$time.array[ni,6, i] <- fit.permute$time.use


    res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = 'polydot')
    timeList$time.array[ni,7, i] <- res.IG$time.use


    res.IG.g <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = "rbfdot")
    timeList$time.array[ni,8, i] <- res.IG.g$time.use
    
    
    apply(timeList$time.array, c(1,2), mean, na.rm=T)
    
  }
  save(timeList, file = 'timeList_simu2_iris_varyn_FS.rds')
}
apply(timeList$time.array, c(1,2), mean, na.rm=T)
save(timeList, file = 'timeList_simu2_iris_varyn_FS.rds')
# Importance score --------------------------------------------------------

data("iris")
str(iris)


dim.PCs <- 2
nsample <- 1000
set.seed(1)
idx_resample <- sample( nrow(iris), nsample,replace = TRUE)
iris.upsample <- iris[idx_resample, ]
n <- nrow(iris.upsample)

methodNames <- c("SKPCA.FS-L", "KPCA-Permute-L",
                  'KPCA-IG-L',"SKPCA.FS-G", "KPCA-Permute-G", "KPCA-IG-G")
q.vec <- c(8, 18, 24, 50, 100)
q.use <- q.vec[3]
N <- 10
p <- 4+ q.use
metricList <- list(imp.score=array(dim=c(p, length(methodNames),N)))
for(ii in 1:length(metricList)) colnames(metricList[[ii]]) <- methodNames
scale2one <- function(x) x/sum(x)

for(i in 1:N){
  #i <- 1
  message('i = ', i)
  Y <- factor(iris.upsample$Species); K <- length(unique(Y))
  set.seed(i)
  X.u <- matrix(rnorm(n*q.use), n, q.use)
  X.f <- as.matrix(iris.upsample[,1:4]) + matrix(rnorm(n*4,  sd=0.1), n, 4)
  X <- cbind(X.f, X.u)
  p <- ncol(X)
  w0 <- rep(0, p); w0[1:4] <- 1
  relist.l <- SKPCA.FS(X, dim.PCs=2, num.fea=4, use.R=T)
  str(relist.l)
  barplot(relist.l$importance.scores)
  metricList$imp.score[,1, i] <- relist.l$importance.scores



  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4)
  barplot(fit.permute$cc.dist)
  metricList$imp.score[,2, i] <- scale2one(fit.permute$cc.dist)


  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = 'polydot')
  str(res.IG)
  barplot(res.IG$means_norms)
  metricList$imp.score[,3, i] <- scale2one(res.IG$means_norms)


  ### Gaussian kernel
  relist.g <- SKPCA.FS(X, dim.PCs=2, num.fea=4,  kern.type = 'gaussian',
                        kernel = rbfdot(sigma = 0.001), use.R=T)
  str(relist.g)
  barplot(relist.g$importance.scores)
  metricList$imp.score[,4, i] <- relist.g$importance.scores

  fit.permute <- KPCA.Permute.fit(X, dim.PCs=2, num.fea = 4, kern.type ="gaussian.radial.basis")
  str(fit.permute)
  barplot(fit.permute$cc.dist)
  metricList$imp.score[,5, i] <- scale2one(fit.permute$cc.dist)

  res.IG <- kpcaIG.fit(X, dim.PCs = 2, num.fea = 4, kern.type = "rbfdot")
  str(res.IG)
  barplot(res.IG$means_norms, ylab='Score')
  metricList$imp.score[,6, i] <- scale2one(res.IG$means_norms)


}
save(metricList, file=paste0("Impscore_q",q.use, '_metricList_iris.rds'))


imp.mat <- apply(metricList$imp.score, c(1,2), mean, na.rm=T)
head(imp.mat)

par(mfrow=c(1,3))
barplot(imp.mat[,1], main='SKPCA.FS')
barplot(imp.mat[,2], main='KPCA-Permute')
barplot(imp.mat[,3], main='KPCA-IG')

par(mfrow=c(1,3))
barplot(imp.mat[,4], main='SKPCA.FS')
barplot(imp.mat[,5], main='KPCA-Permute')
barplot(imp.mat[,6], main='KPCA-IG')




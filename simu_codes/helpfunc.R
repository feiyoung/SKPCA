gendata <- function(n, p = 1000, q=20, seed=1, tau=10){
  set.seed(1)

  Z <- matrix(rnorm(p*q), p, q)
  V <- qr.Q(qr(Z))
  t(V)%*%V
  lambdavec <- 1 - ((1:q)-1)/ p
  Lambda <- diag(lambdavec)
  set.seed(seed)
  U <- matrix(rnorm(n*q), n, q)
  tE <- matrix(rnorm(n*p), n, p) / tau
  X <- U%*%Lambda%*% t(V) + tE
  return(list(X = X, U=U))
}



# Metrics -----------------------------------------------------------------
#-- assist to evaluate reconstruction error
recerr2 <- function(tsdata,trdata,alpha,kern){
  # x <- tsdata
  # data <- trdata
  ns <- nrow(tsdata)
  d <- ncol(alpha)
  if(!inherits(tsdata,'matrix')) tsdata <- as.matrix(tsdata)
  if(!inherits(trdata,'matrix')) trdata <- as.matrix(trdata)
  K <- kernelMatrix(kern, as.matrix(trdata))
  Krow <- apply(K,2, mean)
  Ksum <- mean(Krow)
  sumalpha <- apply(alpha, 2, sum)
  alphaKrow <- Krow %*% alpha
  kmat = kernelMatrix(kern, trdata,tsdata)
  f <- t(kmat) %*% (alpha) - matrix(apply(kmat, 2,mean)- Ksum, ncol=1)%*% t(sumalpha) -
    matrix(alphaKrow, nrow=ns, ncol=d, byrow=TRUE)
  f <- as.matrix(f)
  res <- diag(kernelMatrix(kern, tsdata,tsdata)) - 2*apply(kmat, 2, mean) + Ksum - apply(f*Conj(f),1, sum)
  return(as.vector(res))
}


mat2orth <- function(mat){
  # mat <- R0
  qr1 <- qr(mat)
  Q <- qr.Q(qr1)
  return(Q * sqrt(nrow(mat)))
}
colSD <- function(xmat, na.rm=TRUE){

  apply(xmat, 2, sd, na.rm=na.rm)
}
###
dis_fun <- function(H, H0){
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% ginv(t(H) %*% H) %*% t(H) %*% H0


  val <- 1- tr_fun(mat1)/tr_fun(t(H0) %*% H0)
  return(val)
}
normvec_fun <- function(mu, hmu){
  mean(abs(mu-hmu))
}
project.err <- function(H, H0){
  require(MASS)
  dH <-  H0- H %*% ginv(t(H) %*% H) %*% t(H) %*% H0
  sqrt(sum(dH^2))/ sqrt(sum(H0^2))
}
#### Clustering metric
cluster_metric <- function(hy, y, type=c('ARI', 'NMI') ){

  require(mclust)
  require(aricode)
  type <- match.arg(type)
  switch(type,
         ARI= adjustedRandIndex(hy, y),
         NMI = NMI(as.vector(hy), y))
}
classific_metric <- function(hy, y, type=c('F1score', 'ACC', 'Precision','Recall' )){
  library(Metrics)
  library(MLmetrics)
  type <- match.arg(type)
  ## There are many metrics in this package
  switch(type,
         F1score= F1_Score(as.vector(hy), y),
         ACC = Accuracy(as.vector(hy), y),
         Precision=Precision(as.vector(hy), y),
         Recall = Recall(as.vector(hy), y))
}


get_r2_mcfadden <- function(embeds, y){
  library(nnet)
  library(performance)
  dat <- NULL
  y <- as.numeric(as.factor(y))
  hq <- ncol(embeds)
  dat <- as.data.frame(cbind(y=y, x=embeds))
  dat$y <- factor(dat$y)
  name <-  c('y', paste0('V', 1:hq))
  names(dat) <-name
  formu <- paste0("y~")
  for(i in 1:hq){
    if(i < hq){
      formu <- paste(formu, name[i+1], seq='+')
    }else{
      formu <- paste(formu, name[i+1], seq='')
    }

  }
  model1 <- nnet::multinom(as.formula(formu), data = dat)
  R2 <- r2_mcfadden(model1)
  return(R2$R2_adjusted)
}
# My methods --------------------------------------------------------------

select.bw <- function(X, scale=TRUE, nsample=100, seed=1){

  set.seed(seed)
  if(scale) X <- scale(X)
  n_spots <- nrow(X)
  idx <- sample(n_spots, min(nsample, n_spots))
  dis <-stats::dist(X[idx,])
  sigma <- 1/sqrt(mean(dis^2))
  return(sigma)
}

Tr <- function(x) sum(diag(x))
Diag <- function (vec)
{
  q <- length(vec)
  if (q > 1) {
    y <- diag(vec)
  }
  else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}

linear.prog <- function(y, s= 50){
  library(Rglpk)
  p <- length(y)
  A <- matrix(1, nrow = 1, ncol=p, byrow = TRUE)

  bounds <- list(lower=list(ind=1:p, val=rep(0, p)),
                 upper = list(ind=1:p, val=rep(1, p)))
  result <- Rglpk_solve_LP(obj = y, mat = A, dir = c("<="), rhs = s,
                           bounds = bounds, types = "C",max=TRUE, maxit = 10000, tol.abs = 1e-09, tol.rel = 1e-09)

  optimal_solution <- result$solution
  return(optimal_solution)
}
sketchfun <- function(m, n, type=c("subGaussian", "ROS", "subSampling")){
  library(Matrix)
  Praw <- Diagonal(n, x = 1)
  ind <- sample(n, m)
  R <- Diagonal(n, 1/2 - rbinom(n, 1, 0.5))
  H <- Praw
  S <- switch(type,
              subGaussian = matrix(rnorm(n*m),m,n),
              ROS = sqrt(n/m)*Praw[ind,]%*% H %*% R,
              subSampling = sqrt(n/m)*Praw[ind,])
  return(S)
}
# oKPCA <- function(X,d, kern=rbfdot(sigma=0.5)){
#
#
#   n <- nrow(X)
#   tic <- proc.time()
#   K <- kernelMatrix(kern, as.matrix(X))
#   eigK <- eigen(K, symmetric = TRUE)
#   toc <- proc.time()
#   res <- list()
#   res$alpha <- eigK$vectors[,1:d]
#   res$K <- K
#   res$time.use <- toc[3] - tic[3]
#   return(res)
# }
#
SKPCA.R <- function(X, d, m=ceiling(sqrt(nrow(X))), Stype=c('subGaussian', 'ROS', 'subSampling'),
                  kern=rbfdot(sigma=0.5), rescale=TRUE, seed=2){

  library(irlba)
  Stype <- match.arg(Stype)
  tic <- proc.time()
  n <- nrow(X)
  if(rescale) X <- scale(X) ## scale X
  K <- kernelMatrix(kern, as.matrix(X))
  # Krow <- apply(K,2, mean)
  # Ksum <- mean(Krow)
  #K <- K - matrix(Krow, n,n, byrow=T) - matrix(Krow, n, n) + Ksum
  set.seed(seed) # reproducible
  S <- sketchfun(m,n, Stype)
  SKS <- S %*% K %*% t(S)
  eigSKS <- eigen(SKS, symmetric = TRUE)
  V <- K %*% t(S) %*% (eigSKS$vectors) %*% Diag(1/sqrt(eigSKS$values))
  svdV <- irlba(V, nv=d)
  toc <- proc.time()

  res <- list()
  res$A <- svdV$u
  res$K <- V %*% t(V)
  res$time.use <- toc[3] - tic[3]
  class(res) <- 'sKPCA'
  return(res)
}


# sKPCA.Cpp <- function(X, d, m=ceiling(sqrt(nrow(X))), Stype=c('subGaussian', 'ROS', 'subSampling'),
#            kern=rbfdot(sigma=0.5), rescale=TRUE,seed=2){
#
#   # d=4; m=50;Stype='subGaussian';
#   # kern=rbfdot(sigma=0.5); seed=2;
#   Stype <- match.arg(Stype)
#   if(rescale) X <- scale(X) ## scale X
#   n <- nrow(X)
#   K <- matrix(0, nrow=n, ncol=n); A <- matrix(0, nrow=n, ncol=d)
#   tic <- proc.time()
#   set.seed(seed)
#   res.cpp <- SKPCA_cpp(X, kern, d=d, m=m, Stype=Stype)
#   toc <- proc.time()
#
#   # res.cpp <- list(A=A, K= K)
#   res.cpp$time.use <- toc[3] - tic[3]
#   return(res.cpp)
# }
#
# gendata_LFM <- function(seed=1, n=400, p=200, p_nz=50, q=6, sigma2_e=0.01){
#   set.seed(seed)
#   library(dplyr)
#   H <- rnorm(n*q) %>% matrix(nrow=n, ncol=q) %>% qr %>% qr.Q()
#   Bg <- rnorm(p_nz*q) %>% matrix(nrow=p_nz, ncol=q) %>% qr %>% qr.Q()
#   Bg <- Bg %*% diag(c(10, 8, 7, 5, 3, 1))
#   B <- rbind(Bg, matrix(0, p-p_nz, q))
#   sigma2_e <- 0.01
#   U <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = sigma2_e*diag(rep(1,p)))
#   X <- H %*% t(B) + U
#   nz.index <- 1:p_nz
#   return(list(X, H0=H, nz.index=nz.index))
# }
# SKPCA.FS.R <- function(X, dim.PCs=6, num.fea=50, dim.sketch=round(sqrt(nrow(X))),
#                      kern.type = c("linear", "gaussian"), gaussian.method=c("M1", "M0"),
#                      kernel=NULL, Stype=c('subGaussian', 'ROS', 'subSampling'),
#                      use.sketch=TRUE,  maxIter=20, epsObj = 1e-8){
#   # kern.type = c( "gaussian"); use.sketch=TRUE;  maxIter=20; epsObj = 1e-8
#   #  num.fea=4; dim.PCs=2; m <- 50; kern <- rbfdot(sigma=0.1)
#   require(kernlab)
#   require(irlba)
#   Stype <- match.arg(Stype)
#   X <- scale(X)
#   tic <- proc.time()
#   kern.type <- match.arg(kern.type)
#   gaussian.method <- match.arg(gaussian.method)
#   objfun <- function(X, alpha1, wvec){
#     n <- nrow(X)
#     #alpha1_col <- matrix(alpha1, n, 1)
#     y <- t(alpha1) %*% X %*% t(X * matrix(wvec, n, p, byrow=T)) %*% alpha1
#     return(Tr(y))
#   }
#
#   p <- ncol(X); n <- nrow(X);
#   s <- num.fea; q <- dim.PCs; m <- dim.sketch
#   obj.vec <- rep(NA, maxIter)
#   obj.vec[1] <- -1e10
#
#   if(q >s) stop("dim.PCs must be less than num.fea!")
#   if(q<2) stop("dim.PCs must be greater than 1!")
#   ### Initialize
#   wvec <- rep(1, p)
#   alpha1 <- irlba(X, nv=q)$u[,1:q] ## can be speed up using ilrba
#
#   if(is.null(kernel)){
#     if(kern.type=='linear'){
#       kern <- polydot(degree = 1)
#     }else if(kern.type=='gaussian'){
#       kern <- rbfdot(sigma=0.1)
#     }
#
#   }
#   if(!is.null(kernel)){
#     kern <- kernel
#   }
#
#   if(kern.type == 'gaussian' && gaussian.method == 'M0'){
#     Kcube <- array(dim=c(n, n, p))
#     for(l in 1:p){
#       Kcube[,,l] <-  kernelMatrix(kern,X[,l])
#     }
#
#   }
#   if(kern.type == 'gaussian' && gaussian.method == 'M1'){
#     res.skpca <- sKPCA.Cpp(X=X * matrix(wvec, n, p, byrow=TRUE), d=q, m=m, kern = kern)
#     K <- res.skpca$K
#   }
#
#
#   # message("check point1????")
#   for(iter in 2:maxIter){
#
#     # iter <- 2
#
#     # update wvec
#     if(kern.type=='linear'){
#       y <- rowSums((t(X) %*% alpha1)^2)
#     }else if(kern.type=='gaussian'){
#       y <- rep(NA, p)
#       if(gaussian.method=='M1'){
#         for(l in 1:p){
#           tmpMat <- t(alpha1) %*% (-2*K * (matrix(X[,l], n, n)- matrix(X[,l], n, n, byrow=T))^2*kern@kpar$sigma) %*% alpha1
#           y[l] <- Tr(tmpMat)   # * wvec[l]
#         }
#
#       }else if(gaussian.method == 'M0'){
#         for(l in 1:p){
#           tmpMat <- t(alpha1) %*% (Kcube[,,l]- 0.1) %*% alpha1
#           y[l] <- Tr(tmpMat)   # * wvec[l]
#         }
#       }
#
#
#     }
#     cat(y, '\n')
#     # message("check point2????")
#     #wvec <- linear.prog(y, s=s)
#     wvec <- rep(0,p)
#     idx <- order(y, decreasing = TRUE)
#     wvec[idx[1:s]] <- 1
#     # wy <- update_w( X, alpha1, as.matrix(K), s, kern.type, Kg_h=1/kern@kpar$sigma)
#     # wvec <- wy[,1]
#     # y <- wy[,2]
#     cat(wvec, '\n')
#     cat(y, '\n')
#     ## update alpha: svd
#     # obj1 <- objfun(X, alpha1, wvec)
#     ### speed using sktech
#     if(use.sketch){
#       res.skpca <- sKPCA.Cpp(X=X[,idx[1:s]], d=q, m=m, kern = kern, Stype= Stype)
#       alpha1 <- res.skpca$A
#       K <- res.skpca$K
#     }else{
#       K <- kernelMatrix(kern, X[,idx[1:s]])
#       alpha1 <- eigen(K)$vectors[,1:q]
#     }
#     # obj2 <-  objfun(X, alpha1, wvec)
#     # message("dalpha1 = ", obj2-obj1)
#     # message("check point3????")
#
#     obj3 <-  objfun(X, alpha1, wvec)
#     # message("dpvec = ", obj3-obj2)
#     obj.vec[iter] <- obj3
#     dObj <- abs(obj.vec[iter]-obj.vec[iter-1])/abs(obj.vec[iter-1])
#     message("Obj=", round(obj.vec[iter], 4) , ", dObj = ",  dObj)
#     if(dObj < epsObj){
#       break
#     }
#   }
#   toc <- proc.time()
#   time.use <- toc[3] - tic[3]
#   if(kern.type=='linear'){
#     importance.scores <- y/sum(y)
#   }else if(kern.type=='gaussian'){
#     if(gaussian.method=='M1'){
#       y.scale <- y/max(abs(y))
#       importance.scores <- exp(y.scale/2)/sum(exp(y.scale/2))
#
#     }else if(gaussian.method == 'M0'){
#       for(l in 1:p){
#         tmpMat <- t(alpha1) %*% (Kcube[,,l]- 0.1) %*% alpha1
#         y[l] <- Tr(tmpMat)   # * wvec[l]
#       }
#     }
#   }
#
#   reslist <- list(A = alpha1, w = as.integer(wvec>0), importance.scores=importance.scores, obj.seq = obj.vec[1:iter][-1],
#                   time.use = time.use)
#   return(reslist)
# }
#
#
# SKPCA.FS <- function(X, dim.PCs=6, num.fea=50, dim.sketch=round(sqrt(nrow(X))),
#                      kern.type = c("linear", "gaussian"), gaussian.method=c("M1", "M0"),
#                      kernel.use=NULL, Stype=c('subGaussian', 'ROS', 'subSampling'),
#                      rescale = TRUE,
#                      use.sketch=TRUE,  maxIter=20, epsObj = 1e-8){
#   # kern.type = c( "gaussian"); use.sketch=TRUE;  maxIter=20; epsObj = 1e-8
#   #  num.fea=4; dim.PCs=2; m <- 50; kern <- rbfdot(sigma=0.1)
#   # kern.type <- 'linear'
#   require(kernlab)
#   require(irlba)
#   Stype <- match.arg(Stype)
#   tic <- proc.time()
#   kern.type <- match.arg(kern.type)
#   gaussian.method <- match.arg(gaussian.method)
#
#   if(rescale) X <- scale(X)
#   p <- ncol(X); n <- nrow(X);
#   s <- num.fea; q <- dim.PCs; m <- dim.sketch
#   obj.vec <- rep(NA, maxIter)
#   obj.vec[1] <- -1e10
#
#   if(q >s) stop("dim.PCs must be less than num.fea!")
#   if(q<2) stop("dim.PCs must be greater than 1!")
#   ### Initialize
#   wvec <- rep(1, p)
#   alpha1 <- irlba(X, nv=q)$u[,1:q]
#
#   if(!is.null(kernel.use)){
#     kern <- kernel.use
#   }
#   if(is.null(kernel.use)){
#     if(kern.type=='linear'){
#       kern <- polydot(degree = 1)
#     }else if(kern.type=='gaussian'){
#       kern <- rbfdot(sigma=0.1)
#     }
#   }
#   if(kern.type=='linear'){
#      h <- 1
#   }else if(kern.type=='gaussian'){
#     h <- 1/kern@kpar$sigma
#   }
#   # if(use.sketch){
#   #   res.skpca <- sKPCA.Cpp(X=X * matrix(wvec, n, p, byrow=TRUE), d=q, m=m, kern = kern)
#   #   K <- res.skpca$K
#   # }else{
#   #   K <- kernelMatrix(kern, X)
#   # }
#   K <- kernelMatrix(kern, X)
#   if(kern.type == 'gaussian' && gaussian.method == 'M0'){
#     Kcube <- array(dim=c(n, n, p))
#     for(l in 1:p){
#       Kcube[,,l] <-  kernel.useMatrix(kern,X[,l])
#     }
#
#   }
#
#
#
#   # message("check point1????")
#   for(iter in 2:maxIter){
#
#     # iter <- 2
#     wy <- update_w( X, alpha1, as.matrix(K), s, kern.type, Kg_h=h)
#     wvec <- wy[,1]
#     y <- wy[,2]
#     idx <- order(y, decreasing = TRUE)
#     # cat(wvec, '\n')
#     # cat(y, '\n')
#     ## update alpha: svd
#     # obj1 <- objfun(X, alpha1, wvec)
#     ### speed using sktech
#     if(use.sketch){
#       res.skpca <- sKPCA.Cpp(X=X[,idx[1:s]], d=q, m=m, kern = kern, Stype= Stype)
#       alpha1 <- res.skpca$A
#       K <- res.skpca$K
#     }else{
#       K <- kernel.useMatrix(kern, X[,idx[1:s]])
#       alpha1 <- eigen(K)$vectors[,1:q]
#     }
#     # obj2 <-  objfun(X, alpha1, wvec)
#     # message("dalpha1 = ", obj2-obj1)
#     # message("check point3????")
#
#     # obj3 <-  objfun(X, alpha1, wvec)
#     obj3 <-  objfun_Cpp(X[,idx[1:s]], alpha1)
#     # message("dpvec = ", obj3-obj2)
#     obj.vec[iter] <- obj3
#     dObj <- abs(obj.vec[iter]-obj.vec[iter-1])/abs(obj.vec[iter-1])
#     message("Obj=", round(obj.vec[iter], 4) , ", dObj = ",  dObj)
#     if(dObj < epsObj){
#       break
#     }
#   }
#   toc <- proc.time()
#   time.use <- toc[3] - tic[3]
#   if(kern.type=='linear'){
#     importance.scores <- y/sum(y)
#   }else if(kern.type=='gaussian'){
#     if(gaussian.method=='M1'){
#       y.scale <- y/max(abs(y))
#       importance.scores <- exp(y.scale/2)/sum(exp(y.scale/2))
#
#     }else if(gaussian.method == 'M0'){
#       for(l in 1:p){
#         tmpMat <- t(alpha1) %*% (Kcube[,,l]- 0.1) %*% alpha1
#         y[l] <- Tr(tmpMat)   # * wvec[l]
#       }
#     }
#   }
#
#   reslist <- list(A = alpha1, w = as.integer(wvec>0), importance.scores=importance.scores,
#                   time.use = time.use)
#   return(reslist)
# }
#




# Compared methods --------------------------------------------------------

kpcaIG.fit <- function(X,dim.PCs=6, num.fea = 50, kern.type=c('polydot', "rbfdot")){

  # dim.PCs=6; num.fea = 50; kern.type=c('polydot')
  library(kpcaIG)
  kern.type <- match.arg(kern.type)
  tic <- proc.time()
  # compute one kernel.use for the psychem dataset
  p <- ncol(X)
  colnames(X) <- 1:p
  kpca_tan <- kernelpca(X,  kernel = kern.type, features = dim.PCs)
  kpca_ig_tan <- kpca_igrad(kpca_tan, dim = 1:dim.PCs)
  row.names(kpca_ig_tan) <- kpca_ig_tan$column_names
  kpca_ig_tan_raw <- kpca_ig_tan[as.character(1:p),]
  imp.var.permute <- as.numeric(kpca_ig_tan$column_names)[1:num.fea]
  toc <- proc.time()
  time.use <- toc[3] - tic[3]
  return(list(PCs=kpca_tan@pcv, fea.select=imp.var.permute, means_norms =kpca_ig_tan_raw$means_norms ,time.use=time.use, fit=kpca_ig_tan_raw))



}

KPCA.Permute.fit <- function(X, dim.PCs=6, num.fea = 50, kern.type= c("linear", "gaussian.radial.basis")){
  library(mixKernel)

  kern.type <- match.arg(kern.type)
  tic <- proc.time()
  # compute one kernel.use for the psychem dataset
  p <- ncol(X)
  colnames(X) <- 1:p


  phychem.kernel.use <- compute.kernel(X, kernel.func = kern.type)
  # perform a KPCA
  kernel.use.pca.result <- kernel.pca(phychem.kernel.use, ncomp = dim.PCs)

  # compute importance for all variables in this kernel.use
  kernel.use.pca.result <- kernel.pca.permute(kernel.use.pca.result, phychem=1:p)
  imp.var.permute <- order(kernel.use.pca.result$cc.distances, decreasing = T)[1:num.fea]
  toc <- proc.time()
  time.use <- toc[3] - tic[3]
  return(list(PCs=kernel.use.pca.result$x, fea.select=imp.var.permute, cc.dist =kernel.use.pca.result$cc.distances ,time.use=time.use))

}

SPCA.fit <- function(X, dim.PCs=6, num.fea=50, method=c("SPCA", "RSPCA"), deltaK=2,
                     logalpha = -4, logalpha.start=-5, logalpha.end=-1, log.step=0.05){


  library(sparsepca)
  method <- match.arg(method)
  tic <- proc.time()
  q <- dim.PCs
  message("search the optimal penalty parameters...")
  logalpha.set <- seq(logalpha.start, logalpha.end, by=log.step)
  for(i in seq_along(logalpha.set)){
    # i <- 1
    logalpha <- logalpha.set[i]
    if(method=='SPCA'){
      out.spca <- spca(X, k=q, alpha=10^logalpha, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)

    }else if(method=='RSPCA'){
      out.spca <- rspca(X, k=q, alpha=10^logalpha, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
    }
     num.ftmp <- sum(abs(rowSums(out.spca$loadings)) >0)
    message('i = ', i, ", logalpha=", logalpha, ", num: ", num.ftmp)
    if(num.ftmp <= num.fea+deltaK){
      logalpha.opt <- logalpha
      break
    }
  }
  if(i == length(logalpha.set)){
    logalpha.opt <- logalpha
  }

  if(method=='SPCA'){
    out.spca <- spca(X, k=q, alpha=10^logalpha.opt, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)

  }else if(method=='RSPCA'){
    out.spca <- rspca(X, k=q, alpha=10^logalpha.opt, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
  }
  toc <- proc.time()
  out.spca$time.use <- toc[3] - tic[3]
  return(out.spca)


}







# Generate data -----------------------------------------------------------
gendata_LFM <- function(seed=1, n=400, p=200, p_nz=50, q=6, sigma2_e=0.01){
  set.seed(seed)
  library(dplyr)
  H <- rnorm(n*q) %>% matrix(nrow=n, ncol=q) %>% qr %>% qr.Q()
  Bg <- rnorm(p_nz*q) %>% matrix(nrow=p_nz, ncol=q) %>% qr %>% qr.Q()
  Bg <- Bg %*% diag(c(10, 8, 7, 5, 3, 1))
  B <- rbind(Bg, matrix(0, p-p_nz, q))
  sigma2_e <- 0.01
  U <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = sigma2_e*diag(rep(1,p)))
  X <- H %*% t(B) + U
  nz.index <- 1:p_nz
  return(list(X, H0=H, nz.index=nz.index))
}

# Methods -----------------------------------------------------------------
SKPCA.R <- function(X, d, m=ceiling(sqrt(nrow(X))), Stype=c('subGaussian', 'ROS', 'subSampling'),
                    kern=rbfdot(sigma=0.5), rescale=TRUE, seed=1){

  sketchfun <- function(m, n, type=c("subGaussian", "ROS", "subSampling")){
    if(n < 2) stop("SKPCA.R: n must be gteater than 1!")
    ind <- sample(n, m)
    rvec <-  1/2 - rbinom(n, 1, 0.5)
    S <- switch(type,
                subGaussian = matrix(rnorm(n*m),m,n),
                ROS = sqrt(n/m)*(matrix(rnorm(n*n)*rvec,n,n, byrow=TRUE))[ind,] /sqrt(n),
                subSampling = sqrt(n/m)* diag(n)[ind,])
    return(S)
  }
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
  SKS <- S %*% K %*% base::t(S)

  eigSKS <- eigen(SKS, symmetric = TRUE)
  eigSKS$values[eigSKS$values<1e-12] <- 1e-12 ### avoid the outliears
  # if(Stype != 'subGaussian'){
  #   V <- K %*% Matrix::t(S) %*% (eigSKS$vectors) %*% Diag(1/sqrt(eigSKS$values))
  # }else{
  V <- K %*% base::t(S) %*% (eigSKS$vectors) %*% Diag(1/sqrt(eigSKS$values))
  # }

  svdV <- irlba(V, nv=d)
  toc <- proc.time()

  res <- list()
  res$A <- svdV$u
  res$K <- V %*% base::t(as.matrix(V))
  res$eigvalues <- svdV$d^2
  res$time.use <- toc[3] - tic[3]
  class(res) <- 'SKPCA'
  return(res)
}


### original kernel PCA
oKPCA <- function(X,d, kern=rbfdot(sigma=0.1), rescale=TRUE){
  require(kernlab)
  if(!inherits(X, "matrix")) stop("X must be a matrix object!")
  if(rescale) X <- scale(X)
  n <- nrow(X)
  tic <- proc.time()
  K <- kernelMatrix(kern, as.matrix(X))
  eigK <- eigen(K, symmetric = TRUE)
  toc <- proc.time()
  res <- list()
  res$alpha <- eigK$vectors[,1:d]
  res$eigvalues <- eigK$values
  res$K <- K
  res$time.use <- toc[3] - tic[3]
  return(res)
}

SKPCA <- function(X, d, m=ceiling(sqrt(nrow(X))), Stype=c('subGaussian', 'ROS', 'subSampling'),
                      kern=rbfdot(sigma=1e-4), use.R = FALSE, rescale=TRUE,seed=1){

  # d=4; m=50;Stype='subGaussian';
  # kern=rbfdot(sigma=0.5); seed=2;
  require(kernlab)
  if(!inherits(X, "matrix")) stop("X must be a matrix object!")
  Stype <- match.arg(Stype)
  if(use.R){
    message("Implement SKPCA using R code...")
   res.R <- SKPCA.R(X, d=d, m=m, Stype=Stype,
                        kern=kern, rescale=rescale, seed=seed)
   return(res.R)
  }
  message("Implement SKPCA using C++ code...")
  if(rescale) X <- scale(X) ## scale X
  n <- nrow(X)
  K <- matrix(0, nrow=n, ncol=n); A <- matrix(0, nrow=n, ncol=d)
  tic <- proc.time()
  set.seed(seed)
  res.cpp <- SKPCA_cpp(X, kern, d=d, m=m, Stype=Stype)
  toc <- proc.time()
  res.cpp$eigvalues <- res.cpp$svals^2
  # res.cpp <- list(A=A, K= K)
  res.cpp$time.use <- toc[3] - tic[3]
  res.cpp$svals <- NULL
  return(res.cpp)
}

SKPCA.FS <- function(X, dim.PCs=6, num.fea=50, dim.sketch=round(sqrt(nrow(X))),
                     kern.type = c("linear", "gaussian"),
                     kernel.use=NULL, Stype=c('subGaussian', 'ROS', 'subSampling'),
                     rescale = TRUE,
                     use.sketch=TRUE, use.R = FALSE,  maxIter=20, epsObj = 1e-8){
  # kern.type = c( "gaussian"); use.sketch=TRUE;  maxIter=20; epsObj = 1e-8
  #  num.fea=4; dim.PCs=2; m <- 50; kern <- rbfdot(sigma=0.1)
  # kern.type <- 'linear'

  if(!inherits(X, "matrix")) stop("X must be a matrix object!")
  require(kernlab)
  require(irlba)
  Stype <- match.arg(Stype)
  tic <- proc.time()
  kern.type <- match.arg(kern.type)

  if(rescale) X <- scale(X)
  p <- ncol(X); n <- nrow(X);
  s <- num.fea; q <- dim.PCs; m <- dim.sketch
  obj.vec <- rep(NA, maxIter)
  obj.vec[1] <- -1e10

  if(q >s) stop("dim.PCs must be less than num.fea!")
  if(q<2) stop("dim.PCs must be greater than 1!")
  ### Initialize
  wvec <- rep(1, p)
  alpha1 <- irlba(X, nv=q)$u[,1:q]

  if(!is.null(kernel.use)){
    kern <- kernel.use
  }
  if(is.null(kernel.use)){
    if(kern.type=='linear'){
      kern <- polydot(degree = 1)
    }else if(kern.type=='gaussian'){
      kern <- rbfdot(sigma=0.1)
    }
  }
  if(kern.type=='linear'){
    h <- 1
  }else if(kern.type=='gaussian'){
    h <- 1/kern@kpar$sigma
  }

  K <- kernelMatrix(kern, X)


  # message("check point1????")
  for(iter in 2:maxIter){

    # iter <- 2
    wy <- update_w( X, alpha1, as.matrix(K), s, kern.type, Kg_h=h)
    wvec <- wy[,1]
    y <- wy[,2]
    idx <- order(y, decreasing = TRUE)
    ### speed using sktech
    if(use.sketch){
      res.skpca <- SKPCA(X=X[,idx[1:s]], d=q, m=m, kern = kern, Stype= Stype, use.R=use.R)
      alpha1 <- res.skpca$A
      K <- res.skpca$K
    }else{
      K <- kernel.useMatrix(kern, X[,idx[1:s]])
      alpha1 <- eigen(K)$vectors[,1:q]
    }
    obj3 <-  objfun_Cpp(X[,idx[1:s]], alpha1)
    obj.vec[iter] <- obj3
    dObj <- abs(obj.vec[iter]-obj.vec[iter-1])/abs(obj.vec[iter-1])
    message("Obj=", round(obj.vec[iter], 4) , ", dObj = ",  dObj)
    if(dObj < epsObj){
      break
    }
  }
  toc <- proc.time()
  time.use <- toc[3] - tic[3]
  if(kern.type=='linear'){
    importance.scores <- y/sum(y)
  }else if(kern.type=='gaussian'){
      y.scale <- y/max(abs(y))
      importance.scores <- exp(y.scale/2)/sum(exp(y.scale/2))

  }

  reslist <- list(A = alpha1, w = as.integer(wvec>0), importance.scores=importance.scores,
                  time.use = time.use)
  return(reslist)
}



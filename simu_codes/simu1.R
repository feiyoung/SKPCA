### Scenario 1: SKPCA vs KPCA
dir.file <- "D:\\Working\\Research paper\\cooperateProjects\\KernelPCA\\simu"
setwd(dir.file)
source("helpfunc.R")
library(MASS)
library(SKPCA)
library(kernlab)
library(GFM)

n <- 1000; p <- 100; q <- 20

datlist <- gendata(n=n, p = p, q=q, seed=1)
str(datlist)



# 1. Fix sketch type and dimension of PCs, varing m -----------------------
hsigma <- 1e-4
message('hsigma=', hsigma)
# kern <- rbfdot(sigma = hsigma)
kern <- polydot()
Stype.vec <- c('subGaussian', 'ROS', 'subSampling')
Stype <- Stype.vec[1]

dim.PCs <- 20
N <- 1
m.vec <- round(seq(dim.PCs+1, round(sqrt(n))+30, len=10))
time.kpca.vec <- rep(NA, N)
n_methods <- length(m.vec) + 1
metricList <- list(Err.Tr=matrix(NA,N, n_methods), Err.F=matrix(NA, N, n_methods),
                   Dist.eigs=matrix(NA,N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(i in 1:length(metricList)) colnames(metricList[[i]]) <- c(m.vec, "oKPCA")

for(i in 1:N){
  # i <- 1
  message('i = ', i)
  datlist <- gendata(n=n, p = p, seed=i)
  str(datlist)


  relist.opca <- oKPCA(datlist$X, d=dim.PCs,  kern=kern)
  metricList$Err.Tr[i,n_methods] <- dis_fun( relist.opca$alpha, datlist$U)
  metricList$Err.F[i,n_methods] <- project.err(relist.opca$alpha,  datlist$U)
  metricList$timeMat[i,n_methods] <- relist.opca$time.use
  for(mi in 1:length(m.vec)){
    relist.skpca.sg <- SKPCA(datlist$X, d=dim.PCs, kern=kern, m=m.vec[mi], Stype = Stype, seed=10)
    metricList$Err.Tr[i,mi] <- dis_fun( relist.skpca.sg$A, datlist$U)
    metricList$Err.F[i,mi] <- project.err(relist.skpca.sg$A,  datlist$U)
    metricList$Dist.eigs[i,mi] <- mean(abs(relist.opca$eigvalues[1:dim.PCs]- relist.skpca.sg$eigvalues[1:dim.PCs]))
    metricList$timeMat[i,mi] <- relist.skpca.sg$time.use
  }
}
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD, na.rm=T)

# 2. Compare the running time by varying n------------------------------------------------
hsigma <- 1e-4
message('hsigma=', hsigma)
kern <- rbfdot(sigma = hsigma)
# kern <- polydot()
n.vec <- c(500, 1000, 3000,  6000, 8000, 10000, 15000)
#n.vec <- c(3000, 6000)


N <- 10
m.vec <- round(seq(dim.PCs+1, round(sqrt(n))+30, len=10))
time.kpca.vec <- rep(NA, N)
methods.use <- c("KPCA", "SKPCA-subG", "SKPCA-ROS", "SKPCA-subS")
n_methods <- length(methods.use)
metricList <- list(timeMat = array(NA, dim=c(length(n.vec), n_methods, N)))
for(i in 1:length(metricList)) colnames(metricList[[i]]) <-methods.use

for(i.n in 1:length(n.vec)){
  #i.n <- 1
  message("i.n = ", i.n, "/", length(n.vec))
  m.use <- round(sqrt(n.vec[i.n]))
  for(i in 1:N){
    # i <- 1
    message('i = ', i)
    datlist <- gendata(n=n.vec[i.n], p = p, seed=i)
    str(datlist)
    relist.opca <- oKPCA(datlist$X, d=dim.PCs,  kern=kern)
    metricList$timeMat[i.n,1, i] <- relist.opca$time.use
    relist.skpca.sg <- SKPCA(datlist$X, d=dim.PCs, kern=kern, m=m.use, Stype = 'subGaussian', use.R=T)
    metricList$timeMat[i.n,2, i] <- relist.skpca.sg$time.use
    relist.skpca.ros <- SKPCA(datlist$X, d=dim.PCs, kern=kern, m=m.use, Stype = 'ROS', use.R=T)
    metricList$timeMat[i.n,3, i] <- relist.skpca.ros$time.use
    relist.skpca.ss <- SKPCA(datlist$X, d=dim.PCs, kern=kern, m=m.use, Stype = 'subSampling', use.R=T)
    metricList$timeMat[i.n,4, i] <- relist.skpca.ss$time.use
  }
}
apply(metricList$timeMat, c(1,2), mean)
apply(metricList$timeMat, c(1,2), sd)
save(n.vec, metricList, file=paste0("varyn_Gaussian_kernel_comptime","dimPC",dim.PCs,'S12_case1V2.rds'))


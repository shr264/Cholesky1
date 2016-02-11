library(Matrix)
library(MASS)
library(matrixcalc)
library(clusterGeneration)
library(mvtnorm)
library(graph)
library(igraph)
library(ggm)
library(xtable)

source('threshS.r')
source('graphfiller.r')
source('cholcalc.r')
source('concentrationgraph.r')
source('covargraph2.r')
source('covargraph3.R')
source('choladjest.r')
source('concentrationgraphChol.r')
source('medgeset.r')
source('orderchol.r')
source('covargraphicf.r')

negloglik <- function(S,Sigma){
    (out = matrix.trace(as.matrix(solve(Sigma)%*%S)) + log(det(as.matrix(Sigma))))
    return(out)
}

kldiv <- function(S,Sigma,p){
    (out = (1/2)*(matrix.trace(as.matrix(solve(Sigma)%*%S)) - p + log(det(as.matrix(Sigma))/det(as.matrix(S)))))
    return(out)}

steinsloss<- function(S,Sigma,p){
    (out = matrix.trace(as.matrix(S%*%solve(Sigma))) - log(det(as.matrix(S%*%solve(Sigma)))) - p)
    return(out)}

getrndmsigma <- function(dim,n,thresh){
pdmat = genPositiveDefMat(dim)$Sigma
pdmat = threshS(pdmat,thresh)
while(min(eigen(pdmat)$value)<=0){
    pdmat = pdmat -(min(eigen(pdmat)$value)-0.001)*diag(dim)}
adjcencymatrix = pdmat>10^(-10)
Sigma = pdmat
sum(adjcencymatrix)/(dim^2)
mu = rep(0,dim)
n = ceiling(n)
Y = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
covar = covargraph2(Y,adjctMat=adjcencymatrix)
Sigma = covar$omegahat
adjcencymatrix = covar$adjMat
return(list(Sigma=Sigma,adjcencymatrix=adjcencymatrix))
}

gettdata <- function(p,n,nu,thresh){
Sigma = getrndmsigma(p,p+2,thresh)
mu = rep(0,p)
Y = rmvt(n,sigma=Sigma$Sigma,df=nu,delta=mu)
(S = t(Y)%*%Y/n)
(tvar = (nu/(nu-2))*Sigma$Sigma)
return(list(data = Y,amat=Sigma$adjcencymatrix,S = S,tvar = tvar))
}


icfamatS <- function(adjcencymatrix,S,dim){
(adjct = adjcencymatrix - diag(dim))
dimnames(adjct)[[1]] = as.list(seq(1:dim))
dimnames(adjct)[[2]] = as.list(seq(1:dim))
dimnames(S)[[1]] = as.list(seq(1:dim))
dimnames(S)[[2]] = as.list(seq(1:dim))
return(list(amat = adjct,S = S))}

frobnorm <- function(test5Shat,test2omegahat){
frobsigma = norm(((test5Shat)-(test2omegahat)), type = "F")^2/norm((test2omegahat), type = "F")^2
return(frobsigma)
}

replicationtdist <- function(p,n,maxreplicate,mu,thresh){
(table = matrix(0, nrow = maxreplicate, ncol = 6))
(dimnames(table)[[2]] = c('Frob-Algo1','Frob-Icf','Stein-Algo1','Stein-Icf','Time-Algo1','Time-Icf'))
for (replicate in 1:maxreplicate){
    cat("Replication:",replicate)
Y = gettdata(p,n,mu,thresh)

Yicf = icfamatS(Y$amat,Y$S,p)

time1 = proc.time()
icffit = fitCovGraph(Yicf$amat, Yicf$S,n ,alg = "icf", start.icf = NULL, tol = 1e-06)
time2 = proc.time()
icftime = (time2-time1)[3]
icffrob = frobnorm(icffit$Shat,Y$tvar)
icfstein = steinsloss(icffit$Shat,Y$tvar,p)

algo1fit = covargraph2(Y=Y$data,adjctMat=Y$amat)
algo1frob = frobnorm(algo1fit$omegahat,Y$tvar)
algo1stein = steinsloss(algo1fit$omegahat,Y$tvar,p)
algo1time = algo1fit$time[3]

(table[replicate,1] =algo1frob)
(table[replicate,2] =icffrob)
(table[replicate,3] =algo1stein)
(table[replicate,4] =icfstein)
(table[replicate,5] =algo1time)
(table[replicate,6] =icftime)
mats = list(icffit$Shat,algo1fit$omegahat)
}
return(c(table,mats))
}

replicationtdistalgo1 <- function(p,n,maxreplicate,mu,thresh){
(table = matrix(0, nrow = maxreplicate, ncol = 3))
(dimnames(table)[[2]] = c('Frob-Algo1','Stein-Algo1','Time-Algo1'))
for (replicate in 1:maxreplicate){
    cat("Replication:",replicate)
Y = gettdata(p,n,mu,thresh)

Yicf = icfamatS(Y$amat,Y$S,p)

#time1 = proc.time()
#icffit = fitCovGraph(Yicf$amat, Yicf$S,n ,alg = "icf", start.icf = NULL, tol = 1e-06)
#time2 = proc.time()
#icftime = (time2-time1)[3]
#icffrob = frobnorm(icffit$Shat,Y$tvar)
#icfstein = steinsloss(icffit$Shat,Y$tvar,p)

algo1fit = covargraph2(Y=Y$data,adjctMat=Y$amat)
algo1frob = frobnorm(algo1fit$omegahat,Y$tvar)
algo1stein = steinsloss(algo1fit$omegahat,Y$tvar,p)
algo1time = algo1fit$time[3]

(table[replicate,1] =algo1frob)
(table[replicate,2] =algo1stein)
(table[replicate,3] =algo1time)
mats = algo1fit$omegahat
}
return(list(table=table,mats=mats))
}

replicationtdisticf <- function(p,n,maxreplicate,mu,thresh){
(table = matrix(0, nrow = maxreplicate, ncol = 3))
(dimnames(table)[[2]] = c('Frob-Icf','Stein-Icf','Time-Icf'))
for (replicate in 1:maxreplicate){
    cat("Replication:",replicate)
Y = gettdata(p,n,mu,thresh)

Yicf = icfamatS(Y$amat,Y$S,p)

time1 = proc.time()
icffit = fitCovGraph(Yicf$amat, Yicf$S,n ,alg = "icf", start.icf = NULL, tol = 1e-06)
time2 = proc.time()
icftime = (time2-time1)[3]
icffrob = frobnorm(icffit$Shat,Y$tvar)
icfstein = steinsloss(icffit$Shat,Y$tvar,p)

#algo1fit = covargraph2(Y=Y$data,adjctMat=Y$amat)
#algo1frob = frobnorm(algo1fit$omegahat,Y$tvar)
#algo1stein = steinsloss(algo1fit$omegahat,Y$tvar,p)
#algo1time = algo1fit$time[3]

(table[replicate,1] =icffrob)
(table[replicate,2] =icfstein)
(table[replicate,3] =icftime)
mats = icffit$Shat
}
return(list(table=table,mats=mats))
}

set.seed(12345AA)
table2 = replicationtdistalgo1(100,200,1,3,0.025)
write.table((data.frame((table2$table))),file="table-algo1tdistp100n200AA.txt",row.names = TRUE)
write.table((data.frame((table2$mats))),file="table-algo1mattdistp100n200AA.txt",row.names = TRUE)

set.seed(12345AA)
table3 = replicationtdisticf(100,200,1,3,0.025)
write.table((data.frame((table3$table))),file="table-icftdistp100n200AA.txt",row.names = TRUE)
write.table((data.frame((table3$mats))),file="table-icfmattdistp100n200AA.txt",row.names = TRUE)


#download.packages(c("Matrix","MASS","matrixcalc","clusterGeneration","mvtnorm","ggm"),destdir="/Users/syedrahman/Documents/Spring2015/Research/Cholesky/cholesky03242015")
#download.packages(c("qpgraph"),destdir="/Users/syedrahman/Documents/Spring2015/Research/Cholesky/cholesky03242015")
#download.packages(c("qtl"),destdir="/Users/syedrahman/Documents/Spring2015/Research/Cholesky/cholesky03242015")
#setwd()
#library("gtools")
#getDependencies("ggm")




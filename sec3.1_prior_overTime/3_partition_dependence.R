library(mclust) ## adjustedRandIndex
library(mcclust) ## for variation of information
library(salso) ## for psm - pairwise similarity matrix
###############################
## Source useful functions
###############################
## function that simulate random partitions from the informed partition model
source("utils/randomPartitionGenerator.R")

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition
c0 <- rep(1:4, each=5)
alpha <- matrix(rep(c(0.95,0.75,0.5,0.25), each=5), nrow=length(c0), ncol=ntime, byrow=FALSE)
nMC <- 1000

markov <- t(replicate(n=nMC, rtpartition(N=length(c0), M=1, alpha=alpha, ntime=ntime, FirstPart=c0, model="markovian")))
independent <- t(replicate(n=nMC, rtpartition(N=length(c0), M=1, alpha=alpha, ntime=ntime, FirstPart=c0, model="independent")))

## latent autoregressive 
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)

latent_autoPos <- t(replicate(n=nMC, rtpartition_AR(N=length(c0), M=1, ntime=ntime, FirstPart=c0, 
  zeta0 = zeta0, phi = 0.8, kappa = 0.5, clusterSpecific  = FALSE)))
latent_autoNSPos <- t(replicate(n=nMC, rtpartition_AR(N=length(c0), M=1, ntime=ntime, FirstPart=c0, 
  zeta0 = zeta0, phi = 1.5, kappa = 0.5, clusterSpecific  = FALSE)))
latent_autoNeg <- t(replicate(n=nMC, rtpartition_AR(N=length(c0), M=1,  ntime=ntime, FirstPart=c0, 
  zeta0 = zeta0, phi = -1.5, kappa = 0.5, clusterSpecific  = FALSE)))

ari_lag <- function(ll){
  ci <- ll$ciMat
  ari_lag1 <- numeric()
  ari_lag2 <- matrix(NA, nrow=nrow(ci)-1, ncol=nrow(ci)-1)
  for(i in 2:nrow(ci)){
  	ari_lag1[i-1] <- adjustedRandIndex(ci[1,], ci[i,])
  	for(ii in 2:nrow(ci)){
  	  ari_lag2[i-1,ii-1] <- adjustedRandIndex(ci[i,], ci[ii,])
  	}
  }
  list(ari_lag_c0 = ari_lag1, ari_lag = ari_lag2)
}

ari_mar <- apply(markov,1,ari_lag)
ari_ind <- apply(independent,1,ari_lag)
ari_autoPos <- apply(latent_autoPos,1,ari_lag)
ari_autoNSPos <- apply(latent_autoNSPos,1,ari_lag)
ari_autoNeg <- apply(latent_autoNeg,1,ari_lag)

library(fields)

tmp1 <- apply(sapply(ari_mar,function(x) x[[1]]),1,mean)
tmp2 <- apply(sapply(ari_ind,function(x) x[[1]]),1,mean)
tmpautoPos <- apply(sapply(ari_autoPos,function(x) x[[1]]),1,mean)
tmp4 <- apply(sapply(ari_autoNeg,function(x) x[[1]]),1,mean)
tmpautoNSPos <- apply(sapply(ari_autoNSPos,function(x) x[[1]]),1,mean)

pdf("figures/ARI_NS.pdf", height=5, width=18)
  par(mfrow=c(1,4))
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_mar,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Markovian", xlab="Lag", ylab="Lag" ,cex.lab=2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_ind,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Independent", xlab="Lag", ylab="Lag" ,cex.lab=2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_autoNSPos,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Latent AR, NS, pos", xlab="Lag", ylab="Lag" ,cex.lab=2)
             
  plot(1:9, tmp1, type='b', ylab="", xlab="Lag", cex=1.5, ,cex.lab=2, ylim = c(0, 0.4))
  points(1:9, tmp2, type='b', pch=2, cex=1.5)
  points(1:9, tmpautoNSPos, type='b', pch=4, cex=1.5)
  legend(x=5, y=0.3, legend=c("Markovian","Independent", "AR, NS, pos"), 
    pch=c(1,2, 4, 5), cex=1.5)
  mtext(side=2, "Adjusted Rand Index", line = 2)
dev.off()



pdf("figures/ARI_comparison.pdf", height=5, width=22)
  par(mfrow=c(1,5))
  plot(1:9, tmp1, type='b', ylab="", xlab="Lag", cex=1.5, ,cex.lab=2, ylim = c(0, 0.4))
  points(1:9, tmp2, type='b', pch=2, cex=1.5)
  points(1:9, tmpautoPos, type='b', pch=4, cex=1.5)
    points(1:9, tmpautoNSPos, type='b', pch=5, cex=1.5)
  legend(x=5, y=0.32, legend=c("Markovian","Independent", "AR, S", "AR, NS"), 
    pch=c(1,2, 4, 5), cex=1.5)
  mtext(side=2, "Adjusted Rand Index", line = 2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_mar,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Markovian", xlab="Lag", ylab="Lag" , cex.main =2.5, cex.lab=2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_ind,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Independent", xlab="Lag", ylab="Lag" ,cex.main =2.5 ,cex.lab=2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_autoPos,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Latent AR, S", xlab="Lag", ylab="Lag" ,cex.main =2.5 ,cex.lab=2)
  image.plot(1:9, 1:9, 
             matrix(apply(sapply(ari_autoNSPos,function(x) x[[2]]),1,mean), nrow=ntime-1, ncol=ntime-1, byrow=TRUE),
             zlim=c(0,1), main = "Latent AR, NS", xlab="Lag", ylab="Lag" ,cex.main =2.5 ,cex.lab=2)
             
dev.off()


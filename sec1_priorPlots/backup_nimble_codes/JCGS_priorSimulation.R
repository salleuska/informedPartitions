library(mclust) ## adjustedRandIndex
source("RandomPartitionGenerator_sal.R")

# Sim2 
# n of units to partition
N <- 20;
## DP concentration parameter
M <- 0.5;
## time
ntime <- 10;
## Centering partition
# FirstPart= seq(1:20)
# FirstPart= rep(c(1,2,3,4), each = 5)

# cSimulatePartition(N, M, rho, ntime, FirstPart)
nsim <- 10000
out.ours005 <- list()
out.ours025 <- list()
out.ours050 <- list()
out.ours075 <- list()
out.ours095 <- list()
out.ours100 <- list()
# Generate partitions to see differences between our approach and that of Caron
omat7 <- omat8 <- omat9 <- omat10 <- omat11 <- omat12 <- matrix(NA, ntime, ntime)
set.seed(1)
for(ii in 1:nsim){
	if(ii %% 1000 == 0) cat("ii = ", ii, "\n")
	ours005 <- rtpartition(N=N,M=M,rho=0.05,ntime=ntime, FirstPart=NULL)
	ours025 <- rtpartition(N=N,M=M,rho=0.25,ntime=ntime, FirstPart=NULL)
	ours050 <- rtpartition(N=N,M=M,rho=0.50,ntime=ntime, FirstPart=NULL)
	ours075 <- rtpartition(N=N,M=M,rho=0.75,ntime=ntime, FirstPart=NULL)
	ours095 <- rtpartition(N=N,M=M,rho=0.95,ntime=ntime, FirstPart=NULL)
	ours100 <- rtpartition(N=N,M=M,rho=1,ntime=ntime, FirstPart=NULL)

	for(j in 1:ntime){
		for(jj in 1:ntime){

			omat7[j,jj] <- adjustedRandIndex(ours005$ciMat[j, ], ours005$ciMat[jj, ])			
			omat8[j,jj] <- adjustedRandIndex(ours025$ciMat[j, ], ours025$ciMat[jj, ])			
			omat9[j,jj] <- adjustedRandIndex(ours050$ciMat[j, ], ours050$ciMat[jj, ])			
			omat10[j,jj] <- adjustedRandIndex(ours075$ciMat[j, ], ours075$ciMat[jj, ])			
			omat11[j,jj] <- adjustedRandIndex(ours095$ciMat[j, ], ours095$ciMat[jj, ])			
			omat12[j,jj] <- adjustedRandIndex(ours100$ciMat[j, ], ours100$ciMat[jj, ])			
		}
	}
	out.ours005[[ii]] <- omat7
	out.ours025[[ii]] <- omat8
	out.ours050[[ii]] <- omat9
	out.ours075[[ii]] <- omat10
	out.ours095[[ii]] <- omat11
	out.ours100[[ii]] <- omat12
}


ours005.mn <- Reduce('+', out.ours005)/nsim
ours025.mn <- Reduce('+', out.ours025)/nsim
ours050.mn <- Reduce('+', out.ours050)/nsim
ours075.mn <- Reduce('+', out.ours075)/nsim
ours095.mn <- Reduce('+', out.ours095)/nsim
ours100.mn <- Reduce('+', out.ours100)/nsim

zlim <- range(c(c(ours005.mn),c(ours025.mn),c(ours050.mn),c(ours075.mn),c(ours095.mn),c(ours100.mn)))


out <- list(out.ours005, 
						out.ours025,
						out.ours050,
						out.ours075,
						out.ours095,
						out.ours100)
saveRDS(out, file = "sampledDistances.rds")

# out <- readRDS("sampledDistances.rds")
# nsim <- length(out[[2]])
# ours005.mn <- Reduce('+', out[[1]])/nsim
# ours025.mn <- Reduce('+', out[[2]])/nsim
# ours050.mn <- Reduce('+', out[[3]])/nsim
# ours075.mn <- Reduce('+', out[[4]])/nsim
# ours095.mn <- Reduce('+', out[[5]])/nsim
# ours100.mn <- Reduce('+', out[[6]])/nsim

# zlim <- range(c(c(ours005.mn),c(ours025.mn),c(ours050.mn),c(ours075.mn),c(ours095.mn),c(ours100.mn)))


library(MixSim)
library(fields)
ntime <- 10;
# Reproduce Figure 2 in JCGS paper
# pdf("~/Research/BYU/SpaceTimePPM/latex/plots/ours3.pdf", height=7, width=10)
par(mfrow=c(2,3))
image.plot(ours005.mn, zlim=zlim,  axes=FALSE, main=expression(alpha==0.05))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

image.plot(ours025.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.25))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

image.plot(ours050.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.50))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

image.plot(ours075.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.75))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

image.plot(ours095.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.95))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

image.plot(ours100.mn, zlim=zlim, axes=FALSE, main=expression(alpha==1))
mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
# dev.off()
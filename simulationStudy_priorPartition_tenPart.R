# Aug 17 2023
# N - number of experimental units
# M - DP scale (dispersion) parameter
# ntime - number of time points
# alpha - probability of reallocation 
#         This needs to be a N x ntime matrix
#         If it is constant, then it simply is 
#         matrix of the form alpha*J (J is matrix of 1s)
# tau - sd associated with atom generation
# sig - sd associated with data generation
# c0 - Allows me to provide the initial partition if so desired
# N <- 20;M <- 1;ntime <- 10;c0=NULL; alpha <- matrix(0.5, nrow=N, ncol=ntime);model="markovian"
rtpartition <- function(N, M, alpha, ntime, c0=NULL, model="markovian"){
  library(plyr)
  out <- NULL

  ci <- c0

  Tm <- ntime
  if(is.null(c0)){
    # Create cluster at time period t = 1
    ci <- 1 # cluster label
    mh <- 1 # vector of cluster sizes
    K <- 1 # number of clusters

	for(k in 2:N){
	  p <- c(mh/(k-1+M), M/(k-1+M))
	  ci <- c(ci, sample(1:(K+1), 1, prob=p))
	  mh <- table(ci)
	  K <- length(unique(ci))
	}
  } else {
    mh <- table(ci)
	K <- length(unique(ci))
	N <- length(ci)
  }


# Now to generate 10-1=9 more partitions that are temporally dependent
# Will first try uniform deletion
  ciMat <- matrix(NA, nrow=Tm+1, ncol=N)
  ciMat[1,]  <- ci
  K <-  length(unique(ci))
  gMat <- matrix(0,nrow=Tm+1, ncol=N)

# print(ciMat)
  if(Tm > 1){
    for(t in 2:(Tm+1)){
#	  cat("t = ", t, "\n")
	  r <- rbinom(N, 1, alpha[,t]) #recall if gamma = 1 then not reallocated
	  dnk <- c(1:N)[r==0]
	  gMat[t,] <- r
	  if(sum(r) > 0){
	    ci[dnk]<- 0
	  }

#	  print(ci)
	  mh <- tabulate(ci[ci!=0]);
	  if((K - length(mh)) > 0){
	    mh <- c(mh, rep(0, K-length(mh)))
	  }

      K <- length(unique(ci[ci!=0]))
#	  print(K); print(length(mh))
#	  print(mh)
      # Relabel after removing those that will be reallocated (i.e., gamma = 0)
      if(sum(r) < N){
	    ci[ci!=0]<-as.numeric(as.character(factor(ci[ci!=0], labels=order(unique(ci[ci!=0]),decreasing=FALSE))))
	  }
#	  print(ci)
  	
	  for(k in dnk){
#	    cat("k = ", k, "\n")
#	    print(mh)
		p <- 1
		if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+M), M/(sum(mh[mh!=0])+M))
#		print(p)
#		cat("K + 1 = ", K+1, "\n")

		ci[k] <- sample(1:(K+1), 1, prob=p)
		mh <- table(ci[ci!=0])
		K <- length(unique(ci[ci!=0]))
#		print(ci)
#		print(table(ci))
	  }
#	  print(ci)
	  citmp <- factor(ci,  labels=order(unique(ci), decreasing=FALSE))

#	  print(citmp)
	  ciMat[t,] <- as.numeric(as.character(citmp))
	  if(model=="markovian") ci <- as.numeric(as.character(citmp))
	  if(model=="independent") ci <- c0
    }  	
  }
  out$ciMat <- ciMat
  out$gMat <- gMat
  out
}

relabel <- function(x){
  nc <- length(unique(x))
  tmp <- x
  if(x[1] != 1){
  	tmp <- x
  	tmp[x==1] <- x[1]
  	tmp[x==x[1]] <- 1
  }
  as.numeric(factor(tmp, levels=unique(tmp), labels=c(1:nc)))
}




Tm <- 10
c0 <- rep(1:4, each=5)
alpha <- matrix(rep(c(0.95,0.75,0.5,0.25), each=5), nrow=length(c0), ncol=Tm+1, byrow=FALSE)
nMC <- 1000

markov <- t(replicate(n=nMC, rtpartition(N=length(c0), M=1, alpha=alpha, ntime=Tm, c0=c0, model="markovian")))
independent <- t(replicate(n=nMC, rtpartition(N=length(c0), M=1, alpha=alpha, ntime=Tm, c0=c0, model="independent")))
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

library(fields)

tmp1 <- apply(sapply(ari_mar,function(x) x[[1]]),1,mean)
tmp2 <- apply(sapply(ari_ind,function(x) x[[1]]),1,mean)

pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/mark_indep.pdf", height=5, width=13)
  par(mfrow=c(1,3))
  image.plot(1:10, 1:10, 
             matrix(apply(sapply(ari_mar,function(x) x[[2]]),1,mean), nrow=Tm, ncol=Tm, byrow=TRUE),
             zlim=c(0,1), main = "Markovian", xlab="Lag", ylab="Lag" ,cex.lab=1.5)
  image.plot(1:10, 1:10, 
             matrix(apply(sapply(ari_ind,function(x) x[[2]]),1,mean), nrow=Tm, ncol=Tm, byrow=TRUE),
             zlim=c(0,1), main = "Independdent", xlab="Lag", ylab="Lag" ,cex.lab=1.5)
             
  plot(1:10, tmp1, type='b', ylab="", xlab="Lag", cex=1.5, ,cex.lab=1.5)
  points(1:10, tmp2, type='b', pch=2, cex=1.5)
  legend(x=6.5, y=0.3, legend=c("Markovian","Independent"), pch=c(1,2), cex=1.5)
  mtext(side=2, "Adjusted Rand Index", line = 2)
dev.off()


## Functions needed to generate data

# November 2 2018
# N - number of experimental units
# M - DP scale (dispesrsion) parameter
# rho - temporal dependence parameter
# ntime - number of time points
# tau - sd associated with atom generation
# sig - sd associated with data generation
# Caron - logical indicating if temporal dependence follows Caron
# FirstPart - Allows me to provide the first partition if so desired
# Space - logical indicating of sPPM should be used to generate partition or not
# N <- 3;M <- 1;rho <- 0.5;ntime <- 2;Caron=FALSE;FirstPart=NULL
rtpartition <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL){
	library(plyr)
	out <- NULL

	ci <- FirstPart

	if(is.null(FirstPart)){
		# Create cluster at time period t = 1
		ci <- 1 # cluster label
		mh <- 1 # vector of cluster sizes
		K <- 1 # number of clusters
		Tm <- ntime

		for(k in 2:N){

			p <- c(mh/(k-1+M), M/(k-1+M))
			ci <- c(ci, sample(1:(K+1), 1, prob=p))
			mh <- table(ci)
			K <- length(unique(ci))
		}

	}

#	print(ci)
#	print(mustar)
#	plot(1:N, ci, col=ci, pch=19, cex=1.25)

	# Now to generate Tm-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	ciMat <- matrix(NA, nrow=Tm, ncol=N)
	ciMat[1,]  <- ci
	K <-  length(unique(ci))
	gMat <- matrix(0,nrow=Tm, ncol=N)

#	print(ciMat)
	if(Tm > 1){
		if(Caron){
			for(t in 2:Tm){
				cit <- NULL
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci <- ci[-dnk]
				}	
				mh <- tabulate(ci)
#				K <- length(unique(ci))
				K <- length(mh)

				for(k in 1:N){

					p <- c(mh/(sum(mh)+M), M/(sum(mh)+M))
					cit[k] <- sample(1:(K+1), 1, prob=p)
					ci <- c(ci, cit[k])
					mh <- table(ci)
					K <- length(unique(ci))
				}
	
				ciMat[t,] <- cit
	
			}

		}else{
		
			for(t in 2:Tm){
#				cat("t = ", t, "\n")
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci[dnk]<- 0
					gMat[t,dnk] <- 1

				}	

#				print(ci)

				mh <- tabulate(ci[ci!=0]);
				if((K - length(mh)) > 0){
					mh <- c(mh, rep(0, K-length(mh)))
				}

				K <- length(unique(ci[ci!=0]))


#				print(K); print(length(mh))

#				print(mh)

				if(r < N){
				
					ci[ci!=0]<-as.numeric(as.character(factor(ci[ci!=0], labels=order(unique(ci[ci!=0]),decreasing=FALSE))))

				}
#				print(ci)



				for(k in dnk){
#					cat("k = ", k, "\n")
#					print(mh)
					p <- 1
					if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+M), M/(sum(mh[mh!=0])+M))
#					print(p)
#					cat("K + 1 = ", K+1, "\n")

					ci[k] <- sample(1:(K+1), 1, prob=p)
					mh <- table(ci[ci!=0])
					K <- length(unique(ci[ci!=0]))
#					print(ci)
#					print(table(ci))
				}
#				print(ci)
				citmp <- factor(ci,  labels=order(unique(ci), decreasing=FALSE))

#				print(citmp)
				ciMat[t,] <- as.numeric(as.character(citmp))
#				print(ciMat)

			}
		
		}
	}
	out$ciMat <- ciMat
	out$gMat <- gMat
	out
}

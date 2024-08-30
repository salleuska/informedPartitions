## need to write an autoregressive process that generate
## alphas = probability of reallocation
####################
## N - number of experimental units
## M - DP scale (dispersion) parameter
## alpha - temporal dependence parameter
## ntime - number of time points
## FirstPart - Allows to provide the first partition if desired
## clusterSpecific - use if partition probability depending on clustering 
####################

logit <- function(x) log(x/(1 - x))
expit <- function(x) exp(x)/(1 + exp(x))

N <- 3;M <- 1;ntime <- 10;FirstPart=c(1,2,1)
# alpha <-  matrix(0.1, ntime, N)
# rtpartition(N, M, alpha, ntime, FirstPart )

## Generate a sequence of alphas (one individual - ntime points)

rep <- 1000
zeta0 <- 3

alpha <- numeric(ntime)
zeta  <- numeric(ntime)

## run simulations multiple times
alphaTest <- matrix(NA, ncol = ntime, nrow = 1000)


for(i in 1:rep){
	kappa  <- runif(1, min = 0, max = 0.5)
	phi     <- rnorm(1, mean = 0, sd = 1)

	epsilon <- rnorm(ntime, mean = 0, sd = kappa2)
	zeta[1] <- phi * zeta0 + epsilon[1]
	alpha[1] <- expit(zeta[1])
	for(t in 2:ntime) {
		zeta[t] <-  phi * zeta[t -1] + epsilon[t]
		alpha[t] <- expit(zeta[t])
	}

	alphaTest[i, ] <- alpha

}
hist(alphaTest[, 6])

######################

## zeta_1 = phi x zeta_0 + eps --- eps ~ N(0, kappa2)
kappa2 <- numeric(N)
phi    <- numeric(N)



alpha1 <- expit(zeta0)

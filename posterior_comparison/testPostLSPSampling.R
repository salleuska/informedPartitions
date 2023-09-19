## test
library(Rcpp)
library(RcppArmadillo)

sourceCpp("utils/LSP_dmrpfunctions.cpp")

muTrue <- c(rnorm(50, 0, 1), rnorm(50, 4, 1))
data <- rnorm(100, muTrue, 1)

## initial partition and variance
initialPart<- rep(1:2, each = 50)

## useful constant values
n <- length(data)
tau <-  2
tau_proposal <- 1

# tau <-  0.1/(n*log(n))
# tau_proposal <- 1/(n*log(n))

## fixed variance for model likelihood
likelihood_variance<- 1

## priors for clusters means 
## same prior for each cluster
tau_k <- 10  ## prior variance
mu_k <- 0   ## prior mean


set.seed(878)

## initialize MCMC quantities and create output matrices
niter = 10000

## store niter + 1 values (initial values)
partitionSamples <- matrix(NA,  nrow = niter, ncol = n)
muSamples <- matrix(NA,  nrow = niter, ncol = n)

## sample initial values and relabel  from 1, no gaps
partitionInit <- sample(1:4, n, replace = T)
partitionInit <- as.numeric(as.character(factor(partitionInit, labels=order(unique(partitionInit),decreasing=FALSE))))
# partitionInit <- 1:100
muInit <- rnorm(n, mean = mu_k, sd = sqrt(tau_k))

acceptanceVec <- numeric(niter)
# rlsp(10, partitionInit, tau)

for(i in 1:niter){

	## current values of partition and mu 
	if(i == 1) muCurrent = muInit; partCurrent = partitionInit;

	## Propose a value for the partition and mu
	## propose a partition
	partProposal <- rlsp(1, partCurrent, tau_proposal)
	# table(partProposal, partCurrent)
	## and corresponding mu vector
	KProposal <- length(unique(partProposal[1, ]))
	muProposal <- numeric(n)
	
	for(k in 1:KProposal){
		idx_k <- which(partProposal == k)
		n_k <- length(data[idx_k])

		var_k <- var(data[idx_k])
		posterior_variance_k <- 1 / ((1 / tau_k) + (n_k / likelihood_variance))
	 	posterior_mean_k <- posterior_variance_k * ((mu_k / tau_k) + (sum(data[idx_k]) / likelihood_variance))

		muProposal[idx_k] <- rnorm(n_k, posterior_mean_k, sd = sqrt(posterior_variance_k))
	}

	# Calculate the likelihood, prior, and unnormalized posterior for the current value of mu
	likelihoodCurrent<- sum(dnorm(data, mean = muCurrent, sd = sqrt(likelihood_variance), log = TRUE))
	priorCurrent<- sum(dnorm(muCurrent, mean = mu_k, sd = sqrt(tau_k), log = TRUE)) + log(dlsp(partCurrent,initialPart,tau))
	posteriorCurrent<- likelihoodCurrent + priorCurrent
  
  	## Calculate the likelihood for the proposed value of mu
	likelihoodProposal<- sum(dnorm(data, mean = muProposal, sd = sqrt(likelihood_variance), log = TRUE))
    # Calculate the prior for the proposed value of mu (assuming a normal prior)
	priorProposal<- sum(dnorm(muProposal, mean = mu_k, sd = tau_k, log = TRUE)) + log(dlsp(partProposal,initialPart,tau))
    # Calculate the unnormalized posterior for the proposed value of mu
	posteriorProposal<- likelihoodProposal+ priorProposal


	## Calculate transition probabilities
	## current given proposal
	transCurrentProposal <- sum(dnorm(muCurrent, mean = muProposal, sd = tau_k, log = TRUE)) + log(dlsp(partCurrent,partProposal, tau_proposal))
	## proposal given current 
	transProposalCurrent <- sum(dnorm(muProposal , mean = muCurrent, sd = tau_k, log = TRUE)) + log(dlsp(partProposal, partCurrent, tau_proposal))

	# Calculate the acceptance ratio
	acceptance_ratio <- exp(posteriorProposal + transCurrentProposal - posteriorCurrent - transProposalCurrent)
	
	# Accept or reject the proposed value based on the acceptance ratio
	if (runif(1) < min(acceptance_ratio, 1)) {
	  partCurrent <- partProposal
	  muCurrent <- muProposal  # Accept the proposal
	  acceptanceVec[i] <- 1
	}

	muSamples[i, ] <- muCurrent	  
	partitionSamples[i, ] <- partCurrent

	if(i%%100 == 0) cat(i, " iterations\n")
}


sum(acceptanceVec)

plot(muSamples[, 100], type = "l")
plot(muSamples[, 34], type = "l")

density(apply(muSamples, 2, mean)) |> plot()
density(muInit) |> plot()
plot(apply(muSamples, 2, mean), muInit)

table(partitionSamples[2, ], initialPart)

plot(apply(muSamples, 2, mean), muTrue)


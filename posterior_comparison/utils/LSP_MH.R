###############
## this implement an MCMC for DP model using LSP proposal 
## using a metropolis-hastings algorithm
###############
## libraries to compile function from 
## https://github.com/adam-n-smith/demand-models-random-partitions/tree/master

library(Rcpp)
library(RcppArmadillo)

sourceCpp("utils/LSP_dmrpfunctions.cpp")

## data - vector
## niter - iteration
## likelihood variance
## mu_k s2_k - cluster piors paraemters
## tau - LDP variance
## tau_proposal = lsp proposal variance

LSP_MH <- function(data, niter, likelihood_variance, 
	mu_k = 0, s2_k = 10, tau,  initialPart, partitionInit, tau_proposal){

	n <- length(data)
	# store niter values (wiht first equal to initial values)
	partitionSamples <- matrix(NA,  nrow = niter, ncol = n)
	muSamples <- matrix(NA,  nrow = niter, ncol = n)
	acceptanceVec <- numeric(niter)

	for(i in 1:niter){

	  ########################################
	  ## Initialization - first iteration
	  ########################################
	  ## set current partition to initialization value
	  ## initialize muCurrent conditionally to the partition
	  if(i == 1) {
	    partCurrent = partitionInit;

	    Kcurrent <- length(unique(partCurrent))
	    muCurrent <- numeric(Kcurrent)
	    clustMeansCurrent <- numeric(Kcurrent)
	    clustVarsCurrent  <- numeric(Kcurrent)

	    for(k in 1:Kcurrent){
	      idx_k <- which(partCurrent == k)
	      n_k <- length(data[idx_k])

	      var_k <- var(data[idx_k])
	      posterior_variance_k <- 1 / ((1 / s2_k) + (n_k / likelihood_variance))
	      posterior_mean_k <- posterior_variance_k * ((mu_k / s2_k) + (sum(data[idx_k]) / likelihood_variance))
	  		
	  	  # alternative computation - same results
	      # post_mean_Alt <-  mu_k *(likelihood_variance/(likelihood_variance/n + s2_k) ) +
	      #             mean(data[idx_k]) *(s2_k /(likelihood_variance/n + s2_k))
	  
	      clustVarsCurrent[k] <- posterior_variance_k
	      clustMeansCurrent[k] <- posterior_mean_k
	  
	      muCurrent[k] <-  rnorm(1, posterior_mean_k, sd = sqrt(posterior_variance_k))
	    }
	  }
	  ########################################

	  ########################################
	  ## Proposal 
	  ########################################
	  ## Propose a value for the partition and mu
	  ## propose a partition
	  partProposal <- rlsp(1, partCurrent, tau_proposal)[1, ]
	  # table(partProposal, partCurrent)

	  ## and corresponding mu vector - conditionally on the partition
	  KProposal <- length(unique(partProposal))
	  muProposal <- numeric(KProposal)
	  
	  ## collect values of the posterior distribution for each cluster
	  clustMeansProposal <- numeric(KProposal)
	  clustVarsProposal  <- numeric(KProposal)

	  for(k in 1:KProposal){
	    idx_k <- which(partProposal == k)
	    n_k <- length(data[idx_k])

	    var_k <- var(data[idx_k])
	    posterior_variance_k <- 1 / ((1 / s2_k) + (n_k / likelihood_variance))
	    posterior_mean_k <- posterior_variance_k * ((mu_k / s2_k) + (sum(data[idx_k]) / likelihood_variance))

	    clustMeansProposal[k] <- posterior_mean_k
	    clustVarsProposal[k] <- posterior_variance_k

	    muProposal[k] <- rnorm(1, posterior_mean_k, sd = sqrt(posterior_variance_k))
	  }

	  # plot(1:length(muProposal), muProposal, col = partProposal, pch = 16)
	  ########################################
	  ## Posterior under the two states
	  ########################################

	  # Calculate the likelihood, prior, and unnormalized posterior for the current value of mu
	  likelihoodCurrent<- sum(dnorm(data, mean = muCurrent[partCurrent], sd = sqrt(likelihood_variance), log = TRUE))
	  priorCurrent<- sum(dnorm(muCurrent, mean = mu_k, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partCurrent,initialPart,tau))
	  posteriorCurrent<- likelihoodCurrent + priorCurrent
	  
	  ## Calculate the likelihood for the proposed value of mu
	  likelihoodProposal<- sum(dnorm(data, mean = muProposal[partProposal], sd = sqrt(likelihood_variance), log = TRUE))
	  # Calculate the prior for the proposed value of mu 
	  priorProposal<- sum(dnorm(muProposal, mean = mu_k, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partProposal,initialPart,tau))
  	  # Calculate the unnormalized posterior for the proposed value of mu
	  posteriorProposal<- likelihoodProposal + priorProposal

	  ########################################
	  ## Calculate transition probabilities
	  ########################################
  
 	  ## this should be what is done in Smith and Allenby?
	  ## current given proposal - p(partCurrent | partProposal) * p(muCurrent |partProposal = (clustMeansProposal,  clustVarsProposal))
	  transCurrentProposal <- sum(dnorm(muCurrent, mean = clustMeansCurrent , sd = sqrt(clustVarsCurrent), log = TRUE)) + log(dlsp(partCurrent, partProposal, tau_proposal))
	  ## proposal given current - p(partProposal| partCurrent) * p(muProposal |partCurrent = c(lustMeansCurrent,  clustVarsCurrent))
	  transProposalCurrent <- sum(dnorm(muProposal, mean = clustMeansProposal, sd = sqrt(clustVarsProposal), log = TRUE)) + log(dlsp(partProposal, partCurrent, tau_proposal))  	

	  
	  # Calculate the acceptance ratio
	  acceptance_ratio <- exp(posteriorProposal + transCurrentProposal - posteriorCurrent - transProposalCurrent)
	  
	  # Accept or reject the proposed value based on the acceptance ratio
	  if (runif(1) < min(acceptance_ratio, 1)) {
	    partCurrent <- partProposal
	    muCurrent <- muProposal  # Accept the proposal
	 
	    ## Record clusters parameters
	    clustMeansCurrent <- clustMeansProposal
	    clustVarsCurrent <- clustVarsProposal

	    acceptanceVec[i] <- 1
	  }

	  muSamples[i, ] <- muCurrent[partCurrent]  
	  partitionSamples[i, ] <- partCurrent

	  if(i%%1000 == 0) cat(i, " iterations\n")
	}
	## return list
	list(muSamples = muSamples, partitionSamples = partitionSamples, acceptanceVec = acceptanceVec)	
}



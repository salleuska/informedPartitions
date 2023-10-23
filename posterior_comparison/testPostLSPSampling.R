## test
library(Rcpp)
library(RcppArmadillo)

sourceCpp("utils/LSP_dmrpfunctions.cpp")

##############
## check random sim from rlsp
# part <- rep(1:5, each = 20)
# draw <- rlsp(10000, part , tau = 0.5)

# png(filename = "smithA_05.png", width = 600, height = 600)
# x <- compsm(draw)
# lattice::levelplot(x)
# dev.off()

# part <- rep(1:5, each = 20)
# draw <- rlsp(10000, part , tau = 0.05)

# png(filename = "smithA_005.png", width = 600, height = 600)
# x <- compsm(draw)
# lattice::levelplot(x)
# dev.off()

# # levelplot(x, col.regions = rev(gray(0:100/100)))

# png(filename = "smithA_onePart.png", width = 600, height = 600)
# draw <- rlsp(10000, partitionInit , tau = 0.1)
# x <- compsm(draw)
# lattice::levelplot(x)
# dev.off()

##############
## generate 100 observation from two separate mixtures
muTrue <- c(rep(-2, 50), rep(2, 50))
data <- rnorm(100, muTrue, 1)

## initial partition equal to true partition
initialPart<- rep(1:2, each = 50)

## useful constant values
n <- length(data)
## informative prior (small variance)
## and proposal variance
# tau <-  0.1/(n*log(n))
# tau_proposal <- 1/(n*log(n))

# tau <-  2
# tau_proposal <- 1

## fixed variance for model likelihood
likelihood_variance<- 1

## priors for clusters means 
## same prior for each cluster
s2_k <- 10  ## prior variance
mu_k <- 0   ## prior mean

##################
set.seed(878)

## initialize MCMC quantities and create output matrices
niter = 10000

## store niter + 1 values (initial values)
partitionSamples <- matrix(NA,  nrow = niter, ncol = n)
muSamples <- matrix(NA,  nrow = niter, ncol = n)

# sample initial values and relabel  from 1, no gaps
# partitionInit <- sample(1:5, n, replace = T)
# partitionInit <- as.numeric(as.character(factor(partitionInit, labels=order(unique(partitionInit),decreasing=FALSE))))
partitionInit <- rep(1, n)

# tau <-  1/(n*log(n))
# tau <-  1/n

# draw <- rlsp(10000, initialPart , tau = tau)
# x <- compsm(draw)
# lattice::levelplot(x)

# # tau_proposal <- 1/(n*log(n))
# # tau_proposal <- 10/n
# tau_proposal <- 1/(n*sqrt(log(n))) ## try this!!

# draw <- rlsp(10000, initialPart , tau = tau_proposal)
# x <- compsm(draw)
# lattice::levelplot(x)


tau <-  0.5/(n*log(n))
tau_proposal <- 1/(n*log(n)) 
acceptanceVec <- numeric(niter)

for(i in 1:niter){

	########################################
	## Initialization - first iteration
	########################################
	## set current partition to initialization value
	## initialize muCurrent conditionally to the partition
	if(i == 1) {
		partCurrent = partitionInit;

		muCurrent <- numeric(n)
		Kcurrent <- length(unique(partCurrent))
		clustMeansCurrent <- numeric(n)
		clustVarsCurrent  <- numeric(n)

		for(k in 1:Kcurrent){
			idx_k <- which(partCurrent == k)
			n_k <- length(data[idx_k])

			var_k <- var(data[idx_k])
			posterior_variance_k <- 1 / ((1 / s2_k) + (n_k / likelihood_variance))
		 	posterior_mean_k <- posterior_variance_k * ((mu_k / s2_k) + (sum(data[idx_k]) / likelihood_variance))
	
		 	clustVarsCurrent[idx_k] <- posterior_variance_k
			clustMeansCurrent[idx_k] <- posterior_mean_k
	
			muCurrent[idx_k] <-  rnorm(1, posterior_mean_k, sd = sqrt(posterior_variance_k))
		}
	}
	########################################

	########################################
	## Proposal 
	########################################
	## Propose a value for the partition and mu
	## propose a partition
	partProposal <- rlsp(1, partCurrent, tau_proposal)[1, ]
	#table(partProposal, partCurrent)

 	## and corresponding mu vector - conditionally on the partition
	KProposal <- length(unique(partProposal))
	muProposal <- numeric(n)
	
	## collect values of the posterior distribution for each cluster
	clustMeansProposal <- numeric(n)
	clustVarsProposal  <- numeric(n)

	for(k in 1:KProposal){
		idx_k <- which(partProposal == k)
		n_k <- length(data[idx_k])

		var_k <- var(data[idx_k])
		posterior_variance_k <- 1 / ((1 / s2_k) + (n_k / likelihood_variance))
	 	posterior_mean_k <- posterior_variance_k * ((mu_k / s2_k) + (sum(data[idx_k]) / likelihood_variance))

		clustMeansProposal[idx_k] <- posterior_mean_k
	 	clustVarsProposal[idx_k] <- posterior_variance_k

		muProposal[idx_k] <- rnorm(1, posterior_mean_k, sd = sqrt(posterior_variance_k))
	}

	# plot(1:length(muProposal), muProposal, col = partProposal, pch = 16)
	########################################
	## Posterior under the two states
	########################################

	# Calculate the likelihood, prior, and unnormalized posterior for the current value of mu
	likelihoodCurrent<- sum(dnorm(data, mean = muCurrent, sd = sqrt(likelihood_variance), log = TRUE))
	priorCurrent<- sum(dnorm(muCurrent, mean = mu_k, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partCurrent,initialPart,tau))
	posteriorCurrent<- likelihoodCurrent + priorCurrent
  
  	## Calculate the likelihood for the proposed value of mu
	likelihoodProposal<- sum(dnorm(data, mean = muProposal, sd = sqrt(likelihood_variance), log = TRUE))
    # Calculate the prior for the proposed value of mu 
	priorProposal<- sum(dnorm(muProposal, mean = mu_k, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partProposal,initialPart,tau))
    # Calculate the unnormalized posterior for the proposed value of mu
	posteriorProposal<- likelihoodProposal+ priorProposal

	########################################
	## Calculate transition probabilities
	########################################
	## current given proposal - p(partCurrent | partProposal) * p(muCurrent |partProposal = (clustMeansProposal,  clustVarsProposal))
	transCurrentProposal <- sum(dnorm(muCurrent, mean = clustMeansProposal, sd = sqrt(clustVarsProposal), log = TRUE)) + log(dlsp(partCurrent,partProposal, tau_proposal))
	## proposal given current - - p(partProposal| partCurrent) * p(muProposal |partCurrent = c(lustMeansCurrent,  clustVarsCurrent))
	transProposalCurrent <- sum(dnorm(muProposal , mean = clustMeansCurrent, sd = sqrt(clustVarsCurrent), log = TRUE)) + log(dlsp(partProposal, partCurrent, tau_proposal))

	## SP: OLD (ERROR?) Calculate transition probabilities
	## current given proposal
	## Is this normal distribution the right one?
	# transCurrentProposal <- sum(dnorm(muCurrent, mean = muProposal, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partCurrent,partProposal, tau_proposal))
	# ## proposal given current 
	# transProposalCurrent <- sum(dnorm(muProposal , mean = muCurrent, sd = sqrt(s2_k), log = TRUE)) + log(dlsp(partProposal, partCurrent, tau_proposal))


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

	muSamples[i, ] <- muCurrent	  
	partitionSamples[i, ] <- partCurrent

	if(i%%100 == 0) cat(i, " iterations\n")
}

burn <- niter/2
keep <- seq((burn+1), niter, by = 1)
sum(acceptanceVec)

plot(muSamples[, 100], type = "l")
plot(muSamples[, 34], type = "l")


x <- compsm(partitionSamples[keep,])
lattice::levelplot(x)

cc <- salso(partitionSamples[keep, ])
table(cc, initialPart)
density(apply(muSamples[keep, ], 2, mean)) |> plot()
plot(apply(muSamples[keep, ], 2, mean), muTrue)


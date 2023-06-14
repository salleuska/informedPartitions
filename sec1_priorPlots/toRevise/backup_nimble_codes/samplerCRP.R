

#' @rdname samplers
#' @export
sampler_CRP <- nimbleFunction(
  name = 'sampler_CRP',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    if(!is.null(control$printTruncation))
      printMessage <- control$printTruncation else printMessage <- TRUE

    targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    targetVar <- model$getVarNames(nodes = target)
    n <- length(targetElements)
    
    ## Find nodes indexed by the CRP node.
    if(is.null(control$clusterVarInfo)) {
        clusterVarInfo <- findClusterNodes(model, target)
    } else clusterVarInfo <- control$clusterVarInfo
    tildeVars <- clusterVarInfo$clusterVars

    ##  Various checks that model structure is consistent with our CRP sampler. 
    
    if(!is.null(clusterVarInfo$targetNonIndex))
      stop("sampler_CRP: Detected that the CRP variable is used in some way not as an index: ", clusterVarInfo$targetNonIndex, ". NIMBLE's CRP sampling not set up to handle this case.")
    
    if(any(clusterVarInfo$nTilde == 0)) {
      var <- which(clusterVarInfo$nTilde == 0)
      stop("sampler_CRP: Detected unusual indexing in ", safeDeparse(clusterVarInfo$indexExpr[[var[1]]]),
           " . NIMBLE's CRP MCMC sampling is not designed for this case.")
    }
    
    if(is.null(tildeVars))
      stop('sampler_CRP: The model should have at least one cluster variable.\n')
    
    ## Cases like 'muTilde[xi[n-i+1]]'. sampler_CRP may be ok with this, but when we wrap the cluster node sampling
    ## to avoid sampling empty clusters, this kind of indexing will cause incorrect behavior.
    ## This case is trapped in findClusterNodes.
    
    allTildeNodes <- unlist(clusterVarInfo$clusterNodes)
    dataNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE) # 'data' from the perspective of the clustering model
    stochDepsTildeNodes <- model$getDependencies(allTildeNodes, self = FALSE, stochOnly = TRUE)
    
    ## Make sure tildeNodes as determined from clustering actually are in model.
    if(!all(allTildeNodes %in% model$getNodeNames())) {
      missingNodes <- allTildeNodes[which(!allTildeNodes %in% model$getNodeNames())]
      stop("sampler_CRP: These cluster parameters are not nodes in the model: ", paste(missingNodes, collapse = ','))
    }
    
    ## Check that no other non-data nodes depend on cluster variables. 
    if(!identical(sort(dataNodes), sort(stochDepsTildeNodes)))
      stop("sampler_CRP: Only the variables being clustered can depend on the cluster parameters.")  

    ## Check that nodes in different clusters are distinct.
    ## E.g., this would not be the case if a user specified a joint prior on all cluster nodes
    ## Also check that clustering is done on stochastic nodes.
    for(varIdx in seq_along(clusterVarInfo$clusterVars)) {
        if(length(unique(clusterVarInfo$clusterNodes[[varIdx]])) != length(clusterVarInfo$clusterNodes[[varIdx]]))
            stop("sampler_CRP: cluster parameters in different clusters must be part of conditionally independent nodes.")
        if(any(model$isDeterm(clusterVarInfo$clusterNodes[[varIdx]])))
            stop("findClusterNodes: detected that deterministic nodes are being clustered. Please use the dCRP-distributed node to cluster stochastic nodes.")
    }
    
    ## Check that membership variable is independent of cluster nodes.
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(target %in% stochDepsTildeNodes)
      stop("sampler_CRP: Cluster membership variable has to be independent of cluster parameters.")
    
    ## Check that cluster nodes are independent of membership variable
    ## (dataNodes are the dependents of target and should not contain cluster parameters).
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(length(intersect(dataNodes, allTildeNodes)))
      stop("sampler_CRP: Cluster parameters have to be independent of cluster membership variable.")
    
    ## Check that observations are independent of each other.
    ## In non-conjugate case, this could potentially be relaxed within each cluster, provided we figure
    ## out correct ordering of dataNodes plus intermNodes in calculate().
    dataNodesIDs <- model$getDependencies(target, stochOnly = TRUE, self = FALSE, returnType = 'ids')
    sapply(dataNodes, function(x) {
      if(any(dataNodesIDs %in% model$getDependencies(x, self = FALSE, stochOnly = TRUE, returnType = 'ids')))
        stop("sampler_CRP: Variables being clustered must be conditionally independent. To model dependent variables being clustered jointly, you may use a multivariate distribution.")
    })

    ## Check that if same clusterNodes used twice in a declaration that they are used identically
    ## e.g., dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]])) is ok
    ## but dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]+1])) is not as can't properly add a new cluster
    ## (or account for when not to sample an empty cluster).
    if(length(unique(tildeVars)) != length(tildeVars)) 
        for(idx1 in 2:length(tildeVars))
            for(idx2 in 1:(length(tildeVars)-1))
                if(!identical(clusterVarInfo$clusterNodes[[idx1]], clusterVarInfo$clusterNodes[[idx2]]) &&
                   any(clusterVarInfo$clusterNodes[[idx1]] %in% clusterVarInfo$clusterNodes[[idx2]]))
                    stop("sampler_CRP: Inconsistent indexing in or inconsistent dependencies of ",
                         safeDeparse(clusterVarInfo$indexExpr[[idx1]]), " and ",
                         safeDeparse(clusterVarInfo$indexExpr[[idx2]]), ".")
    
    nData <- length(dataNodes)

    ## Check that no use of multiple clustering variables, such as 'thetaTilde[xi[i], eta[j]]'.
    ## It's likely that if we set the non-conjugate sampler and turn off wrapping omits sampling of empty clusters
    ## (which is not set up correctly for this case), that the existing code would give correct sampling.
    if(any(clusterVarInfo$multipleStochIndexes))
        stop("sampler_CRP: Detected use of multiple stochastic indexes of a variable: ", safeDeparse(clusterVarInfo$indexExpr[[1]]), ". NIMBLE's CRP sampling is not yet set up to handle this case. Please contact the NIMBLE development team if you are interested in this functionality.")

    ## Check there is at least one or more "observation" per random index.
    ## Note that cases like mu[xi[i],xi[j]] are being trapped in findClusterNodes().
    if(n > nData)
       stop("sampler_CRP: At least one variable has to be clustered for each cluster membership ID.")
    
   
    p <- length(tildeVars) 
    nTilde <- clusterVarInfo$nTilde / clusterVarInfo$numNodesPerCluster  ## implied number of potential clusters
    if(length(unique(nTilde)) != 1)
        stop('sampler_CRP: In a model with multiple cluster parameters, the number of those parameters must all be the same.\n')
    min_nTilde <- nTilde[1]
    if(min_nTilde < n)
      messageIfVerbose('  [Warning] sampler_CRP: The number of clusters based on the cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if it ever proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.')
    
    ## Determine if concentration parameter is fixed or random (code similar to the one in sampleDPmeasure function).
    ## This is used in truncated case to tell user if model is proper or not.
    fixedConc <- TRUE
    parentNodesTarget <- model$getParents(target, stochOnly = TRUE)
    if(length(parentNodesTarget)) {
      fixedConc <- FALSE
    }

    ## Here we try to set up some data structures that allow us to do observation-specific
    ## computation, to save us from computing for all observations when a single cluster membership is being proposed.
    ## At the moment, this is the only way we can easily handle dependencies for multiple node elements in a
    ## 'vectorized' way.
    nObsPerClusID <- nData / n # equal to one in standard CRP model (former J)
    dataNodes <- rep(targetElements[1], nData) ## this serves as dummy nodes that may be replaced below
    ## needs to be legitimate nodes because run code sets up calculate even if if() would never cause it to be used
    for(i in seq_len(n)) { # dataNodes are always needed so only create them before creating  intermNodes
      stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self = FALSE)
      if(length(stochDeps) != nObsPerClusID)
          stop("sampler_CRP: The number of nodes that are jointly clustered must be the same for each group. For example if 'mu[i,j]' are clustered such that all nodes for a given 'i' must be in the same cluster, then the number of nodes for each 'i' must be the same. Group ", safeDeparse(i), " has ", safeDeparse(length(stochDeps)), " nodes being clustered where ", safeDeparse(nObsPerClusID), " are expected.")
      dataNodes[((i-1)*nObsPerClusID + 1) : (i*nObsPerClusID)] <- stochDeps
    }
    nIntermClusNodesPerClusID <- length(model$getDependencies(targetElements[1], determOnly = TRUE))  #nInterm
    intermNodes <- dataNodes # initialize nodes in case not used (needed for compilation to go through)
    if(nIntermClusNodesPerClusID > 0) {
      intermNodes <- rep(as.character(NA), nIntermClusNodesPerClusID * n)
      for(i in seq_len(n)) {
        detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
        if(length(detDeps) != nIntermClusNodesPerClusID)
            stop("sampler_CRP: The number of intermediate deterministic nodes that are jointly clustered must be the same for each group. For example if 'mu[i,j]' are clustered such that all nodes for a given 'i' must be in the same cluster, then the number of nodes for each 'i' must be the same. Group ", safeDeparse(i), " has ", safeDeparse(length(detDeps)), " intermediate nodes being clustered where ", safeDeparse(nIntermClusNodesPerClusID), " are expected.")
        intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)] <- detDeps
      }
    }
    
    helperFunctions <- nimbleFunctionList(CRP_helper)
        
    ## use conjugacy to determine which helper functions to use
    if(control$checkConjugacy) {
        conjugacyResult <- checkCRPconjugacy(model, target) 
    } else conjugacyResult <- NULL
    
    if(is.null(conjugacyResult)) {
      sampler <- 'CRP_nonconjugate'
    } else 
      sampler <- switch(conjugacyResult,
                        conjugate_dnorm_dnorm_identity = 'CRP_conjugate_dnorm_dnorm',
                        conjugate_dmnorm_dmnorm_identity = 'CRP_conjugate_dmnorm_dmnorm',
                        conjugate_dnorm_dnorm_additive_nonidentity = 'CRP_conjugate_dnorm_dnorm_nonidentity',
                        conjugate_dnorm_dnorm_multiplicative_nonidentity = 'CRP_conjugate_dnorm_dnorm_nonidentity',
                        conjugate_dnorm_dnorm_linear_nonidentity = 'CRP_conjugate_dnorm_dnorm_nonidentity',
                        conjugate_dinvgamma_dnorm_identity = 'CRP_conjugate_dinvgamma_dnorm',
                        conjugate_dinvwish_dmnorm_identity = 'CRP_conjugate_dinvwish_dmnorm',
                        conjugate_dwish_dmnorm_identity = 'CRP_conjugate_dwish_dmnorm',
                        conjugate_dnorm_invgamma_dnorm = 'CRP_conjugate_dnorm_invgamma_dnorm',
                        conjugate_dnorm_gamma_dnorm = 'CRP_conjugate_dnorm_gamma_dnorm',
                        conjugate_dmnorm_invwish_dmnorm = 'CRP_conjugate_dmnorm_invwish_dmnorm',
                        conjugate_dmnorm_wish_dmnorm = 'CRP_conjugate_dmnorm_wish_dmnorm',
                        conjugate_dbeta_dbern_identity  = 'CRP_conjugate_dbeta_dbern',
                        conjugate_dbeta_dbin_identity = 'CRP_conjugate_dbeta_dbin',
                        conjugate_dbeta_dnegbin_identity = 'CRP_conjugate_dbeta_dnegbin',
                        conjugate_dgamma_dpois_identity = 'CRP_conjugate_dgamma_dpois',
                        conjugate_dgamma_dexp_identity = 'CRP_conjugate_dgamma_dexp',
                        conjugate_dgamma_dgamma_identity = 'CRP_conjugate_dgamma_dgamma',
                        conjugate_dgamma_dnorm_identity = 'CRP_conjugate_dgamma_dnorm',
                        conjugate_dgamma_dweib_identity = 'CRP_conjugate_dgamma_dweib',
                        conjugate_dgamma_dinvgamma_identity = 'CRP_conjugate_dgamma_dinvgamma',
                        conjugate_ddirch_dmulti_identity = 'CRP_conjugate_ddirch_dmulti',
                        'CRP_nonconjugate')  ## default if we don't have sampler set up for a conjugacy


    clusterIDs <- unique(clusterVarInfo$clusterIDs[[1]])
    nClusters <- length(clusterVarInfo$clusterIDs)
    if(nClusters > 1) {
        ## Check that each set of tildeNodes indicate same number of clusters.        
        sapply(2:nClusters, function(i) {
            if(!identical(clusterIDs, unique(clusterVarInfo$clusterIDs[[i]])))
                stop("sampler_CRP: differing number of clusters indicated by ", paste0(clusterVarInfo$clusterNodes[[1]], collapse = ', '), " and ", paste0(clusterVarInfo$clusterNodes[[i]], collapse = ', '), ".")
        })
    }

    nClusNodesPerClusID <- sum(clusterVarInfo$numNodesPerCluster)

    ## Determine correct order of clusterNodes, including any intermediate nodes
    ## standing between different clusterNodes. This set of nodes are the
    ## 'marginalized' nodes.
    ## Also check independence of cluster parameters across clusters.
    ## This block seems to increase buildMCMC time by about 30%.
    allNodes <- unlist(clusterVarInfo$clusterNodes)
    marginalizedNodes <- allNodes
    dataIntermNodes <- c(dataNodes, intermNodes)
    ids <- unlist(clusterVarInfo$clusterIDs)
    if(nClusNodesPerClusID > 1) {
        for(id in seq_along(clusterIDs)) {
            origNodes <- allNodes[ids == id]
            nodes <- model$getDependencies(origNodes)
            nodes <- nodes[!nodes %in% dataIntermNodes]        
            if(id == min(ids)) {
                nClusNodesPerClusID <- length(nodes)  # now includes determ intermediates between cluster nodes
                marginalizedNodes <- rep('', nClusNodesPerClusID * min_nTilde)
            }
            if(length(nodes) != nClusNodesPerClusID)
                stop("sampler_CRP: detected differing number of cluster parameters across clusters or dependence of parameters across clusters.") # unequal number could occur because of dependence across clusters
            marginalizedNodes[((id-1)*nClusNodesPerClusID+1):(id*nClusNodesPerClusID)] <- nodes
            if(any(nodes %in% allNodes[!allNodes %in% origNodes]))
                stop("sampler_CRP: cluster parameters must be independent across clusters.")
        }
        
    } else {
        marginalizedNodes <- clusterVarInfo$clusterNodes[[1]]
        deps <- unique(unlist(lapply(marginalizedNodes, function(x)
            model$getDependencies(x, self = FALSE, stochOnly = TRUE))))
        if(any(marginalizedNodes %in% deps))
            stop("sampler_CRP: cluster parameters must be independent across clusters.")
        if(!identical(sort(deps), sort(dataNodes)))
            messageIfVerbose("  [Warning] sampler_CRP: dependencies of cluster parameters include unexpected nodes: ",
                    paste0(deps[!deps %in% dataNodes], collapse = ', '))
    }
    
    identityLink <- TRUE

    if(p == 2 && sampler %in% c("CRP_conjugate_dnorm_invgamma_dnorm", "CRP_conjugate_dnorm_gamma_dnorm", 
                                "CRP_conjugate_dmnorm_invwish_dmnorm", "CRP_conjugate_dmnorm_wish_dmnorm")) {
        if(sampler == "CRP_conjugate_dnorm_invgamma_dnorm") {
            dist1 <- 'dnorm'
            dist2 <- 'dinvgamma'
        }
        if(sampler == "CRP_conjugate_dnorm_gamma_dnorm") {
            dist1 <- 'dnorm'
            dist2 <- 'dgamma'
        }
        if(sampler == "CRP_conjugate_dmnorm_invwish_dmnorm") {
            dist1 <- 'dmnorm'
            dist2 <- 'dinvwish'
        }
        if(sampler == "CRP_conjugate_dmnorm_wish_dmnorm") {
            dist1 <- 'dmnorm'
            dist2 <- 'dwish'
        }
        for(i in seq_along(tildeVars)) {
            if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == dist1) {
                marginalizedNodes1 <- clusterVarInfo$clusterNodes[[i]]
            } 
            if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == dist2) {
                marginalizedNodes2 <- clusterVarInfo$clusterNodes[[i]]
            }
        }  
      helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes1, marginalizedNodes2, dataNodes, nObsPerClusID, nClusNodesPerClusID)
      calcNodes <- model$getDependencies(c(target, marginalizedNodes1, marginalizedNodes2))
    } else {
      calcNodes <- model$getDependencies(c(target, marginalizedNodes))
      if(sampler == "CRP_conjugate_dnorm_dnorm_nonidentity") {
          identityLink <- FALSE
          helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes, dataNodes, intermNodes, nIntermClusNodesPerClusID, nObsPerClusID, nClusNodesPerClusID)
      } else {
          if(sampler != "CRP_nonconjugate" && nObsPerClusID != nClusNodesPerClusID)
              ## Our conjugate samplers use J (nObsPerClusID) instead of M (nClusNodesPerClusID) because they assume J=M.
              ## We shouldn't get to this point (based on conjugacy checking), but catch this if we do.
              stop("Number of observations per group does not equal number of cluster parameters per group. NIMBLE's CRP sampling is not set up to handle this except for the non-conjugate sampler.")
          helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes, dataNodes, nObsPerClusID, nClusNodesPerClusID)
      }
    }

    curLogProb <- numeric(n)
  },
  
  
  run = function() {
    
    conc <- model$getParam(target, 'conc')
    helperFunctions[[1]]$storeParams()
    
    xi <- model[[target]]
    
    ## Find unique values in model[[target]].
    ## We don't relabel the unique values, but we do create each new cluster as the lowest unused positive integer.
    ## k denotes the number of unique labels in xi
    
    xiUniques <- numeric(min_nTilde)
    xiCounts <- numeric(n)
    
    aux <- min(xi):max(xi) 
    k <- 1
    for(i in seq_along(aux)) { 
      nMembers <- sum(aux[i] == xi)
      if(nMembers > 0) {
        xiCounts[aux[i]] <- nMembers
        xiUniques[k] <- aux[i]
        k <- k + 1
      }
    }
    k <- k-1 # number of unique labels in xi
    
    kNew <- 1 # kNew is the new label that can be sampled
    while(xiCounts[kNew] > 0 & kNew < n) { 
      kNew <- kNew + 1
    }
    if( kNew == n & xiCounts[kNew] > 0 ) {  # all clusters are filled (with singletons)
      kNew <- 0
    }
    if(kNew > min_nTilde & min_nTilde < n) {
      if(printMessage) {
        if(fixedConc) {
          nimCat('  [Warning] CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
        } else {
          nimCat('  [Warning] CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
        }
      }
      kNew <- 0 
      printMessage <<- FALSE 
    }
    
    
    for(i in 1:n) { # updates one cluster membership at the time , i=1,...,n
      
      xi <- model[[target]]
      xiCounts[xi[i]] <- xiCounts[xi[i]] - 1
      
      # Computing sampling probabilities and sampling an index.
      if( xiCounts[xi[i]] == 0 ) { # cluster is a singleton.
        
        ## First, compute probability of sampling an existing label.
        reorderXiUniques <- numeric(min_nTilde) # here we save reordered version of xiUniques when there is a singleton. This is used later for updating xiUniques if a component is deleted.
        iprob <- 1
        for(j in 1:k) {
          if( xiCounts[xiUniques[j]] >= 1 ) { 
            model[[target]][i] <<- xiUniques[j] # <<-
            if(nIntermClusNodesPerClusID > 0) {
              model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
            }
            curLogProb[iprob] <<- log(xiCounts[xiUniques[j]]) +
                model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])
            reorderXiUniques[iprob] <- xiUniques[j]
            iprob <- iprob + 1
          }
        }
        
        ## Second, compute probability of sampling a new cluster, here, new cluster is the current cluster!
        model[[target]][i] <<- xi[i] # <<- label of new component
        if(sampler == 'CRP_nonconjugate'){ # simulate tildeVars[xi[i]] # do this everytime there is a singleton so we ensure this comes always from the prior
          helperFunctions[[1]]$sample(i, model[[target]][i])
          if(nIntermClusNodesPerClusID > 0) {
            model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
          }
          model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
        }
        if(!identityLink) 
            helperFunctions[[1]]$calculate_offset_coeff(i, model[[target]][i])
        curLogProb[k] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # <<- probability of sampling a new label, only k components because xi_i is a singleton
        
        ## Sample new cluster.
        index <- rcat( n=1, exp(curLogProb[1:k]-max(curLogProb[1:k])) )
        if(index == k) {
          newLab <- xi[i] 
          newLabCond <- TRUE
        } else {
          newLab <- reorderXiUniques[index]
          newLabCond <- FALSE
        }
        
      } else { # cluster is not a singleton.
        ## First, compute probability of sampling an existing label.
        for(j in 1:k) { 
          model[[target]][i] <<- xiUniques[j]  
          if(nIntermClusNodesPerClusID > 0) {
            model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
          }
          curLogProb[j] <<- log(xiCounts[xiUniques[j]]) +
              model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
        }
        ## Second, compute probability of sampling a new cluster depending on the value of kNew.       
        if(kNew == 0) { # no new cluster can be created 
          curLogProb[k+1] <<- log(0)  # <<- k+1 <= n always because k==n requires all singletons, handled above
        } else { # a new cluster can be created
          model[[target]][i] <<- kNew 
          if(sampler == 'CRP_nonconjugate'){
            helperFunctions[[1]]$sample(i, model[[target]][i])
            if(nIntermClusNodesPerClusID > 0) {
              model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
            }
            model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
          }
          if(!identityLink) 
            helperFunctions[[1]]$calculate_offset_coeff(i, model[[target]][i])
          curLogProb[k+1] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # <<- probability of sampling a new label
        }
        
        # sample an index from 1 to (k+1)
        index <- rcat( n=1, exp(curLogProb[1:(k+1)]-max(curLogProb[1:(k+1)])) )
        if(index == (k+1)) {
          newLab <- kNew
          newLabCond <- TRUE
        } else {
          newLab <- xiUniques[index]
          newLabCond <- FALSE
        }
      }
      
      ## Update metadata about clustering.
      model[[target]][i] <<- newLab 
      
      if( newLabCond ) { # a component is created. It can really create a new component or keep the current label if xi_i is a singleton
        if(sampler != 'CRP_nonconjugate') { # updating the cluster parameters of the new cluster
          helperFunctions[[1]]$sample(i, model[[target]][i])
        }
        if( xiCounts[xi[i]] != 0) { # a component is really created
          k <- k + 1
          xiUniques[k] <- newLab 
          kNew <- kNew + 1
          mySum <- sum(xi == kNew) 
          while(mySum > 0 & kNew < (n+1)) { # need to make sure don't go beyond length of vector
            kNew <- kNew+1
            mySum <- sum(xi == kNew)
          }
          if(kNew > min_nTilde & min_nTilde < n) {
            if( printMessage ) {
              if(fixedConc) {
                nimCat('CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
              } else {
                nimCat('CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
              }
            }
            kNew <- 0
            printMessage <<- FALSE 
          }
        }
        xiCounts[model[[target]][i]] <- 1
      } else { # an existing label is sampled
        ## Reset to previous marginalized node value; we choose to store information on what elements to be restored in sample()
        ## but an alternative would be to have i=0 determine reset and pass j=kNew here.
        if(sampler == 'CRP_nonconjugate')   
          helperFunctions[[1]]$sample(i, 0)
        if( xiCounts[xi[i]] == 0 ) { # xi_i is a singleton, a component was deleted
          k <- k - 1
          xiUniques <- reorderXiUniques
          if( kNew == 0 ) { # the sampler was not nonparametric or xi=1:n
            kNew <- xi[i] 
          } else { # the sampler was and remains nonparametric.
            if( kNew > xi[i] ) {
              kNew <- xi[i]
            }
          }
        }
        xiCounts[model[[target]][i]] <- xiCounts[model[[target]][i]] + 1
      }
    }
    
    ## We have updated cluster variables but not all logProb values are up-to-date.
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( 
    reset = function () {
      printMessage <<- TRUE
    }
  )
)

findClusterNodes <- function(model, target) {
  ## Determine which model nodes are the cluster parameters by processing expressions to look
  ## for what is indexed by the dCRP clusterID nodes. This also determine which clusterID
  ## each cluster parameter is associated with.
  targetVar <- model$getVarNames(nodes = target)
    if(model$getVarInfo(targetVar)$nDim > 1)
        stop("findClusterNodes: CRP variable, '", targetVar, "' cannot be a matrix or array.")
  targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  deps <- model$getDependencies(target, self = FALSE, returnType = 'ids')
  declIDs <- model$modelDef$maps$graphID_2_declID[deps] ## declaration IDs of the nodeIDs
  uniqueIDs <- unique(declIDs)
  depsByDecl <- lapply(uniqueIDs, function(x) deps[which(x == declIDs)])

  ## Find one example dependency per BUGS declaration for more efficient processing
  exampleDeps <- model$modelDef$maps$graphID_2_nodeName[sapply(depsByDecl, `[`, 1)]
    
  ## Once we find the cluster parameter variables below, we want to evaluate the cluster membership
  ## values (e.g., xi[1],...,xi[n]) for all possible values they could take, this will
  ## allow us to determine all possible cluster nodes in the model (though some may
  ## not actually be specified in the model, if there is truncation).
  ## Therefore, set up an evaluation environment in which (xi[1],...,xi[n]) = (1,2,...,n)
  ## first try was: e[[targetVar]] <- seq_along(targetElements)
  ## However in first try, that wouldn't handle xi[3:10] ~ dCRP(), but next construction does.
  e <- list()
  idxExpr <- model$getDeclInfo(target)[[1]]$indexExpr[[1]]
  eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(targetElements)), list(VAR = targetVar, IDX = idxExpr)))
  ## For cases of cross clustering (e.g., mu[xi[i],eta[j]]) we need the other dcrp node(s)
  nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  dists <- model$getDistribution(nodes)
  if(length(dists == 'dCRP') > 1) { 
      dcrpNodes <- nodes[dists == 'dCRP' & nodes != target]
      for(i in seq_along(dcrpNodes)) {
          dcrpElements <- model$expandNodeNames(dcrpNodes[i], returnScalarComponents = TRUE)
          dcrpVar <- model$getVarNames(nodes = dcrpNodes[i])
          idxExpr <- model$getDeclInfo(dcrpNodes[i])[[1]]$indexExpr[[1]]
          eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(dcrpElements)), list(VAR = dcrpVar, IDX = idxExpr)))

      }
  }
  clusterNodes <- indexExpr <- clusterIDs <- list()
  clusterVars <- indexPosition <- numIndexes <- targetIsIndex <- targetIndexedByFunction <-
      loopIndex <- NULL
  varIdx <- 0

  targetNonIndex <- NULL
  multipleStochIndexes <- NULL

  modelVars <- model$getVarNames()
  modelVars <- modelVars[!modelVars == targetVar]

  ## Process model declaration expressions to find stochastic indexing and the indexed variable.
  for(idx in seq_along(exampleDeps)) {
    ## Pull out expressions, either as RHS of deterministic or parameters of stochastic
    fullExpr <- cc_getNodesInExpr(model$getValueExpr(exampleDeps[idx]))
    for(j in seq_along(fullExpr)) {
      subExpr <- parse(text = fullExpr[j])[[1]]  # individual parameter of stochastic or RHS of deterministic
      len <- length(subExpr)
      ## Look for target variable within expression, but only when used within index
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' &&
        sum(all.vars(subExpr) == targetVar) && subExpr[[2]] != targetVar) {
        varIdx <- varIdx + 1
        multipleStochIndexes <- c(multipleStochIndexes, FALSE)
          
        clusterVars <- c(clusterVars, safeDeparse(subExpr[[2]], warn = TRUE))
        
        ## Determine which index the target variable occurs in.
        k <- whichIndex <- 3
        foundTarget <- FALSE
        while(k <= len) {
            if(sum(all.vars(subExpr[[k]]) == targetVar)) {
                if(foundTarget) {
                    stop("findClusterNodes: CRP variable used multiple times in ", safeDeparse(subExpr),
                           ". NIMBLE's CRP MCMC sampling not designed for this situation.")
                } else {
                    foundTarget <- TRUE
                    whichIndex <- k
                }
            }
            ## We will need to relax this when allow crossed clustering.
            if(sum(all.vars(subExpr[[k]]) %in% modelVars)) { ## cases like mu[xi[i],eta[j]]
                ## We are adding support for this case.
                ## warning("findClusterNodes: multiple indexing variables in '", deparse(subExpr),
                ##          "'. NIMBLE's CRP MCMC sampling not designed for this situation.", call. = FALSE)
                multipleStochIndexes[varIdx] <- TRUE
            }                
            k <- k+1
        }
        if(!foundTarget) stop("findClusterNodes: conflicting information about presence of CRP variable in expression.")

        declInfo <-  model$getDeclInfo(exampleDeps[idx])[[1]]
        
        ## Determine how target variable enters into cluster node definition
        indexPosition[varIdx] <- whichIndex-2
        numIndexes[varIdx] <- len - 2
        indexExpr[[varIdx]] <- subExpr
        ## Is target used directly as index, e.g., "mu[xi[.]]" as opposed to something like "mu[1+xi[.]]".
        targetIsIndex[varIdx] <- length(subExpr[[whichIndex]]) == 3 &&
          subExpr[[whichIndex]][[1]] == '[' &&
          subExpr[[whichIndex]][[2]] == targetVar
        ## Is indexing of target a simple index, e.g. xi[i], as opposed to something like "xi[n-i+1]".
        targetIndexedByFunction[varIdx] <- any(sapply(declInfo$symbolicParentNodes,
                                                    function(x) 
                                                        length(x) >= 3 && x[[1]] == '[' &&
                                                        x[[2]] == targetVar && length(x[[3]]) > 1))
        ## Determine all sets of index values so they can be evaluated in context of possible values of target element values.
        unrolledIndices <- declInfo$unrolledIndicesMatrix

        if(targetIndexedByFunction[varIdx] && ncol(unrolledIndices) > 1)  ## Now that we allow cluster parameters with multiple indexes, this is very hard to handle in terms of identifying what column of unrolledIndices to use for sorting clusterNodes.
            stop("findClusterNodes: Detected that a cluster parameter is indexed by a function such as 'mu[xi[n-i+1]]' rather than simple indexing such as 'mu[xi[i]]'. NIMBLE's CRP MCMC sampling not designed for this case.")
        loopIndexes <- unlist(sapply(declInfo$symbolicParentNodes,
                                    function(x) {
                                        if(length(x) >= 3 && x[[1]] == '[' &&
                                           x[[2]] == targetVar) return(safeDeparse(x[[3]], warn = TRUE))
                                        else return(NULL) }))
        if(length(loopIndexes) != 1)
            stop("findClusterNodes: found cluster membership parameters that use different indexing variables; NIMBLE's CRP sampling not designed for this case.")
        ## Note not clear when NULL would be the result...
        loopIndex[varIdx] <- loopIndexes

        ## Determine potential cluster nodes by substituting all possible clusterID values into the indexing expression. 
        n <- nrow(unrolledIndices)
        if(n > 0 && loopIndex[varIdx] %in% dimnames(unrolledIndices)[[2]]) {  # catch cases like use of xi[2] rather than xi[i]
            ## Order so that loop over index of cluster ID in order of cluster ID so that
            ## clusterNodes will be grouped in chunks of unique cluster IDs for correct
            ## sampling of new clusters when have multiple obs per cluster.
            ord <- order(unrolledIndices[ , loopIndex[varIdx]])
            unrolledIndices <- unrolledIndices[ord, , drop = FALSE]
            clusterIDs[[varIdx]] <- unrolledIndices[ , loopIndex[varIdx]]

            clusterNodes[[varIdx]] <- rep(NA, n)
            
            ## Determine unevaluated expression, e.g., muTilde[xi[i],j] not muTilde[xi[1],2]
            expr <- declInfo$valueExprReplaced
            expr <- parse(text = cc_getNodesInExpr(expr)[[j]])[[1]]
            templateExpr <- expr   
            
            ## Now evaluate index values for all possible target element values, e.g.,
            ## xi[i] for all 'i' values with xi taking values 1,...,n
            for(i in seq_len(n)) { 
                for(k in 3:len) # this will deal with muTilde[xi[i], j] type cases
                    if(length(all.vars(expr[[k]])))  # prevents eval of things like 1:3, which the as.numeric would change to c(1,3)
                        templateExpr[[k]] <- as.numeric(eval(substitute(EXPR, list(EXPR = expr[[k]])),
                                                             c(as.list(unrolledIndices[i,]), e)))  # as.numeric avoids 1L, 2L, etc.
                clusterNodes[[varIdx]][i] <- safeDeparse(templateExpr, warn = TRUE)  # convert to node names
            }
        } else {
            clusterNodes[[varIdx]] <- character(0)
            clusterIDs[[varIdx]] <- numeric(0)
        }
      } 
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' && subExpr[[2]] == targetVar)
          targetNonIndex <- safeDeparse(model$getDeclInfo(exampleDeps[idx])[[1]]$codeReplaced, warn = TRUE)
    }
  }
  ## Find the potential cluster nodes that are actually model nodes,
  ## making sure that what we decide are real cluster nodes are the full potential set
  ## or a truncated set that starts with the first cluster node, e.g., muTilde[1], ..., muTilde[3] is ok;
  ## muTilde[2], ..., muTilde[4] is not (unless the model nodes are muTilde[2], ...., muTilde[4]).
  nTilde <- sapply(clusterNodes, length)
  modelNodes <- model$getNodeNames()
  for(varIdx in seq_along(clusterVars)) {
    if(nTilde[varIdx]) {
        if(any(is.na(clusterNodes[[varIdx]])))  
            stop("findClusterNodes: fewer cluster IDs in ", target, " than elements being clustered.")
        
        ## Handle cases where indexing of variables in dynamic indexing does not correspond to actual
        ## stochastic model nodes.
        if(any(!clusterNodes[[varIdx]] %in% modelNodes)) {
            clusterNodes[[varIdx]] <- lapply(clusterNodes[[varIdx]], function(x) model$expandNodeNames(x))
            clusterIDs[[varIdx]] <- rep(clusterIDs[[varIdx]], times = sapply(clusterNodes[[varIdx]], length))
            clusterNodes[[varIdx]] <- unlist(clusterNodes[[varIdx]])
        }
        ## Now remove duplicates when indexed variables correspond to same model node,
        ## but only for duplicates within a cluster.
        groups <- split(clusterNodes[[varIdx]], clusterIDs[[varIdx]])
        dups <- unlist(lapply(groups, duplicated))
        clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][!dups]
        clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][!dups]

        ## Formerly we were checking that we had a contiguous set of cluster nodes
        ## starting with the first one, but for clusterNodes with more than one index and
        ## truncation this is hard to do, so just fall back to returning the clusterNodes
        ## that are actually part of the model.
        validNodes <- clusterNodes[[varIdx]] %in% modelNodes
        
        if(!all(validNodes)) {  # i.e., truncated representation
            clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][validNodes]
            clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][validNodes]
        }
    }
  }

  nTilde <- sapply(clusterNodes, length)
  numNodesPerCluster <- sapply(clusterIDs, function(x) {
      tbl <- table(x)
      num <- unique(tbl)
      if(length(num) > 1) stop("findClusterNodes: detected differing numbers of nodes (i.e., parameters) per cluster. NIMBLE's CRP sampling not designed for this case.")
      return(num)})
      
  return(list(clusterNodes = clusterNodes, clusterVars = clusterVars, nTilde = nTilde,
              numNodesPerCluster = numNodesPerCluster, clusterIDs = clusterIDs, loopIndex = loopIndex,
              targetIsIndex = targetIsIndex, indexPosition = indexPosition, indexExpr = indexExpr,
              numIndexes = numIndexes, targetIndexedByFunction = targetIndexedByFunction,
              targetNonIndex = targetNonIndex, multipleStochIndexes = multipleStochIndexes))
}
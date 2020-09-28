# this is the modelfile for separate meta-analysis 
# for humerus and femur data 
# the two meta-analyses are run in the same file to allow easy 
# construction of results 
#

model {
    
#### Humerus data
    for(k1 in 1:K.hum) {
        
        # observational part of the model
        x.down.1dot[k1] ~ dbin(p.down.1dot[k1], N.down.hum[k1])
        x.nodown.1dot[k1] ~ dbin(p.nodown.1dot[k1], N.nodown.hum[k1])
        
        p.down.1dot[k1] <- ilogit(eta.xi.hum[k1, 1]) 
        p.nodown.1dot[k1] <- ilogit(eta.xi.hum[k1, 2]) 
        
        # structural part of the model 
        eta.xi.hum[k1, 1:2] ~dmnorm(Eta.Xi.hum[1:2], Omega.hum)
        
    }
    
    # prior for mean -- independent vague priors suffice 
    #                -- otherwise, mvnormal with vague Wishart 
    Eta.Xi.hum[1] ~ dnorm(mu0.hum[1], prec.mu0.hum[1])
    Eta.Xi.hum[2] ~ dnorm(mu0.hum[2], prec.mu0.hum[2])


    # S.hum is a diagonal matrix of sd's
    S.hum[1,1] ~ dunif(s.lower.hum[1], s.upper.hum[1])  # or use inverse gamma
    S.hum[2,2] ~ dunif(s.lower.hum[2], s.upper.hum[2])
    S.hum[1,2] <- 0
    S.hum[2,1] <- S.hum[1,2]
    
    # R.hum is a correlation matrix
    R.hum[1,1] <- 1
    R.hum[2,2] <- 1
    R.hum[1,2] ~ dunif(r.lower.hum[1], r.upper.hum[1])  # or use atanh
    R.hum[2,1] <- R.hum[1,2]
    
    # Factorize
    Tau.hum <- S.hum %*% R.hum %*% S.hum
    Omega.hum <- inverse(Tau.hum) 



#### Femur data
    for(k2 in 1:K.fem) {
    
        # observational part of the model
        x.down.dot1[k2] ~ dbin(p.down.dot1[k2], N.down.fem[k2])
        x.nodown.dot1[k2] ~ dbin(p.nodown.dot1[k2], N.nodown.fem[k2])
    
        p.down.dot1[k2] <- ilogit(eta.xi.fem[k2, 1]) 
        p.nodown.dot1[k2] <- ilogit(eta.xi.fem[k2, 2]) 

        # structural part of the model 
        eta.xi.fem[k2, 1:2] ~dmnorm(Eta.Xi.fem[1:2], Omega.fem)
    
    }

    # prior for mean -- independent vague priors suffice 
    #                -- otherwise, mvnormal with vague Wishart 
    Eta.Xi.fem[1] ~ dnorm(mu0.fem[1], prec.mu0.fem[1])
    Eta.Xi.fem[2] ~ dnorm(mu0.fem[2], prec.mu0.fem[2])


    # S.fem is a diagonal matrix of sd's
    S.fem[1,1] ~ dunif(s.lower.fem[1], s.upper.fem[1])  # or use inverse gamma
    S.fem[2,2] ~ dunif(s.lower.fem[2], s.upper.fem[2])
    S.fem[1,2] <- 0
    S.fem[2,1] <- S.fem[1,2]

    # R.fem is a correlation matrix
    R.fem[1,1] <- 1
    R.fem[2,2] <- 1
    R.fem[1,2] ~ dunif(r.lower.fem[1], r.upper.fem[1])  # or use atanh
    R.fem[2,1] <- R.fem[1,2]

    # Factorize
    Tau.fem <- S.fem %*% R.fem %*% S.fem
    Omega.fem <- inverse(Tau.fem) 
        
        
    # Data to report 
    tpr.hum <- ilogit(Eta.Xi.hum[1])
    fpr.hum <- ilogit(Eta.Xi.hum[2])

    tpr.fem <- ilogit(Eta.Xi.fem[1])
    fpr.fem <- ilogit(Eta.Xi.fem[2])
    
    diff.tpr <- tpr.hum - tpr.fem
    diff.fpr <- fpr.hum - fpr.fem
    
    or.tpr <- exp(Eta.Xi.hum[1]-Eta.Xi.fem[1])
    or.fpr <- exp(Eta.Xi.hum[2]-Eta.Xi.fem[2])

    diff.eta <- Eta.Xi.hum[1]-Eta.Xi.fem[1]
    diff.xi <- Eta.Xi.hum[2]-Eta.Xi.fem[2]

}

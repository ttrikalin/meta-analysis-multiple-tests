# load libraries 
source("./00_functions_library.r")

library(coda)
library(rjags)



# paths 
path.stata <- file.path('data', 'database_counts_all.dta', sep='')
path.model.separate <- file.path('model_separate.r', sep='')
res.path <- file.path('bayes_results', sep='')
stem <- 'separate'

# jags fitting 
model.file <- path.model.separate
n.adapt <- 10000
n.burn.in <- 50000
n.record <- 100000 


#########################################
# DATA
#########################################

data <- read.dta(path.stata)


#### humerus 
    valid.hum <- !is.na(data[,"down_n_1dot"])  # only 11 studies have humerus 
    data.hum <- data[valid.hum,]

    ## data 
    K.hum <- dim(data.hum)[1]
    
    x.down.1dot <- as.vector(data.hum[,'down_n_1dot'])
    N.down.hum <- as.vector(data.hum[,'down_n_1dot'] + data.hum[,'down_n_0dot'])

    x.nodown.1dot <- as.vector(data.hum[,'nodown_n_1dot']) 
    N.nodown.hum <- as.vector(data.hum[,'nodown_n_1dot'] + data.hum[,'nodown_n_0dot'])
    
    ## priors -- flat 
    mu0.hum <- c(0, 0)
    prec.mu0.hum <- rep(10^-6, 2)
    s.lower.hum <- rep(10^-4, 2)
    s.upper.hum <- rep(5, 2)
    r.lower.hum <- rep(-.999, 2)
    r.upper.hum <- rep(+.999, 2)
    
    
#### femur 
    valid.fem <- !is.na(data[,"down_n_dot1"]) # all studies have femur  - superfluous
    data.fem <- data[valid.fem,]

    ## data 
    K.fem <- dim(data.fem)[1]

    x.down.dot1 <- as.vector(data.fem[,'down_n_dot1']) 
    N.down.fem <- as.vector(data.fem[,'down_n_dot1'] + data.fem[,'down_n_dot0'])

    x.nodown.dot1 <- as.vector(data.fem[,'nodown_n_dot1']) 
    N.nodown.fem <- as.vector(data.fem[,'nodown_n_dot1'] + data.fem[,'nodown_n_dot0'])

    ## priors -- flat 
    mu0.fem <- c(0, 0)
    prec.mu0.fem <- rep(10^-6, 2)
    s.lower.fem <- rep(10^-4, 2)
    s.upper.fem <- rep(5, 2)
    r.lower.fem <- rep(-.999, 2)
    r.upper.fem <- rep(+.999, 2)


jags.data <- list('K.hum' = K.hum, 
                  'x.down.1dot' = x.down.1dot, 'N.down.hum' = N.down.hum, 
                  'x.nodown.1dot' = x.nodown.1dot, 'N.nodown.hum' = N.nodown.hum, 
                  #'alpha.down.hum' = alpha.down.hum, 'beta.down.hum' = beta.down.hum, 
                  'mu0.hum' = mu0.hum, 'prec.mu0.hum' = prec.mu0.hum, 
                  's.lower.hum' = s.lower.hum, 's.upper.hum' = s.upper.hum, 
                  'r.lower.hum' = r.lower.hum, 'r.upper.hum' = r.upper.hum,
                  
                  'K.fem' = K.fem, 
                  'x.down.dot1' = x.down.dot1, 'N.down.fem' = N.down.fem, 
                  'x.nodown.dot1' = x.nodown.dot1, 'N.nodown.fem' = N.nodown.fem, 
                  #'alpha.down.fem' = alpha.down.fem, 'beta.down.fem' = beta.down.fem, 
                  'mu0.fem' = mu0.fem, 'prec.mu0.fem' = prec.mu0.fem, 
                  's.lower.fem' = s.lower.fem, 's.upper.fem' = s.upper.fem, 
                  'r.lower.fem' = r.lower.fem, 'r.upper.fem' = r.upper.fem )

#########################################
# RUN THE MODEL AND MONITOR RESULTS 
#########################################

the.jags.model <- jags.model(file=model.file, 
                  data=jags.data, 
                  n.chains=3, n.adapt=n.adapt)


update(the.jags.model, n.iter=n.burn.in)

the.mcmc.res <- my.coda.samples(the.jags.model, 
                              c('deviance', 'pD','tpr.hum', 'fpr.hum', 
                              'tpr.fem', 'fpr.fem', 'diff.tpr', 'diff.fpr',
                              'or.tpr', 'or.fpr', 'Eta.Xi.hum', 'Eta.Xi.fem', 'diff.eta', 'diff.xi',
                              'Tau.hum[1,1]', 'Tau.hum[1,2]', 'Tau.hum[2,2]', 'R.hum[1,2]',
                              'Tau.fem[1,1]', 'Tau.fem[1,2]', 'Tau.fem[2,2]', 'R.fem[1,2]'),
                              n.iter=n.record, thin=10)


#########################################
# SPIT RESULTS TO STATA 
#########################################

## collect stats -- the first part in the list is the non-pD chains
stats <- summary(the.mcmc.res[[1]])

pD <- c("pD", summary(the.mcmc.res[[2]], mean)$stat, NA, NA, NA) 
stats.m <- cbind(rownames(stats[[1]]), stats[[1]])
stats.m <- rbind(stats.m, pD)
colnames(stats.m)[1:5]<-c("node", "mean", "sd", "naivese", "timeseriesse")



problist <- c(0.025, 0.25, 0.50, 0.75, 0.975)
pD <- c("pD", quantile(the.mcmc.res[[2]], probs=problist))
stats.q <- cbind(rownames(stats[[2]]), stats[[2]])
stats.q <- rbind(stats.q, pD)
colnames(stats.q)[1:6] <- c("node","p025","p250", "p500", "p750", "p975")

write.table(stats.m, file=paste(res.path,"/","stats_m_", stem, ".txt", sep=''), sep='\t',row.names = FALSE)
write.table(stats.q, file=paste(res.path,"/","stats_q_", stem, ".txt", sep=''), sep='\t',row.names = FALSE)

# save the trace plots
png(paste(res.path, "/", stem, "_%02d.png", sep=""))
plot(the.mcmc.res[[1]])
dev.off()

# collect diagnostics 
gelman <- gelman.diag(the.mcmc.res[[1]])

gelman.row.names <- rownames(gelman[[1]])
gelman.table <-cbind(gelman.row.names, gelman[[1]], rep(gelman[[2]], length(gelman.row.names)))
colnames(gelman.table)<- c("node", "prsf500", "prsf975", "prsf975_multi")  
write.table(gelman.table, file=paste(res.path,"/","gelman_", stem, ".txt", sep=''), sep='\t',row.names = FALSE)


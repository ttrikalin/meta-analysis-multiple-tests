# load libraries 
source("./00_functions_library.r")

library(foreign)
library(coda)
library(rjags)

load.module('dic')


# paths 
path.stata <- file.path('data', 'database_counts_all.dta', sep='')
path.model.joint <- file.path('model_joint_nocorrelations_all.r', sep='')
res.path <- file.path('bayes_results', sep='')
stem <- 'joint_nocorrelations'

# jags fitting 
model.file <- path.model.joint
n.adapt <- 50000
n.burn.in <- 100000
n.record <- 400000 


#########################################
# DATA
#########################################

data <- read.dta(path.stata)

# Organize the data into 4 categories:
# 1. Both tests, crosstabs for TPR AND FPR
# 2. Both tests, crosstabs for TPR BUT NOT FOR FPR
# 3. Both tests, no crosstabs
# 4. Only femur


both.tests.tpr <- (!is.na(data["down_n_1dot"]) & !is.na(data[,"down_n_dot1"]))
miss.crosstab.tpr <- (both.tests.tpr==TRUE & is.na(data[,"down_n_00"]))
with.crosstab.tpr <- (both.tests.tpr==TRUE & !is.na(data[,"down_n_00"]))

both.tests.fpr <- (!is.na(data["nodown_n_1dot"]) & !is.na(data[,"nodown_n_dot1"]))
miss.crosstab.fpr <- (both.tests.fpr==TRUE & is.na(data[,"nodown_n_00"]))
with.crosstab.fpr <- (both.tests.fpr==TRUE & !is.na(data[,"nodown_n_00"]))

only.femur <- (both.tests.tpr==FALSE) # true also for fpr in the dataset


## Group 1
# observational part of the model
# x order is 11 10 01 00 
# first humerus, second femur 
#
x.down.11.G1 <- as.vector(data[with.crosstab.fpr, 'down_n_11'])
x.down.1dot.G1 <- as.vector(data[with.crosstab.fpr, 'down_n_1dot'])
x.down.dot1.G1 <- as.vector(data[with.crosstab.fpr, 'down_n_dot1'])

x.nodown.11.G1 <- as.vector(data[with.crosstab.fpr, 'nodown_n_11'])
x.nodown.1dot.G1 <- as.vector(data[with.crosstab.fpr, 'nodown_n_1dot'])
x.nodown.dot1.G1 <- as.vector(data[with.crosstab.fpr, 'nodown_n_dot1'])


N.down.G1 <- as.vector(data[with.crosstab.fpr, 'down_n'])
N.nodown.G1 <- as.vector(data[with.crosstab.fpr, 'nodown_n'])

K.G1 <- length(N.down.G1)


## Group 4
# observational part of the model
#
x.down.dot1.G4 <- as.vector(data[only.femur, 'down_n_dot1'])
x.nodown.dot1.G4 <- as.vector(data[only.femur, 'nodown_n_dot1'])

N.down.G4 <- as.vector(data[only.femur, 'down_n'])
N.nodown.G4 <- as.vector(data[only.femur, 'nodown_n'])

K.G4 <- length(N.down.G4)

I <- diag(rep(1,6))
P.G4 <- I[,c(2, 5, 1, 3, 4, 6)]



# Priors
## priors -- flat 
mu0 <- rep(0, 6)
prec.mu0 <- rep(10^-6, 6)
s.lower <- rep(10^-3, 6)
s.upper <- rep(5, 6)
chol.r.lower <- rep(10^-3, 15)
chol.r.upper <- rep(pi-10^-3, 15)



######################

jags.data <- list('K.G1' = K.G1, 'K.G4' = K.G4, 'P.G4' = P.G4, 
                  'x.down.11.G1' = x.down.11.G1, 'x.down.1dot.G1' = x.down.1dot.G1, 
                  'x.down.dot1.G1' = x.down.dot1.G1, 'N.down.G1' = N.down.G1, 
                  'x.nodown.11.G1' = x.nodown.11.G1, 'x.nodown.1dot.G1' = x.nodown.1dot.G1, 
                  'x.nodown.dot1.G1' = x.nodown.dot1.G1, 'N.nodown.G1' = N.nodown.G1, 
                  'x.down.dot1.G4' = x.down.dot1.G4, 'N.down.G4' = N.down.G4, 
                  'x.nodown.dot1.G4' = x.nodown.dot1.G4, 'N.nodown.G4' = N.nodown.G4, 

                  'mu0' = mu0, 'prec.mu0' = prec.mu0, 
                  's.lower' = s.lower, 's.upper' = s.upper, 
                  'chol.r.lower' = chol.r.lower, 'chol.r.upper' = chol.r.upper)



#########################################
# RUN THE MODEL AND MONITOR RESULTS 
#########################################

the.jags.model <- jags.model(file=model.file, 
                data=jags.data, 
                n.chains=3, n.adapt=n.adapt)


update(the.jags.model, n.iter=n.burn.in)



the.mcmc.res <- my.coda.samples(the.jags.model, 
                              c('deviance', 'pD','tpr.hum', 'fpr.hum', 'jtpr', 'jfpr',
                              'tpr.fem', 'fpr.fem', 'diff.tpr', 'diff.fpr',
                              'or.tpr', 'or.fpr', 'Eta.Xi', 'diff.eta', 'diff.xi',
                              'S[1,1]', 'S[2,2]', 'S[3,3]', 'S[4,4]', 'S[5,5]', 'S[6,6]', 
                              'R[1,2]', 'R[1,3]', 'R[1,4]', 'R[1,5]', 'R[1,6]',
                                        'R[2,3]', 'R[2,4]', 'R[2,5]', 'R[2,6]',
                                                  'R[3,4]', 'R[3,5]', 'R[3,6]', 
                                                            'R[4,5]', 'R[4,6]',
                                                                      'R[5,6]'),
                              n.iter=n.record, thin=20)





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

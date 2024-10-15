#===============
# Code kindly provided by Qing Zhao (Hierarchical Modeling in Ecology Google Group)
#===============


#===============
# Simulate data
#===============
library(boot)

n_spp <- 200  # Number of genera
n_bins <- 30   # Number of time bins
n_surveys <- 10  # Number of surveys (geological formations) per time bin

# Initial occupancy probability
initial_psi <- 0.3

# Detection probability
p <- runif(n_bins,0,1)  

# covariate effect on origination
intercept_gamma <- logit(0.2)
beta_gamma1 <- 0.5
beta_gamma2 <- 1

# covariate effect on persistence
intercept_phi <- logit(0.3)  
beta_phi1 <- 1
beta_phi2 <- -1

# covariates (generate just once and save to be used in the other two simulation sets)
X1 <- runif (n_bins-1, -2, 2)
X2 <- runif (n_bins-1, -2, 2)  

phi <- gamma <- array(NA, dim = (n_bins-1)) # persistence, colonisation
         
# produce transition probs
for(t in 1:(n_bins-1)){
         
  # origination
  gamma[t] <- inv.logit(intercept_gamma + beta_gamma1 * X1[t] + beta_gamma2 * X2[t])
           
  # persistence probability
  phi[t] <- inv.logit(intercept_phi + beta_phi1 * X1[t] + beta_phi2 * X2[t])
           
} # t
       
# Simulate occupancy states
z <- array(NA, dim = c(n_spp, n_bins))
muZ <- array(NA, dim = c(n_spp, n_bins-1))
z[,1] <- rbinom(n_spp, 1, initial_psi)
       
# true incidence
for (i in 1:n_spp) {
  for (t in 2:n_bins) {
    muZ[i,t-1] <- z[i,t-1] * phi[t-1] + (1 - z[i,t-1]) * gamma[t-1]
    z[i,t] <- rbinom(1, 1, muZ[i,t-1])
  } # t
} # i
       
# Simulate detections
y <- array(NA, dim = c(n_spp, n_bins))
for (i in 1:n_spp) {
  for (t in 1:n_bins) {
           
    y[i,t] <- rbinom(1, n_surveys, z[i,t] * p[t])
           
  } # t
} # i

#========================
# Define model in Nimble
#========================
library(nimble)

code <- nimbleCode({
     
  # Priors
  initial_psi ~ dunif(0,1)
  intercept_gamma ~ dnorm(0, sd=10)
  beta_gamma1 ~ dnorm(0, sd=10)
  beta_gamma2 ~ dnorm(0, sd=10)
  intercept_phi ~ dnorm(0, sd=10)
  beta_phi1 ~ dnorm(0, sd=10)
  beta_phi2 ~ dnorm(0, sd=10)
  for (t in 1:n_bins) {
    p[t] ~ dunif(0,1)
  }

  # Process
  for (t in 1:(n_bins-1)){
     
    # speciation
    logit(gamma[t]) <-  intercept_gamma +
                        beta_gamma1 * X1[t]+
                        beta_gamma2 * X2[t]
                                   
    # persistence
    logit(phi[t]) <-  intercept_phi +
                      beta_phi1 * X1[t]+
                      beta_phi2 * X2[t]

  } # t

  for (g in 1:n_spp) {
         
    z[g,1] ~ dbern(initial_psi) # occupancy status initialization
     
    for (t in 2:n_bins){
             
      muZ[g,t-1] <- z[g,t-1] * phi[t-1] + (1 - z[g,t-1]) * gamma[t-1]
      z[g,t] ~ dbern(muZ[g,t-1])
             
    }#t
       
  } #g

  # Observation
  for (g in 1:n_spp) { ## loop over genera
     
    for (t in 1:n_bins) { ## loop over time bins
   
      y[g,t] ~ dbin(z[g,t] * p[t], n_surveys)
                 
    } # t
     
  } # g

}) # nimbleCode

#============
# Run Nimble
#============
# Data
constants <- list(
  n_spp=n_spp, n_bins=n_bins, n_surveys=n_surveys,
  X1=X1, X2=X2
)

data <- list(
  y=y
)

# Initial values
inits <- list(initial_psi=0.5,
              intercept_gamma=0, beta_gamma1=0, beta_gamma2=0,
              intercept_phi=0, beta_phi1=0, beta_phi2=0,
              p=rep(0.5,n_bins), z=ifelse(y==0,0,1)
)

start_time <- Sys.time() # start time of computing

model <- nimbleModel(code, constants=constants, data=data, inits=inits)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

fit <- runMCMC(compiled$mcmc, nchains = 3, niter = 5000, nburnin = 0)

end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

#===============
# Check results
#===============
parms <- c('initial_psi', 'intercept_gamma', 'beta_gamma1', 'beta_gamma2', 'intercept_phi', 'beta_phi1', 'beta_phi2')
parm_values <- c(initial_psi, intercept_gamma, beta_gamma1, beta_gamma2, intercept_phi, beta_phi1, beta_phi2)
fitnms <- colnames(fit$chain1)

library(rstan)

par(mfrow=c(3,3))
par(mar=c(2,2,2,1))
for (i in 1:length(parms)) {
  tt <- cbind(fit$chain1[,which(fitnms == parms[i])],
              fit$chain2[,which(fitnms == parms[i])],
              fit$chain3[,which(fitnms == parms[i])])
  diff <- max(abs(tt[2001:5000,] - parm_values[i])) * 1.2
  plot(1, xlim=c(1,5000), ylim=c(parm_values[i]-diff,parm_values[i]+diff),
       type='n', xlab='', ylab='', main=parms[i])
  for (k in 1:3) {
    lines(tt[,k], type='l', col=k+1)
  } # k
  abline(h=parm_values[i], col='grey20')
  text(x=3000, y=parm_values[i]+diff, pos=1,
       labels=paste(c('Rhat = ', format(round(Rhat(tt[2001:5000,]), digits=2), nsmall=2)), collapse=''))
} # i
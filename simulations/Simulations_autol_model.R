# ----------------------------------------


# Simulation study

# using a simplified autologistic model in which dynamics on psi_gt accounts for the realized incidence in the previous time t-1

# ----------------------------------------
rm(list=ls())
require(here)

# create dir
dir.create (here ("simulations","output_autolog"))

# Write the JAGS model
model_string <- "
model {
  
      #  Occupancy Dynamics (Priors)
       
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        
        # intercepts
        psi.u ~ dunif(0,1) # range origination
        intercept.psi <- logit(psi.u) # intercept origination
        
        # regression coeff
        beta1 ~ dunif(-20,20)
        beta2 ~ dunif(-20,20)
        for (t in 2:n_bins){
          phi[t] ~ dunif(-20,20)
        }
        
        ## set initial conditions for occupancy of each genus
        initial_psi ~ dunif(0,1)
      
       ############      Model       #############
        
        # Occupancy dynamics ---------------------
       
        for (g in 1:n_spp) {
         
            z[g,1]~dbern(initial_psi) # occupancy status initialization
      
                for (t in 2:n_bins){
              
                 # model likelihood
                 logit(muZ[g,t]) <-  intercept.psi + 
                                     beta1*X1[t]+
                                     beta2*X2[t]+
                                     phi[t]*z[g,t-1]
                                       
                  
                 # realized occurrence
                 z[g,t] ~ dbern(muZ[g,t])
              
          }#t
        
        } #g
  
  
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    ###  detection intercept
    # intercept    
    for (t in 1:n_bins) {
      p[t] ~ dunif(0,1) 
    }
    
    # observation submodel
    for (g in 1:n_spp) { ## loop over genera 
      
      for (t in 1:n_bins) { ## loop over time bins 
    
          # observation
          # Specify the binomial observation model conditional on occupancy state
          y[g,t] ~ dbin(muY[g,t], n_surveys[t])
          muY[g,t] <- z[g,t]*p[t]
                  
        }
      
    }
      

}

"

# Load necessary libraries
library(jagsUI)
library(here)

# Set parameters
set.seed(42)
n_spp <- 200  # Number of genera
n_bins <- 30   # Number of time bins
n_surveys <- 10  # Number of surveys (geological formations) per time bin

# Initial occupancy probability
initial_psi <- 0.3

# Detection probability
p <- runif(n_bins,0,1)  

# Autologistic term
phi <- runif (n_bins, -5,5)
plot(phi,type="b")

# covariate effect on origination
intercept.psi <- qlogis(0.2) 
beta1 <- 0.5
beta2 <- -1

# covariates (generate just once and save to be used in the other two simulation sets)
X1 <- runif (n_bins, -2, 2) 
X2 <- runif (n_bins, -2, 2)  
# save(X1,X2, file=here("covariates.RData"))
# load(file=here("simulations","covariates.RData")) # activate after the creation
cor(cbind(X1,X2))

# scale covariates
X1<-scale(X1)[,1]
X2<-scale(X2)[,1]


# start simulations ---------------------------------

n.sims <- 20
my.seeds <- floor(runif (n.sims,0,5000))

# run
lapply (seq(1,n.sims), function (s) {
  
  # set seed
  set.seed(my.seeds[s])
  
  # Simulate occupancy states
  z <- array(NA, dim = c(n_spp, n_bins))
  muZ <- array(NA, dim = c(n_spp, n_bins))
  z[, 1] <- rbinom(n_spp, 1, initial_psi)
  
  # true incidence
  for (i in 1:n_spp) {
    for (t in 2:n_bins) {
      # dynamics
      muZ[i,t] <- plogis(intercept.psi+
                           beta1+X1[t]+
                           beta2+X2[t]+ 
                           phi[t]*z[i,t-1])
      z[i, t] <- rbinom(1, 1, muZ[i,t])
    }
  }
  
  # Simulate detections
  y <- array(NA, dim = c(n_spp, n_bins))
  for (i in 1:n_spp) {
    for (t in 1:n_bins) {
      
      y[i, t] <- rbinom(1, rep(n_surveys,n_bins), z[i, t] * p[t])
      
    }
  }
  
  #save data
  save (y,z ,muZ,p, phi,gamma, file = here ("simulations","output_autolog", paste0("data_", s,".RData")))
  
  # Prepare data for JAGS
  jags_data <- list(
    y = y,
    n_spp = n_spp,
    n_bins = n_bins,
    n_surveys = rep(n_surveys,n_bins),
    X1=X1,
    X2=X2
  )
  
  # Initial values for the latent states
  init_values <- function() {
    list(z = ifelse(y>0,1,0),
         psi.u = 0.2
         
    )
  }
  
  # Parameters to monitor
  parameters <- c("intercept.psi",
                  "beta1",
                  "beta2",
                  "psi.u",
                  "phi",
                  "initial_psi", 
                  "p", 
                  "muZ")
  
  # Run JAGS model
  ## MCMC runs
  samples <-jags (data = jags_data, 
                  parameters.to.save = parameters, 
                  model.file = textConnection(model_string), 
                  inits = init_values, 
                  n.chains = 3, 
                  n.thin = 1, 
                  n.iter = 500,
                  n.adapt = 200,
                  n.burnin = 400, 
                  DIC = T,  
                  parallel=F
  )
  
  # extract summary
  point_estimates <- samples$summary
  
  #save
  save (point_estimates, file = here ("simulations","output_autolog", paste0("sims_run", s,".RData")))
  
  
})

# end

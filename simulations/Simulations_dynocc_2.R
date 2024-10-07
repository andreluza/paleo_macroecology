# ----------------------------------------


# Simulation study

# Can the same set of covarites be included in the origination and extinction models??
# global - scale analysis ( so the genera are in the rows of the detection table )

# ----------------------------------------
rm(list=ls())
require(here)
dir.create (here ("simulations", "output2"))

# Write the JAGS model
model_string <- "
model {
  
      # Site Occupancy Dynamics (Priors)
       
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        
        # intercepts
        gamma.u ~ dunif(0,1) # range origination
        intercept.gamma <- logit(gamma.u) # intercept origination
        # regression coeff
        beta_gamma1 ~ dunif(-20,20)
        beta_gamma2 ~ dunif(-20,20)
        
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta_phi2 ~ dunif(-20,20)
        
        
        ## set initial conditions for occupancy of each genus
        initial_psi ~ dunif(0,1)
        
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        for (t in 1:(n_bins-1)){
      
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   beta_gamma1*X1[t]+
                                   beta_gamma2*X2[t]
                                   
              # persistence
              logit(phi[t]) <-  intercept.phi + 
                                  beta_phi2*X2[t]
                                  
                                  
        }
        
       
                       
      # occupancy dynamics
       
        for (g in 1:n_spp) {
         
            z[g,1]~dbern(initial_psi) # occupancy status initialization
      
                for (t in 2:n_bins){
              
                  # model likelihood
                  ### modeling dynamics conditional on previous time realized occurrence z
                  muZ[g,t] <- z[g,t-1] *  phi[t-1] + ### if occupied, p of not getting extinct/persist in the next time
                                (1-z[g,t-1]) *  gamma[t-1] ###  if not occupied, p of originate in the next time
                  
                 # realized occurrence
                 muZW[g,t] <- muZ[g,t]
          		   z[g,t] ~ dbern(muZW[g,t])
              
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
      p[t] ~ dunif(0,1) # constant detection
    }
    
    # observation submodel
    for (g in 1:n_spp) { ## loop over observations 
      
      for (t in 1:n_bins) { ## loop over observations 
    
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
n_spp <- 200  # Number of sites
n_bins <- 30   # Number of time bins
n_surveys <- 10  # Number of surveys (geological formations) per time bin

# Initial occupancy probability
initial_psi <- 0.3

# Detection probability
p <- runif(n_bins,0,1)  

# covariate effect on origination
intercept_gamma <- qlogis(0.2) 
beta_gamma1 <- 0.5
beta_gamma2 <- 1

# covariate effect on persistence
intercept_phi <- qlogis(0.3)  
#beta_phi1 <- 1
beta_phi2 <- -1

# covariates
load(file=here("simulations","covariates.RData"))
cor(cbind(X1,X2))

# scale
X1<-scale(X1)[,1]
X2<-scale(X2)[,1]


# start simulations ---------------------------------

n.sims <- 20
my.seeds <- floor(runif (n.sims,0,5000))

# run
lapply (seq(1,n.sims), function (s) {
  
        # set seed
        set.seed(my.seeds[s])
        
        # Origination probability
        phi <- gamma <- array(NA, dim = (n_bins-1)) # persistence, colonisation
         
        # produce transition probs
        for(t in 1:(n_bins-1)){
          
            # origination 
            gamma[t] <- plogis(intercept_gamma+beta_gamma1+X1[t]+
                               beta_gamma2+X2[t])
            
            # persistence probability
            phi[t] <- plogis(intercept_phi+
                             beta_phi2+X2[t])# back to prob scale
            
            
        }
        
        # Simulate occupancy states
        z <- array(NA, dim = c(n_spp, n_bins))
        muZ <- array(NA, dim = c(n_spp, n_bins))
        z[, 1] <- rbinom(n_spp, 1, initial_psi)
        
        # true incidence
        for (i in 1:n_spp) {
            for (t in 2:n_bins) {
              muZ[i,t] <- z[i, t - 1] * phi[t-1] + (1 - z[i, t - 1]) * gamma[t-1]
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
        save (y,z ,muZ,p, phi,gamma, file = here ("simulations", "output2", paste0("data_", s,".RData")))
        
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
               phi.u = 0.3,
               gamma.u=0.2
               
               )
        }
        
        # Parameters to monitor
        parameters <- c("intercept.phi",
                        "beta_phi2",
                        "intercept.gamma",
                        "beta_gamma1",
                        "beta_gamma2",
                        "initial_psi", 
                        "gamma",
                        "phi",
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
        save (point_estimates, file = here ("simulations", "output2", paste0("sims_run", s,".RData")))
        
        
})

# end

# ----------------------------------------


# Simulation study

# Can the same set of covarites be included in the origination and extinction models??
# global - scale analysis ( so the species are in the rows of the detection table )
# around 3 hours per run

# ----------------------------------------


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
        beta_phi1 ~ dunif(-20,20)
        beta_phi2 ~ dunif(-20,20)
        
        
        ## set initial conditions for occupancy of each genus
        initial_psi ~ dunif(0,1)
        
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        for (t in 1:(n_years-1)){
      
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   beta_gamma1*X1[t]+
                                   beta_gamma2*X2[t]
                                   
              # persistence
              logit(phi[t]) <-  intercept.phi + 
                                  beta_phi1*X1[t]+
                                  beta_phi2*X2[t]
                                  
                                  
        }
        
       
                       
      # occupancy dynamics
       
        for (g in 1:n_spp) {
         
            z[g,1]~dbern(initial_psi) # occupancy status initialization
      
                for (t in 2:n_years){
              
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
    p ~ dunif(0,1) # constant detection
     
    # observation submodel
    for (g in 1:n_spp) { ## loop over observations 
      
      for (t in 1:n_years) { ## loop over observations 
    
          # observation
          # Specify the binomial observation model conditional on occupancy state
          y[g,t] ~ dbin(muY[g,t], n_surveys[t])
          muY[g,t] <- z[g,t]*p
                  
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
n_years <- 30   # Number of time bins
n_surveys <- 6  # Number of surveys (geological formations) per time bin
initial_psi <- 0.6  # Initial occupancy probability

# covariate effect on origination
beta_gamma1 <- 0.5
beta_gamma2 <- 1

# covariate effect on extinction
beta_phi1 <- 0.5
beta_phi2 <- -1

# covariates
# PS: we repeated the generation of variables up to a point in which one correlation was of intermediate strength rho ~ 0.5
#X1 <- rnorm (n_years-1, 0, 1) 
#X2 <- rnorm (n_years-1, 0, 1) 
#X3 <- rnorm (n_years-1, 0, 1) 

#save(X1,X2,X3, file=here("simulations","covariates.RData"))
load(file=here("simulations","covariates.RData"))
cor(cbind(X1,X2,X3))

# create dir
dir.create (here ("simulations", "output"))

# start simulations ---------------------------------

n.sims <- 20
my.seeds <- floor(runif (n.sims,0,5000))

# run
lapply (seq(2,n.sims), function (s) {
  
        # set seed
        set.seed(my.seeds[s])
        
        # Origination probability
        intercept_gamma <- qlogis(0.5)  
        gamma <- intercept_gamma+beta_gamma1+X1+
                                  beta_gamma2+X2
        gamma <- plogis(gamma) # back to prob scale
        
        # Extinction probability
        intercept_phi <- qlogis(0.5)  
        phi <- intercept_phi+beta_phi1+X1+
                              beta_phi2+X2
        phi <- plogis(phi) # back to prob scale
        
        # Detection probability
        p <- 0.2  
        
        # Simulate occupancy states
        z <- array(0, dim = c(n_spp, n_years))
        z[, 1] <- rbinom(n_spp, 1, initial_psi)
        
        for (t in 2:n_years) {
          for (i in 1:n_spp) {
            z[i, t] <- rbinom(1, 1, z[i, t - 1] * (1 - phi[t-1]) + (1 - z[i, t - 1]) * gamma[t-1])
          }
        }
        
        # Simulate detections
        y <- array(0, dim = c(n_spp, n_years))
        for (i in 1:n_spp) {
          for (t in 1:n_years) {
            
              y[i, t] <- rbinom(1, n_surveys, z[i, t] * p)
            
          }
        }
        
        # Prepare data for JAGS
        jags_data <- list(
          y = y,
          n_spp = n_spp,
          n_years = n_years,
          n_surveys = rep(n_surveys,n_years),
          X1=X1,
          X2=X2
        )
        
        # Initial values for the latent states
        init_values <- function() {
          list(z = matrix(1, n_spp, n_years))
        }
        
        # Parameters to monitor
        parameters <- c("intercept.phi","beta_phi1",
                          "beta_phi2",
                        "intercept.gamma",
                        "beta_gamma2",
                        "initial_psi", "gamma", "phi", "p", "muZ")
        
        # Run JAGS model
        ## MCMC runs
        samples <-jags (data = jags_data, 
                                                  parameters.to.save = parameters, 
                                                  model.file = textConnection(model_string), 
                                                  inits = init_values, 
                                                  n.chains = 3, 
                                                  n.thin = 5, 
                                                  n.iter = 5000,
                                                  n.adapt = 2000,
                                                  n.burnin = 3000, 
                                                  DIC = T,  
                                                  n.cores=3,
                                                  parallel=T
        )
        
        # extract summary
        point_estimates <- samples$summary
        #save
        save (point_estimates, file = here ("simulations", "output", paste0("sims_run", s,".RData")))
        
        
        
        
        
})

# end of the simulations
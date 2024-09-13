# ----------------------------------------


# Simulation study

# Can the same set of covarites be included in the origination and extinction models??

# global - scale analysis ( so the species are in the rows of the detection table )

# ----------------------------------------


# Load necessary libraries
library(rjags)

# Set parameters
set.seed(42)
n_spp <- 500  # Number of sites
n_years <- 30   # Number of years (time bins)
n_surveys <- 6  # Number of surveys (geological formations) per year
initial_psi <- 0.6  # Initial occupancy probability

# covariate effect on origination
beta_gamma1 <- 0.5
beta_gamma2 <- 2
beta_gamma3 <- -0.5

# covariate effect on extinction
beta_epsilon1 <- 0.5
beta_epsilon2 <- 2
beta_epsilon3 <- -0.5

# covariates
X1 <- rnorm (n_years, 0, 1) 
X2 <- rnorm (n_years, 0, 1) 
X3 <- rnorm (n_years, 0, 1) 
cor(cbind(X1,X2,X3))

# Origination probability
intercept_gamma <- qlogis(0.1)  
gamma <- intercept_gamma+beta_gamma1+X1+
                          beta_gamma2+X2+
                          beta_gamma3+X3
gamma <- plogis(gamma) # back to prob scale

# Extinction probability
intercept_epsilon <- qlogis(0.7)  
epsilon <- intercept_epsilon+beta_epsilon1+X1+
                            beta_epsilon2+X2+
                            beta_epsilon3+X3
epsilon <- plogis(epsilon) # back to prob scale

# Detection probability
p <- 0.7  

# Simulate occupancy states
z <- array(0, dim = c(n_spp, n_years))
z[, 1] <- rbinom(n_spp, 1, initial_psi)

for (t in 2:n_years) {
  for (i in 1:n_spp) {
    z[i, t] <- rbinom(1, 1, z[i, t - 1] * (1 - epsilon[t]) + (1 - z[i, t - 1]) * gamma[t])
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
  X2=X2,
  X3=X3
)

# Initial values for the latent states
init_values <- function() {
  list(z = matrix(1, n_spp, n_years))
}

# Parameters to monitor
parameters <- c("psi", "gamma", "epsilon", "p")

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
        beta_gamma1 ~ dnorm(0,1/3^2)# precipitation
        beta_gamma2 ~ dnorm(0,1/3^2)# temperature
        beta_gamma3 ~ dnorm(0,1/3^2) # area
        
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta_epsilon1 ~ dnorm(0,1/3^2)# precipitation
        beta_epsilon2 ~ dnorm(0,1/3^2)# precipitation
        beta_epsilon3 ~ dnorm(0,1/3^2)# precipitation
        
        
        ## set initial conditions for occupancy of each genus
        initial_psi ~ dunif(0,1)
        
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        for (t in 2:n_years){
      
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   beta_gamma1*X1[t]+
                                   beta_gamma2*X2[t]+
                                   beta_gamma3*X3[t]
                                   
              # persistence
              logit(phi[t]) <-  intercept.phi + 
                                  beta_epsilon1*X1[t]+
                                  beta_epsilon2*X2[t]+
                                  beta_epsilon3*X3[t]
                                  
                                  
        }
        
       
                       
      # occupancy dynamics
       
        for (g in 1:n_spp) {
         
            z[g,1]~dbern(initial_psi) # occupancy status initialization
      
                for (t in 2:n_years){
              
                  # model likelihood
                  ### modeling dynamics conditional on previous time realized occurrence z
                  muZ[g,t] <- z[g,t-1] *  phi[t] + ### if occupied, p of not getting extinct/persist in the next time
                                (1-z[g,t-1]) *  gamma[t] ###  if not occupied, p of originate in the next time
                  
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

# Run JAGS model
jags_model <- jags.model(textConnection(model_string), data = jags_data, inits = init_values, n.chains = 3, n.adapt = 500)
update(jags_model, n.iter = 500)
samples <- coda.samples(jags_model, variable.names = parameters, n.iter = 1000, thin = 2)

# Print results
print(summary(samples))

samples <- do.call(rbind,samples)

plot(gamma,
     colMeans(samples[grep("gamma\\[",colnames(samples)),]))

plot(epsilon,
     colMeans(samples[grep("epsilon\\[",colnames(samples)),]))

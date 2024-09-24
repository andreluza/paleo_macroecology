# --------------------------------------

# 2 . modeling (GLOBAL SCALE)
# conditional on the first and last detection of the taxon


# saveJAGS help
# https://web.archive.org/web/20220408094124/http://mmeredith.net/blog/2018/Intro_saveJAGS.htm

## Global MCMC settings
######################
require(jagsUI)

# number of interactions for a global-scale analyses
na <- 40000; nb <- 60000; ni <- 80000; nc <- 3; nt <- 20

## short form for tests
na <- 6000; nb <- 7000; ni <- 10000; nc <- 3; nt <- 10
# na <- 40; nb <- 70; ni <- 100; nc <- 3; nt <- 1

# load packages
source("R/packages.R")

# load data
load(here ("processed_data","CMR_data_observation_cov.RData"))


# find the period of first and last detection
function_stages <- function (x) {
  
  # the first stage with the taxon,    
  sel_cols <- if  (max(which(colSums (x)>0)) == 33) {
    
    a <- seq(
      ifelse (min(which(colSums (x)>0)) ==1,
              min(which(colSums (x)>0)),
              min(which(colSums (x)>0))-1),
      
      
      max(which(colSums (x)>0)),
      
      1)
    
  } else {
    
    a <- seq(ifelse (min(which(colSums (x)>0)) ==1,
                     min(which(colSums (x)>0)),
                     min(which(colSums (x)>0))-1),
             max(which(colSums (x)>0))+1,
             1)
  }
  
  return(a)
  
}


# ------------------------------
# model for CMR data style
# model with no covariates

sink("global_CMRmodel_matrix_no_covariates.txt")
cat("
   
    
     model {
    
    
      #############################################################
      #                                                           #
      #                  Biological process                       #
      #                                                           #
      #############################################################
      
      # prior to describe the full community size
      #omega~dunif(0,1)
      omega ~ dbeta(0.001,1)
      
      # superpopulation process
      for (g in 1:ngen) {
        
        w[g]~dbern(omega)
        
      }
      
       # Site Occupancy Dynamics (Priors)
       # intercepts
       # gamma
       for(t in 1:(nint-1)){
       
         phi[t]~dunif(0,1)
         gamma[t]~dunif(0,1)
         
       }
       
       ## set initial conditions for each genus
       psi1 ~ dunif(0,1)
     
   ############      Model       #############
   
   # model for phi and gamma
   ### model dynamic parameters
    
    for (g in 1:ngen) {
         
            z[g,1]~dbern(psi1) # occupancy status initialization
      
                for (t in 2:nint){
      
                  # model likelihood
                  ### modeling dynamics conditional on previous time realized occurrence z
                  muZ[g,t] <- z[g,t-1] *  phi[t-1] + ### if occupied, p of not getting extinct/persist in the next time
                              (1-z[g,t-1]) *  gamma[t-1] ###  if not occupied, p of originate in the next time
                  
                 # realized occurrence
                 muZW[g,t] <- muZ[g,t]*w[g]
          		   z[g,t] ~ dbern(muZW[g,t])
              
          }#t
        
         
        
        } #g
    
  
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    # Priors for detection probability
    ###  detection intercept and coefficients
    
    for (t in 1:nint){
      p[t] ~ dunif(0,1)
    }
    
    # observation submodel
     for (g in 1:ngen) {
      for (t in 1:nint){
    
            # observation
            # Specify the binomial observation model conditional on occupancy state
            y [g,t] ~ dbin(muY[g,t],nform[t])
            muY[g,t] <- z[g,t]*p[t]      
                  
                }
        }
      
     
      
    # derived parameters
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
    # total metacommunity size
    Ntotal <- sum(w[])
    
    
    }## end of the model
    
    ",fill = TRUE)


sink()

# dir create
dir.create (here ("output", "No_covariate_model"))

# run model
no_cov_samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                           "Non-mammalian Mammaliaformes",
                                                           "Mammalia" ), function (i) {
                                                             
                                                             
                # augment the dataset
                A0_input <- array (0, dim = c(1000,
                                              length(time_bins)))
                
                
                # bundle the data                                                       
                str(jags.data <- list(
                  y =rbind (
                    
                    data.matrix(array_genus_bin [which(clades %in% i),]),
                    A0_input
                  ),
                  ngen = dim(array_genus_bin [which(clades %in% i),])[1] + nrow(A0_input),
                  nint = dim(array_genus_bin [which(clades %in% i),])[2],
                  nform = formations_per_interval
                  
                )
                )
                
                # condition the dataset to one stage before the first appearance in the fossil record, and one stage after its last appearance (detection)    
                sel_cols <- function_stages(jags.data$y)        
                
                # bundle conditional data
                str(jags.data <- list(
                  y = jags.data$y[,sel_cols],
                  ngen = nrow (jags.data$y[,sel_cols]),
                  nint = ncol(jags.data$y[,sel_cols]),
                  nform = formations_per_interval[function_stages(jags.data$y)]
                  
                )
                )
                
                
                # Set initial values
                # function inits (to be placed in each chain)
                inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                                          w = rep(1,nrow(jags.data$y)),
                                          psi1 = runif (1)
                                          
                                          
                )}
                
                
                ## Parameters to monitor
                params <- c("psi1","gamma", "omega","w","phi","p","z", "muZW","muY", "SRexp","Ntotal")
                
                # MCMC runs
                # models
                ## MCMC runs
                samples_paleo_cynodontia_binomial <-jags (data = jags.data, 
                                                          parameters.to.save = params, 
                                                          model.file = "global_CMRmodel_matrix_no_covariates.txt", 
                                                          inits = inits, 
                                                          n.chains = nc, 
                                                          n.thin = nt, 
                                                          n.iter = ni, 
                                                          n.burnin = nb, 
                                                          DIC = T,  
                                                          n.cores=nc,
                                                          parallel=F
                )
                
                
                
                save(samples_paleo_cynodontia_binomial,
                     file = paste0 ("output/No_covariate_model/CMR_global_binomial_no_cov", i,".RData"))
                
                
                
              })


rm(list=ls())
# end
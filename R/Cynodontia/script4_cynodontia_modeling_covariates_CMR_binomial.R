# --------------------------------------

# 2 . modeling
# sensitivity analysis conditioning detection history in the first detection

## Global MCMC settings
######################

# number of interactions for a global-scale analyses
na <- 10000; nb <- 15000; ni <- 20000; nc <- 3; nt <- 10

## short form for tests
#na <- 400; nb <- 700; ni <- 1000; nc <- 3; nt <- 1
# na <- 40; nb <- 70; ni <- 100; nc <- 3; nt <- 1

# load packages
source("R/packages.R")

# load data
load(here ("processed_data","CMR_data_observation_cov.RData"))

# load environmental data
load (here ("processed_data", "site_covs.RData"))


# conditionate the data to the first detection
# genus ids
ids_genus <- unique(vectorized_occ_global$genusID)
cond_det <- lapply (ids_genus, function (i) {
  
  rows_taxon <- which(vectorized_occ_global$genusID %in% i)
  rows_det <- vectorized_occ_global[rows_taxon,]
  cond_det <- rows_det [which(rows_det$detection == 1):nrow(rows_det),]
  cond_det
  
})

cond_det <- do.call(rbind, cond_det)
# adjust stage
cond_det$stage <- cond_det$stage-6
# order
cond_det <- cond_det[order(cond_det$stage),]

# ------------------------------

# model for CMR data style
# model with no covariates

sink("global_CMRmodel_vectorized_no_covariates.txt")
cat("
   
    model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    
    
   ## set initial conditions for each taxon
   psi1 ~ dunif(0,1)
    
    # dynamic params   
    for (t in 2:nint) {
        
          gamma [t]~dunif(0,1)
          phi [t]~dunif(0,1)
        
        
   }
   
   ############      Model       #############
    
    for (g in 1:ngen) {
    
          z[g,1]~dbern(psi1) # occupancy status initialization
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[g,t] <- z[g,t-1] *  phi[t] + ### if occupied, p of not getting extinct/persist in the next time
                              (1-z[g,t-1]) *  gamma[t] ###  if not occupied, p of originate in the next time
                
          # realized occurrence
  		    z[g,t] ~ dbern(muZ[g,t])
      
        } # t
        
      } # g
      

   
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
  
    # observation submodel
    for (n in 1:nobs) { ## loop over observations 
    
        p[n] ~ dunif (0,1)
        # observation
        # Specify the binomial observation model conditional on occupancy state
        y[n] ~ dbin(z[genus[n],stage[n]]*p[n], nform[n])
                  
       }     
   
      
   # derived parameters
   for (t in 1:nint){
     
       SRexp[t] <- sum (z[,t])  # true richness estimate
       
     }
  

  
   # descriptors of change
   for (t in 2:nint){
    
      FSS[t] <- sum (muZ[,t])  # finite sample size
       CH[t] <- (sum(z[,t]) - sum(z[,t-1])) ### turnover (proportional gain or loss)
       CH1[t]<-sum(z[,t-1])
       RER[t] <- (1-phi[t])/gamma[t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t] <- (1-phi[t])-gamma[t]  ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
       psi.eq[t]<-gamma[t]/((1-phi[t])+gamma[t]) # equilibrium occupancy
       turnover[t]<-(gamma[t]*(1-phi[t]))/(gamma[t]+(1-phi[t]))
      
   }
   
    
        
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()


# --------------------------------------

dir.create (here ("output", "res_det_conditional"))

## bundle data
# run
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                    "Non-mammalian Mammaliaformes",
                                                    "Mammalia" ), function (i) {
                                                      
                                                      
                                                      # select clade                                                             
                                                      cond_det <- cond_det %>% 
                                                        filter (cladeID == i) %>%
                                                        mutate (stageMod = stage-(min(stage)-1),
                                                                genus = as.numeric(as.factor(genus)))
                                                      
                                                      
                                                      # bundle data
                                                      str(jags.data <- list(
                                                        y = cond_det$detection,
                                                        ngen = length(unique(cond_det$genusID)),
                                                        nint = length(unique(cond_det$stageMod)),
                                                        nobs = nrow(cond_det),
                                                        genus = cond_det$genus,
                                                        stage = cond_det$stageMod,
                                                        nform = cond_det$formations
                                                      )
                                                      )
                                                      
                                                      
                                                      # Set initial values
                                                      # function inits (to be placed in each chain)
                                                      inits <- function(){ list(z = tapply (jags.data$y,
                                                                                            list(jags.data$genus,
                                                                                                 jags.data$stage),
                                                                                            max)
                                                                                
                                                      )}
                                                      
                                                      ## Parameters to monitor
                                                      ## long form
                                                      ## long form
                                                      params <- c(
                                                        "gamma",
                                                        "phi",
                                                        "p",
                                                        "SRexp",
                                                        "FSS",
                                                        "CH",
                                                        "CH1",
                                                        "psi1",
                                                        "R0",
                                                        "RER"
                                                      )
                                                      
                                                      
                                                      
                                                      # MCMC runs
                                                      samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                                                                                                    modelFile="global_CMRmodel_vectorized_no_covariates.txt",
                                                                                                    chains=nc, 
                                                                                                    sample2save=((ni-nb)/nt), 
                                                                                                    nSaves=6, 
                                                                                                    burnin=nb, 
                                                                                                    thin=nt,
                                                                                                    fileStub=paste ("output/res_det_conditional/CMR_global_binomial_no_cov", i,sep="_"))
                                                      
                                                    })


# -------------------------------------------------------

# covariates per interval using binomial distribution in P


# model with covariates

sink("global_CMRmodel_vectorized_covariates.txt")
cat("
   
     model {
    
    
      #############################################################
      #                                                           #
      #                  Biological process                       #
      #                                                           #
      #############################################################
      
       # Site Occupancy Dynamics (Priors)
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        # intercepts
        gamma.u ~ dunif(0,1) # range origination
        intercept.gamma <- logit(gamma.u) # intercept origination
        # regression coeff
        beta.gamma.prec ~ dnorm(0, 0.01)# precipitation
        beta.gamma.temp ~ dnorm(0, 0.01)# temperature
        beta.gamma.area ~ dnorm(0,0.01) # latitude
        
     
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta.phi.prec ~ dnorm(0, 0.01)# precipitation
        beta.phi.temp ~ dnorm(0, 0.01)# temperature
        beta.phi.area ~ dnorm(0,0.01) # area
       
        ## set initial conditions for each genus
        psi1 ~ dunif(0,1)
           
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        for (t in 2:nint){
      
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   beta.gamma.prec*precipitation[t]+
                                   beta.gamma.temp*temperature[t]+
                                   beta.gamma.area*area[t]
                                   
                                         
              # persistence
              logit(phi[t]) <-  intercept.phi + 
                                  beta.phi.prec*precipitation[t]+
                                  beta.phi.temp*temperature[t]+
                                  beta.phi.area*area[t]
                                                  
        }
       
                     
        # occupancy dynamics
         
          for (g in 1:ngen) {
           
              z[g,1]~dbern(psi1) # occupancy status initialization
        
                  for (t in 2:nint){
                
                    # model likelihood
                    ### modeling dynamics conditional on previous time realized occurrence z
                    muZ[g,t] <- z[g,t-1] *  phi[t] + ### if occupied, p of not getting extinct/persist in the next time
                                  (1-z[g,t-1]) *  gamma[t] ###  if not occupied, p of originate in the next time
                    
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
        p.u ~ dunif(0,1) # range origination
        intercept.p <- logit(p.u) # intercept origination
        beta.p.lat ~ dnorm(0,0.01)
        beta.p.range ~ dnorm(0,0.01)
        beta.p.time ~ dnorm(0,0.01)
        beta.p.temp ~ dnorm(0,0.01)
        
    
     
    # observation submodel
    for (n in 1:nobs) { ## loop over observations 
      
          # p independent on the number of formations
          logit(p[n])<- intercept.p+
                          beta.p.lat*lat_obs[n]+
                          beta.p.range*range[n]+
                          beta.p.time*time[n]+
                          beta.p.temp*temp_obs[n]
            
            
          # observation
          # Specify the binomial observation model conditional on occupancy state
          y[n] ~ dbin(z[genus[n],stage[n]]*p[n], nform[n])
                  
                  
        }
      
      
    # derived parameters
    
      for (t in 1:nint){
      
        SRexp[t] <- sum (z[,t])  # true richness estimate
        
      }
  
    
   
    
   # descriptors of change
    for (t in 2:nint){
     
       FSS[t] <- sum (muZ[,t])  # finite sample size
       CH[t] <- (sum(z[,t]) - sum(z[,t-1])) ### turnover (proportional gain or loss)
       CH1[t]<-sum(z[,t-1])
       RER[t] <- (1-phi[t])/gamma[t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t] <- (1-phi[t])-gamma[t]  ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
       
    }
   
    
    }## end of the model
    
    ",fill = TRUE)
sink()


## bundle data
# run
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                          "Non-mammalian Mammaliaformes",
                                                          "Mammalia" ), function (i) {
                                                            
                                                            
                                                            
                          # select clade                                                             
                          cond_det <- cond_det %>% 
                            filter (cladeID == i) %>%
                            mutate (stageMod = stage-(min(stage)-1),
                                    genus = as.numeric(as.factor(genus)))
                          
                          # subset time covariates
                          time_covariates <- time_covariates [which(rownames(time_covariates) %in% unique (cond_det$stage)),]
                          
                          # bundle data 
                          str(jags.data <- list(
                            
                                  y = cond_det$detection,
                                  ngen = length(unique(cond_det$genusID)),
                                  nint = length(unique(cond_det$stageMod)),
                                  nobs = nrow(cond_det),
                                  genus = cond_det$genus,
                                  stage = cond_det$stageMod,
                                  nform = cond_det$formations,
                                  # covariates
                                  # time = vectorized_occ$
                                  range = (cond_det$range_size-mean(cond_det$range_size))/sd(cond_det$range_size),
                                  lat_obs = (cond_det$latitude-mean(cond_det$latitude))/sd(cond_det$latitude),
                                  time = (cond_det$time-mean(cond_det$time))/sd(cond_det$time),
                                  temp_obs = (cond_det$temperature-mean(cond_det$temperature))/sd(cond_det$temperature),
                                  # time covs 
                                  # time covs 
                                  temperature = as.vector(scale(time_covariates$temperature)),
                                  precipitation = as.vector(scale(time_covariates$precipitation)),
                                  area = as.vector(scale(time_covariates$area))
                                  
                                                      
                                )
                                )
                                
                          
                          # Set initial values
                          # function inits (to be placed in each chain)
                          inits <- function(){ list(z = tapply (jags.data$y,
                                                                list(jags.data$genus,
                                                                     jags.data$stage),
                                                                max),
                                                    
                                                    gamma.u = runif (1),
                                                    phi.u = runif (1),
                                                    p.u = runif (1),
                                                    beta.gamma.prec = rnorm (1),
                                                    beta.gamma.temp = rnorm (1),
                                                    beta.gamma.area = rnorm (1),
                                                    beta.phi.prec = rnorm (1),
                                                    beta.phi.temp = rnorm (1),
                                                    beta.phi.area = rnorm (1),
                                                    beta.p.time = rnorm (1),
                                                    beta.p.range = rnorm (1),
                                                    beta.p.lat = rnorm (1),
                                                    beta.p.temp = rnorm (1),
                                                    psi1 = runif (1)
                                                    
                                                    )}
                                
                          ## Parameters to monitor
                          ## long form
                          ## long form
                          params <- c(
                            "gamma.u",
                            "phi.u",
                            "intercept.gamma",
                            "beta.gamma.prec",
                            "beta.gamma.temp",
                            "beta.gamma.area",
                            "intercept.phi",
                            "beta.phi.prec",
                            "beta.phi.temp",
                            "beta.phi.area",
                            "intercept.p",
                            "beta.p.time",
                            "beta.p.range",
                            "beta.p.lat",
                            "beta.p.temp",
                            "gamma",
                            "phi",
                            "p",
                            "SRexp",
                            "FSS",
                            "CH",
                            "CH1",
                            "psi1",
                            "R0",
                            "RER",
                            "turnover",
                            "psi.eq"
                            
                          )
                          
                      
          
              # MCMC runs
              samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                                                                  modelFile="global_CMRmodel_vectorized_covariates.txt",
                                                                  chains=nc, 
                                                                  sample2save=((ni-nb)/nt), 
                                                                  nSaves=6, 
                                                                  burnin=nb, 
                                                                  thin=nt,
                                                                  fileStub=paste ("output/res_det_conditional/CMR_global_binomial", i,sep="_"))

})

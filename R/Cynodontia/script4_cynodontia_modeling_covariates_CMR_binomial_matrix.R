# --------------------------------------

# 2 . modeling


## Global MCMC settings
######################

# number of interactions for a global-scale analyses
na <- 15000; nb <- 18000; ni <- 20000; nc <- 3; nt <- 5

## short form for tests
#na <- 400; nb <- 700; ni <- 1000; nc <- 3; nt <- 1
na <- 40; nb <- 70; ni <- 100; nc <- 3; nt <- 1

# load packages
source("R/packages.R")

# load data
load(here ("processed_data","CMR_data_observation_cov.RData"))

# load environmental data
load (here ("processed_data", "site_covs.RData"))


# ------------------------------
# model with covariates at time level (global scale analyses)

#scaled_temperature <- scale (time_covariates$temperature,T,T)
#scaled_temperature_sd <- scale (time_covariates$temperature_sd,T,T)
#scaled_precipitation <-  scale (time_covariates$precipitation,T,T)
#scaled_precipitation_sd <-  scale (time_covariates$precipitation_sd,T,T)
#
## time ( preservation )
#scaled_time_stage <- scale (bins[time_bins,"mid_ma"],T,T)
#
## genus spatial range
#range_taxon_analysis <- (range_area_taxon[which(range_area_taxon$taxon %in% cynodontia_data),])
#scaled_range <- scale(range_taxon_analysis$range_area,T,T)

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
      
       # Site Occupancy Dynamics (Priors)
       
       # intercepts
       # gamma
       for (t in 1:nint) {
       
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
            y [g,t] ~ dbin(z[g,t]*p[t], nform[t])
                  
                  
                }
        }
      
    # derived parameters
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
    
        
   # descriptors of change
   for (t in 2:nint) {  
     
     
       FSS[t] <- sum (muZ[,t])  # finite sample size
       CH[t] <- (sum(z[,t]) - sum(z[,t-1])) ### turnover (proportional gain or loss)
       CH1[t]<-sum(z[,t-1])
       RER[t] <- (1-phi[t])/gamma[t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t] <- (1-phi[t])-gamma[t]  ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
       
    }
    
    
    }## end of the model
    
    

    
    
    ",fill = TRUE)
sink()



# run model
no_cov_samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                    "Non-mammalian Mammaliaformes",
                                                    "Mammalia" ), function (i) {
        
        
        str(jags.data <- list(
          y = data.matrix(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)]),
          ngen = dim(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)])[1],
          nint = dim(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)])[2],
          nform = formations_per_interval$formations_per_interval
          
        )
        )
        
        # Set initial values
        # function inits (to be placed in each chain)
        inits <- function(){ list(z = data.matrix(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)]),
                                  
                                  psi1 = runif (1)
                                  
                                  
        )}
        
        
        ## Parameters to monitor
        params <- c("psi1","gamma", "phi","p",
          "FSS","SRexp","RER", "R0",
          "CH","CH1"
          
        )
        
        # MCMC runs
        # models
        samples_paleo_cynodontia_no_covariate <- jags (data = jags.data, 
                                          parameters.to.save = params, 
                                          model.file = "global_CMRmodel_matrix_no_covariates.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.thin = nt, 
                                          n.iter = ni, 
                                          n.burnin = nb, 
                                          DIC = T,  
                                          n.cores=nc,
                                          parallel=T
        )

                                                    
        })


# check

plot(no_cov_samples_paleo_cynodontia_binomial_run[[3]]$mean$SRexp,col="red",type="b")
lines(no_cov_samples_paleo_cynodontia_binomial_run[[2]]$mean$SRexp, col = "blue")
lines(no_cov_samples_paleo_cynodontia_binomial_run[[1]]$mean$SRexp, col = "green")

# save
save (no_cov_samples_paleo_cynodontia_binomial_run,
      file = here ("output","no_covariates_CMR_global_binomial.RData"))


# -------------------------------------------------------

# covariates per interval using binomial distribution in P


# model with covariates

sink("global_CMRmodel_matrix_covariates.txt")
cat("
   
     model {
    
    
      #############################################################
      #                                                           #
      #                  Biological process                       #
      #                                                           #
      #############################################################
      
       # Site Occupancy Dynamics (Priors)
       
             # intercepts
             # gamma
             intercept.gamma ~ dnorm(0, 0.001)
             # phi
             intercept.phi ~ dnorm(0, 0.001)

            # ----------------------
            #     Gamma (origination)
            # ----------------------
            # regression coefficients
            beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
            beta.gamma.temp ~ dnorm(0, 0.001)# temperature
            
            # ----------------------
            #     Phi (persistence)
            # ----------------------
            beta.phi.prec ~ dnorm(0, 0.001)# precipitation
            beta.phi.temp ~ dnorm(0, 0.001)# temperature
            
            
            ## set initial conditions for each genus
            psi1 ~ dunif(0,1)
           
       
       
   ############      Model       #############
   
   # model for phi and gamma
   ### model dynamic parameters
    
    for (t in 2:nint){
  
         # speciation
         logit(gamma[t]) <-  intercept.gamma + 
                               beta.gamma.prec*precipitation[t]+
                               beta.gamma.temp*temperature[t]
                                     
          # persistence
          logit(phi[t]) <-  intercept.phi + 
                              beta.phi.prec*precipitation[t]+
                              beta.phi.temp*temperature[t]
                                              
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
        
    intercept.p ~ dnorm(0,0.001)
    beta.p.lat ~ dnorm(0,0.001)
    beta.p.range ~ dnorm(0,0.001)
    beta.p.time ~ dnorm(0,0.001)
    beta.p.temp ~ dnorm(0,0.001)
     
    # observation submodel
    for (g in 1:ngen) { ## loop over observations 
      
      for (t in 1:nint) { ## loop over observations 
    
          # p independent on the number of formations
          logit(p[g,t])<- intercept.p+
                          beta.p.lat*lat_obs[g]+
                          beta.p.range*range[g]+
                          beta.p.time*time[t]+
                          beta.p.temp*temp_obs[g]
            
            
          # observation
          # Specify the binomial observation model conditional on occupancy state
          y[g,t] ~ dbin(z[g,t]*p[g,t], nform[t])
                  
        }
      
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
                                                            
                                                            
                                                                            
                str(jags.data <- list(
                        y = data.matrix(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)]),
                        ngen = dim(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)])[1],
                        nint = dim(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)])[2],
                        nform = formations_per_interval$formations_per_interval,
                        # covariates
                        # time = vectorized_occ$
                        range = as.vector(scale(range_area_taxon[which(range_area_taxon$taxon %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"range_area"])),
                        #range_area_taxon[which(range_area_taxon$taxon %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"taxon"] == array_genus_bin [which(array_genus_bin$clade %in% i), 2]
                        lat_obs = as.vector(scale(observation_covariates[which(observation_covariates$unique_name %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"latitude"])),
                        time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)])),
                        temp_obs = as.vector(scale(observation_covariates[which(observation_covariates$unique_name %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"temperature"])),
                        # time covs 
                        temperature = as.vector(scale(time_covariates$temperature)),
                        precipitation = as.vector(scale(time_covariates$precipitation))
                        
                                            
                      )
                      )
                      
                
                
                # Set initial values
                # function inits (to be placed in each chain)
                inits <- function(){ list(z = data.matrix(array_genus_bin [which(array_genus_bin$clade %in% i), -c(1,2)]),
                                          
                                          intercept.gamma = rnorm (1),
                                          beta.gamma.prec = rnorm (1),
                                          beta.gamma.temp = rnorm (1),
                                          intercept.phi = rnorm (1),
                                          beta.phi.prec = rnorm (1),
                                          beta.phi.temp = rnorm (1),
                                          intercept.p = rnorm (1),
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
                  "intercept.gamma",
                  "beta.gamma.prec",
                  "beta.gamma.temp",
                  "intercept.phi",
                  "beta.phi.prec",
                  "beta.phi.temp",
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
                  "RER"
                )
                
                # usual jags to test
                samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                         modelFile="global_CMRmodel_matrix_covariates.txt",
                         chains=nc, 
                         sample2save=((ni-nb)/nt), 
                         nSaves=5, 
                         burnin=nb, 
                         thin=nt,
                         fileStub=paste ("output/CMR_global_binomial", i,sep="_"))
})


## MCMC runs
#samples_paleo_cynodontia_binomial <-jags (data = jags.data, 
#                                          parameters.to.save = params, 
#                                          model.file = "global_CMRmodel_matrix_covariates.txt", 
#                                          inits = inits, 
#                                          n.chains = nc, 
#                                          n.thin = nt, 
#                                          n.iter = ni, 
#                                          n.burnin = nb, 
#                                          DIC = T,  
#                                          n.cores=nc,
#                                          parallel=T
#)



#--------------------------------------------------------------------- #

# - binomial model with covariates at the site/region level

# try nimble to speed up the process
#library (nimble)

#Section6p4_code <- nimbleCode( {

sink("regional_CMRmodel_matrix_covariates.txt")
cat("
   
     
     model {
    
    
      #############################################################
      #                                                           #
      #                  Biological process                       #
      #                                                           #
      #############################################################
      
       # Site Occupancy Dynamics (Priors)
       # one parameter per group
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        # intercept
        intercept.gamma ~ dnorm(0, 0.001)
        # regression coefficients
        beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
        beta.gamma.temp ~ dnorm(0, 0.001)# temperature
        beta.gamma.lat ~ dnorm(0,0.001) # latitude
        
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercept
        intercept.phi ~ dnorm(0, 0.001)
        # regression coeff 
        beta.phi.prec ~ dnorm(0, 0.001)# precipitation
        beta.phi.temp ~ dnorm(0, 0.001)# temperature
        beta.phi.lat ~ dnorm(0,0.001)
       
        ## set initial conditions for each genus
        psi1 ~ dunif(0,1)
        
      # initial state
      for (g in 1:ngen) {
        for (i in 1:nsites) {
        
           z[g,1,i]~dbern(psi1) # occupancy status initialization
                
        }
      }
    
    ############      Model       #############
   
   # occupancy dynamics
    
     for (i in 1:nsites) {
       for (t in 2:nint){
    
            # origination
            logit(gamma[t,i]) <-  intercept.gamma + 
                                    beta.gamma.prec*precipitation[t,i]+
                                    beta.gamma.temp*temperature[t,i]+
                                    beta.gamma.lat*lat[i]
            # persistence
            logit(phi[t,i]) <-    intercept.phi + 
                                    beta.phi.prec*precipitation[t,i]+
                                    beta.phi.temp*temperature[t,i]+
                                    beta.phi.lat*lat[i]
    
            # dynamics across genus
              
            for (g in 1:ngen) {
                    
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[g,t,i] <- z[g,t-1,i] *  phi[t,i] + ### if occupied, p of not getting extinct/persist in the next time
                              (1-z[g,t-1,i]) *  gamma[t,i] ###  if not occupied, p of originate in the next time
                
               # realized occurrence
        		   z[g,t,i] ~ dbern(muZ[g,t,i])
            
             }#gen
            } #stages
          }#i sites
      
    
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    # Priors for detection probability
   
   
    ###  detection intercept
    intercept.p ~ dnorm(0,0.001)
    beta.p.lat ~ dnorm(0,0.001)
    beta.p.range ~ dnorm(0,0.001)
    beta.p.time ~ dnorm(0,0.001)
    beta.p.temp ~ dnorm(0,0.001)
        
     
    # observation submodel
    for (g in 1:ngen) { ## loop over observations 
      
      for (t in 1:nint) { ## loop over observations 
    
          # p independent on the number of formations
          logit(p[g,t])<- intercept.p+
                          beta.p.lat*lat_obs[g]+
                          beta.p.range*range[g]+
                          beta.p.time*time[t]+
                          beta.p.temp*temp_obs[g]
            
          for (i in 1:nsites) {  
      
              # observation process
              # Specify the binomial observation model conditional on occupancy state
              y[g,t,i] ~ dbin(z[g,t,i]*p[g,t], nform[t,i])
                      
        }
      
      }
    
    }
    
    # derived parameters
    # true richness
    for (t in 1:nint){
       for (i in 1:nsites) {
       
         SRexp[t,i] <- sum (z[,t,i]) # true richness estimate
         
    
    }
   }
 
    
 #  
 # # change
 
 for (t in 2:nint) {  
      
        for (i in 1:nsites) {  
         
           FSS[t,i] <- sum (muZ[,t,i]) # finite sample size 
           CH[t,i] <- (sum(z[,t,i]) - sum(z[,t-1,i])) ### turnover (proportional gain or loss)
           CH1[t,i]<-sum(z[,t-1,i])
           RER[t,i] <- (1-phi[t,i])/gamma[t,i] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
           R0[t,i] <- (1-phi[t,i])-gamma[t,i] ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
        
  
     }
  
    }
  
  
    
        
    } ## end of the model
    
    ",fill = TRUE)
sink()

#  ----------------------------------------------------

# bundle data
# replace NAs by a very small number
region_temperature[is.na(region_temperature)] <- -999
region_precipitation[is.na(region_precipitation)] <- -999


# run
samples_paleo_cynodontia_sites_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                          "Non-mammalian Mammaliaformes",
                                                          "Mammalia" ), function (i) {
                        
                        
                        
                        str(jags.data <- list(y = (array_genus_bin_site [which(clades %in% i),,]),
                                              ngen = dim(array_genus_bin_site [which(clades %in% i),,])[1],
                                              nint = dim(array_genus_bin_site [which(clades %in% i),,])[2],
                                              nsites = dim(array_genus_bin_site [which(clades %in% i),,])[3],
                                              nform = as.matrix(t(formations_per_site_interval)),
                                              # covariates
                                              # time = vectorized_occ$
                                              range = as.vector(scale(range_area_taxon[which(range_area_taxon$taxon %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"range_area"])),
                                              #range_area_taxon[which(range_area_taxon$taxon %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"taxon"] == array_genus_bin [which(array_genus_bin$clade %in% i), 2]
                                              lat_obs = as.vector(scale(observation_covariates[which(observation_covariates$unique_name %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"latitude"])),
                                              time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)])),
                                              temp_obs = as.vector(scale(observation_covariates[which(observation_covariates$unique_name %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"temperature"])),
                                              # lat = site_covs
                                              temperature = (region_temperature-mean (region_temperature,na.rm=T))/sd(region_temperature,na.rm=T),
                                              precipitation = (region_precipitation-mean (region_precipitation,na.rm=T))/sd(region_precipitation,na.rm=T),
                                              lat = as.vector(scale (region_latitude))
                        )
                        )
                        
                        
                        # Set initial values
                        z_inits<-array_genus_bin_site [which(clades %in% i),,]
                        
                        #z_inits<-ifelse (is.na(z_inits),1,z_inits)
                        # function inits (to be placed in each chain)
                        inits <- function(){ list(z=z_inits,
                                                  intercept.gamma = rnorm (1),
                                                  beta.gamma.prec = rnorm (1),
                                                  beta.gamma.temp = rnorm (1),
                                                  intercept.phi = rnorm (1),
                                                  beta.phi.prec = rnorm (1),
                                                  beta.phi.temp = rnorm (1),
                                                  intercept.p = rnorm (1),
                                                  beta.p.time = rnorm (1),
                                                  beta.p.range = rnorm (1),
                                                  beta.p.lat = rnorm (1),
                                                  beta.p.temp = rnorm (1),
                                                  psi1 = runif (1)
                        )}
                        
                        
                        ## Parameters to monitor
                        ## long form
                        params <- c(
                          
                          "intercept.gamma",
                          "beta.gamma.prec",
                          "beta.gamma.temp",
                          "beta.gamma.lat",
                          "intercept.phi",
                          "beta.phi.prec",
                          "beta.phi.temp",
                          "beta.phi.lat",
                          "intercept.p",
                          "beta.p.time",
                          "beta.p.range",
                          "beta.p.lat",
                          "beta.p.temp",
                          "p",
                          "SRexp",
                          "FSS",
                          "CH",
                          "CH1",
                          "psi1",
                          "gamma",
                          "phi"
                          
                          
                        )
                        
                        
                        # MCMC runs
                        # run save jags
                        samples_paleo_cynodontia_sites_binomial <- saveJAGS(jags.data, inits, params, 
                                                                            modelFile="regional_CMRmodel_matrix_covariates.txt",
                                                                            chains=nc, 
                                                                            sample2save=((ni-nb)/nt), 
                                                                            nSaves=5, 
                                                                            burnin=nb, 
                                                                            thin=nt,
                                                                            fileStub=paste ("output/CMR_regional_binomial", i,sep="_"))
                        
                        
})


# usual jags to test
samples_paleo_cynodontia_sites_binomial <- jags (data = jags.data, 
                                                 parameters.to.save = params, 
                                                 model.file = "regional_CMRmodel_matrix_covariates.txt", 
                                                 inits = inits, 
                                                 n.chains = nc, 
                                                 n.thin = nt, 
                                                 n.iter = ni, 
                                                 n.burnin = nb, 
                                                 DIC = T,  
                                                 n.cores=nc,
                                                 parallel=F
)

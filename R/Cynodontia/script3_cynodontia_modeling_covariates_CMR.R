# --------------------------------------

# 3 . modeling


## Global MCMC settings
######################

na <- 8000; nb <-10000; ni <- 20000; nc <- 3; nt <- 25

## short form
#na <- 25; nb <- 30; ni <- 50; nc <- 2; nt <- 1

# load packages
source("R/packages.R")

# load data
load(here ("processed_data","CMR_data.RData"))

# load environmental data
load (here ("processed_data", "site_covs.RData"))


# ------------------------------------------------------------

# model for CMR data style
# model with no covariates

sink("dyn_model_vectorized_no_covariates_CMR.txt")
cat("
   
    model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    
   for (t in 1:nint) {
    
      gamma [t]~dunif(0,1)
      phi [t]~dunif(0,1)
    
    
    }
   
   ## set initial conditions for each genus
   for (g in 1:ngen) {
      psi1[g] ~ dunif(0,1)#dbeta(a,b)
   }
   
    # Specify the hyperparameters for the Beta distribution
    #a <- 10   # probability of success
    #b <- 90   # probability of failures
    
   ############      Model       #############
    
    for (g in 1:ngen) {
    
          z[g,1]~dbern(psi1[g]) # occupancy status initialization
    
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
    # constant detection probability
    
    ###  detection intercept
    intercept.p ~ dnorm(0,0.001)
    beta.p.form ~ dnorm(0,0.001)   
    
    # observation submodel
    for (t in 1:nint){
        
          # model
          # p dependent on the number of formations
          logit(p[t])<- intercept.p+beta.p.form*nform[t]
          
    
            for (g in 1:ngen) {
            
                  # observation
                  y [g,t] ~ dbern(muY[g,t])
                  # Specify the observation model conditional on occupancy state
                  muY [g,t] <- z[g,t] * p[t]
                  
                  
                }
        
    }
      
      
    # derived parameters
    for (t in 1:nint){
      SRexp[t] <- sum (z[,t])
    }
        
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()



# --------------------------------------


## bundle data
str(jags.data <- list(y = array_genus_bin,
                      nint= ncol(array_genus_bin), 
                      ngen = nrow(array_genus_bin),
                      nform = as.vector(unname(scale(formations_per_interval)))
                      
                      )
)

# Set initial values
# inits for psi1
# psi1 <- apply (array_genus_bin,1,sum)/nrow(array_genus_bin)

# function inits (to be placed in each chain)
inits <- function(){ list(z = array_genus_bin)}

## Parameters to monitor
## long form
params <- c(
  
  "gamma", "phi","p","z","beta.p.form", "intercept.p","SRexp"
)

## MCMC settings
######################
## short form
na <- 6000; nb <-10000; ni <- 20000; nc <- 3; nt <- 50
#na <- 250; nb <- 500; ni <- 1000; nc <- 3; nt <- 1 # test

# MCMC runs
# models

samples_paleo_cynodontia <- jags (data = jags.data, 
                       parameters.to.save = params, 
                       model.file = "dyn_model_vectorized_no_covariates_CMR.txt", 
                       inits = inits, 
                       n.chains = nc, 
                       n.thin = nt, 
                       n.iter = ni, 
                       n.burnin = nb, 
                       DIC = T,  
                       n.cores=nc,
                       parallel=T
)

# save
save (samples_paleo_cynodontia,file = here ("output","samples_paleo_cynodontia_no_covariates_CMR.RData"))

# --------------------------------------------------------------------

# model with covariates at time level

sink("dyn_model_vectorized_covariates_CMR.txt")
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
      beta.gamma.elev ~ dnorm(0, 0.001)# elevation
      beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
      beta.gamma.temp ~ dnorm(0, 0.001)# temperature
      
      # ----------------------
      #     Phi (persistence)
      # ----------------------
      
      beta.phi.prec ~ dnorm(0, 0.001)# precipitation
      beta.phi.temp ~ dnorm(0, 0.001)# temperature
      beta.phi.elev ~ dnorm(0, 0.001)# elevation
      
      
     ## set initial conditions for each genus
     psi1 ~ dunif(0,1)
     
    
    ############      Model       #############
   
   # model for phi and gamma
   ### model dynamic parameters
  
    for (t in 2:nint){
  
         # speciation
        logit(gamma[t]) <-  intercept.gamma + 
                                     beta.gamma.elev*elevation[t]+
                                     beta.gamma.prec*precipitation[t]+
                                     beta.gamma.temp*temperature[t]
          # persistence
          logit(phi[t]) <-    intercept.phi + 
                                     beta.phi.elev*elevation[t]+
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
    
    # Priors for detection probability
    # constant detection probability
    
    ###  detection intercept
    intercept.p ~ dnorm(0,0.001)
    beta.p.form ~ dnorm(0,0.001)   
    beta.p.time ~ dnorm(0,0.001)   
    beta.p.range ~ dnorm(0,0.001)   
    
    # observation submodel
    for (g in 1:ngen) {
      for (t in 1:nint){
      
          # model
          # p dependent on the number of formations
          logit(p[g,t])<- intercept.p+beta.p.form*nform[t]+beta.p.time*time[t]+beta.p.range*range[g]
          
          # observation
          y [g,t] ~ dbern(muY[g,t])
          # Specify the observation model conditional on occupancy state
          muY [g,t] <- z[g,t] * p[g,t]
                  
                  
                }
        }
      
    # derived parameters
    # true richness
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
   # descriptors of change
   for (t in 2:nint) {  
     
       propcH[t] <-(sum(z[,t])-sum (z[,t-1]))/sum(z[,t-1]) ### turnover (proportional gain or loss)
       RER[t] <- (1-phi[t])/gamma[t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t] <- (1-phi[t])-gamma[t]  ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
      
   }
    
    
        
    }## end of the model
    
    ",fill = TRUE)
sink()


# ------------------------------
# model with covariates at time level

# bins in data
elevation <- site_covs$elevation [which(names(site_covs$elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
temperature <- site_covs$temp [which(names(site_covs$temp) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
precipitation <- site_covs$prec [which(names(site_covs$prec) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]

# cells among the studied sites
elevation<-elevation[which(rownames (elevation) %in% cells),] 
temperature<-temperature[which(rownames (temperature) %in% cells),] 
precipitation<-precipitation[which(rownames (precipitation) %in% cells),] 

# scale variables
# for elevation, we will consider only cells above and at the sea level in a given stage
list_elevation <- as.vector(elevation)
mean_elevation_stage <- lapply (list_elevation, function (i)
  

  mean (i[which(i >=0)]) # mean elevation for cells above and at the sea level

)
#unlist
mean_elevation_stage <- unlist(mean_elevation_stage)
scaled_elevation <- scale (mean_elevation_stage)
# temperature, precipitation
scaled_temperature <- (scale (colMeans(temperature ),T,T))
scaled_precipitation <-  (scale (colMeans(precipitation),T,T))

# reverse the order of columns
scaled_elevation<-scaled_elevation[rev(rownames(scaled_elevation)),]
scaled_temperature<-scaled_temperature[rev(rownames(scaled_temperature)),]
scaled_precipitation<-scaled_precipitation[rev(rownames(scaled_precipitation)),]

# time ( preservation )
scaled_time_stage <- scale (bins[time_bins,"mid_ma"],T,T)

# genus spatial range
range_taxon_analysis <- (range_area_taxon[which(range_area_taxon$taxon %in% cynodontia_data),])
scaled_range <- scale(range_taxon_analysis$range_area,T,T)

## bundle data
# bundle one dataset per site
### run in parallel processing
str(jags.data <- list(y = array_genus_bin,
                      nint= ncol(array_genus_bin), 
                      ngen = nrow(array_genus_bin),
                      nform = as.vector(unname(scale(formations_per_interval))), 
                      elevation =  as.vector (scaled_elevation),
                      temperature = as.vector(scaled_temperature),
                      precipitation = as.vector(scaled_precipitation),
                      time = as.vector(scaled_time_stage),
                      range = as.vector(scaled_range)
                      
                      
)
)

# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = array_genus_bin)}

## Parameters to monitor
## long form
params <- c(
  "intercept.gamma",
  "beta.gamma.elev",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "beta.phi.elev",
  "intercept.p",
  "beta.p.form",
  "beta.p.time",
  "beta.p.range",
  "gamma",
  "phi",
  "p",
  "SRexp",
  "propcH",
  "R0",
  "RER"
)

# MCMC runs
# models
samples_paleo_cynodontia_stage <- jags (data = jags.data, 
                                  parameters.to.save = params, 
                                  model.file = "dyn_model_vectorized_covariates_CMR.txt", 
                                  inits = inits, 
                                  n.chains = nc, 
                                  n.thin = nt, 
                                  n.iter = ni, 
                                  n.burnin = nb, 
                                  DIC = T,  
                                  n.cores=nc,
                                  parallel=T
)

save (samples_paleo_cynodontia_stage,file = here ("output","samples_paleo_cynodontia_covariates_CMR_stage.RData"))

#--------------------------------------------------------------------- #

# - model with covariates at the site level


sink("dyn_model_vectorized_covariates_CMR_sites.txt")
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
      beta.gamma.elev ~ dnorm(0, 0.001)# elevation
      beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
      beta.gamma.temp ~ dnorm(0, 0.001)# temperature
      beta.gamma.lat ~ dnorm(0,0.001) # latitude
      
      # ----------------------
      #     Phi (persistence)
      # ----------------------
      
      beta.phi.prec ~ dnorm(0, 0.001)# precipitation
      beta.phi.temp ~ dnorm(0, 0.001)# temperature
      beta.phi.elev ~ dnorm(0, 0.001)# elevation
      beta.phi.lat ~ dnorm(0,0.001)
    
      ## set initial conditions for each genus
      for (i in 1:nsites) {
        psi1[i]~dunif(0,1)
      }
     
     # initial state
     
     for (g in 1:ngen) {
     
        for (i in 1:nsites) {
        
        z[g,1,i]~dbern(psi1[i]) # occupancy status initialization
                
        }
        
    }
    
    ############      Model       #############
   
   # occupancy dynamics
    
    
      for (i in 1:nsites) {
     
          for (t in 2:nint){
    
            # origination
            logit(gamma[t,i]) <-  intercept.gamma + 
                              beta.gamma.elev*elevation[t,i]+
                              beta.gamma.prec*precipitation[t,i]+
                              beta.gamma.temp*temperature[t,i]+
                              beta.gamma.lat*lat[i]
            # persistence
            logit(phi[t,i]) <-    intercept.phi + 
                              beta.phi.elev*elevation[t,i]+
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
            
        }#t
      
      
      } #g
      
  
    }#i
    
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    # Priors for detection probability
    # constant detection probability
    
    ###  detection intercept
    intercept.p ~ dnorm(0,0.001)
    beta.p.form ~ dnorm(0,0.001)   
    beta.p.lat ~ dnorm(0,0.001)
    beta.p.range ~ dnorm(0,0.001)
    beta.p.time ~ dnorm(0,0.001)
    
    # observation submodel
    for (g in 1:ngen) {
        for (t in 1:nint){
          for (i in 1:nsites){
    
            # model
            # p dependent on the number of formations
            logit(p[g,t,i])<- intercept.p+beta.p.form*nform[t,i]+
                                        beta.p.lat*lat[i]+
                                        beta.p.range*range[g]+
                                        beta.p.time*time[t]
            
      
            # observation
            y [g,t,i] ~ dbern(muY[g,t,i])
            # Specify the observation model conditional on occupancy state
            muY [g,t,i] <- z[g,t,i] * p[g,t,i]
                    
                  
            }
        }
    }
    
    
    # derived parameters
    # true richness
    for (t in 1:nint){
      for (i in 1:nsites) {
        SRexp[t,i] <- sum (z[,t,i])
      }
    }
    
    
    
   # change
   for (t in 2:nint) {  
    for (i in 1:nsites) {  
       
       #propcH[t,i] <-(sum(z[,t,i])-sum(z[,t-1,i]))/sum(z[,t-1,i]) ### turnover (proportional gain or loss)
       RER[t,i] <- (1-phi[t,i])/gamma[t,i] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t,i] <- (1-phi[t,i])-gamma[t,i] ## net diversification rate (r= μ - λ; Rabosky 2018) of each tim
     
      }
   }
    
        
    }## end of the model
    
    ",fill = TRUE)
sink()


# --------------------

# bins in data
elevation <- site_covs$elevation [which(names(site_covs$elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
temperature <- site_covs$temp [which(names(site_covs$temp) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
precipitation <- site_covs$prec [which(names(site_covs$prec) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
lat <- ( (site_covs$paleolat))

# cells among the studied sites
elevation<-elevation[which(rownames (elevation) %in% cells),] 
temperature<-temperature[which(rownames (temperature) %in% cells),] 
precipitation<-precipitation[which(rownames (precipitation) %in% cells),] 
lat<-lat[which(rownames (precipitation) %in% cells)] 

# scale variables
scaled_elevation <- (elevation - mean(as.matrix(elevation))) / sd(as.matrix(elevation))
scaled_temperature <- (temperature - mean(as.matrix(temperature))) / sd(as.matrix(temperature))
scaled_precipitation <- (precipitation - mean(as.matrix(precipitation))) / sd(as.matrix(precipitation))
scaled_lat <- (lat - mean(lat)) / sd(lat)

# reverse the order of columns
scaled_elevation<-scaled_elevation[,rev(colnames(scaled_elevation))]
scaled_temperature<-scaled_temperature[,rev(colnames(scaled_temperature))]
scaled_precipitation<-scaled_precipitation[,rev(colnames(scaled_precipitation))]

# time ( preservation )
scaled_time_stage <- scale (bins[time_bins,"mid_ma"])

# genus spatial range
range_taxon_analysis <- (range_area_taxon[which(range_area_taxon$taxon %in% cynodontia_data),])
scaled_range <- scale(range_taxon_analysis$range_area)

# formations per intervals and site
scaled_formations <- (formations_per_site_interval - mean(as.matrix(formations_per_site_interval))) / sd(as.matrix(formations_per_site_interval))

# sel sites with 5
sel_sites <- which(colSums(apply(array_genus_bin_site, c(1,3), sum))>=5)


## bundle data
# bundle one dataset per site
str(jags.data <- list(y = array_genus_bin_site[,,sel_sites],
                          ngen = dim(array_genus_bin_site[,,sel_sites])[1],
                          nint= dim(array_genus_bin_site[,,sel_sites])[2], 
                          nsites=dim(array_genus_bin_site[,,sel_sites])[3],
                          elevation = t((scaled_elevation))[,sel_sites],
                          temperature = t((scaled_temperature))[,sel_sites],
                          precipitation = t((scaled_precipitation))[,sel_sites],
                          nform = t(scaled_formations)[,sel_sites],
                          range = scaled_range[,1],
                          time = scaled_time_stage[,1],
                          lat = scaled_lat[sel_sites]
                          
    )
    )
    
# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = array_genus_bin_site[,,sel_sites])}

## Parameters to monitor
## long form
params <- c(
 "intercept.gamma",
 "beta.gamma.elev",
 "beta.gamma.prec",
 "beta.gamma.temp",
 "beta.gamma.lat",
 "intercept.phi",
 "beta.phi.prec",
 "beta.phi.temp",
 "beta.phi.elev",
 "beta.phi.lat",
 "intercept.p",
 "beta.p.form",
 "beta.p.range",
 "beta.p.time",
 "beta.p.lat",
 "gamma",
 "phi",
 "p",
 "SRexp",
 #"z",
 "R0",
 "RER"
)

# MCMC runs
# models
samples_cynodontia_CMR_sites <- jags (data = jags.data, 
                                  parameters.to.save = params, 
                                  model.file = "dyn_model_vectorized_covariates_CMR_sites.txt", 
                                  inits = inits, 
                                  n.chains = nc, 
                                  n.thin = nt, 
                                  n.iter = ni, 
                                  n.burnin = nb, 
                                  DIC = T,  
                                  n.cores=nc,
                                  parallel=T
)
# save
#save (samples_cynodontia_CMR_sites,file = here ("output","samples_cynodontia_CMR_sites.RData"))

require(saveJAGS)
samples_cynodontia_CMR_sites <- saveJAGS(jags.data, inits, params, 
         modelFile="dyn_model_vectorized_covariates_CMR_sites.txt",
         chains=nc, 
         sample2save=((ni-nb)/nt), 
         nSaves=3, 
         burnin=nb, 
         thin=nt,
         fileStub="output/CMR_sites_bernoulli")


# if need to run longer chains
# newRes_samples_cynodontia_CMR_sites <- resumeJAGS(fileStub="output/CMR_sites_bernoulli", nSaves=1)
# combine saves
# newRes_samples_cynodontia_CMR_sites <- combineSaves(recoverSaves("output/CMR_sites_bernoulli"))

# ------------------------------------------------------------ #
# interesting parameters
# temporal trends
# MUZ_ALLB <- combineSaves(res_ALLB, params = "Nweekinf",thin = 10)
# MUZ_ALLB_TEMPRAIN <- combineSaves(res_ALLB_TEMPRAIN, params = "Nweekinf", thin = 10)
# MUZ_ALLB_NoSpace <- combineSaves(res_ALLB_NoSPace, params = "Nweekinf", thin = 10)

# -------------------------------------------------------

# covariates per interval using binomial distribution in P


# model with covariates

sink("dyn_model_vectorized_covariates_CMR_binomial.txt")
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
      beta.gamma.elev ~ dnorm(0, 0.001)# elevation
      beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
      beta.gamma.temp ~ dnorm(0, 0.001)# temperature
      
      # ----------------------
      #     Phi (persistence)
      # ----------------------
      
      
      beta.phi.prec ~ dnorm(0, 0.001)# precipitation
      beta.phi.temp ~ dnorm(0, 0.001)# temperature
      beta.phi.elev ~ dnorm(0, 0.001)# elevation
      
      
     ## set initial conditions for each genus
     psi1 ~ dunif(0,1)
     
    
    ############      Model       #############
   
   # model for phi and gamma
   ### model dynamic parameters
  
    for (t in 2:nint){
  
         # speciation
         logit(gamma[t]) <-  intercept.gamma + 
                                     beta.gamma.elev*elevation[t]+
                                     beta.gamma.prec*precipitation[t]+
                                     beta.gamma.temp*temperature[t]
          # persistence
          logit(phi[t]) <-    intercept.phi + 
                                     beta.phi.elev*elevation[t]+
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
    
    # Priors for detection probability
    ###  detection intercept and coefficients
    intercept.p ~ dnorm(0,0.001)
    beta.p.time ~ dnorm(0,0.001)   
    beta.p.range ~ dnorm(0,0.001)   
    
    # observation submodel
     for (g in 1:ngen) {
      for (t in 1:nint){
    
            # observation
            # p dependent on the number of formations
            logit(p[g,t])<- intercept.p+beta.p.time*time[t]+beta.p.range*range[g]
    
            # Specify the binomial observation model conditional on occupancy state
            y [g,t] ~ dbin(z[g,t]*p[g,t], nform[t])
                  
                  
                }
        }
      
    # derived parameters
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
    
        
    }## end of the model
    
    ",fill = TRUE)
sink()


## bundle data
# bundle one dataset per site
### run in parallel processing
str(jags.data <- list(y = array_genus_bin,
                      nint= ncol(array_genus_bin), 
                      ngen = nrow(array_genus_bin),
                      nform = as.vector(unname((formations_per_interval))), 
                      elevation =  as.vector (scaled_elevation),
                      temperature = as.vector(scaled_temperature),
                      precipitation = as.vector(scaled_precipitation),
                      time = as.vector(scaled_time_stage),
                      range = as.vector(scaled_range)
                      
)
)

# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = array_genus_bin)}

## Parameters to monitor
## long form
## long form
params <- c(
  "intercept.gamma",
  "beta.gamma.elev",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "beta.phi.elev",
  "intercept.p",
  "beta.p.time",
  "beta.p.range",
  "gamma",
  "phi",
  "p",
  "SRexp",
  "propcH",
  "R0",
  "RER"
)

# MCMC runs
# models
samples_paleo_cynodontia <- jags (data = jags.data, 
                                  parameters.to.save = params, 
                                  model.file = "dyn_model_vectorized_covariates_CMR_binomial.txt", 
                                  inits = inits, 
                                  n.chains = nc, 
                                  n.thin = nt, 
                                  n.iter = ni, 
                                  n.burnin = nb, 
                                  DIC = T,  
                                  n.cores=nc,
                                  parallel=T
)



save (samples_paleo_cynodontia,file = here ("output","samples_paleo_cynodontia_covariates_CMR_stages_binomial.RData"))




#--------------------------------------------------------------------- #

# - binomial model with covariates at the site level

# try nimble to speed up the process
#library (nimble)

#Section6p4_code <- nimbleCode( {
  
sink("dyn_model_vectorized_covariates_CMR_sites_binomial.txt")
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
      beta.gamma.elev ~ dnorm(0, 0.001)# elevation
      beta.gamma.prec ~ dnorm(0, 0.001)# precipitation
      beta.gamma.temp ~ dnorm(0, 0.001)# temperature
      beta.gamma.lat ~ dnorm(0,0.001) # latitude
      
      # ----------------------
      #     Phi (persistence)
      # ----------------------
      
      beta.phi.prec ~ dnorm(0, 0.001)# precipitation
      beta.phi.temp ~ dnorm(0, 0.001)# temperature
      beta.phi.elev ~ dnorm(0, 0.001)# elevation
      beta.phi.lat ~ dnorm(0,0.001)
    
      ## set initial conditions for each genus
      for (i in 1:nsites) {
        psi1[i]~dunif(0,1)
      }
     
     # initial state
     
     for (g in 1:ngen) {
     
        for (i in 1:nsites) {
        
        z[g,1,i]~dbern(psi1[i]) # occupancy status initialization
                
        }
        
    }
    
    ############      Model       #############
   
   # occupancy dynamics
    
    
      for (i in 1:nsites) {
     
          for (t in 2:nint){
    
            # origination
            logit(gamma[t,i]) <-  intercept.gamma + 
                              beta.gamma.elev*elevation[t,i]+
                              beta.gamma.prec*precipitation[t,i]+
                              beta.gamma.temp*temperature[t,i]+
                              beta.gamma.lat*lat[i]
            # persistence
            logit(phi[t,i]) <-    intercept.phi + 
                              beta.phi.elev*elevation[t,i]+
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
            
        }#t
      
      
      } #g
      
  
    }#i
    
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
    
     
    # observation submodel
    for (g in 1:ngen) {
    
      for (i in 1:nsites) {
    
        for (t in 1:nint){
    
            
          #prior 
          # model
          # p dependent on the number of formations
          logit(p[g,t,i])<- intercept.p+beta.p.lat*lat[i]+beta.p.time*time[t]+beta.p.range*range[g]
          
          
                  # observation
                  # Specify the binomial observation model conditional on occupancy state
                  y [g,t,i] ~ dbin(z[g,t,i]*p[g,t,i], nform[t,i])
                  
                  
                }
        }
    
    }
    # derived parameters
    # true richness
    for (t in 1:nint){
      for (i in 1:nsites) {
        SRexp[t,i] <- sum (z[,t,i])
      }
    }
    
    
    
   # change
   for (t in 2:nint) {  
    for (i in 1:nsites) {  
      
       #propcH[t,i] <-(sum(z[,t,i])-sum (z[,t-1,i]))/sum(z[,t-1,i]) ### turnover (proportional gain or loss)
       RER[t,i] <- (1-phi[t,i])/gamma[t,i] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[t,i] <- (1-phi[t,i])-gamma[t,i] ## net diversification rate (r= μ - λ; Rabosky 2018) of each tim
     
      }
   }
    
        
    } ## end of the model
    
    ",fill = TRUE)
sink()

#  ----------------------------------------------------


# load data
load (here ("processed_data", "site_covs.RData"))

# scale variables
elevation<-(scale (site_covs$elevation)) # scaled elevation
temperature <- (scale (site_covs$temp)) # scaled temp
precipitation <-  (scale (site_covs$prec))
lat <- ( (site_covs$paleolat))

# reverse the order of columns
elevation<-elevation[,rev(colnames(elevation))]
temperature<-temperature[,rev(colnames(temperature))]
precipitation<-precipitation[,rev(colnames(precipitation))]


# bins in data
elevation <- elevation [, which(colnames(elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(elevation)==bins [which(bins$bin %in% time_bins),"interval_name"]
temperature <- temperature [, which(colnames(temperature) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(temperature)==bins [which(bins$bin %in% time_bins),"interval_name"]
precipitation <- precipitation [, which(colnames(precipitation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(precipitation)==bins [which(bins$bin %in% time_bins),"interval_name"]

# subsetting of cells
precipitation<-(precipitation[which(rownames(precipitation) %in% cells),])
temperature<-(temperature[which(rownames(temperature) %in% cells),])
elevation<-(elevation[which(rownames(elevation) %in% cells),])
lat<-(lat[which(rownames(elevation) %in% cells)])


# sel sites with 5
sel_sites <- which(colSums(apply(array_genus_bin_site, c(1,3), sum))>=5)


## bundle data
# bundle one dataset per site
str(jags.data <- list(y = array_genus_bin_site[,,sel_sites],
                      ngen = dim(array_genus_bin_site[,,sel_sites])[1],
                      nint= dim(array_genus_bin_site[,,sel_sites])[2], 
                      nsites=dim(array_genus_bin_site[,,sel_sites])[3],
                      elevation = t((scaled_elevation))[,sel_sites],
                      temperature = t((scaled_temperature))[,sel_sites],
                      precipitation = t((scaled_precipitation))[,sel_sites],
                      nform = t(formations_per_site_interval)[,sel_sites],
                      range = scaled_range[,1],
                      time = scaled_time_stage[,1],
                      lat = scaled_lat[sel_sites]
                      
)
)

# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = array_genus_bin_site[,,sel_sites])}


## Parameters to monitor
## long form
params <- c(
  "intercept.gamma",
  "beta.gamma.elev",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "beta.phi.elev",
  "intercept.p",
  "beta.p.time",
  "beta.p.range",
  "gamma",
  "phi",
  "p",
  "SRexp",
  #"propcH",
  "R0",
  "RER"
)

  
# MCMC runs
# models
samples_paleo_cynodontia_sites_binomial <- jags (data = jags.data, 
                                    parameters.to.save = params, 
                                    model.file = "dyn_model_vectorized_covariates_CMR_sites_binomial.txt", 
                                    inits = inits, 
                                    n.chains = nc, 
                                    n.thin = nt, 
                                    n.iter = ni, 
                                    n.burnin = nb, 
                                    DIC = T,  
                                    n.cores=nc,
                                    parallel=T
  )

#devtools::install_github("mikemeredith/saveJAGS")
require(saveJAGS)

#samples_paleo_cynodontia_sites_binomial <- saveJAGS(jags.data, inits, params, 
#                                       modelFile="dyn_model_vectorized_covariates_CMR_sites_binomial.txt",
#                                       chains=nc, 
#                                       sample2save=((ni-nb)/nt), 
#                                       nSaves=5, 
#                                       burnin=nb, 
#                                       thin=nt,
#                                       fileStub="output/CMR_sites_binomial")
#
#
# run save jags
samples_paleo_cynodontia_sites_binomial <- saveJAGS(jags.data, inits, params, 
                                         modelFile="dyn_model_vectorized_covariates_CMR_sites_binomial.txt",
                                         chains=nc, 
                                         sample2save=((ni-nb)/nt), 
                                         nSaves=3, 
                                         burnin=nb, 
                                         thin=nt,
                                         fileStub="output/CMR_sites_binomial")

newRes_samples_cynodontia_CMR_sites_binomial <- combineSaves(recoverSaves("output/CMR_sites_binomial"))

# save each step
# save (samples_paleo_cynodontia_sites_binomial,file = here ("output","samples_paleo_cynodontia_sites_binomial.RData"))
  


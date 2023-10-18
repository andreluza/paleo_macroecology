# --------------------------------------

# 2 . modeling


## Global MCMC settings
######################

# number of interactions for a global-scale analyses
na <- 30000; nb <- 40000; ni <- 50000; nc <- 3; nt <- 20

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

sink("global_CMRmodel_vectorized_no_covariates.txt")
cat("
   
    model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    
    
   for (k in 1:ntaxa) {  
       
      ## set initial conditions for each taxon
      psi1[k] ~ dunif(0,1)
      #psi1[k] ~ dbeta(a[k],b[k])
      # Specify the hyperparameters for the Beta distribution
      #a[k] ~ dunif (0,100)   # probability of success
      #b[k] ~ dunif (0,100)   # probability of failures
      
       
       for (t in 1:nint) {
        
          gamma [k,t]~dunif(0,1)
          phi [k,t]~dunif(0,1)
        
        }
   }
   
   ############      Model       #############
    
    for (k in 1:ntaxa){
    
      for (g in 1:ngen) {
    
          z[k,g,1]~dbern(psi1[k]) # occupancy status initialization
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[k,g,t] <- z[k,g,t-1] *  phi[k,t] + ### if occupied, p of not getting extinct/persist in the next time
                              (1-z[k,g,t-1]) *  gamma[k,t] ###  if not occupied, p of originate in the next time
                
        # realized occurrence
		    z[k,g,t] ~ dbern(muZ[k,g,t])
    
        } # t
        
      } # g
      
    } # k
    
   
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
   
    ###  detection intercept
    for (k in 1:ntaxa) {
        
      # Specify the hyperparameters for the Beta distribution
      #a.d[k] ~ dunif (0,100)   # probability of success
      #b.d[k] ~ dunif (0,100)   # probability of failures
        
        for (t in 1:nint) {
        
            p[k,t] ~ dunif (0,1)
        
            #p[k,t] ~ dbeta(a.d[k],b.d[k])
              
        }
        
    }
    
    
    
     
    # observation submodel
    for (n in 1:nobs) { ## loop over observations 
    
         
        # observation
        # Specify the binomial observation model conditional on occupancy state
        y[n] ~ dbin(z[taxon[n],genus[n],stage[n]]*p[taxon[n],stage[n]], nform[n])
                  
       }     
   
      
      
      
    # derived parameters
    for (k in 1:ntaxa) {
      for (t in 1:nint){
      
        SRexp[k,t] <- sum (z[k,,t]) # true richness estimate
        
      }
    }
    
   # descriptors of change
   for (k in 1:ntaxa) {
    
    for (t in 2:nint){
       FSS[k,t] <- sum (muZ[k,,t]) # finite sample size
       CH[k,t] <- (sum(z[k,,t]) - sum(z[k,,t-1])) ### turnover (proportional gain or loss)
       CH1[k,t]<-sum(z[k,,t-1])
       RER[k,t] <- (1-phi[k,t])/gamma[k,t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
       R0[k,t] <- (1-phi[k,t])-gamma[k,t]  ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
       
    }
   
   } 
        
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()



# initial vals
#z_init <- tapply (vectorized_occ_global$detection,
#                  list(vectorized_occ_global$clade,
#                       vectorized_occ_global$genus,
#                       vectorized_occ_global$stage),
#                  max,na.rm=T)

# conditionate in the first appearance of each genus

#vectorized_occ_global<- lapply (unique (vectorized_occ_global$genusID), function (i) {
#  lines_sel <- which(vectorized_occ_global$genusID == i)
#  t1 <- vectorized_occ_global [lines_sel,"detection"]
#  output <- vectorized_occ_global[lines_sel [which(t1 == 1):length(lines_sel)],]
#  output
#})
# melt
#vectorized_occ_global<-do.call(rbind , vectorized_occ_global)

# --------------------------------------

str(jags.data <- list(
  y = vectorized_occ_global$detection,
  ntaxa = length(unique(vectorized_occ_global$clade)),
  ngen = length(unique(vectorized_occ_global$genusID)),
  nint = length(unique(vectorized_occ_global$stage)),
  nobs = nrow(vectorized_occ_global),
  taxon = vectorized_occ_global$clade,
  genus = vectorized_occ_global$genus,
  stage = vectorized_occ_global$stage-6,
  nform = vectorized_occ_global$formations
  
)
)

table(rowSums(tapply (jags.data$y,
                      list(jags.data$genus,
                           jags.data$stage),
                      max)) == rowSums(array_genus_bin [order(array_genus_bin$unique_name),-c(1,2)]))

# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = tapply (jags.data$y,
                                      list(jags.data$taxon,
                                           jags.data$genus,
                                           jags.data$stage),
                                      max),
                          
                          psi1 = runif (3)
                          
                          
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
                                  model.file = "global_CMRmodel_vectorized_no_covariates.txt", 
                                  inits = inits, 
                                  n.chains = nc, 
                                  n.thin = nt, 
                                  n.iter = ni, 
                                  n.burnin = nb, 
                                  DIC = T,  
                                  n.cores=nc,
                                  parallel=F
)

# check
plot(samples_paleo_cynodontia_no_covariate$mean$FSS[1,],type="b",col = "red")
lines(samples_paleo_cynodontia_no_covariate$mean$SRexp[1,],col = "black",lty=2)


plot(samples_paleo_cynodontia_no_covariate$mean$FSS[2,],type="b",col = "green")
lines(samples_paleo_cynodontia_no_covariate$mean$SRexp[2,],col = "black",lty=2)

plot(samples_paleo_cynodontia_no_covariate$mean$FSS[3,],type="b",col = "blue")
lines(samples_paleo_cynodontia_no_covariate$mean$SRexp[3,],col = "black",lty=2)

# save
save (samples_paleo_cynodontia_no_covariate,
      file = here ("output","no_covariates_CMR_global_binomial.RData"))


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
                                                            vectorized_occ_global <- vectorized_occ_global %>% 
                                                              filter (cladeID == i)
                                                                            
                str(jags.data <- list(
                        y = vectorized_occ_global$detection,
                        #ntaxa = length(unique(vectorized_occ_global$clade)),
                        ngen = length(unique(vectorized_occ_global$genusID)),
                        nint = length(unique(vectorized_occ_global$stage)),
                        nobs = nrow(vectorized_occ_global),
                        #taxon = vectorized_occ_global$clade,
                        genus = vectorized_occ_global$genus,
                        stage = vectorized_occ_global$stage-6,
                        nform = vectorized_occ_global$formations,
                        # covariates
                        # time = vectorized_occ$
                        range = (vectorized_occ_global$range_size-mean(vectorized_occ_global$range_size))/sd(vectorized_occ_global$range_size),
                        lat_obs = (vectorized_occ_global$latitude-mean(vectorized_occ_global$latitude))/sd(vectorized_occ_global$latitude),
                        time = (vectorized_occ_global$time-mean(vectorized_occ_global$time))/sd(vectorized_occ_global$time),
                        temp_obs = (vectorized_occ_global$temperature-mean(vectorized_occ_global$temperature))/sd(vectorized_occ_global$temperature),
                        # time covs 
                        temperature = (time_covariates$temperature-mean (time_covariates$temperature,na.rm=T))/sd(time_covariates$temperature,na.rm=T),
                        precipitation = (time_covariates$precipitation-mean (time_covariates$precipitation,na.rm=T))/sd(time_covariates$precipitation,na.rm=T)
                        
                                            
                      )
                      )
                      
                
                table(rowSums(tapply (jags.data$y,
                        list(jags.data$genus,
                             jags.data$stage),
                        max)) == rowSums(array_genus_bin [order(array_genus_bin$unique_name),-c(1,2)]))
                
                # Set initial values
                # function inits (to be placed in each chain)
                inits <- function(){ list(z = tapply (jags.data$y,
                                                      list(jags.data$taxon,
                                                           jags.data$genus,
                                                           jags.data$stage),
                                                      max),
                                          
                                          intercept.gamma = rnorm (3),
                                          beta.gamma.prec = rnorm (3),
                                          beta.gamma.temp = rnorm (3),
                                          intercept.phi = rnorm (3),
                                          beta.phi.prec = rnorm (3),
                                          beta.phi.temp = rnorm (3),
                                          intercept.p = rnorm (3),
                                          beta.p.time = rnorm (3),
                                          beta.p.range = rnorm (3),
                                          beta.p.lat = rnorm (3),
                                          beta.p.temp = rnorm (3),
                                          psi1 = runif (3)
                                          
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
                samples_paleo_cynodontia_binomial <- jags (data = jags.data, 
                                                               parameters.to.save = params, 
                                                               model.file = "global_CMRmodel_vectorized_covariates.txt", 
                                                               inits = inits, 
                                                               n.chains = nc, 
                                                               n.thin = nt, 
                                                               n.iter = ni, 
                                                               n.burnin = nb, 
                                                               DIC = T,  
                                                               n.cores=nc,
                                                               parallel=F
                )


#plot(tapply (samples_paleo_cynodontia_binomial$mean$p,
#        list(jags.data$taxon,
#             jags.data$stage),
#        mean)[1,],type="b")
#plot(tapply (samples_paleo_cynodontia_binomial$mean$p,
#             list(jags.data$taxon,
#                  jags.data$stage),
#             mean)[2,],type="b")


    # MCMC runs
    samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                                                        modelFile="global_CMRmodel_vectorized_covariates.txt",
                                                        chains=nc, 
                                                        sample2save=((ni-nb)/nt), 
                                                        nSaves=5, 
                                                        burnin=nb, 
                                                        thin=nt,
                                                        fileStub="output/CMR_global_binomial")

})

#--------------------------------------------------------------------- #

# - binomial model with covariates at the site/region level

# try nimble to speed up the process
#library (nimble)

#Section6p4_code <- nimbleCode( {

sink("regional_CMRmodel_vectorized_covariates.txt")
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
    for (n in 1:nobs) { ## loop over observations 
    
        # p independent on the number of formations
        logit(p[n])<- intercept.p+
                        beta.p.lat*lat_obs[n]+
                        beta.p.range*range[n]+
                        beta.p.time*time[n]+
                        beta.p.temp*temp_obs[n]
          
          
        # observation
        # Specify the binomial observation model conditional on occupancy state
        y[n] ~ dbin(z[genus[n],stage[n],site[n]]*p[n], nform[n])
                  
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
                                                            
                                                            
                                                            # select clade                                                             
                                                            vectorized_occ <- vectorized_occ %>% 
                                                              filter (cladeID == i)
                                                            
                                                            str(jags.data <- list(y = vectorized_occ$detection,
                                                                                  #ntaxa = length(unique(vectorized_occ$clade)),
                                                                                  nsites = length(unique(vectorized_occ$site)),
                                                                                  ngen = length(unique(vectorized_occ$genusID)),
                                                                                  nint = length(unique(vectorized_occ$stage)),
                                                                                  nobs = nrow(vectorized_occ),
                                                                                  # taxon = vectorized_occ$clade,
                                                                                  genus = vectorized_occ$genus,
                                                                                  stage = vectorized_occ$stage,
                                                                                  site = vectorized_occ$site,
                                                                                  nform = vectorized_occ$formations,
                                                                                  # covariates
                                                                                  # time = vectorized_occ$
                                                                                  range = (vectorized_occ$range_size-mean(vectorized_occ$range_size))/sd(vectorized_occ$range_size),
                                                                                  lat_obs = (vectorized_occ$latitude-mean(vectorized_occ$latitude))/sd(vectorized_occ$latitude),
                                                                                  time = (vectorized_occ$time-mean(vectorized_occ$time))/sd(vectorized_occ$time),
                                                                                  temp_obs = (vectorized_occ$temperature-mean(vectorized_occ$temperature))/sd(vectorized_occ$temperature),
                                                                                  # lat = site_covs
                                                                                  temperature = (region_temperature-mean (region_temperature,na.rm=T))/sd(region_temperature,na.rm=T),
                                                                                  precipitation = (region_precipitation-mean (region_precipitation,na.rm=T))/sd(region_precipitation,na.rm=T),
                                                                                  lat = as.vector(scale (region_latitude))
                                                            )
                                                            )
                                                            
                                                            
                                                            # Set initial values
                                                            z_inits<-tapply (vectorized_occ$detection,
                                                                             list(vectorized_occ$genus,
                                                                                  vectorized_occ$stage,
                                                                                  vectorized_occ$site),
                                                                             max)
                                                            
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
                                                                                                                modelFile="regional_CMRmodel_vectorized_covariates.txt",
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
                                                 model.file = "regional_CMRmodel_vectorized_covariates.txt", 
                                                 inits = inits, 
                                                 n.chains = nc, 
                                                 n.thin = nt, 
                                                 n.iter = ni, 
                                                 n.burnin = nb, 
                                                 DIC = T,  
                                                 n.cores=nc,
                                                 parallel=F
)

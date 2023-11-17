# --------------------------------------

# 2 . modeling


# conditional on the first and last detection of the taxon


# saveJAGS help
# https://web.archive.org/web/20220408094124/http://mmeredith.net/blog/2018/Intro_saveJAGS.htm

## Global MCMC settings
######################

# number of interactions for a global-scale analyses
na <- 15000; nb <- 30000; ni <- 50000; nc <- 3; nt <- 20

## short form for tests
#na <- 400; nb <- 700; ni <- 1000; nc <- 3; nt <- 1
#na <- 40; nb <- 70; ni <- 100; nc <- 3; nt <- 1

# load packages
source("R/packages.R")

# load data
load(here ("processed_data","CMR_data_observation_cov.RData"))

# load environmental data
load (here ("processed_data", "site_covs.RData"))




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
      
       # Site Occupancy Dynamics (Priors)
       
       # intercepts
       # gamma
       for (t in 2:nint) {
       
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
            y [g,t] ~ dbin(muY[g,t],nform[t])
            muY[g,t] <- z[g,t]*p[t]      
                  
                }
        }
      
     
      
    # derived parameters
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
    
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
                                                      A0_input <- array (0, dim = c(ifelse (i == "Non-mammalian Mammaliaformes",
                                                                                            200,
                                                                                            100),
                                                                                    length(time_bins)))
            
            
            # bundle the data                                                       
            str(jags.data <- list(
              y =rbind (
                
                data.matrix(array_genus_bin [which(clades %in% i),]),
                A0_input
              ),
              ngen = dim(array_genus_bin [which(clades %in% i),])[1] + nrow(A0_input),
              nint = dim(array_genus_bin [which(clades %in% i),])[2],
              nform = formations_per_interval$formations_per_interval
          
        )
        )
        
        # condition the dataset to one stage before the first appearance in the fossil record, and one stage after its last appearance (detection)    
        sel_cols <- function_stages(jags.data$y)        
        
        # bundle conditional data
        str(jags.data <- list(
          y = jags.data$y[,sel_cols],
          ngen = nrow (jags.data$y[,sel_cols]),
          nint = ncol(jags.data$y[,sel_cols]),
          nform = formations_per_interval$formations_per_interval[function_stages(jags.data$y)]
          
        )
        )
        
            
        # Set initial values
        # function inits (to be placed in each chain)
        inits <- function(){ list(z = jags.data$y,
                                  
                                  psi1 = runif (1)
                                  
                                  
        )}
        
        
        ## Parameters to monitor
        params <- c("psi1","gamma", "phi","p","z", "muY", "SRexp")
        
        # MCMC runs
        # models
        # usual jags to test
        samples_paleo_cynodontia_binomial_no_cov <- saveJAGS(jags.data, inits, params, 
                                                      modelFile="global_CMRmodel_matrix_no_covariates.txt",
                                                      chains=nc, 
                                                      sample2save=((ni-nb)/nt), 
                                                      nSaves=6, 
                                                      burnin=nb, 
                                                      thin=nt,
                                                      fileStub=paste ("output/No_covariate_model/CMR_global_binomial_no_cov", i,sep="_"))

                                                    
        })


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
       
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        # intercepts
        gamma.u ~ dunif(0,1) # range origination
        intercept.gamma <- logit(gamma.u) # intercept origination
        # regression coeff
        beta.gamma.prec ~ dnorm(0, 0.01)# precipitation
        beta.gamma.temp ~ dnorm(0, 0.01)# temperature
        beta.gamma.lat ~ dnorm(0,0.01) # latitude
        beta.gamma.area ~ dnorm(0,0.01) # area
     
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta.phi.prec ~ dnorm(0, 0.01)# precipitation
        beta.phi.temp ~ dnorm(0, 0.01)# temperature
        beta.phi.lat ~ dnorm(0,0.01) # lat
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
    # regression coef
    beta.p.lat ~ dnorm(0,0.01)
    beta.p.range ~ dnorm(0,0.01)
    beta.p.time ~ dnorm(0,0.01)
    beta.p.temp ~ dnorm(0,0.01)
     
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
          y[g,t] ~ dbin(muY[g,t], nform[t])
          muY[g,t] <- z[g,t]*p[g,t]
                  
        }
      
    }
      
    # derived parameters
    for (t in 1:nint){
        SRexp[t] <- sum (z[,t])
      }
    
    
    }## end of the model
    
    ",fill = TRUE)
sink()



# dir create
dir.create (here ("output", "global"))


# correlation between variables
require(corrplot)

corrplot (cor(time_covariates [,c("temperature", "precipitation", "area")]),
          method = "square",
          diag=F,
          col = COL2('RdYlBu', 8),
          mar = rep(4,4),
          title = "Global-scale covariates")


# region scale
cor (as.vector(region_area), as.vector(region_precipitation))
cor (as.vector(region_temperature), as.vector(region_precipitation))
cor (as.vector(region_temperature), as.vector(region_area))
cor (as.vector(region_temperature), as.vector(matrix(region_latitude,nrow=33,ncol=9,byrow=T)))
cor (as.vector(region_precipitation), as.vector(matrix(region_latitude,nrow=33,ncol=9,byrow=T)))
cor (as.vector(region_area), as.vector(matrix(region_latitude,nrow=33,ncol=9,byrow=T)))


corrplot (cor(data.frame (
  area = as.vector(region_area),
  precipitation = as.vector(region_precipitation),
  temperature = as.vector(region_temperature),
  latitude = as.vector(matrix(region_latitude,nrow=33,ncol=9,byrow=T)))),
          method = "square",
          diag=F,
          col = COL2('RdYlBu', 8),
  mar = rep(4,4),
  title = "Region-scale covariates")




## bundle data
# run
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                          "Non-mammalian Mammaliaformes",
                                                          "Mammalia" ), function (i) {
                                                            
                      
                                                            
                   # augment the dataset
                   A0_input <- array (0, dim = c(ifelse (i == "Non-mammalian Mammaliaformes",
                                      100, # 200
                                      100),
                                      length(time_bins)))
                   range_input<- array(0,dim = nrow(A0_input))                                            
                   temp_input<- array(0,dim = nrow(A0_input))                                            
                   lat_input<- array(0,dim = nrow(A0_input))                                            
                                                            
                                                                                                  
                 # bundle the data                                                       
                 str(jags.data <- list(
                        y =rbind (
                          
                          data.matrix(array_genus_bin [which(clades %in% i), ]),
                          A0_input
                        ),
                        ngen = dim(array_genus_bin [which(clades %in% i), ])[1] + nrow(A0_input),
                        nint = dim(array_genus_bin [which(clades %in% i), ])[2],
                        nform = formations_per_interval$formations_per_interval,
                        # covariates
                        # time = vectorized_occ$
                        range = c(
                          
                            as.vector(scale(range_area_taxon[which(clades %in% i),"range_area"])),
                            
                            range_input),
                        
                        
                        #range_area_taxon[which(range_area_taxon$taxon %in% array_genus_bin [which(array_genus_bin$clade %in% i), 2]),"taxon"] == array_genus_bin [which(array_genus_bin$clade %in% i), 2]
                        lat_obs = c( 
                          
                          as.vector(scale(observation_covariates[which(clades %in% i),"latitude"])),
                          
                          lat_input
                        ),
                        
                        temp_obs = c (
                          
                          as.vector(scale(observation_covariates[which(clades %in% i),"temperature"])),
                          
                          temp_input
                        )
                                            
                      )
                      )
                      
                
                 
                 # condition the dataset to one stage before the first appearance in the fossil record, and one stage after its last appearance (detection)    
                 
                 
                 # bundle conditional data
                 str(jags.data <- list(
                   y = jags.data$y[,function_stages(jags.data$y)],
                   ngen = nrow (jags.data$y[,function_stages(jags.data$y)]),
                   nint = ncol(jags.data$y[,function_stages(jags.data$y)]),
                   
                   # observation covs
                   nform = formations_per_interval$formations_per_interval[function_stages(jags.data$y)],
                   range = jags.data$range,
                   lat_obs = jags.data$lat_obs,
                   temp_obs = jags.data$temp_obs,
                   # add time
                   time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
                   # time covs 
                   temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
                   precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
                   area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)]))
                   
                 ))
                 
                 
                 
                                                            
                # Set initial values
                # function inits (to be placed in each chain)
                inits <- function(){ list(z = jags.data$y,
                                          
                                          gamma.u = runif (1),
                                          beta.gamma.prec = rnorm (1),
                                          beta.gamma.temp = rnorm (1),
                                          beta.gamma.area = rnorm (1),
                                          phi.u = runif (1),
                                          beta.phi.prec = rnorm (1),
                                          beta.phi.temp = rnorm (1),
                                          beta.phi.area = rnorm (1),
                                          p.u = runif (1),
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
                  "p.u",
                  "intercept.p",
                  "beta.p.time",
                  "beta.p.range",
                  "beta.p.lat",
                  "beta.p.temp",
                  "gamma",
                  "phi",
                  "p",
                  "muY",
                  "psi1",
                  "muY",
                  "z",
                  "SRexp"
                  
                )
                
                # usual jags to test
                samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                         modelFile="global_CMRmodel_matrix_covariates.txt",
                         chains=nc, 
                         sample2save=((ni-nb)/nt), 
                         nSaves=5, 
                         burnin=nb, 
                         thin=nt,
                         fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
})


# extend chains
newRes <-lapply (c ("Non-mammaliaform cynodonts",
           "Non-mammalian Mammaliaformes",
           "Mammalia" ), function (i) {
             
             newRes <- resumeJAGS(fileStub=
                                    paste ("output/global/CMR_global_binomial", i,sep="_"),
                        nSaves=3)
           })

## MCMC runs
samples_paleo_cynodontia_binomial <-jags (data = jags.data, 
                                          parameters.to.save = params, 
                                          model.file = "global_CMRmodel_matrix_covariates.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.thin = nt, 
                                          n.iter = ni, 
                                          n.burnin = nb, 
                                          DIC = T,  
                                          n.cores=nc,
                                          parallel=T
)


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
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        # intercepts
        gamma.u ~ dunif(0,1) # range origination
        intercept.gamma <- logit(gamma.u) # intercept origination
        # regression coeff
        beta.gamma.prec ~ dnorm(0, 0.01)# precipitation
        beta.gamma.temp ~ dnorm(0, 0.01)# temperature
        #beta.gamma.lat ~ dnorm(0,0.01) # latitude
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
        #beta.phi.lat ~ dnorm(0,0.01)
        beta.phi.area ~ dnorm(0,0.01)
        
       
        ## set initial conditions for the whole system
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
                                    #beta.gamma.lat*lat[i]+
                                    beta.gamma.area*area[t,i]
            # persistence
            logit(phi[t,i]) <-    intercept.phi + 
                                    beta.phi.prec*precipitation[t,i]+
                                    beta.phi.temp*temperature[t,i]+
                                    #beta.phi.lat*lat[i]+
                                    beta.phi.area*area[t,i]
    
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
    # intercept    
    p.u ~ dunif(0,1) # range origination
    intercept.p <- logit(p.u) # intercept origination
    
    # regression coef
    beta.p.lat ~ dnorm(0,0.01)
    beta.p.range ~ dnorm(0,0.01)
    beta.p.time ~ dnorm(0,0.01)
    beta.p.temp ~ dnorm(0,0.01)
        
     
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
              y[g,t,i] ~ dbin(muY[g,t,i], nform[t,i])
              muY[g,t,i] <- z[g,t,i]*p[g,t]
                      
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
 
    
        
    } ## end of the model
    
    ",fill = TRUE)
sink()

#  ----------------------------------------------------

# bundle data
# replace NAs by a very small number
# create dir
dir.create (here ("output", "region"))

# subset the data
range_area_taxon <- range_area_taxon[match(genus ,range_area_taxon$taxon),]

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
                                              range = as.vector(scale(range_area_taxon[which(clades %in% i),"range_area"])),
                                              lat_obs = as.vector(scale(observation_covariates[which(clades %in% i),"latitude"])),
                                              temp_obs = as.vector(scale(observation_covariates[which(clades %in% i),"temperature"]))
                        )
                        )
                        sel_int <- function_stages(apply (jags.data$y,c(1,2), max))
                        # condition on the first and last appearance
                         str(jags.data <- list(y = jags.data$y [,sel_int,],
                                               ngen = jags.data$ngen,
                                               nint = length(sel_int),
                                               nsites = dim(jags.data$y)[3],
                                               nform = as.matrix(t(formations_per_site_interval[,sel_int])),
                                               # add time
                                               time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][sel_int])),
                                               
                                               # covariates
                                               # time = vectorized_occ$
                                               range = jags.data$range,
                                               lat_obs = jags.data$lat_obs,
                                               temp_obs = jags.data$temp_obs,
                                               # lat = site_covs
                                               temperature = (region_temperature[sel_int,]-mean (region_temperature[sel_int,],na.rm=T))/sd(region_temperature[sel_int,],na.rm=T),
                                               precipitation = (region_precipitation[sel_int,]-mean (region_precipitation[sel_int,],na.rm=T))/sd(region_precipitation[sel_int,],na.rm=T),
                                               area = (region_area[sel_int,]-mean (region_area[sel_int,],na.rm=T))/sd(region_area[sel_int,],na.rm=T)
                                               #lat = as.vector(scale (region_latitude))
                         )
                         )
                                                            
                                                                                                
                                                            
                        # Set initial values
                        z_inits<-jags.data$y 
                        
                        #z_inits<-ifelse (is.na(z_inits),1,z_inits)
                        # function inits (to be placed in each chain)
                        inits <- function(){ list(z=z_inits,
                                                  gamma.u = runif (1),
                                                  beta.gamma.prec = rnorm (1),
                                                  beta.gamma.temp = rnorm (1),
                                                  #beta.gamma.lat = rnorm (1),
                                                  beta.gamma.area = rnorm (1),
                                                  phi.u = runif (1),
                                                  beta.phi.prec = rnorm (1),
                                                  beta.phi.temp = rnorm (1),
                                                  #beta.phi.lat = rnorm (1),
                                                  beta.phi.area = rnorm (1),
                                                  p.u = runif (1),
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
                          #"beta.gamma.lat",
                          "beta.gamma.area",
                          "intercept.phi",
                          "beta.phi.prec",
                          "beta.phi.temp",
                          #"beta.phi.lat",
                          "beta.phi.area",
                          "intercept.p",
                          "beta.p.time",
                          "beta.p.range",
                          "beta.p.lat",
                          "beta.p.temp",
                          "gamma",
                          "phi",
                          "p",
                          "muY",
                          "z"
                          
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
                                                                            fileStub=paste ("output/region/CMR_regional_binomial", i,sep="_"))
                        
                        
})

# extend chains
newRes <-lapply (c ("Non-mammaliaform cynodonts",
                    "Non-mammalian Mammaliaformes",
                    "Mammalia" ), function (i) {
                      
                      newRes <- resumeJAGS(fileStub=
                                             paste ("output/region/CMR_regional_binomial", i,sep="_"),
                                           nSaves=7)
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



# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# region with area effect

#--------------------------------------------------------------------- #

# - binomial model with covariates at the site/region level (random effect intercept)


sink("regional_CMRmodel_matrix_covariates_rdm_area.txt")
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
        for (i in 1:nsites) {
        
          intercept.phi[i]~dnorm(mu.phi[i], tau.phi[i]) # intercept persistence
          intercept.gamma[i]~dnorm(mu.gamma[i], tau.gamma[i]) # intercept origination
        
        }
        
        # hyperpriors
        for (i in 1:nsites) {
          # phi
          m.phi[i]~dunif(0,1)
          mu.phi[i]<-logit(m.phi[i])
          sigma.phi[i]~dunif(0,4)
          tau.phi[i]<-1/sigma.phi[i]^2
          
          # gamma
          m.gamma[i]~dunif(0,1)
          mu.gamma[i]<-logit(m.gamma[i])
          sigma.gamma[i]~dunif(0,4)
          tau.gamma[i]<-1/sigma.gamma[i]^2
          
        }        
        
        # regression coeff
        beta.gamma.prec ~ dnorm(0, 0.01)# precipitation
        beta.gamma.temp ~ dnorm(0, 0.01)# temperature
        beta.gamma.area ~ dnorm(0, 0.01)# temperature
        #beta.gamma.lat ~ dnorm(0,0.01) # latitude
     
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        
        # regression coeff
        beta.phi.prec ~ dnorm(0, 0.01)# precipitation
        beta.phi.temp ~ dnorm(0, 0.01)# temperature
        beta.phi.area ~ dnorm(0, 0.01)# temperature
        #beta.phi.lat ~ dnorm(0,0.01)
       
        ## set initial conditions for the whole system
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
            logit(gamma[t,i]) <-  intercept.gamma[i] + 
                                    beta.gamma.prec*precipitation[t,i]+
                                    beta.gamma.temp*temperature[t,i]+
                                    #beta.gamma.lat*lat[i]+
                                    beta.gamma.area*area[t,i]
                                    
            # persistence
            logit(phi[t,i]) <-    intercept.phi[i] + 
                                    beta.phi.prec*precipitation[t,i]+
                                    beta.phi.temp*temperature[t,i]+
                                    #beta.phi.lat*lat[i]+
                                    beta.phi.area*area[t,i]
    
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
    # intercept    
    p.u ~ dunif(0,1) # range origination
    intercept.p <- logit(p.u) # intercept origination
    
    # regression coef
    beta.p.lat ~ dnorm(0,0.01)
    beta.p.range ~ dnorm(0,0.01)
    beta.p.time ~ dnorm(0,0.01)
    beta.p.temp ~ dnorm(0,0.01)
        
     
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
              y[g,t,i] ~ dbin(muY[g,t,i], nform[t,i])
              muY[g,t,i] <- z[g,t,i]*p[g,t]
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
 
#for (t in 2:nint) {  
#     
#       for (i in 1:nsites) {  
#        
#          FSS[t,i] <- sum (muZ[,t,i]) # finite sample size 
#          CH[t,i] <- (sum(z[,t,i]) - sum(z[,t-1,i])) ### turnover (proportional gain or loss)
#          CH1[t,i]<-sum(z[,t-1,i])
#          RER[t,i] <- (1-phi[t,i])/gamma[t,i] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
#          R0[t,i] <- (1-phi[t,i])-gamma[t,i] ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
#          psi.eq[t,i]<-gamma[t,i]/((1-phi[t,i])+gamma[t,i]) # equilibrium occupancy
#          turnover[t,i]<-(gamma[t,i]*(1-phi[t,i]))/(gamma[t,i]+(1-phi[t,i]))
#      
# 
#    }
# 
#   }
  
  
    
        
    } ## end of the model
    
    ",fill = TRUE)
sink()

#  ----------------------------------------------------

# bundle data
# replace NAs by a very small number

# create dir
dir.create (here ("output", "region_rdm_area"))

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
                                                                                  range = as.vector(scale(range_area_taxon[which(clades %in% i),"range_area"])),
                                                                                  lat_obs = as.vector(scale(observation_covariates[which(clades %in% i),"latitude"])),
                                                                                  temp_obs = as.vector(scale(observation_covariates[which(clades %in% i),"temperature"]))
                                                            )
                                                            )
                                                            sel_int <- function_stages(apply (jags.data$y,c(1,2), max))
                                                            # condition on the first and last appearance
                                                            str(jags.data <- list(y = jags.data$y [,sel_int,],
                                                                                  ngen = jags.data$ngen,
                                                                                  nint = length(sel_int),
                                                                                  nsites = dim(jags.data$y)[3],
                                                                                  nform = as.matrix(t(formations_per_site_interval[,sel_int])),
                                                                                  # add time
                                                                                  time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][sel_int])),
                                                                                  
                                                                                  # covariates
                                                                                  # time = vectorized_occ$
                                                                                  range = jags.data$range,
                                                                                  lat_obs = jags.data$lat_obs,
                                                                                  temp_obs = jags.data$temp_obs,
                                                                                  # lat = site_covs
                                                                                  temperature = (region_temperature[sel_int,]-mean (region_temperature[sel_int,],na.rm=T))/sd(region_temperature[sel_int,],na.rm=T),
                                                                                  precipitation = (region_precipitation[sel_int,]-mean (region_precipitation[sel_int,],na.rm=T))/sd(region_precipitation[sel_int,],na.rm=T),
                                                                                  area = (region_area[sel_int,]-mean (region_area[sel_int,],na.rm=T))/sd(region_area[sel_int,],na.rm=T)
                                                                                  #lat = as.vector(scale (region_latitude))
                                                            )
                                                            )
                                                            
                                                            
                                                            
                                                            # Set initial values
                                                            z_inits<-jags.data$y 
                                                            
                                                            #z_inits<-ifelse (is.na(z_inits),1,z_inits)
                                                            # function inits (to be placed in each chain)
                                                            inits <- function(){ list(z=z_inits,
                                                                                      gamma.u = runif (1),
                                                                                      beta.gamma.prec = rnorm (1),
                                                                                      beta.gamma.temp = rnorm (1),
                                                                                      #beta.gamma.lat = rnorm (1),
                                                                                      beta.gamma.area = rnorm (1),
                                                                                      phi.u = runif (1),
                                                                                      beta.phi.prec = rnorm (1),
                                                                                      beta.phi.temp = rnorm (1),
                                                                                      #beta.phi.lat = rnorm (1),
                                                                                      beta.phi.area = rnorm (1),
                                                                                      p.u = runif (1),
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
                                                              # "beta.gamma.lat",
                                                              "beta.gamma.area",
                                                              "intercept.phi",
                                                              "beta.phi.prec",
                                                              "beta.phi.temp",
                                                              #"beta.phi.lat",
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
                                                              "muY",
                                                              "z"
                                                              
                                                            )
                                                            
                                                            # MCMC runs
                                                            # run save jags
                                                            samples_paleo_cynodontia_sites_binomial_rdm <- saveJAGS(jags.data, inits, params, 
                                                                                                                    modelFile="regional_CMRmodel_matrix_covariates_rdm_area.txt",
                                                                                                                    chains=nc, 
                                                                                                                    sample2save=((ni-nb)/nt), 
                                                                                                                    nSaves=6, 
                                                                                                                    burnin=nb, 
                                                                                                                    thin=nt,
                                                                                                                    fileStub=paste ("output/region_rdm_area/CMR_regional_binomial_rdm", i,sep="_"))
                                                            
                                                            
                                                          })


# extend chains
newRes <-lapply (c ("Non-mammaliaform cynodonts",
                    "Non-mammalian Mammaliaformes",
                    "Mammalia" ), function (i) {
                      
                      newRes <- resumeJAGS(fileStub=
                                             paste ("output/region_rdm_area/CMR_regional_binomial_rdm", i,sep="_"),
                                           nSaves=3)
                    })

# usual jags to test
samples_paleo_cynodontia_sites_binomial_rdm_area <- jags (data = jags.data, 
                                                          parameters.to.save = params, 
                                                          model.file = "regional_CMRmodel_matrix_covariates_rdm_area.txt", 
                                                          inits = inits, 
                                                          n.chains = nc, 
                                                          n.thin = nt, 
                                                          n.iter = ni, 
                                                          n.burnin = nb, 
                                                          DIC = T,  
                                                          n.cores=nc,
                                                          parallel=F
)



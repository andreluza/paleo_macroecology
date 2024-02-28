# --------------------------------------

# 2 . modeling (GLOBAL SCALE)


# conditional on the first and last detection of the taxon


# saveJAGS help
# https://web.archive.org/web/20220408094124/http://mmeredith.net/blog/2018/Intro_saveJAGS.htm

## Global MCMC settings
######################

# number of interactions for a global-scale analyses
na <- 15000; nb <- 30000; ni <- 50000; nc <- 3; nt <- 20

## short form for tests
#na <- 600; nb <- 700; ni <- 1000; nc <- 3; nt <- 1
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
      
      # prior to describe the full community size
      omega~dunif(0,1)
      
      # superpopulation process
      for (g in 1:ngen) {
        
        w[g]~dbern(omega)
        
      }
      
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
                                                             A0_input <- array (0, dim = c(400,
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

# option save jags
samples_paleo_cynodontia_binomial_no_cov <- saveJAGS(jags.data, inits, params, 
                                                     modelFile="global_CMRmodel_matrix_no_covariates.txt",
                                                     chains=nc, 
                                                     sample2save=((ni-nb)/nt), 
                                                     nSaves=6, 
                                                     burnin=nb, 
                                                     thin=nt,
                                                     fileStub=paste ("output/No_covariate_model/CMR_global_binomial_no_cov", i,sep="_"))



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
      
      # prior to describe the full community size
      omega~dunif(0,1)
      
      # superpopulation process
      for (g in 1:ngen) {
        
        w[g]~dbern(omega)
        
      }   
      
      # Site Occupancy Dynamics (Priors)
       
       
        # ----------------------
        #     Gamma (origination)
        # ----------------------
        # intercepts
        gamma.u ~ dunif(0,1) # range origination
        intercept.gamma <- logit(gamma.u) # intercept origination
        # regression coeff
        beta.gamma.prec ~ dnorm(0,1/3^2)# precipitation
        beta.gamma.temp ~ dnorm(0,1/3^2)# temperature
        beta.gamma.area ~ dnorm(0,1/3^2) # area
        beta.gamma.coast ~ dnorm(0,1/3^2) # area
        
        # tropics
        beta.gamma.area.t ~ dnorm(0,1/3^2) # area
        beta.gamma.coast.t ~ dnorm(0,1/3^2) # area
     
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta.phi.prec ~ dnorm(0,1/3^2)# precipitation
        beta.phi.temp ~ dnorm(0,1/3^2)# temperature
        beta.phi.area ~ dnorm(0,1/3^2) # area
        beta.phi.coast ~ dnorm(0,1/3^2) # area
        
        # tropics
        beta.phi.area.t ~ dnorm(0,1/3^2) # area
        beta.phi.coast.t ~ dnorm(0,1/3^2) # area
        
        ## set initial conditions for occupancy of each genus
        psi1 ~ dunif(0,1)
        
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        for (t in 2:nint){
      
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   beta.gamma.prec*precipitation[t]+
                                   beta.gamma.temp*temperature[t]+
                                   beta.gamma.area*area[t]+
                                   beta.gamma.coast*coast[t]+
                                   beta.gamma.area.t*area.t[t]+
                                   beta.gamma.coast.t*coast.t[t]
                                   
                                   
                                         
              # persistence
              logit(phi[t]) <-  intercept.phi + 
                                  beta.phi.prec*precipitation[t]+
                                  beta.phi.temp*temperature[t]+
                                  beta.phi.area*area[t]+
                                  beta.phi.coast*coast[t]+
                                  beta.phi.area.t*area.t[t]+
                                  beta.phi.coast.t*coast.t[t]
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
                 muZW[g,t] <- muZ[g,t]*w[g]
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
    p.u ~ dunif(0,1) # range origination
    intercept.p <- logit(p.u) # intercept origination
    # regression coef
    beta.p.lat ~ dnorm(0,1/3^2)
    beta.p.range ~ dnorm(0,1/3^2)
    beta.p.time ~ dnorm(0,1/3^2)
    beta.p.temp ~ dnorm(0,1/3^2)
     
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
    
    # total metacommunity size
    Ntotal <- sum(w[])
    
    
    }## end of the model
    
    ",fill = TRUE)
sink()


# dir create
dir.create (here ("output", "global"))


# correlation between variables
require(corrplot)

corrplot (cor(
  
  cbind (
      time_covariates [,c("temperature", "precipitation", "area", "coastalLine")],
      data.frame (region_temperature, region_precipitation, region_area, region_coastalLine)
)),
          method = "square",
          diag=F,
          col = COL2('RdYlBu', 8),
          mar = rep(4,4),
          title = "Global-scale covariates")



## bundle data
# run

samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                    "Non-mammalian Mammaliaformes",
                                                    "Mammalia" ), function (i) {
                                                      
                                                      
                                                      # augment the dataset
                                                      A0_input <- array (0,dim=c(400,
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
                                                        nform = formations_per_interval,
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
                                                        nform = formations_per_interval[function_stages(jags.data$y)],
                                                        range = jags.data$range,
                                                        lat_obs = jags.data$lat_obs,
                                                        temp_obs = jags.data$temp_obs,
                                                        # add time
                                                        time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
                                                        # time covs 
                                                        temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
                                                        precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
                                                        area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
                                                        coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
                                                        area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
                                                        coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
                                                        
                                                      ))
                                                      
                                                      
                                                      
                                                      
                                                      # Set initial values
                                                      # function inits (to be placed in each chain)
                                                      inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                                                                                w = rep(1,nrow(jags.data$y)),
                                                                                gamma.u = runif (1),
                                                                                phi.u = runif (1),
                                                                                p.u = runif (1),
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
                                                        "beta.gamma.area.t",
                                                        "beta.gamma.coast",
                                                        "beta.gamma.coast.t",
                                                        "intercept.phi",
                                                        "beta.phi.prec",
                                                        "beta.phi.temp",
                                                        "beta.phi.area",
                                                        "beta.phi.area.t",
                                                        "beta.phi.coast",
                                                        "beta.phi.coast.t",
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
                                                        "z",
                                                        "SRexp",
                                                        "w",
                                                        "omega",
                                                        "muZW",
                                                        "Ntotal"
                                                        
                                                      )
                                                      
                                                      # usual jags to test
                                                      #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
                                                      #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
                                                      #                                              chains=nc, 
                                                      #                                              sample2save=((ni-nb)/nt), 
                                                      #                                              nSaves=5, 
                                                      #                                              burnin=nb, 
                                                      #                                              thin=nt,
                                                      #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
                                                      
                                                      
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
                                                                                                parallel=F
                                                      )
                                                      
                                                      save(samples_paleo_cynodontia_binomial,
                                                           file = paste0 ("output/global/CMR_global_binomial_", i,".RData"))
                                                      
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
                                          parallel=F
)

# ===================================================

# sensitivity

# more iterations
# only cynodonts
na <- 50000; nb <- 80000; ni <- 100000; nc <- 3; nt <- 20 # only cynodonts # put in the other script

# run
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                    #"Non-mammalian Mammaliaformes",
                                                    #"Mammalia" 
                                                    ), function (i) {
    
    
    # augment the dataset
    A0_input <- array (0,dim=c(400,
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
      nform = formations_per_interval,
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
      nform = formations_per_interval[function_stages(jags.data$y)],
      range = jags.data$range,
      lat_obs = jags.data$lat_obs,
      temp_obs = jags.data$temp_obs,
      # add time
      time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
      # time covs 
      temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
      precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
      area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
      coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
      area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
      coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
      
    ))
    
    
    
    
    # Set initial values
    # function inits (to be placed in each chain)
    inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                              w = rep(1,nrow(jags.data$y)),
                              gamma.u = runif (1),
                              phi.u = runif (1),
                              p.u = runif (1),
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
      "beta.gamma.area.t",
      "beta.gamma.coast",
      "beta.gamma.coast.t",
      "intercept.phi",
      "beta.phi.prec",
      "beta.phi.temp",
      "beta.phi.area",
      "beta.phi.area.t",
      "beta.phi.coast",
      "beta.phi.coast.t",
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
      "z",
      "SRexp",
      "w",
      "omega",
      "muZW",
      "Ntotal"
      
    )
    
    # usual jags to test
    #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
    #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
    #                                              chains=nc, 
    #                                              sample2save=((ni-nb)/nt), 
    #                                              nSaves=5, 
    #                                              burnin=nb, 
    #                                              thin=nt,
    #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
    
    
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
                                              parallel=F
    )
    
    save(samples_paleo_cynodontia_binomial,
         file = paste0 ("output/global/CMR_global_binomial_", i,"_more_iterations.RData"))
    
  })





# -------------------------------


# 500 genus imputed (only mammals)

# less   interactions
na <- 9000; nb <- 10000; ni <- 20000; nc <- 3; nt <- 10

# run
samples_paleo_cynodontia_binomial_run <- lapply (c (#"Non-mammaliaform cynodonts",
  #"Non-mammalian Mammaliaformes",
  "Mammalia" ), function (i) {
    
    
    # augment the dataset
    A0_input <- array (0,dim=c(500,
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
      nform = formations_per_interval,
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
      nform = formations_per_interval[function_stages(jags.data$y)],
      range = jags.data$range,
      lat_obs = jags.data$lat_obs,
      temp_obs = jags.data$temp_obs,
      # add time
      time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
      # time covs 
      temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
      precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
      area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
      coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
      area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
      coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
      
    ))
    
    
    
    
    # Set initial values
    # function inits (to be placed in each chain)
    inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                              w = rep(1,nrow(jags.data$y)),
                              gamma.u = runif (1),
                              phi.u = runif (1),
                              p.u = runif (1),
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
      "beta.gamma.area.t",
      "beta.gamma.coast",
      "beta.gamma.coast.t",
      "intercept.phi",
      "beta.phi.prec",
      "beta.phi.temp",
      "beta.phi.area",
      "beta.phi.area.t",
      "beta.phi.coast",
      "beta.phi.coast.t",
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
      "z",
      "SRexp",
      "w",
      "omega",
      "muZW",
      "Ntotal"
      
    )
    
    # usual jags to test
    #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
    #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
    #                                              chains=nc, 
    #                                              sample2save=((ni-nb)/nt), 
    #                                              nSaves=5, 
    #                                              burnin=nb, 
    #                                              thin=nt,
    #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
    
    
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
                                              parallel=F
    )
    
    save(samples_paleo_cynodontia_binomial,
         file = paste0 ("output/global/CMR_global_binomial_500spp_", i,".RData"))
    
  })


# ---------------------------------------------

# 800 genus imputed


# less   interactions
na <- 9000; nb <- 10000; ni <- 20000; nc <- 3; nt <- 10



# run
samples_paleo_cynodontia_binomial_run <- lapply (c (#"Non-mammaliaform cynodonts",
  #"Non-mammalian Mammaliaformes",
  "Mammalia" ), function (i) {
    
    
    # augment the dataset
    A0_input <- array (0,dim=c(800,
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
      nform = formations_per_interval,
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
      nform = formations_per_interval[function_stages(jags.data$y)],
      range = jags.data$range,
      lat_obs = jags.data$lat_obs,
      temp_obs = jags.data$temp_obs,
      # add time
      time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
      # time covs 
      temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
      precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
      area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
      coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
      area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
      coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
      
    ))
    
    
    
    
    # Set initial values
    # function inits (to be placed in each chain)
    inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                              w = rep(1,nrow(jags.data$y)),
                              gamma.u = runif (1),
                              phi.u = runif (1),
                              p.u = runif (1),
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
      "beta.gamma.area.t",
      "beta.gamma.coast",
      "beta.gamma.coast.t",
      "intercept.phi",
      "beta.phi.prec",
      "beta.phi.temp",
      "beta.phi.area",
      "beta.phi.area.t",
      "beta.phi.coast",
      "beta.phi.coast.t",
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
      "z",
      "SRexp",
      "w",
      "omega",
      "muZW",
      "Ntotal"
      
    )
    
    # usual jags to test
    #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
    #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
    #                                              chains=nc, 
    #                                              sample2save=((ni-nb)/nt), 
    #                                              nSaves=5, 
    #                                              burnin=nb, 
    #                                              thin=nt,
    #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
    
    
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
                                              parallel=F
    )
    
    save(samples_paleo_cynodontia_binomial,
         file = paste0 ("output/global/CMR_global_binomial_800spp_", i,".RData"))
    
  })


# =================================

# 1000 imputed genera

# less   interactions
na <- 9000; nb <- 10000; ni <- 20000; nc <- 3; nt <- 10


# run
samples_paleo_cynodontia_binomial_run <- lapply (c (#"Non-mammaliaform cynodonts",
  #"Non-mammalian Mammaliaformes",
  "Mammalia" ), function (i) {
    
    
    # augment the dataset
    A0_input <- array (0,dim=c(1000,
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
      nform = formations_per_interval,
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
      nform = formations_per_interval[function_stages(jags.data$y)],
      range = jags.data$range,
      lat_obs = jags.data$lat_obs,
      temp_obs = jags.data$temp_obs,
      # add time
      time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
      # time covs 
      temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
      precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
      area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
      coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
      area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
      coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
      
    ))
    
    
    
    
    # Set initial values
    # function inits (to be placed in each chain)
    inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                              w = rep(1,nrow(jags.data$y)),
                              gamma.u = runif (1),
                              phi.u = runif (1),
                              p.u = runif (1),
                              psi1 = runif (1)
                              
    )}
    
    ## Parameters to monitor
    ## long form
    ## long form
    params <- c(
      "w",
      "omega",
      "muZW",
      "Ntotal"
      
    )
    
    # usual jags to test
    #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
    #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
    #                                              chains=nc, 
    #                                              sample2save=((ni-nb)/nt), 
    #                                              nSaves=5, 
    #                                              burnin=nb, 
    #                                              thin=nt,
    #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
    
    
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
                                              parallel=F
    )
    
    save(samples_paleo_cynodontia_binomial,
         file = paste0 ("output/global/CMR_global_binomial_1000spp_", i,".RData"))
    
  })





# cynodonts, 500 imputed genera


# number of interactions for a global-scale analyses
na <- 15000; nb <- 30000; ni <- 50000; nc <- 3; nt <- 20


# run
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts"#,
  #"Non-mammalian Mammaliaformes",
  #"Mammalia" 
  ), function (i) {
    
    
    # augment the dataset
    A0_input <- array (0,dim=c(500,
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
      nform = formations_per_interval,
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
      nform = formations_per_interval[function_stages(jags.data$y)],
      range = jags.data$range,
      lat_obs = jags.data$lat_obs,
      temp_obs = jags.data$temp_obs,
      # add time
      time = as.vector(scale(bins$mid_ma[7:length(bins$mid_ma)][function_stages(jags.data$y)])),
      # time covs 
      temperature = as.vector(scale(time_covariates$temperature[function_stages(jags.data$y)])),
      precipitation = as.vector(scale(time_covariates$precipitation[function_stages(jags.data$y)])),
      area = as.vector(scale(time_covariates$area[function_stages(jags.data$y)])),
      coast = as.vector(scale(time_covariates$coastalLine[function_stages(jags.data$y)])),
      area.t = as.vector(scale(region_area[,1][function_stages(jags.data$y)])),
      coast.t = as.vector(scale(region_coastalLine[,1][function_stages(jags.data$y)]))
      
    ))
    
    
    
    
    # Set initial values
    # function inits (to be placed in each chain)
    inits <- function(){ list(z = ifelse (jags.data$y>1,1,1),
                              w = rep(1,nrow(jags.data$y)),
                              gamma.u = runif (1),
                              phi.u = runif (1),
                              p.u = runif (1),
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
      "beta.gamma.area.t",
      "beta.gamma.coast",
      "beta.gamma.coast.t",
      "intercept.phi",
      "beta.phi.prec",
      "beta.phi.temp",
      "beta.phi.area",
      "beta.phi.area.t",
      "beta.phi.coast",
      "beta.phi.coast.t",
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
      "z",
      "SRexp",
      "w",
      "omega",
      "muZW",
      "Ntotal"
      
    )
    
    # usual jags to test
    #samples_paleo_cynodontia_binomial <- saveJAGS(jags.data, inits, params, 
    #                                              modelFile="global_CMRmodel_matrix_covariates.txt",
    #                                              chains=nc, 
    #                                              sample2save=((ni-nb)/nt), 
    #                                              nSaves=5, 
    #                                              burnin=nb, 
    #                                              thin=nt,
    #                                              fileStub=paste ("output/global/CMR_global_binomial", i,sep="_"))
    
    
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
                                              parallel=F
    )
    
    save(samples_paleo_cynodontia_binomial,
         file = paste0 ("output/global/CMR_global_binomial_500spp_", i,".RData"))
    
  })




# end of the analysis





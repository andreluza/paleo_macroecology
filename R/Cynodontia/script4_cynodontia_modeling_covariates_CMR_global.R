# --------------------------------------

# 2 . modeling (GLOBAL SCALE)


# conditional on the first and last detection of the taxon


# saveJAGS help
# https://web.archive.org/web/20220408094124/http://mmeredith.net/blog/2018/Intro_saveJAGS.htm

## Global MCMC settings
######################
require(jagsUI)

# number of interactions for a global-scale analyses
na <- 15000; nb <- 30000; ni <- 50000; nc <- 3; nt <- 20

## short form for tests
na <- 600; nb <- 700; ni <- 1000; nc <- 3; nt <- 1
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
      #omega~dunif(0,1)
      omega ~ dbeta(0.001,1)
      
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
        beta.gamma.prec ~ dunif(-20,20)# dnorm(0,1/3^2)# precipitation
        beta.gamma.temp ~ dunif(-20,20)#dnorm(0,1/3^2)# temperature
        beta.gamma.area ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        beta.gamma.coast ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        
        # tropics
        beta.gamma.area.t ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        beta.gamma.coast.t ~ dunif(-20,20)#dnorm(0,1/3^2) # area
     
        # ----------------------
        #     Phi (persistence)
        # ----------------------
        # intercepts
        phi.u ~ dunif(0,1) # range persistence
        intercept.phi <- logit(phi.u) # intercept persistence
        # regression coeff
        beta.phi.prec ~ dunif(-20,20)#dnorm(0,1/3^2)# precipitation
        beta.phi.temp ~ dunif(-20,20)#dnorm(0,1/3^2)# temperature
        beta.phi.area ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        beta.phi.coast ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        
        # tropics
        beta.phi.area.t ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        beta.phi.coast.t ~ dunif(-20,20)#dnorm(0,1/3^2) # area
        
        ## set initial conditions for occupancy of each genus
        psi1 ~ dunif(0,1)
        
           
       ############      Model       #############
       
       # model for phi and gamma
       ### model dynamic parameters
        
        #for (t in 2:nint){
        for(t in 1:(nint-1)){
        
             # speciation
             logit(gamma[t]) <-  intercept.gamma + 
                                   #beta.gamma.prec*precipitation[t]+
                                   beta.gamma.temp*temperature[t]+
                                   beta.gamma.area*area[t]+
                                   #beta.gamma.coast*coast[t]+
                                   beta.gamma.area.t*area.t[t]#+
                                   #beta.gamma.coast.t*coast.t[t]
                                   
                                   
                                         
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
    
    ###  detection intercept
    # intercept    
    p.u ~ dunif(0,1) # range origination
    intercept.p <- logit(p.u) # intercept origination
    # regression coef
    beta.p.lat ~ dunif(-20,20)#dnorm(0,1/3^2)
    beta.p.range ~ dunif(-20,20)#dnorm(0,1/3^2)
    beta.p.time ~ dunif(-20,20)#dnorm(0,1/3^2)
    beta.p.temp ~ dunif(-20,20)#dnorm(0,1/3^2)
     
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
par(mfrow=c(1,1))
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
require(jagsUI)
samples_paleo_cynodontia_binomial_run <- lapply (c ("Non-mammaliaform cynodonts",
                                                    "Non-mammalian Mammaliaformes",
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
                        file = paste0 ("output/global/CMR_global_binomial1000sp_", i,".RData"))
                   
                 })



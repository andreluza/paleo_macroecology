# ------------------------------------------------------------

# design the model


sink("dyn_model_vectorized_covariates.txt")
cat("
   
    model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    
   
       for (g in 1:ngen) {
    
	    # intercept
          # gamma
          gamma.d[g] ~ dunif(0,1)
          intercept.gamma[g] <- logit(gamma.d[g])
          
          # phi
          phi.d[g] ~ dunif(0,1)
          intercept.phi[g] <- logit(phi.d[g])
          
          # regression coefficients
          
          # ----------------------
          #     Gamma (origination)
          # ----------------------
          
          # elevation
          beta.gamma.elev[g] ~ dnorm (mu.int.gamma[g],tau.mu.gamma[g])  
          mu.int.gamma[g] ~ dnorm(0, 0.001)
          tau.mu.gamma[g] <- 1/(sigma.int.gamma[g]*sigma.int.gamma[g])
          sigma.int.gamma[g] ~ dunif(0,10)
          
          # precipitation
          beta.gamma.prec[g] ~ dnorm (mu.int.gamma.prec[g],tau.mu.gamma.prec[g])  
          mu.int.gamma.prec[g] ~ dnorm(0, 0.001)
          tau.mu.gamma.prec[g] <- 1/(sigma.int.gamma.prec[g]*sigma.int.gamma.prec[g])
          sigma.int.gamma.prec[g] ~ dunif(0,10)
          
          # temperature
          beta.gamma.temp[g] ~ dnorm (mu.int.gamma.temp[g],tau.mu.gamma.temp[g])  
          mu.int.gamma.temp[g] ~ dnorm(0, 0.001)
          tau.mu.gamma.temp[g] <- 1/(sigma.int.gamma.temp[g]*sigma.int.gamma.temp[g])
          sigma.int.gamma.temp[g] ~ dunif(0,10)
          
          
          # ----------------------
          #     Phi (persistence)
          # ----------------------
          
          # precipitation
          beta.phi.prec[g] ~ dnorm (mu.int[g],tau.mu[g])  
          mu.int[g] ~ dnorm(0, 0.001)
          tau.mu[g] <- 1/(sigma.int[g]*sigma.int[g])
          sigma.int[g] ~ dunif(0,10)
          
          # temperature
          beta.phi.temp[g] ~ dnorm (mu.int.temp[g],tau.mu.temp[g])  
          mu.int.temp[g] ~ dnorm(0, 0.001)
          tau.mu.temp[g] <- 1/(sigma.int.temp[g]*sigma.int.temp[g])
          sigma.int.temp[g] ~ dunif(0,10)
          
          
        }
       
       
  
      ## set initial conditions
        ## priors for occupancy in time 1
        for (i in 1:nsites) {
        
            for (g in 1:ngen) {
         
              psi1[i,g] ~ dbeta(a,b)
             
            }
         }
    
    # Specify the hyperparameters for the Beta distribution
    a <- 10   # probability of success
    b <- 90   # probability of failures
    
    
    
    ############      Model       #############
    
    
    for (i in 1:nsites) {
    
      for (g in 1:ngen) {
    
          z[i,1,g]~dbern(psi1[i,g]) # occupancy status initialization
    
              for (t in 2:nint){
            
               ### model dynamic parameters
                logit(gamma[i,t,g]) <-  intercept.gamma[g] + 
                                        beta.gamma.elev[g]*elevation[i,t]+
                                        beta.gamma.prec[g]*precipitation[i,t]+
                                        beta.gamma.temp[g]*temperature[i,t]
        
                logit(phi[i,t,g]) <-  intercept.phi[g] + 
                                      beta.phi.prec[g]*precipitation[i,t]+
                                      beta.phi.temp[g]*temperature[i,t]
        
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[i,t,g] <- z[i,t-1,g] *  phi[i,t,g] + ### if occupied, p of not getting extinct/persist in the next time
                          (1-z[i,t-1,g]) *  gamma[i,t,g] ###  if not occupied, p of originate in the next time
                
                # realized occurrence
		    z[i,t,g] ~ dbern(muZ[i,t,g])
    
        }#t
        
      } #g
    
    }#i
    
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    
    
        # Priors for detection probability
        
        ### temperature effect on detection
        for (g in 1:ngen) {
           
             alpha.p[g] ~ dunif(0,1)
             intercept.p[g] <- logit(alpha.p[g])
             alpha1.temp[g] ~ dnorm(0, 0.001)
           
           }
        

     ############      Model       #############
     
     
     
     # observation submodel: detection probability based on depth and videos
     for (k in 1:nobs) { ## loop over observations
                            
               y [k] ~ dbern(muY[site[k],form[k],int[k],gen[k]])
               muY [site[k],form[k],int[k],gen[k]] <- z[site[k],int[k],gen[k]] * p[k]
                         
               # model
               logit(p[k])<- intercept.p[gen[k]]+ 
                                      alpha1.temp[gen[k]]*tempObs[k] 
                
             }
      
    
    # -----------------------------------------------
    
    ## derived parameters
    # number of gen per interval
    for (t in 1:nint) {
        Ngen[t]<-sum(z[,t,])
    }
    
    ## number of genera per site
    for (i in 1:nsites) {
        Ngen_site[i]<-sum(z[i,,])
    }
    

    ## average persistence and origination
   #for (g in 1:ngen) {
   #  avphi[g] <- mean(phi[,2:nint,g])
   #  avgamma[g]<- mean(gamma[,2:nint,g])
   #}
   #
   ### turnover (proportional gain or loss)
   #for (t in 2:nint) {  
   #  
   #    propcH [t] <-(sum(z[,t,])-sum (z[,t-1,]))/sum(z[,t-1,]) 
   #  
   #}
    
    ## equilibrium occupancy (which gen decline or increase over time)
    #for (g in 1:ngen) {
    #
    #    psi.eq[g] <- mean(gamma[2:nint,g])/(mean(gamma[2:nint,g])+mean(1-phi[2:nint,g])) # Equilibrium occupancy
    #    #psi.eq[g] <- gamma[g]/(gamma[g]+(1-phi[g])) # Equilibrium occupancy
    #}
    #
    ## relative extinction rate (μ/λ; Rabosky 2018) of each time
    #for (t in 2:nint) { 
    # for (g in 1:ngen) {
    #
    #    #RER[g] <- (1-phi[g])/gamma[g]
    #    RER[t,g] <- (1-phi[t,g])/gamma[t,g]
    #
    #  }
	  #}
    #
    ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
    #for (t in 2:nint) { 
    # for (g in 1:ngen) {
    #
    #    #R0[g] <- (1-phi[g])-gamma[g]
    #    R0[t,g] <- (1-phi[t,g])-gamma[t,g]
    #
     #}
    #}
    
        
    }## end of the model
    
    
    

    
    ",fill = TRUE)
sink()


# --------------------------------------


# load packages
source("R/packages.R")


# --------------------------------------


# synapsida level
# load data
load(here ("processed_data","table_data_array_synapsida.RData"))
load(here ("processed_data","table_space.RData"))
load(here ("processed_data", "table_naive_synapsida.RData"))
load (here ("processed_data", "site_covs.RData"))


# scale variables
elevation<-(scale (site_covs$elevation)) # scaled elevation
temperature <- (scale (site_covs$temp)) # scaled temp
precipitation <-  (scale (site_covs$prec))


# bins in data
elevation <- elevation [, which(colnames(elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(elevation)==bins [which(bins$bin %in% time_bins),"interval_name"]
temperature <- temperature [, which(colnames(temperature) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(temperature)==bins [which(bins$bin %in% time_bins),"interval_name"]
precipitation <- precipitation [, which(colnames(precipitation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(precipitation)==bins [which(bins$bin %in% time_bins),"interval_name"]

# adjust dimensions (rm sites missnig in occ data --- see how to estimate data for them )
temperature <- (as.matrix(temperature [match( rownames(table_data_basis[[1]]),rownames(temperature)),]))
precipitation <- as.matrix(precipitation [match( rownames(table_data_basis[[1]]),rownames(precipitation)),])
elevation <- as.matrix(elevation [match( rownames(table_data_basis[[1]]),rownames(elevation)),])

# data onto 0,1
table_data_long_df$det[which(table_data_long_df$det >0)]<-1



# subset
# table_data_long_df <- table_data_long_df[which(table_data_long_df$taxon %in% seq(1,10)),]

## bundle data
str(jags.data <- list(y = table_data_long_df$det,
                      gen = table_data_long_df$taxon,
                      int=table_data_long_df$int,
                      form = table_data_long_df$form,
                      #lith = table_data_long_df$lith2,
                      site = table_data_long_df$site,
                      nsites = nrow(temperature),
                      nint= length(unique(table_data_long_df$int)), 
                      ngen = max(unique(table_data_long_df$taxon)),
                      #nlith= length(unique(table_data_long_df$lith2)),
                      tempObs = as.vector(scale (table_data_long_df$temp)),
                      nobs=nrow(table_data_long_df),
                      elevation = unname (elevation),
                      temperature = unname(temperature),
                      precipitation = unname(precipitation)
                      
                      )
)

# Set initial values

require(reshape)
table_inits <- lapply (seq (1,max(table_data_long_df$taxon)), function (i)
  
        cast(table_data_long_df[which (table_data_long_df$taxon == i),], site ~ int, 
            max, 
            value = "det",
            fill=0,
            margins="taxon")[,-1]
        )




# melt the list
zst<-array( unlist( table_inits ) , dim = c(nrow(table_inits[[1]]),
                                            ncol(table_inits[[1]]), 
                                            length((table_inits)))
)


# inits for psi1
psi1 <- (apply (zst, c(1,3),mean))

# function inits (to be placed in each chain)
inits <- function(){ list(z = zst,
                          psi1 = psi1)}

## Parameters to monitor
## long form
params <- c(
  
  "gamma", "phi","p",
  "intercept.gamma", 
  "beta.gamma.elev",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "intercept.phi", 
  "beta.phi.prec",
  "beta.phi.temp", 
 "alpha.p",
 "intercept.p",
 "alpha1.temp",
 "psi1",	
 "Ngen",
 "Ngen_site"#, 
 #"avphi",
 #"avgamma",
 #"propcH"
  
  
)

## MCMC settings
######################
## short form
na <- 6000; nb <- 7000; ni <- 15000; nc <- 3; nt <- 8
na <- 250; nb <- 500; ni <- 1000; nc <- 3; nt <- 1

# MCMC runs
# models

samples_paleo_synapsida <- jags (data = jags.data, 
                       parameters.to.save = params, 
                       model.file = "dyn_model_vectorized_covariates.txt", 
                       inits = inits, 
                       n.chains = nc, 
                       n.thin = nt, 
                       n.iter = ni, 
                       n.burnin = nb, 
                       DIC = T,  
                       #n.cores=nc,
                       parallel=F
)



# save
save (samples_paleo_synapsida,file = here ("output","samples_paleo_synapsida.RData"))


# try WinBUGS
samples_paleo_synapsida <- bugs(data = jags.data, 
     parameters.to.save = params, 
     model.file = "dyn_model_vectorized_covariates.txt", 
     inits = inits, 
     n.chains = nc, 
     n.thin = nt, 
     n.iter = ni, 
     n.burnin = nb, 
     codaPkg=F, 
     DIC=TRUE, 
     debug=T,
     bugs.directory="C:/Program Files/WinBUGS14/", program= "WinBUGS")


# save jags output
#dir.create ("output")
save (samples_paleo_synapsida,file = here ("output","samples_paleo_synapsida.RData"))




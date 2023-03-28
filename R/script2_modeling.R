
# ------------------------------------------------------------

# design the model

sink("dyn_model_vectorized.txt")
cat("
    model {
    
    #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Priors
        
    for (t in 2:nint) {
    
       for (g in 1:ngenus) {
    
        
          ## colonization (origination)
          gamma [t,g] ~ dunif (0,1)
          ## extinction
          phi [t,g]~ dunif (0,1)
        
        }
       
       }
  
    
    ## set initial conditions
    ## priors for occupancy in time 1
     for (g in 1:ngenus) {
     
      psi1 [g] ~ dunif (0,1)
     
     }
    
    ############ Model ########################
    for (g in 1:ngenus) {
    
        z[1,g]~dbern(psi1[g]) # occupancy status initialization
    
            for (t in 2:nint){
            
              # model likelihood
              ### modeling dynamics conditional on previous time realized occurrence z
              muZ[t,g] <- z[t-1,g] * (1-phi[t,g]) + ### if occupied, p of not getting extinct in the next time
                          (1-z[t-1,g]) * gamma[t,g] ###  if not occupied, p of getting colonized in the next time
                        
            z[t,g] ~ dbern(muZ[t,g])
    
        }#t
    }#i
    
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    # priors
    ## detection
    #for (t in 1:nint) {
    
      for (g in 1:ngenus) {
    
            p[g] ~ dunif (0,1)

      }
    #}

    ### model
    for (k in 1:nobs){
        
        y [k] ~ dbern (muY[form[k],int[k],genus[k]])
        muY [form[k],int[k],genus[k]] <- z [int[k],genus[k]] * p[genus[k]]
        
    }
    
    
    # -----------------------------------------------
    
    ## derived parameters
    # number of genus per interval
    for (t in 1:nint) {
        Ngen[t]<-sum(z[t,])
    }
    
    # average extinction and origination
    for (g in 1:ngenus) {
      avphi[g] <- mean(phi[2:nint,g])
      avgamma[g]<- mean(gamma[2:nint,g])
    }
      
    # turnover (proportional gain or loss)
    for (t in 2:nint) {  
      
        propcH [t] <-(sum (z[t-1,]) - sum(z[t,]))/sum(z[t-1,]) 
      
    }
    
    # equilibrium occupancy (which genus decline or increase over time)
    for (g in 1:ngenus) {
    
        psi.eq[g] <- mean(gamma[2:nint,g])/(mean(gamma[2:nint,g])+mean(1-phi[2:nint,g])) # Equilibrium occupancy
    
    }
    
    
    }## end of the model
    
    ",fill = TRUE)
sink()


# require packs
require(here); require(jagsUI)

# load data
load(here ("output","table_data_array.RData"))

## bundle data
str(jags.data <- list(y = table_data_long_df$det,
                      genus = table_data_long_df$taxon,
                      int=table_data_long_df$int,
                      form = table_data_long_df$form, 
                      nint= ncol(table_data_basis), 
                      ngenus = max(unique(table_data_long_df$taxon)),
                      nobs=nrow(table_data_long_df)
))

# Set initial values
zst <- matrix(1,
              nrow= 17,
              ncol=max(unique(table_data_long_df$taxon)))
inits <- function(){ list(z = zst)}

## Parameters to monitor
## long form
params <- c(
  
  "gamma", "phi","p","muZ",
  "propcH", "avgamma", "avphi",
  "Ngen",'psi.eq'
  
)

## MCMC settings
######################
## short form
na <- 6000; nb <- 7000; ni <- 15000; nc <- 3; nt <- 8
#na <- 50; nb <- 60; ni <- 100; nc <- 3; nt <- 1

# MCMC runs
samples_paleo <- jags (data = jags.data, 
                       parameters.to.save = params, 
                       model.file = "dyn_model_vectorized.txt", 
                       inits = inits, 
                       n.chains = nc, 
                       n.thin = nt, 
                       n.iter = ni, 
                       n.burnin = nb, 
                       DIC = T,  
                       n.cores=nc,
                       parallel=T
)

# save jags output
#dir.create ("output")
save (samples_paleo,file = here ("output","samples_paleo.RData"))

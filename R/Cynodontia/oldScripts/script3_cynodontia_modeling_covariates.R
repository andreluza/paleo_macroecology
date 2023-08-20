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
    
       # intercept
       # gamma
       gamma.d ~ dunif(0,1)
       intercept.gamma <- logit(gamma.d)
          
       # phi
       phi.d ~ dunif(0,1)
       intercept.phi <- logit(phi.d)
          
       # regression coefficients
          
      # ----------------------
      #     Gamma (origination)
      # ----------------------
      
      # elevation
      beta.gamma.elev ~ dnorm (mu.int.gamma,tau.mu.gamma)  
      mu.int.gamma ~ dnorm(0, 0.001)
      tau.mu.gamma <- 1/(sigma.int.gamma*sigma.int.gamma)
      sigma.int.gamma ~ dunif(0,10)
      
      # precipitation
      beta.gamma.prec ~ dnorm (mu.int.gamma.prec,tau.mu.gamma.prec)  
      mu.int.gamma.prec ~ dnorm(0, 0.001)
      tau.mu.gamma.prec <- 1/(sigma.int.gamma.prec*sigma.int.gamma.prec)
      sigma.int.gamma.prec ~ dunif(0,10)
      
      # temperature
      beta.gamma.temp ~ dnorm (mu.int.gamma.temp,tau.mu.gamma.temp)  
      mu.int.gamma.temp ~ dnorm(0, 0.001)
      tau.mu.gamma.temp <- 1/(sigma.int.gamma.temp*sigma.int.gamma.temp)
      sigma.int.gamma.temp ~ dunif(0,10)
      
      # latitude
      beta.gamma.lat ~ dnorm (mu.int.gamma.lat,tau.mu.gamma.lat)  
      mu.int.gamma.lat ~ dnorm(0, 0.001)
      tau.mu.gamma.lat <- 1/(sigma.int.gamma.lat*sigma.int.gamma.lat)
      sigma.int.gamma.lat ~ dunif(0,10)
      
      # ----------------------
      #     Phi (persistence)
      # ----------------------
      
      # precipitation
      beta.phi.prec ~ dnorm (mu.int,tau.mu)  
      mu.int ~ dnorm(0, 0.001)
      tau.mu <- 1/(sigma.int*sigma.int)
      sigma.int ~ dunif(0,10)
      
      # temperature
      beta.phi.temp ~ dnorm (mu.int.temp,tau.mu.temp)  
      mu.int.temp ~ dnorm(0, 0.001)
      tau.mu.temp <- 1/(sigma.int.temp*sigma.int.temp)
      sigma.int.temp ~ dunif(0,10)
      
      # latitude
      beta.phi.lat ~ dnorm (mu.int.lat,tau.mu.lat)  
      mu.int.lat ~ dnorm(0, 0.001)
      tau.mu.lat <- 1/(sigma.int.lat*sigma.int.lat)
      sigma.int.lat ~ dunif(0,10)
      
    
      ## set initial conditions
        ## priors for occupancy in time 1
        for (i in 1:nsites) {
        
            for (g in 1:ngen) {
         
              psi1[i,g]~dunif(0,1)#dbeta(a,b)
             
            }
         }
    
    # Specify the hyperparameters for the Beta distribution
    #a <- 10   # probability of success
    #b <- 90   # probability of failures
    
    
    ############      Model       #############
    
    
    for (i in 1:nsites) {
    
      for (g in 1:ngen) {
    
          z[i,1,g]~dbern(psi1[i,g]) # occupancy status initialization
    
              for (t in 2:nint){
            
               ### model dynamic parameters
                logit(gamma[i,t,g]) <-  intercept.gamma + 
                                        beta.gamma.elev*elevation[i,t]+
                                        beta.gamma.prec*precipitation[i,t]+
                                        beta.gamma.temp*temperature[i,t]+
                                        beta.gamma.lat*lat[i]
        
                logit(phi[i,t,g]) <-  intercept.phi + 
                                      beta.phi.prec*precipitation[i,t]+
                                      beta.phi.temp*temperature[i,t]+
                                      beta.phi.lat*lat[i]
        
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
        
        ###  detection intercept
        alpha.p ~ dunif(0,1)
        intercept.p <- logit(alpha.p)
        
        # temp effect on detection
        beta.p.temp ~ dnorm (mu.p.temp,tau.mu.p.temp)  
        mu.p.temp ~ dnorm(0, 0.001)
        tau.mu.p.temp <- 1/(sigma.int.p.temp*sigma.int.p.temp)
        sigma.int.p.temp ~ dunif(0,10)
      
        # latitude effect on detection
        beta.p.lat ~ dnorm (mu.p.lat,tau.mu.p.lat)  
        mu.p.lat ~ dnorm(0, 0.001)
        tau.mu.p.lat <- 1/(sigma.int.p.lat*sigma.int.p.lat)
        sigma.int.p.lat ~ dunif(0,10)
      
      
        
        

     ############      Model       #############
     
     
     
     # observation submodel: detection probability based on depth and videos
     for (k in 1:nobs) { ## loop over observations
                            
               y [k] ~ dbern(muY[site[k],form[k],int[k],gen[k]])
               # detection conditional on true occurrence 
               muY [site[k],form[k],int[k],gen[k]] <- z[site[k],int[k],gen[k]] * p[k]
                         
               # model
               logit(p[k])<- intercept.p+
                              beta.p.temp*tempObs[k]+
                              beta.p.lat*latObs[k]
                
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
load(here ("processed_data","table_data_array_cynodontia.RData"))
load(here ("processed_data","table_space.RData"))
load(here ("processed_data", "table_naive_cynodontia.RData"))
load (here ("processed_data", "site_covs.RData"))


#table_data_long_df_original <- table_data_long_df
#table_data_long_df <- table_data_long_df_original

# scale variables
elevation<-(scale (site_covs$elevation)) # scaled elevation
temperature <- (scale (site_covs$temp)) # scaled temp
precipitation <-  (scale (site_covs$prec))
lat <- (scale (site_covs$paleolat))

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
lat <- as.matrix(lat [match( rownames(table_data_basis[[1]]),site_covs$cell_ID),])



## bundle data
str(jags.data <- list(y = table_data_long_df$det,
                      gen = table_data_long_df$taxon_code,
                      int= table_data_long_df$int,
                      form = table_data_long_df$form,
                      #lith = table_data_long_df$lith2,
                      site = table_data_long_df$site,
                      nsites = nrow(temperature),
                      nint= length(unique(table_data_long_df$int)), 
                      ngen = max(unique(table_data_long_df$taxon_code)), # taxon
                      #nlith= length(unique(table_data_long_df$lith2)),
                      tempObs = as.vector(scale (table_data_long_df$temp)),
                      latObs=as.vector(scale (table_data_long_df$paleolat)),
                      nobs=nrow(table_data_long_df),
                      elevation = unname (elevation),
                      temperature = unname(temperature),
                      precipitation = unname(precipitation),
                      lat = as.vector(unname(lat))
                      
                      )
)

# Set initial values

max_row <- lapply (table_naive_cynodontia, function(i)
  
  lapply (i, function (k)
    
    apply (k,1,max,na.rm=T)
    
  )
)


max_rowB<- lapply (max_row, function (i)
  
  do.call(cbind,i)
  
)

# melt the list
zst<-array( unlist( max_rowB ) , dim = c(nrow(max_rowB[[1]]),
                                         ncol(max_rowB[[1]]), 
                                         length((max_rowB)))
)

zst[(zst>=1)] <- 1
zst[is.infinite(zst)] <- 0

# inits for psi1
#psi1 <- (apply (zst, c(1,3),mean,na.rm=T))

# function inits (to be placed in each chain)
inits <- function(){ list(z = zst)}

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
"beta.phi.lat",
"intercept.p",
"beta.p.temp",
"beta.p.lat",

  "gamma", "phi","p",

 "Ngen",
 "Ngen_site",
 "z"


  
  
)

## MCMC settings
######################
## short form
na <- 6000; nb <-10000; ni <- 20000; nc <- 3; nt <- 50
#na <- 250; nb <- 500; ni <- 1000; nc <- 3; nt <- 1

# MCMC runs
# models

samples_paleo_cynodontia <- jags (data = jags.data, 
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
save (samples_paleo_cynodontia,file = here ("output","samples_paleo_cynodontia_covariates.RData"))


# time trend in SR
data.frame (ntaxa=rowSums(apply (zst, c(2,3),max)),
            est=rowSums(apply (samples_paleo_cynodontia$mean$z,c(2,3), max)),
            int = bins$mid_ma [match (time_bins,bins$bin)]) %>%
  
  melt (id.vars="int") %>%
  
  ggplot (aes(x=int,y=value, fill=variable))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  scale_x_reverse("Age (Ma)") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 70,270), 
    ylim = c(0, 200),
    pos = list("b", "b"),
    size = list(4, 6),
    abbrv = list(TRUE, FALSE)
  ) 


# time trend in gamma and phi
data.frame (cbind (apply (samples_paleo_cynodontia$mean$gamma,2, mean),
                   int =  bins$mid_ma [match (time_bins,bins$bin)],
                   var = "gamma")) %>%
  
  melt (id.vars=c("int","var")) %>%
  
  rbind(data.frame (cbind (apply (samples_paleo_cynodontia$mean$phi,2, mean),
                           int =  bins$mid_ma [match (time_bins,bins$bin)],
                           var = "phi")) %>%
          
          melt (id.vars=c("int","var")))%>%
  mutate(int=as.numeric(int),
         value = as.numeric(value),
         var = as.factor (var),
         variable = as.factor (variable)) %>%
  filter (is.na(value) !=T) %>%
  
  ggplot (aes(x=int,
              y=value, 
              fill=var,
              colour=var,
              group=var,
              shape = as.factor (as.numeric(var))))+
  geom_point(aes (col= variable))+
  geom_smooth(alpha=0.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_reverse("Age (Ma)") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 70,270), 
    ylim = c(0, 1),
    pos = list("b", "b"),
    size = list(4, 6),
    abbrv = list(TRUE, FALSE)
  ) 


# plot 

# poison smooth 
binomial_smooth <- function(...) {
  geom_smooth(method = "gam", method.args = list(family = "binomial"), ...)
  # geom_smooth(method = "glm", method.args = list(family = "poisson"), ...) # glm option
}


# persistence probability

samples_paleo_cynodontia$mean$gamma %>%
  melt() %>%
  ggplot (aes(x=X2,y=value,group=X1,fill=X1))+
  geom_point(alpha=0.4)+
  binomial_smooth(alpha=0.01,
                  formula = y ~ splines::ns(x, 2),
                  #formula = y ~ splines::ns(x, 3),
                  se=T,
                  size=1,
                  colour = "black")+
  ylab ("Origination probability")+
  xlab ("Interval")

# persistence probability

samples_paleo_cynodontia$mean$phi %>%
  melt() %>%
  ggplot (aes(x=X2,y=value,group=X1,fill=X1,colour=X1))+
  geom_point(alpha=0.4)+
  binomial_smooth(alpha=0.01,
                  formula = y ~ splines::ns(x, 2),
                  #formula = y ~ splines::ns(x, 3),
                  se=T,
                  size=1,
                  colour = "black")+
  ylab ("Persistence probability")+
  xlab ("Interval")+
  theme_bw()




# -------------------------------------------

# simulations of the dynamic model

# could the dynamic model recover good and bad sites in terms of origination?

# simular dados á mão

# processo: 
# imagine varias regioes no mundo, e entao tem lugares que o animal nunca se extingue
# e antes dele se estinguir por completo do planeta, ele deve começar a reduzir em alguns locais (sites)
# ou seja se extingue em um local, aparece de novo, se extingue de novo, se extingue e demora pra recolonizar , e não aparece mais,
# e isso acontece em um outro lugar
# cada vez menos 1s com o passar do tempo
# assm, fica cada vez mais dificil de ser encontrado 

# encontrar uma espécie comum que vai desaparecendo aos poucos do mundo para comparar

require("AHMbook")
require(dplyr)
require(ggplot2)

set.seed(123)
nyears <- 30 # number of time bins (stages)

# values of gamma and phi
# range of values for gamma (same for phi)
# constant probabilities over time (all rates reach zero)
gamma <- list(c(0.02,0), # suboptimal site
              c(0.05,0), # suboptimal site
              c(0.10,0), # suboptimal site
              c(0.15,0), # suboptimal site
              c(0.40,0), # intermediate origination site
              c(0.45,0), # intermediate origination site
              c(0.5,0), # intermediate origination site
              c(0.55,0), # intermediate origination site
              c(0.75,0), # optimal site 
              c(0.82,0), # optimal site
              c(0.9,0), # optimal site
              c(0.98,0)) # optimal site
phi <- c(0.1,0.05) # small persistence probability
psi1<-0.15


# simulate dataset
sim_data <- lapply (gamma, function (i) 
  
  # function of the AHM Book
  simDynocc(nsites= 1, # 1 site each time
            nyears=nyears,
            nsurveys=1,
            mean.psi1=psi1,
            range.phi = phi,
            range.gamma =i,
            range.p=c(1,1), # perfect detection
            #beta.Xp=2,
            #range.beta1.survey = c(-2,2),
            show.plots = F)
)


# melt data
sim_data <- t(sapply (sim_data, "[[", "z"))
sum(rowSums(sim_data)>0)/length(gamma)

#number of sites occupied over time
plot (colSums(sim_data)/length(gamma), type="b")

# design the dynamic model to detect probabilities
sink("dyn_model_vectorized_no_detection.txt")
cat("
   
     model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    for(i in 1:nsites) {
      gamma[i]~dunif(0,1)
      phi[i]~dunif(0,1)
    }
    
    psi1~dunif(0,1)     
         
    # Initial values
    for (i in 1:nsites) {
      p[i,1] <- psinit    # Initial occupancy probability
      }
    

   ############      Model       #############
   for (i in 1:nsites) {
      
          y[i,1] ~ dbern(psi1) # Initial occupancy probability
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time observed occurrence y
                
                muZ[i,t] <- y[i,t-1] *  phi[i] + ### if occupied, p of not getting extinct/persist in the next time
                          (1-y[i,t-1]) *  gamma[i] ###  if not occupied, p of originate in the next time
                
        # realized occurrence
		    p[i,t] ~ dbern(muZ[i,t])      # Occupancy model
		    y[i,t] ~ dbern(p[i,t])        # Observation model
    
        }#t
      
     }#i
  
    
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()

## bundle data
str(jags.data <- list(y = sim_data,
                      nsites = dim(sim_data)[1],
                      nint= dim(sim_data)[2],
                      psinit=psi1)
)


## Parameters to monitor
## long form
params <- c(
  
  "psi1",
  "gamma", 
  "phi",
  "p",
  "y"
  
)

# Set initial values
# function inits (to be placed in each chain)
inits <- function(){ list(z = sim_data)}

## MCMC settings
######################
## short form
#na <- 6000; nb <- 7000; ni <- 15000; nc <- 3; nt <- 8
na <- 1000; nb <- 5000; ni <- 10000; nc <- 3; nt <- 10

# MCMC runs
# models
require(jagsUI)
samples_sims_no_det <- jags (data = jags.data, 
                             parameters.to.save = params, 
                             model.file = "dyn_model_vectorized_no_detection.txt", 
                             inits = NULL, 
                             n.chains = nc, 
                             n.thin = nt, 
                             n.iter = ni, 
                             n.burnin = nb, 
                             DIC = T,  
                             #n.cores=nc,
                             parallel=F
)

#number of sites occupied over time
plot (seq(1,nyears), colSums(sim_data)/length(gamma), type="b",ylim=c(0,1)) # observed
lines (seq(1,nyears), # estimated
       (samples_sims_no_det$mean$ISS/length(gamma)), type="b",col="red")
lines (seq(1,nyears), # estimated
       (apply (samples_sims_no_det$sims.list$ISS,2,quantile,0.025)/length(gamma)), type="b",col="red")
lines (seq(1,nyears), # estimated
       (apply (samples_sims_no_det$sims.list$ISS,2,quantile,0.975)/length(gamma)), type="b",col="black")



# plot the data
df_res <- do.call(rbind,gamma)
# bind estimates
df_res <- cbind (df_res,
                 sites = seq(1,nrow(df_res)) ,
                 gamma = samples_sims_no_det$mean$gamma,
                 uci = apply(samples_sims_no_det$sims.list$gamma,2,quantile,0.975),
                 lci = apply(samples_sims_no_det$sims.list$gamma,2,quantile,0.025))
colnames(df_res)[1:2]<-c("inf_range","sup_range")

# plot res simulations

ggplot (data=data.frame (df_res),
        aes (x=sites,y=gamma)) + 
  
  geom_point(col="red", size=2)+
  
  geom_errorbar(aes (ymin=lci,ymax= uci),col="red",width=0.2,linewidth=1.3)+
  
  geom_errorbar(aes (ymax=inf_range,ymin= sup_range),width=0.2)+
  
  theme_bw() + 
  
  ylab ("Origination probability") + 
  
  xlab ("Sites") + 
  
  scale_x_continuous(breaks = seq(1, 12, by = 1))


# plot res simulations
# bind estimates

df_res_phi <- cbind (sites = seq(1,nrow(df_res)) ,
                     phi = samples_sims_no_det$mean$phi,
                     uci = apply(samples_sims_no_det$sims.list$phi,2,quantile,0.975),
                     lci = apply(samples_sims_no_det$sims.list$phi,2,quantile,0.025))

# plot
ggplot (data=data.frame (df_res_phi),
        aes (x=sites,y=phi)) + 
  
  geom_point(col="red", size=2)+
  
  geom_errorbar(aes (ymin=lci,ymax= uci),col="red",width=0.2,linewidth=1.3)+
  
  #geom_errorbar(aes (ymin= min(phi), ymax=max(phi)),width=0.2)+
  geom_rect(aes(xmin = 1, xmax = 12, ymin = 0.05, ymax = 0.1),alpha=0.1) +
  
  theme_bw() + 
  
  ylab ("Persistence probability") + 
  
  xlab ("Sites") + 
  
  scale_x_continuous(breaks = seq(1, 12, by = 1))


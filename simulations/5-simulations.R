# simulations of the dynamic model
# mixing Royle & Kery 2021 and Kocsics et al. MEE (paleoTree)
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00223.x

# simular dados á mão
# processo
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
# constant probabilities over time
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


# other simulations

library(paleotree) # load the library
library (here)
library (palaeoverse)
library(reshape)
require(jagsUI)

# set seed
set.seed(1)

# sample size and system characteristics
nsites <- 10 # n sites or cells
nintervals <- 33 # geological intervals
nsurveys <- 10 # formations

# set parameter values to simulate data
psi1 <- 0.05 # probability of initial occupancy (prop occ cells) low value at the initial time
phi <- 0.1 # persistence probability
gamma <- 0.1 # origination probability
p <-0.05 # detection probability
(psi.eq <- gamma/(gamma+(1-phi)))

# binning time intervals
bins <- time_bins(interval = c("Permian", "Cretaceous"), rank = "stage")

# other params of paleoTree
startTaxaN <- 1 # number of taxa to start sim
nruns <- 1 # simulation runs
timeRange <- c (min(bins$min_ma),
                max(bins$max_ma)) # range of time
totalTaxa <- c(1, 50) # min and max n taxa
nExtant <- 0 # number of extant taxa


## Simulate some fossil ranges (stratigraphic range) with simFossilRecord
# it will comprise the true occurrence
record_taxa <- lapply (seq_len (nsites), function (i) { 
  
  record<-simFossilRecord(
    p = gamma, q = phi, 
    startTaxa = startTaxaN,
    nruns = nruns,
    totalTime =timeRange ,
    nTotalTaxa = totalTaxa,
    nExtant = nExtant
  )
  
  # list into dataframe
  taxa <- fossilRecord2fossilTaxa(record)
  
  # bind detection (range)
  taxa <- cbind (taxa, 
                 det=1)

  # discreticing intervals
  # create empty cols to receive these values
  taxa <- cbind (taxa, orig.int=NA)
  taxa <- cbind (taxa, ext.int=NA)
  
  # bind intervals (run the function for each taxon (row))
  taxa <- lapply (seq(1,nrow(taxa)), function (i) {
  
    # bind origination interval
    taxa[i,"orig.int"] <- bins [bins$max_ma >= taxa[i,"orig.time"] & bins$min_ma <= taxa[i,"orig.time"],"bin"] # origination time
    
    # bind extinction interval
    taxa[i,"ext.int"] <- bins [bins$max_ma >= taxa[i,"ext.time"] & bins$min_ma <= taxa[i,"ext.time"],"bin"] # origination time
    
    
    
    ; # retun
    taxa[i,]
    
    })
  # melt
  taxa<- do.call(rbind,taxa)
  taxa <-   melt (data.frame(taxa),id.vars = c("taxon.id",
                         "ancestor.id",
                         "orig.time",
                         "ext.time",
                         "still.alive",
                         "looks.like",
                         "det"))
  
  # subset by taxon
  # apply to all them
  taxon_range <- lapply (unique(taxa$taxon.id), function (i) {
    # subset 
    taxon_subset <- taxa [taxa$taxon.id == i,]
    
    # add missing intervals
    to_add <- seq (min(taxon_subset [,"value"]),
                   max(taxon_subset [,"value"]))
    to_add<-to_add[which(to_add %in% taxon_subset[,"value"] ==F)]
    rep_row <- taxon_subset[rep(1, each = length(to_add)), ]
    rep_row[,"value"] <- to_add
    # replicate row to add
    taxon_subset<-rbind (taxon_subset,
                        rep_row) 
    # order()
    taxon_subset <- taxon_subset[order(taxon_subset$value),]
    ; #return
    taxon_subset
    
  })
  # melt
  taxon_range <-do.call(rbind,taxon_range)
  
  
  # table of true occupancy
  table_z <- cast (data=taxon_range,
        formula = taxon.id~value)
  rownames(table_z) <- table_z$taxon.id; table_z<-table_z[,-1]
  table_z[table_z>1]<-1
  
  # input missing intervals & taxa
  missing_intervals <- matrix (0,
                               nrow = nrow (table_z),
                               ncol = length(seq(1,nintervals) [which(seq(1,nintervals) %in% as.numeric(colnames(table_z)) ==F)]),
                               dimnames = list (rownames(table_z),
                                                seq(1,nintervals) [which(seq(1,nintervals) %in% as.numeric(colnames(table_z)) ==F)])
  )
  # bind
  table_z <-cbind (table_z,
                   missing_intervals)
  # order
  table_z<-table_z[,order(as.numeric(colnames(table_z)))]
  
  
  
  ; # return our simulated stratigraphic ranges
  
  table_z
  
  })


# different number of species were simulated across sites
# bind taxa
max_spp<-max((unlist(lapply (record_taxa,nrow))))

# bind
record_taxa_z <- lapply (record_taxa, function (i) {
  
  # what to bind?
   missing_taxa <- matrix (0, 
                              ncol=ncol (i),
                              nrow = length(seq (1,max_spp)[which(seq (1,max_spp) %in% as.numeric(rownames(i)) ==F)]),
                              dimnames=list(seq (1,max_spp)[which(seq (1,max_spp) %in% as.numeric(rownames(i)) ==F)],
                                            colnames(i)
                                            ))
      
      
      i <- rbind (i,
                  missing_taxa)
      
      ;
      
      i
})

# z matrix for each taxon
z_matrix <- lapply (seq (1,max_spp), function (k)


              do.call(rbind,
        
                    lapply (record_taxa_z, function (i)
          
          
                      # get the matrix for 1 taxon
          
                      i[k,]
  
)))


apply (z_matrix[[1]],2,sum,na.rm=T)/nsites # true occupancy


plot(apply (z_matrix[[1]],2,sum,na.rm=T)/nsites,
     type="l",
     col=rgb(0.1,0.1,0.5,alpha=0.5)) # true occupancy)

# plot 
lapply (z_matrix, function (i){

  points (apply (i,2,sum,na.rm=T)/nsites,col=rgb(0.1,0.1,0.5,alpha=0.5),pch=19)
  lines (apply (i,2,sum,na.rm=T)/nsites,col=rgb(0.1,0.1,0.5,alpha=0.5))
}
)


# Convert the list to an array
z_array <- array(unlist(z_matrix), 
                 dim = c(nsites, nintervals, max_spp))

# simulate the observation process
y <- array (NA,dim=c(nsites,nsurveys, nintervals,max_spp)) # observed data

# simulate sampling

for (t in 1:nintervals){
  
  for (j in 1:nsurveys) {
    
    for (s in 1:max_spp) {
    
    y[,j,t,s] <- rbinom (n=nsites, size=1,prob=z_array[,t,s]*p)
    
    }
  }
}
str(y)  


# long format
long_y <- melt(y)
colnames(long_y) <- c("site","obs", "int","spp","det")
long_y<- long_y[is.na(long_y$det) !=T,]
apply (long_y,2,unique)


# design the model
sink("dyn_model_vectorized_no_covariates.txt")
cat("
   
   model {
    
     #############################################################
    #                                                           #
    #                  Biological process                       #
    #                                                           #
    #############################################################
    
    
    # Site Occupancy Dynamics (Priors)
    
   gamma~dunif(0,1)
   phi~dunif(0,1)
    
   ## set initial conditions
   psi1~dunif(0,1)
       
   ############      Model       #############
   
   
   for (i in 1:nsites) {

    for (s in 1:nspp) {
    
          z[i,1,s]~dbern(psi1) # occupancy status initialization
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[i,t,s] <- z[i,t-1,s] *  phi + ### if occupied, p of not getting extinct/persist in the next time
                          (1-z[i,t-1,s]) *  gamma ###  if not occupied, p of originate in the next time
                
        # realized occurrence
		    z[i,t,s] ~ dbern(muZ[i,t,s])
    
        }#t
      } #spp
     }#i
    
    #############################################################
    #                                                           #
    #         Observation process across formations             #
    #                                                           #
    #############################################################
    
    
    
        # Priors for detection probability
        p~dunif(0,1)
        

     ############      Model       #############
     
     
     # observation submodel
     
     for (k in 1:nobs) { ## loop over observations
                            
               y [k] ~ dbern(muY[site[k],survey[k],int[k],spp[k]])
               # Specify the observation model conditional on occupancy state
               muY [site[k],survey[k],int[k],spp[k]] <- z[site[k],int[k],spp[k]] * p
                
                }
      
        
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()

## bundle data
str(jags.data <- list(y = long_y$det,
                      int=long_y$int,
                      site = long_y$site,
                      survey = long_y$obs,
                      spp = long_y$spp,
                      nsites = nsites,
                      nint= nintervals, 
                      nspp=max_spp,
                      nobs=nrow(long_y)))
#inits
inits <- function(){list(z = (z_array))}

## Parameters to monitor
## long form
params <- c(
  "z",
  "gamma", 
  "phi",
  "p"
  
)

## MCMC settings
######################
## short form
#na <- 6000; nb <- 7000; ni <- 15000; nc <- 3; nt <- 8
na <- 250; nb <- 500; ni <- 1000; nc <- 3; nt <- 1

# MCMC runs
# models
samples_sims <- jags (data = jags.data, 
                      parameters.to.save = params, 
                      model.file = "dyn_model_vectorized_no_covariates.txt", 
                      inits = inits, 
                      n.chains = nc, 
                      n.thin = nt, 
                      n.iter = ni, 
                      n.burnin = nb, 
                      DIC = T,  
                      #n.cores=nc,
                      parallel=F
)


round(mean((apply (samples_sims$mean$z,c(2,3),mean))[1,])/bins$duration_myr[1],2)






# summarize results
# data frame with observed data
obs_data <- data.frame (phi=phi,
                        gamma=gamma,
                        p=p,
                        dat="obs",
                        lic=NA,
                        uic=NA)
obs_data <- melt (obs_data,id.vars=c("dat","dat","uic", "lic"))


# simulated data (with detection)
sim_data_phi <- data.frame (phi=samples_sims$mean$phi,
            uic=quantile (samples_sims$sims.list$phi,0.975),
            lic=quantile (samples_sims$sims.list$phi,0.025),
            dat="sim")
sim_data_gamma <- data.frame (gamma=samples_sims$mean$gamma,
                            uic=quantile (samples_sims$sims.list$gamma,0.975),
                            lic=quantile (samples_sims$sims.list$gamma,0.025),
                            dat="sim")
sim_data_p <- data.frame (p=samples_sims$mean$p,
                            uic=quantile (samples_sims$sims.list$p,0.975),
                            lic=quantile (samples_sims$sims.list$p,0.025),
                            dat="sim")
# bind output sim data
sim_data <- rbind (
  
  melt (sim_data_phi, id.vars=c("dat","uic", "lic")),
  melt (sim_data_gamma, id.vars=c("dat","uic", "lic")),
  melt (sim_data_p, id.vars=c("dat","uic", "lic"))
  
)


# simulated data (no detection)
sim_data_phi_ND <- data.frame (phi=samples_sims_no_det$mean$phi,
                            uic=quantile (samples_sims_no_det$sims.list$phi,0.975),
                            lic=quantile (samples_sims_no_det$sims.list$phi,0.025),
                            dat="simND")
sim_data_gamma_ND <- data.frame (gamma=samples_sims_no_det$mean$gamma,
                              uic=quantile (samples_sims_no_det$sims.list$gamma,0.975),
                              lic=quantile (samples_sims_no_det$sims.list$gamma,0.025),
                              dat="simND")
# bind output sim data
sim_data_ND <- rbind (
  
  melt (sim_data_phi_ND, id.vars=c("dat","uic", "lic")),
  melt (sim_data_gamma_ND, id.vars=c("dat","uic", "lic"))
  
)

# bind simulated without and with detection
sim_data <- rbind (sim_data,sim_data_ND)


require(ggplot2)

ggplot (sim_data, aes(x=value,
                      y=variable,
                      fill=dat,
                      col=dat))+
  geom_point(size=2,alpha=0.5)+
  geom_errorbar(aes(xmin = lic, xmax = uic),width=.1)+
  geom_point(data=obs_data, aes (x=value,y=variable),
             col="red",
             alpha=0.5,
             size=4)+
  theme_bw()



# design the multispp model
sink("dyn_model_vectorized_no_detection_multispp.txt")
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
      for (s in 1:nspp) {
    
        p[i,1,s] <- psinit    # Initial occupancy probability
     
      }
    }
    

   ############      Model       #############
   for (i in 1:nsites) {
      for (s in 1:nspp) {
      
          y[i,1,s] ~ dbern(psi1) # Initial occupancy probability
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time observed occurrence y
                
                muZ[i,t,s] <- y[i,t-1,s] *  phi[i] + ### if occupied, p of not getting extinct/persist in the next time
                          (1-y[i,t-1,s]) *  gamma[i] ###  if not occupied, p of originate in the next time
                
        # realized occurrence
		    p[i,t,s] ~ dbern(muZ[i,t,s])      # Occupancy model
		    y[i,t,s] ~ dbern(p[i,t,s])        # Observation model
    
        }#t
      }# spp
     }#i
    
    
    }## end of the model
    
    
    
    ",fill = TRUE)
sink()


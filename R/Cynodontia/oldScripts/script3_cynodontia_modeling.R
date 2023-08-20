# ------------------------------------------------------------

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
    
   for (i in 1:nsites) {
    
    for (t in 1:nint) {
    
      gamma [i,t]~dunif(0,1)
      phi [i,t]~dunif(0,1)
    
    
    }
   
   }
   
   ## set initial conditions
   for (g in 1:ngen) {
   
      psi1[g] ~ dunif(0,1)
       
   }
    ############      Model       #############
    
    
    for (i in 1:nsites) {
    
      for (g in 1:ngen) {
    
          z[i,1,g]~dbern(psi1[g]) # occupancy status initialization
    
              for (t in 2:nint){
            
                # model likelihood
                ### modeling dynamics conditional on previous time realized occurrence z
                muZ[i,t,g] <- z[i,t-1,g] *  phi[i,t] + ### if occupied, p of not getting extinct/persist in the next time
                          (1-z[i,t-1,g]) *  gamma[i,t] ###  if not occupied, p of originate in the next time
                
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
        for (g in 1:ngen) {
   
            p[g] ~ dunif(0,1)
        
        }

     ############      Model       #############
     
     
     # observation submodel
     
     for (k in 1:nobs) { ## loop over observations
                            
               y [k] ~ dbern(muY[site[k],form[k],int[k],gen[k]])
               # Specify the observation model conditional on occupancy state
               muY [site[k],form[k],int[k],gen[k]] <- z[site[k],int[k],gen[k]] * p[gen[k]]
                
                
             }
      
        
    }## end of the model
    
    
    

    
    ",fill = TRUE)
sink()


# --------------------------------------


# load packages
source("R/packages.R")


# --------------------------------------


# synapsida level
# load data
#load(here ("processed_data","table_data_array_cynodontia.RData"))
load(here ("processed_data","table_data_array_cynodontia_smaller_bins.RData"))
#load(here ("processed_data","table_space.RData"))
load(here ("processed_data","table_space_smaller_bins.RData"))
#load(here ("processed_data", "table_naive_cynodontia.RData"))
load(here ("processed_data", "table_naive_cynodontia_small_bins.RData"))
#load (here ("processed_data", "site_covs.RData"))

## bundle data
str(jags.data <- list(y = table_data_long_df$det,
                      gen = table_data_long_df$taxon_code,
                      int= table_data_long_df$int,
                      form = table_data_long_df$form,
                      site = table_data_long_df$site,
                      nsites = length(unique(table_data_long_df$site)),
                      nint= length(unique(table_data_long_df$int)), 
                      ngen = max(unique(table_data_long_df$taxon_code)), 
                      nobs=nrow(table_data_long_df)
                      
                      )
)

# naive table
naive<-cast (formula=site ~ int,
             data = table_data_long_df,
             value = "det",
             fill=0,
             fun.aggregate = max)[,-1]

#  number of sites per interval
plot(rev(colSums(naive>0)/40),type="l",ylim=c(0,1))

# list
naive_list <- replicate (n=max(unique(table_data_long_df$taxon_code)),
                         naive,
                         simplify = F)
zst <- array(unlist(naive), 
             dim = c(nrow(naive),
                    ncol (naive), 
                    length(naive_list)))
lines(rev(colSums(zst[,,1]>0)/40),type="l",col="red",lty=2)

# inits for psi1
#psi1 <- (apply (zst, c(1,3),mean,na.rm=T))

# function inits (to be placed in each chain)
inits <- function(){ list(z = zst)}

## Parameters to monitor
## long form
params <- c(
  
  "gamma", "phi","p",
   "psi1",	
   "z"
)

## MCMC settings
######################
## short form
#na <- 6000; nb <- 7000; ni <- 15000; nc <- 3; nt <- 20
na <- 250; nb <- 500; ni <- 1000; nc <- 3; nt <- 1

# MCMC runs
# models
samples_paleo_cynodontia <- jags (data = jags.data, 
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



# save
save (samples_paleo_cynodontia,file = here ("output","samples_paleo_cynodontia_smaller_bins.RData"))

# observed data
obs_data <- (cast (data = table_data_long_df,
                   formula = int ~taxon_code,
                   value = "det",
                   fun.aggregate = max))



lapply (seq (1,dim(samples_paleo_cynodontia$mean$z)[3]), function (i)

  lines((colSums(samples_paleo_cynodontia$mean$z[,,i])/40),type="l",col="green",lty=2)

)


# time trend in SR
data.frame (ntaxa=(rowSums(obs_data)),
            est=(rowSums(apply (samples_paleo_cynodontia$mean$z,c(2,3), max))),
            int =names((rowSums(obs_data)))) %>%
  
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



#  number of sites per interval
plot(NA,ylim=c(0,1),xlim = c(1,66))

lapply (seq(1,nrow(samples_paleo_cynodontia$mean$gamma)), function (i)
  lines (samples_paleo_cynodontia$mean$gamma[i,])
)


#  number of sites per interval
plot(NA,ylim=c(0,1),xlim = c(1,66))

lapply (seq(1,nrow(samples_paleo_cynodontia$mean$phi)), function (i)
  lines (samples_paleo_cynodontia$mean$phi[i,])
)


# time trend in gamma and phi
data.frame (cbind (apply (samples_paleo_cynodontia$mean$gamma,2, mean),
            int =  seq(1,66),
            var = "gamma")) %>%
  
  melt (id.vars=c("int","var")) %>%
  
  rbind(data.frame (cbind (apply (samples_paleo_cynodontia$mean$phi,2, mean),
                                int =  seq(1,66),
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





# https://www.canva.com/design/DAFdRXwxzQU/DHsyzZsWq7F2U8gA6g8PFA/edit?utm_source=shareButton&utm_medium=email&utm_campaign=designshare
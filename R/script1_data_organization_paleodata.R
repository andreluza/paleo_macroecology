# ======================================================

# organizing paleodata

# ======================================================

# =====================================================
# load packages
require(here);require(rgeos);require(rgdal);require(sp);require(raster); require(openxlsx)
# data processing & plot
library(ggplot2); library(rasterVis); library(rgdal);library(viridis); require(dplyr)

# =====================================================
# load data
# collections
pbdb_data_collections <- read.csv (here ("data", "PaleoData", "pbdb_data_collections.csv"),
                       skip=16)

# strata
pbdb_data_strata <- read.csv (here ("data", "PaleoData", "pbdb_data_strata.csv"),
                              skip=15)

#taxa
pbdb_data_taxa <- read.csv (here ("data", "PaleoData", "pbdb_data_taxa.csv"),
                              skip=19)

# occ
pbdb_data_occ <- read.csv (here ("data", "PaleoData", "pbdb_data_occ.csv"),
                           skip=16)

# matching occurrence and taxa
colnames(pbdb_data_taxa) [which(colnames(pbdb_data_taxa) %in% colnames(pbdb_data_occ))]
# matching strata and collections
colnames(pbdb_data_collections) [which(colnames(pbdb_data_collections) %in% colnames(pbdb_data_strata))]
# matching collections and occ 
colnames(pbdb_data_collections) [which(colnames(pbdb_data_collections) %in% colnames(pbdb_data_occ))]


# the unifier ID is collection_no
table(unique(pbdb_data_collections$collection_no) %in% unique(pbdb_data_occ$collection_no))

# merging datasets
coll_occ <- merge (pbdb_data_occ,pbdb_data_collections, by="collection_no", all=T)
# aggregate taxa
coll_occ_taxa <- merge (coll_occ,pbdb_data_taxa, by="accepted_name", all=T)

# list of intervals of interest
triassic <- c("Olenekian", "Anisian", "Ladinian", "Carnian", "Norian", "Rhaetian") # Triassic
# triassic %in% coll_occ_taxa$early_interval.x
jurassic <- c("Hettangian", "Sinemurian", "Pliensbachian", "Toarcian", "Aalenian","Bajocian", "Bathonian", "Callovian", "Oxfordian", "Kimmeridgian", "Tithonian")
jurassic %in% coll_occ_taxa$early_interval.x

# subset of intervals within Jurassic and Triassic
coll_occ_taxa_jur_tri <- coll_occ_taxa [which(coll_occ_taxa$early_interval.x %in% 
                                                c(triassic, jurassic)),]

# level of genus to lower levels
coll_occ_taxa_jur_tri <- coll_occ_taxa_jur_tri[which(coll_occ_taxa_jur_tri$taxon_rank %in% 
                              c("genus", "species")),]
# collapsing species
coll_occ_taxa_jur_tri$genus <- sapply (strsplit(coll_occ_taxa_jur_tri$accepted_name, " "), "[",1)
  
# formations (replicates)
all_formations <- unique(coll_occ_taxa_jur_tri$formation)[order(unique(coll_occ_taxa_jur_tri$formation))]
coll_occ_taxa_jur_tri$formation <- as.factor (coll_occ_taxa_jur_tri$formation)

# recode
coll_occ_taxa_jur_tri <- coll_occ_taxa_jur_tri %>% 
  mutate (formation = recode (formation,
                              "Baños del Flaco" = "Banos del Flaco",
                              "Asiento" = "Asientos",
                              "Cantos del Agua" = "Canto del Agua",
                              "Ischigalasto"  = "Ischigualasto", 
                              "Quehita"="Quehuita"
                              ))

# --------------------------------------

# organize a table for each genus in the data
require(reshape)
genus_data <- unique(coll_occ_taxa_jur_tri$genus)

# table 
table_data <- lapply (genus_data, function (i) {

  tryCatch({
            
        # core table with existing data
            table_data <- (cast (formula = formation ~ early_interval.y,
                  data = coll_occ_taxa_jur_tri[which(coll_occ_taxa_jur_tri$genus %in% i),],
                  value = "min_ma.y", 
                  fun.aggregate = mean,
                  drop=F,
                  fill=0
                  ))
            
            
            if (nrow (table_data) ==1) {
              name <- colnames (table_data)[2]
              table_data<-data.frame (name = table_data[,2],
                                      row.names = table_data$formation)
              colnames (table_data) <- name
              
            } else {
            
              rownames (table_data)<- table_data$formation; table_data<-table_data[,-1]
            
            }
            
            # bind missing formations (zeros; true absences)
            input_formations <- matrix (0, 
                                        
                                      nrow = sum (all_formations %in% rownames(table_data) != T),
                                      ncol = ncol (table_data),
                                      byrow=T,
                                      
                                      dimnames = list (
                                          
                                          all_formations[which(all_formations %in% rownames (table_data) == F)],
                                          colnames(table_data)
                                           
                                          )
                              )
            # bind 
            table_data_formations <- rbind (table_data,
                                            input_formations)
            
            # input intervals (NAs, not assessed)
            input_intervals <- matrix (NA, 
                                        
                                        nrow = nrow (table_data_formations),
                                        ncol = sum(c(triassic,jurassic) %in% colnames (table_data_formations) != T),
                                        byrow=T,
                                        
                                        dimnames = list (
                                          
                                          rownames (table_data_formations),
                                          c(triassic,jurassic) [which(c(triassic,jurassic) %in% colnames (table_data_formations) != T)]
                                          
                                        )
            )
            
            # bind 
            table_data_formations_intervals <- cbind (table_data_formations,
                                                      input_intervals)
            
            # order cols
            table_data_formations_intervals <- table_data_formations_intervals[,match(c(triassic,jurassic),
                                                   colnames (table_data_formations_intervals))]
            # remove empty rownames
            table_data_formations_intervals <- table_data_formations_intervals[which(nchar (rownames(table_data_formations_intervals))>0),]
            # order of formations
            table_data_formations_intervals <- table_data_formations_intervals[order (rownames(table_data_formations_intervals)),]
            
            ; # return
            table_data_formations_intervals
        
        
        }, # close tryCatch
    
    
        error = function(e) return ("NULL"))

})

# list into array (dim1= formations, dim2= intervals, dim3= species)
table_data_array <- array(as.numeric(unlist(table_data)), dim = c(nrow(table_data[[1]]),ncol(table_data[[1]]),length(table_data)))
dimnames (table_data_array) <- list (rownames (table_data[[1]]),
                                 colnames (table_data[[1]]),
                                 genus_data)
# changing time (workaround) with 1s
table_data_array <- ifelse (table_data_array >1,1,0)

# remova NA
table_data_array<-table_data_array[,,-1]

# find the number of detections per genus

table (unlist(
  
  lapply (seq (1,dim(table_data_array)[3]), function (i)

  sum(table_data_array[,,i]>0,na.rm=T)
 )
)
)

# select detection in at least one period and formation
sel_genus <- unlist(
  
  lapply (seq (1,dim(table_data_array)[3]), function (i)
    
    sum(table_data_array[,,i]>0,na.rm=T)
  )
)>0

table_data_array <- table_data_array [,,sel_genus]


# ------------------------------------------------------------

# design the model

sink("dyn_model.txt")
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
          epslon [t,g]~ dunif (0,1)
        
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
              muZ[t,g] <- z[t-1,g] * (1-epslon[t,g]) + ### if occupied, p of not getting extinct
                          (1-z[t-1,g]) * gamma[t,g] ###  if not occupied, p of getting colonized
                        
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
    for (j in 1:nform){
      for (t in 1:nint) {
         for (g in 1:ngenus) {
            y [j,t,g] ~ dbern (muY[j,t,g])
            muY [j,t,g] <- z [t,g] * p [g]
          }
       }
    }    

    ## derived parameters
    # number of genus per interval
    for (t in 1:nint) {
        Ngen[t]<-sum(z[t,])
    }
    
    # average extinction and origination
    for (g in 1:ngenus) {
      avepslon[g] <- mean(epslon[2:nint,g])
      avgamma[g]<- mean(gamma[2:nint,g])
    }
      
    # turnover (proportional gain or loss)
    for (t in 2:nint) {  
      
        propcH [t] <-(sum (z[t-1,]) - sum(z[t,]))/sum(z[t-1,]) 
      
    }
    
    # equilibrium occupancy (which genus decline or increase over time)
    for (g in 1:ngenus) {
    
        psi.eq[g] <- mean(gamma[2:nint,g])/(mean(gamma[2:nint,g])+mean(epslon[2:nint,g])) # Equilibrium occupancy
    
    }
    
    
    }## end of the model
    
    ",fill = TRUE)
sink()


## bundle data

str(jags.data <- list(y = unname (table_data_array), 
                      nform = dim(table_data_array)[1], 
                      nint= dim(table_data_array)[2], 
                      ngenus = dim(table_data_array)[3]
                      ))

# Set initial values
zst <- apply(unname (table_data_array), c(2,3), max, na.rm = TRUE)	# Observed occurrence as inits for z
zst[zst == '-Inf'] <- 1 # max of c(NA,NA,NA) with na.rm = TRUE returns -Inf, change to 1
inits <- function(){ list(z = zst)}

## Parameters to monitor
## long form
params <- c(
  
  "gamma", "epslon","p","muZ",
  "propcH", "avgamma", "avepslon",
  "Ngen",'psi.eq'
  
)

## MCMC settings
######################
## short form
na <- 5000; nb <- 6000; ni <- 10000; nc <- 3; nt <- 2
#na <- 50; nb <- 60; ni <- 100; nc <- 3; nt <- 1

require(jagsUI)
samples_paleo <- jags (data = jags.data, 
                      parameters.to.save = params, 
                      model.file = "dyn_model.txt", 
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
dir.create ("output")
save (samples_paleo,file = here ("output","samples_paleo.RData"))
save (table_data_array,file = here ("output","table_data_array.RData"))

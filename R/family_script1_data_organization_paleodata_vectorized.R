# ======================================================

# organizing paleodata

# ======================================================

# =====================================================
# load packages
require(here);require(rgeos);require(rgdal);require(sp);require(raster); require(openxlsx)
# data processing & plot
library(ggplot2); library(rasterVis); library(viridis); require(dplyr); require(magick); require(reshape)

# paleomaps
library(mapast)

# =====================================================
# load data
# collections
pbdb_data_collections <- read.csv (here ("data", "PaleoData", "pbdb_data_collections.csv"),
                       skip=17)

# strata
pbdb_data_strata <- read.csv (here ("data", "PaleoData", "pbdb_data_strata.csv"),
                              skip=16)

#taxa
pbdb_data_taxa <- read.csv (here ("data", "PaleoData", "pbdb_data_taxa.csv"),
                              skip=20)

# occ
pbdb_data_occ <- read.csv (here ("data", "PaleoData", "pbdb_data_occ.csv"),
                           skip=18)


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

# get the taxonomic classification of taxa
# https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c
require(rgbif)
unique_accepted_names <- unique (coll_occ_taxa_jur_tri$accepted_name)
test_lookup<-lapply (unique_accepted_names, name_lookup)


# formations (replicates)
all_formations <- unique(coll_occ_taxa_jur_tri$formation.x)[order(unique(coll_occ_taxa_jur_tri$formation.x))]
coll_occ_taxa_jur_tri$formation.x <- as.factor (coll_occ_taxa_jur_tri$formation.x)

# recode
coll_occ_taxa_jur_tri <- coll_occ_taxa_jur_tri %>% 
  mutate (formation.x = recode (formation.x,
                              "Baños del Flaco" = "Banos del Flaco",
                              "Asiento" = "Asientos",
                              "Cantos del Agua" = "Canto del Agua",
                              "Ischigalasto"  = "Ischigualasto", 
                              "Quehita"="Quehuita",
                              "Agardfjellet" = "Agardhfjellet",
                              "Alcobaca" = "Alcobaça",
                              "Bakkar" = "Bakhar",
                              "Bazanovo and Bokhto Suites" = "Bazanov and Bokhtin Suite",
                              "Bihin Limestone" = "Bihen Limestone",
                              "Birdlip Limestone " = "Birdlip Limestone",
                              "Buchenstein " = "Buchenstein",
                              "Cabacos" = "Cabaços",
                              "Fullers Earth" = "Fuller's Earth",
                              "Halstatt" = "Hallstatt",
                              "Hallstätt" = "Hallstatt",
                              "Inklin " = "Inklin",
                              "Kimmeridge-clay" = "Kimmeridge Clay",
                              "Kossen"="Kössen",
                              "Limanangcong" = "Liminangcong",
                              "Lower Calcaerous Grit"="Lower Calcareous Grit",
                              "luning" = "Luning",
                              "Manual Creek"="Manuel Creek",
                              "Nordenskjold" = "Nordenskjöld",
                              "Parson Bay" = "Parsons Bay",
                              "Ronne" = "Rønne",
                              "Schöll"= "Schnöll",
                              "Sewa" = "Sêwa",
                              "steinalm Limestones" = "Steinalm Limestones",
                              "Timezgadiwine" = "Timezgadiouine",
                              "Tschermakjfellet" = "Tschermakfjellet",
                              "Ulughey" = "Ulugei",
                              "Vikinghogda"="Vikinghøgda"
                              
  ))
# corrected names
all_formations <- unique(coll_occ_taxa_jur_tri$formation.x)

# bind period in data
coll_occ_taxa_jur_tri$Period <- ifelse (coll_occ_taxa_jur_tri$early_interval.x %in%
                                          jurassic, "Jurassic","Triassic")
# create index as value in the table
coll_occ_taxa_jur_tri$detection <- 1

# get the age of each interval do have a plot.

ages <- tapply (coll_occ_taxa_jur_tri$min_ma.y, 
        list (coll_occ_taxa_jur_tri$early_interval.x),
        mean,
        na.rm=T)
ages<- ages[order(ages)]
# doubt 
# agoudim e agoudim marls
# Augusta e Augusta Mountain
# badamu e badamu limestone
# Dolomite
# Dolomitic
# Dolomitique

# paleomaps of formations - take a time ...
#maps_models <- getmap(ma = unname(ages), 
#                         model = "GOLONKA")
#maps_models <- maps_models [order(ages)] # ordering
#save (maps_models, file = here ("output", "paleomaps.RData")) 
load(here ("output", "paleomaps.RData"))
names(maps_models) <- names(ages)

# plot points
lapply (names (ages), function (i) {

    png (here ("output", paste ("maps_", (i), ".png",sep ="")),
               width = 30,height = 15, res=300, 
           units = "cm")
        par (mar=rep (1,4))
        # triassic
        plot (maps_models[[which(names(maps_models) == (i))]],
              main = paste (i, " (",
                            round(mean(coll_occ_taxa_jur_tri$min_ma.y[which(coll_occ_taxa_jur_tri$early_interval.x == i)]),2),
                            " Ma)",sep=""),
                            
              border = "gray",
              col = rgb (0.1,0.1,0.1,alpha=0.1))
        
        with (coll_occ_taxa_jur_tri[which(coll_occ_taxa_jur_tri$early_interval.x == i),], 
              
              points (paleolng,
                      paleolat,
                      pch=19,
                      col = rgb(0.1,0.5,0,alpha=0.05))
        )
        

        dev.off()
})

#anime
list_img <- list.files(path = here ("output"), full.names = T, pattern = "maps_")
# ordering
list_img <- list_img[match (names(ages),
                    gsub (".png", "",gsub ("D:/Pos_Doc_Paleonto_Macroecology/modeling/paleo_macroecology/output/maps_","",list_img))
)]
##https://cran.r-project.org/web/packages/magick/vignettes/intro.html
a_image<-image_read(list_img)
animation <-  image_animate(a_image, fps = 1)
image_write(animation, here ("output","animation_data.gif"))

# --------------------------------------

# organize a table for each genus in the data
genus_data <- unique(coll_occ_taxa_jur_tri$genus)


# the basic table summarizing all detections across intervals and formations
table_data_basis <- (cast (formula = formation.x ~ early_interval.y,
                     data = coll_occ_taxa_jur_tri,
                     value = "detection", 
                     fun.aggregate = sum,
                     drop=F,
                     fill=NA
))
# order names
table_data_basis <- table_data_basis[order(table_data_basis$formation.x),]
rownames(table_data_basis) <- table_data_basis$formation.x; table_data_basis<- table_data_basis[-1,-1]
# order cols
table_data_basis <- table_data_basis[,match(c(triassic,jurassic),
                                      colnames (table_data_basis))]

# there is no observation in 88% of the combinations between intervals and formations
table(is.na(table_data_basis))/sum(table(is.na(table_data_basis)))

# summation of detections (each entry in the database is a detection)
range(table_data_basis,na.rm=T)

# detections
hist(as.numeric(data.matrix(table_data_basis[-1,-1]))[is.na(as.numeric(data.matrix(table_data_basis[-1,-1])))!=T])

# since 1 to 1775 detections
range(as.numeric(data.matrix(table_data_basis[-1,-1]))[is.na(as.numeric(data.matrix(table_data_basis[-1,-1])))!=T])



# table 
table_data <- lapply (genus_data, function (i) {

  tryCatch({
            
        # core table with existing data
            table_data <- (cast (formula = formation.x ~ early_interval.y,
                  data = coll_occ_taxa_jur_tri[which(coll_occ_taxa_jur_tri$genus %in% i),],
                  value = "detection", 
                  fun.aggregate = max,
                  drop=F,
                  fill=NA
                  ))
            
            
            if (nrow (table_data) ==1) {
              name <- colnames (table_data)[2]
              table_data<-data.frame (name = table_data[,2],
                                      row.names = table_data$formation.x)
              colnames (table_data) <- name
              
            } else {
            
              rownames (table_data)<- table_data$formation.x; table_data<-table_data[,-1]
            
            }
            
            # bind missing formations (zeros; true absences)
            input_formations <- matrix (NA, 
                                        
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



# matching with basis table
# if no taxa was detected in any time, set NA
# if at least one taxa was detected in a time, set zero

table_data_long <- lapply (seq(1,length(table_data)), function (i) {
  
  # df frame with detections of each taxon and detections in the complete basis dataset 
  test <- data.frame (det_taxon=as.numeric(unlist(table_data[[i]])),
            det_all_taxon=(as.numeric(unlist(table_data_basis))))
  # matching
  test$test <- ifelse (is.na (test$det_all_taxon) &  is.na(test$det_taxon), # if both are NA
          NA, # NA is the output
          ifelse ((is.na (test$det_all_taxon) != T) &  is.na(test$det_taxon), # if at least one taxon was detected in any time and formation
                  0, # zero is output (non-detection)
                  test$det_taxon # otherwise the data will be filled with ones
  ))

  # matrix
  matrix_data <- matrix(as.numeric(test$test),
                         nrow=nrow(table_data_basis),
                         ncol=ncol(table_data_basis),
                         byrow = F,
                         dimnames = list(seq(1,length(rownames(table_data_basis))),
                                         seq(1,length(colnames(table_data_basis)))
                                         )
                        )

  # check if tables match
  # table (is.na(matrix_data) == is.na(table_data_basis))
  
  # melt to long format
  table_data_long <- melt(matrix_data,as.is=T)
  # range(table((table_data_long$X1))) # check if it's ok
  colnames(table_data_long) <- c("form", "int", "det")
  # removing NAs
  table_data_long<- table_data_long[is.na(table_data_long$det)!=T,]
  
  ; # return
  table_data_long
  
})

# detection/non-detection data relative to the total N of possible observations
mean((unlist(lapply(table_data_long,nrow))))/(nrow (table_data_basis)*ncol(table_data_basis))
# detection/non-detection data
range(unlist(lapply(table_data_long,nrow)))

# detections of each genus
hist(unlist(
lapply (table_data_long, function (i)
  
  sum(i$det)
)
)
)


# find the number of detections per genus
par(mar = rep(4,4))
# number of genus with 5 or more detections
table(unlist(lapply (
  
  table_data_long, function (i)
    
    sum(i$det)
)
) >=5
)

# select detection in at least one period and formation
sel_genus <- unlist(lapply (
  
  table_data_long, function (i)
    
    sum(i$det)
)
) >=5
barplot(table(sel_genus))



# table  for naive occupancy
table_naive <- table_data [which(sel_genus==T)]
save (table_naive,  file = here ("output", "table_naive.RData"))


# filter
table_data_long <- table_data_long [which(sel_genus==T)]

# bind code for genus
table_data_long <- lapply (seq(1,length(table_data_long)), function (i) {
  
  table_data_long [[i]] <- cbind (table_data_long[[i]],
                                  taxon = i)
  
})
# melt
table_data_long_df <- do.call(rbind,table_data_long) # [1:2]

# save to use after
save (genus_data,
      sel_genus,
      table_data_basis,
      table_data_long,
      table_data_long_df,
      all_formations,
      file = here ("output","table_data_array.RData"))
rm(list=ls())

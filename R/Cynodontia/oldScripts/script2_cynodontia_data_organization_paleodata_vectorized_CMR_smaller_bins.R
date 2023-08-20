

# ======================================================


# CMR design (species in the rows and intervals/time bin in the cols)
# covariates (file "site_covs.RData") will gathered from the code "script2_cynodontia_data_organization_paleodata_vectorized.R"
# with alternative bins 

# organizing paleodata
# the metadata for all this
# https://taphonomy.doit.wisc.edu/data1.2/colls/single_doc.html
# geoplates
# https://www.gplates.org/docs/user-manual/Interacting_Features/#30-assign-plate-ids
# https://www.gplates.org/docs/user-manual/Reconstructions/#plate-ids

# ======================================================

# =====================================================

# load packages
source("R/packages.R")
source("R/functions.R")

# =====================================================
# load data
# collections
pbdb_data_collections <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_collections.csv"),
                                   skip=18)

# strata
pbdb_data_strata <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_strata.csv"),
                              skip=17)

#taxa
pbdb_data_taxa <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_taxa.csv"),
                            skip=15)

# occ
pbdb_data_occ <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_occ.csv"),
                           skip=18)




# --------------------------------------------------------------------
# run a cleaning of coordinates
# help here : https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_PBDB_fossils_with_CoordinateCleaner.html

library(CoordinateCleaner)

pbdb_data_occ <- cc_val(pbdb_data_occ, lat = "lat", lon = "lng") # valid coords
flag_equal <- cc_equ(pbdb_data_occ, value = "flagged",lat = "lat", lon = "lng") # equal coordinates
# extract and check the flagged records
pbdb_data_occ <- pbdb_data_occ[flag_equal==T,] 
flag_pol <- cc_cen(pbdb_data_occ, lat = "lat", lon = "lng", value = "flagged") # political coords
pbdb_data_occ <- pbdb_data_occ[flag_pol==T,]  # filter 
flag_inst <- cc_inst(pbdb_data_occ, lat = "lat", lon = "lng",value = "flagged") # flag institutions
pbdb_data_occ <- pbdb_data_occ[flag_inst==T,]  # filter 
flag_gbif <- cc_gbif(pbdb_data_occ, lat = "lat", lon = "lng",value = "flagged") # flag gbif headquarters
flag_zero <- cc_zero(pbdb_data_occ, lat = "lat", lon = "lng",value = "flagged") # flag zero coordinates
pbdb_data_occ <- pbdb_data_occ[flag_zero==T,]  # filter 


# matching all datasets
colnames(pbdb_data_taxa) [which(colnames(pbdb_data_taxa) %in% colnames(pbdb_data_occ))]
# matching strata and collections
colnames(pbdb_data_collections) [which(colnames(pbdb_data_collections) %in% colnames(pbdb_data_strata))]
# matching collections and occ 
colnames(pbdb_data_collections) [which(colnames(pbdb_data_collections) %in% colnames(pbdb_data_occ))]
# the unifier ID is collection_no
table(unique(pbdb_data_collections$collection_no) %in% unique(pbdb_data_occ$collection_no))

# merging datasets
# previous coll_occ
coll_occ_taxa <- merge (pbdb_data_occ,
                        pbdb_data_collections, 
                        by="collection_no", 
                        all=T)
coll_occ_taxa <- cbind (coll_occ_taxa, 
                        # bind rank
                        taxon_rank=pbdb_data_taxa [match (coll_occ_taxa$accepted_name,
                                                          pbdb_data_taxa$accepted_name),
                                                   "taxon_rank"]
)



# --------------------------------------------------------------------
# obtain ages
# all this done with the help of paleoverse package
# https://cran.r-project.org/web/packages/palaeoverse/vignettes/tetrapod-biodiversity.html
# Rename columns


# binning time intervals
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)


# mid point
coll_occ_taxa$interval_mid_ma <- (coll_occ_taxa$min_ma.x + coll_occ_taxa$max_ma.x)/2

# check intervals
coll_occ_taxa %>%
  mutate(age_range = abs(max_ma.x - min_ma.x)) %>%
  ggplot(aes(x = age_range)) +
  geom_histogram() +
  labs(x = 'Age range (My)', y = 'Count') +
  theme_bw()


# rename
#colnames(coll_occ_taxa)[which(colnames (coll_occ_taxa) %in% c("max_ma", "min_ma"))] <- c("old_max_ma", "old_min_ma")
#colnames(coll_occ_taxa)[which(colnames (coll_occ_taxa) %in% c("interval_max_ma", "interval_min_ma"))] <- c("max_ma", "min_ma")

# we're interested in permian to triassic data
# Remove occurrences that are younger than the time intervals we're interested in
coll_occ_taxa_perm_cret <- subset(coll_occ_taxa, max_ma.x <= max(bins$max_ma) & 
                                    min_ma.x >= min(bins$min_ma))


# rename cols to feed the binning function
# the max_ma.y and min_ma.y will be the original ages
colnames(coll_occ_taxa_perm_cret)[which(colnames(coll_occ_taxa_perm_cret) == "max_ma.x")] <- "max_ma"
colnames(coll_occ_taxa_perm_cret)[which(colnames(coll_occ_taxa_perm_cret) == "min_ma.x")] <- "min_ma"


# Generate time bins
# Make a table of potential number of bins
table(coll_occ_taxa_perm_cret$bin) # raw counts
table(coll_occ_taxa_perm_cret$bin) / nrow(coll_occ_taxa_perm_cret) # proportions


## Testing age validity
age_validity <- cf_equal(coll_occ_taxa_perm_cret, min_age = "min_ma", max_age = "max_na")
rang <- coll_occ_taxa_perm_cret$max_ma - coll_occ_taxa_perm_cret$min_ma # age range (= max age - min age) of each record
# the precision is ok
hist(rang, breaks = 40, xlab = "Date range [max age - min age]", main = "Age Range (my)")



## Testing spatio-temporal outliers on dataset level
# Outlier dataset
cl <- cf_outl(coll_occ_taxa_perm_cret, taxon = "", lat = "lat.x", lon = "lng.x",
              min_age = "min_ma", max_age = "max_ma")
# Outlier taxon
cf_out <- cf_outl(coll_occ_taxa_perm_cret, taxon = "accepted_name", lat = "lat.x", lon = "lng.x",
                  min_age = "min_ma", max_age = "max_ma", value = "flagged")
coll_occ_taxa_perm_cret<-coll_occ_taxa_perm_cret[cf_out==T,] # filter


# --------------------------------------------------------------------




# adjust formation and lithology names

# remove quotes and other characters

coll_occ_taxa_perm_cret$formation.x <- gsub ("\\? ", "", coll_occ_taxa_perm_cret$formation.x)
coll_occ_taxa_perm_cret$formation.x <- gsub ("\"", "", coll_occ_taxa_perm_cret$formation.x)
coll_occ_taxa_perm_cret$formation.x <- gsub ("^[:alnum:]", "", coll_occ_taxa_perm_cret$formation.x)
coll_occ_taxa_perm_cret$formation.x <- gsub ("=", "", coll_occ_taxa_perm_cret$formation.x)
coll_occ_taxa_perm_cret$formation.x <-  (noquote(coll_occ_taxa_perm_cret$formation.x))
coll_occ_taxa_perm_cret$formation.x <- gsub("\"", '', coll_occ_taxa_perm_cret$formation.x)
coll_occ_taxa_perm_cret$formation.x <- gsub("'", '', coll_occ_taxa_perm_cret$formation.x)

# remove weird names
coll_occ_taxa_perm_cret<- coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$formation.x %in% c("", "=", "5","6") == F),]
unique(coll_occ_taxa_perm_cret$formation.x ) [order(unique(coll_occ_taxa_perm_cret$formation.x ))]

#as factor
coll_occ_taxa_perm_cret$formation.x <- as.factor (coll_occ_taxa_perm_cret$formation.x)

# recode
coll_occ_taxa_perm_cret <- data.frame (coll_occ_taxa_perm_cret) %>% 
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
                                "Vikinghogda"="Vikinghøgda",
                                "AG 4" = "AG4"
                                
  ))


# corrected names
# all_formations <- unique(coll_occ_taxa_perm_cret$formation.x)

# correct names of lithology
coll_occ_taxa_perm_cret$lithology1.x <- gsub ("\\? ", "", coll_occ_taxa_perm_cret$lithology1.x)
coll_occ_taxa_perm_cret$lithology1.x <- gsub ("\"", "", coll_occ_taxa_perm_cret$lithology1.x)
coll_occ_taxa_perm_cret$lithology1.x[which(coll_occ_taxa_perm_cret$lithology1.x == "not reported")] <- NA
coll_occ_taxa_perm_cret$lithology1.x[which(coll_occ_taxa_perm_cret$lithology1.x == "")] <- NA


# --------------------------------------------------------------------

# bin occurrences spatially into spatial bins


coll_occ_taxa_perm_cret <- bin_space(occdf = coll_occ_taxa_perm_cret %>%
                                       filter (is.na(paleolat2.x) != T),
                                     lng = 'paleolng2.x',
                                     lat = 'paleolat2.x',
                                     spacing = 2000,
                                     #sub_grid = 500,
                                     plot=T,return=T)


# object with the grid
grid_info <- coll_occ_taxa_perm_cret [2:4]


# spatial data
plot(grid_info$grid_base)
plot(grid_info$grid,add=T, col = rgb(0.1,0.1,0.5,alpha=0.2))
#plot(grid_info$sub_grid,add=T, col = rgb(0.5,0.1,0.1,alpha=0.2))


# occurrence data
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret$occdf

# map
points(
  coll_occ_taxa_perm_cret$paleolng2.x,
  coll_occ_taxa_perm_cret$paleolat2.x,
  col = rgb (0.8,0.1,0.1,alpha=0.5),pch=19)

# cell coordinates 
cell_coordinates <- coll_occ_taxa_perm_cret %>% 
  group_by (cell_ID) %>%
  summarise (centroid_lat = mean (cell_centroid_lat),
             centroid_lng = mean(cell_centroid_lng),
             sd_lat = sd (cell_centroid_lat),
             sd_lng = sd (cell_centroid_lng))



#--------------------------------------------------------------

# taxonomic resolution 
# species level : Occurrences identified to a coarser taxonomic resolution than the desired level
# are retained if they belong to a clade which is not otherwise represented in the dataset 
coll_occ_taxa_perm_cret <- tax_unique(occdf = coll_occ_taxa_perm_cret, 
                                      genus = "genus", 
                                      family = "family",
                                      order = "order", 
                                      class = "class", 
                                      name = "accepted_name",
                                      resolution="genus",
                                      append=T)


coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$unique_name != " sp."),]

# -----------------------------------------------

# determine detection and intervals

# create index as value in the table
coll_occ_taxa_perm_cret$detection <- 1

# check preservation quality
table(coll_occ_taxa_perm_cret$preservation_quality.x)
# filter(coll_occ_taxa_perm_cret, preservation_quality.x != "very poor")


# --------------------------------------------------------------------

# select interesting cols (to downsize the dataset)

coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[,  c("collection_no",
                                                        "unique_name", 
                                                        "genus",
                                                        "formation.x",
                                                        "lithology1.x",
                                                        #"bin_assignment",
                                                        #"Interval",
                                                        #"Period",
                                                        "min_ma",
                                                        "max_ma",
                                                        "cell_ID",
                                                        "cell_centroid_lat",
                                                        "paleolng2.x",
                                                        "paleolat2.x",
                                                        "detection"
)]




# alternative bin
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret %>%
  filter(abs(max_ma - min_ma) < 25) %>% # exclude uncertain data
  mutate(mid_ma = (max_ma + min_ma) / 2,
         bin = bin_ages(mid_ma, by = 2))

# Generate time bins
coll_occ_taxa_perm_cret <- bin_time(occdf = coll_occ_taxa_perm_cret,
                                    bins = bins,
                                    method = 'all')

# check intervals
coll_occ_taxa_perm_cret %>%
  mutate(age_range = abs(max_ma - min_ma)) %>%
  ggplot(aes(x = age_range)) +
  geom_histogram() +
  labs(x = 'Age range (My)', y = 'Count') +
  theme_bw()


# load intervals
# find the interval of each raster
intervals <- openxlsx::read.xlsx (here ("data", "periods_intervals_ages.xlsx"))

# match bins to have interval name
coll_occ_taxa_perm_cret <- cbind (coll_occ_taxa_perm_cret,
                                  "Interval" = bins$interval_name [match (coll_occ_taxa_perm_cret$bin_assignment,
                                                                          bins$bin)]
)


# bind
coll_occ_taxa_perm_cret<-cbind (coll_occ_taxa_perm_cret,
                                "Period" = intervals$Period [match (coll_occ_taxa_perm_cret$Interval,
                                                                    intervals$Interval)]
)





# ---------------------------------------------------------------------
# check



# binned intervals per cell
hist (
  
  
  (rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                  coll_occ_taxa_perm_cret$bin
  )>0)
  ),
  main = "intervals per cell"
)


# formations per cell
hist (
  
  
  (rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                  coll_occ_taxa_perm_cret$formation.x
  )>0)
  ),
  main= "formations per cell"
)


# formations per interval
hist (
  
  
  (colSums(
    table (coll_occ_taxa_perm_cret$formation.x,
           coll_occ_taxa_perm_cret$bin
    )>0)
  ),
  main = "formations per interval"
)


# pull of the recent?
plot((colSums(
  table (coll_occ_taxa_perm_cret$formation.x,
         coll_occ_taxa_perm_cret$bin_assignment
  )>0)
),type="l",xlab= "intervals (older (0 - Capitanian) to the most recent (33 - Maastrichian))", ylab= "formations")



# --------------------------------------



# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)

cynodontia_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin)[order(unique(coll_occ_taxa_perm_cret$bin))]
cells <- unique(coll_occ_taxa_perm_cret$cell_ID)[order( unique(coll_occ_taxa_perm_cret$cell_ID))]
formations <- unique(coll_occ_taxa_perm_cret$formation.x)[order(unique(coll_occ_taxa_perm_cret$formation.x))]

# the basic table summarizing all detections across intervals and formations

table_data_basis <- lapply (cells, function (i) {
  
  
  tab_basis <- cast (formula = unique_name ~ bin,
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$cell_ID %in% i),],
                     value = "detection", 
                     fun.aggregate = sum,
                     drop=F,
                     fill=0)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$unique_name

    
  # input formations/subcells not in this data
  tab_basis <- cbind (tab_basis, 
                      matrix (0, 
                              nrow = nrow (tab_basis),
                              ncol= length(time_bins [which(time_bins  %in% colnames (tab_basis) == F)]), # subcells or formations
                              
                              dimnames = list (rownames(tab_basis),
                                               time_bins [which(time_bins  %in% colnames (tab_basis) == F)]
                                               
                              )
                      )
  )
  tab_basis <- tab_basis[,-which(colnames(tab_basis) == "unique_name")]
  
  # input cells not in this data
  tab_basis <- rbind (tab_basis, 
                      matrix (0, 
                              nrow = sum(cynodontia_data %in% rownames (tab_basis) == F),
                              ncol= ncol (tab_basis),
                              
                              dimnames = list (cynodontia_data [which(cynodontia_data  %in% rownames (tab_basis) == F)],
                                               colnames(tab_basis)
                                               
                              )
                      )
  )
  
  # order columns and rows
  tab_basis <- tab_basis[order (rownames(tab_basis)),
                         order(as.numeric(colnames(tab_basis)))]
  
  ; # return
  tab_basis
  
})

# checknames
colnames(table_data_basis[[1]]) == colnames(table_data_basis[[33]])
rownames(table_data_basis[[1]]) == rownames(table_data_basis[[33]])


# data to array
array_genus_bin_site <- array (unlist (table_data_basis),
                               dim = c(length(cynodontia_data),
                                       length(time_bins),
                                       length(table_data_basis)))

array_genus_bin_site[array_genus_bin_site>0]<-1


# array of genus per interval
array_genus_bin <- (apply (array_genus_bin_site, c(1,2), sum))
array_genus_bin[array_genus_bin>0]<-1

# plot observed SR per site
plot(NA,ylim=c(range(apply (array_genus_bin_site,c(2,3),sum))[1],
               150),
     xlim=c(1,33))

lapply (seq(1,dim(array_genus_bin_site)[3]), function (i)
  
  lines (apply (array_genus_bin_site,c(2,3),sum)[,i])
  
)
lines (colSums(array_genus_bin), col="red", lwd=2)




# number of formations per cell and interval
# the basic table summarizing all detections across intervals and formations

table_data_formations <- lapply (cells, function (i) {
  
  
  tab_basis <- cast (formula = bin ~ formation.x,
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$cell_ID %in% i),],
                     value = "detection", 
                     fun.aggregate = max,
                     drop=F,
                     fill=0)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$bin
  
  
  # input formations/subcells not in this data
  tab_basis <- cbind (tab_basis, 
                      matrix (0, 
                              nrow = nrow (tab_basis),
                              ncol= length(formations [which(formations  %in% colnames (tab_basis) == F)]), # subcells or formations
                              
                              dimnames = list (rownames(tab_basis),
                                               formations [which(formations  %in% colnames (tab_basis) == F)]
                                               
                              )
                      )
  )
  tab_basis <- tab_basis[,-which(colnames(tab_basis) == "bin")]
  
  # input cells not in this data
  tab_basis <- rbind (tab_basis, 
                      matrix (0, 
                              nrow = sum(time_bins %in% rownames (tab_basis) == F),
                              ncol= ncol (tab_basis),
                              
                              dimnames = list (time_bins [which(time_bins  %in% rownames (tab_basis) == F)],
                                               colnames(tab_basis)
                                               
                              )
                      )
  )
  
  # order columns and rows
  tab_basis <- tab_basis[order (as.numeric(rownames(tab_basis))),
                         order((colnames(tab_basis)))]
  
  # sum the number of formations per interval
  N_formation_interval <- rowSums(tab_basis)
  
  ; # return
  
  N_formation_interval
  
})
# melt
# formations per site and interval
formations_per_site_interval<-do.call(rbind,table_data_formations)


# get the number of formations per interval

formations_per_interval <- colSums(formations_per_site_interval)

# pull of the recent?
plot(formations_per_interval,type= "b")


# save 
save (coll_occ_taxa_perm_cret,
      array_genus_bin, 
      array_genus_bin_site,
      formations_per_interval,
      formations_per_site_interval,
      cynodontia_data,
      time_bins ,
      cells ,
      formations,
      file= here("processed_data", "CMR_data_smaller_bins.RData")
      )
rm(list=ls())

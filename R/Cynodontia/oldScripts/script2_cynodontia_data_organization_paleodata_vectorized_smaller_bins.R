# ======================================================

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


# check preservation quality
table(coll_occ_taxa_perm_cret$preservation_quality.x)
# filter(coll_occ_taxa_perm_cret, preservation_quality.x != "very poor")

# remove poor quality
#coll_occ_taxa_perm_cret<-coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$preservation_quality.x != "poor"),]


# change the order of bin assignment (more recent first)
# coll_occ_taxa_perm_cret$bin_assignment<-40-(coll_occ_taxa_perm_cret$bin_assignment)



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


# bind the interaction between formation and lithology
coll_occ_taxa_perm_cret$formation_lithology <-  paste (coll_occ_taxa_perm_cret$formation.x,
                                                       coll_occ_taxa_perm_cret$lithology1.x,
                                                       sep="_")





# ------------------------------------------

# save output to PyRate

PyRate_data <- coll_occ_taxa_perm_cret

# replace  blanks in taxon names
PyRate_data$unique_name <- gsub("[[:blank:]]{1,}","_", PyRate_data$unique_name)


#simulated current status, soley for demonstration purposes, replace with your own data
mock_status <- data.frame(unique_name = unique(PyRate_data$unique_name),
                          status = "extinct")


#add current status to fossils
PyRate_data2 <- inner_join(PyRate_data, mock_status, by = "unique_name")



# Write PyRate input to disk
# https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_PBDB_fossils_with_CoordinateCleaner.html
write_pyrate(PyRate_data, fname = "processed_data/pyrate_data_small_bins", status = PyRate_data2$status,
             taxon = "unique_name", min_age = "min_ma", max_age = "max_ma")


# data to divDyn
divDyn_data <- coll_occ_taxa_perm_cret %>%
  select (unique_name, bin,#bin_assignment, 
          Interval,
          min_ma, max_ma, collection_no, paleolat2.x,paleolng2.x)

# write
write.csv(divDyn_data, file = here ("processed_data", "divDyn_data_small_bins.csv"))


# --------------------------------------




# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)



cynodontia_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin)[order(unique(coll_occ_taxa_perm_cret$bin))]
cells <- unique(coll_occ_taxa_perm_cret$cell_ID)[order( unique(coll_occ_taxa_perm_cret$cell_ID))]
formations <- unique(coll_occ_taxa_perm_cret$formation.x)[order(unique(coll_occ_taxa_perm_cret$formation.x))]
formations_lithologies <- unique(coll_occ_taxa_perm_cret$formation_lithology )[order( unique(coll_occ_taxa_perm_cret$formation_lithology ))]
#subcells <- unique(coll_occ_taxa_perm_cret$cell_ID_sub)[order(unique(coll_occ_taxa_perm_cret$cell_ID_sub))]


# the basic table summarizing all detections across intervals and formations

table_data_basis <- lapply (time_bins, function (i) {
  
    
  tab_basis <- cast (formula = cell_ID ~ formation.x,#cell_ID_sub,#
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$bin %in% i),],
                     value = "detection", 
                     fun.aggregate = sum,
                     drop=F,
                     fill=NA)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$cell_ID
  
  
  # input formations/subcells not in this data
  tab_basis <- cbind (tab_basis, 
              matrix (NA, 
                      nrow = nrow (tab_basis),
                      ncol= sum(formations %in% colnames (tab_basis) == F), # subcells or formations
                      
                      dimnames = list (rownames(tab_basis),
                                       formations [which(formations  %in% colnames (tab_basis) == F)]
                                       
                      )
              )
  )
  tab_basis <- tab_basis[,-which(colnames(tab_basis) == "cell_ID")]
  # input cells not in this data
  tab_basis <- rbind (tab_basis, 
              matrix (NA, 
                      nrow = sum(cells %in% rownames (tab_basis) == F),
                      ncol= ncol (tab_basis),
                      
                      dimnames = list (cells [which(cells  %in% rownames (tab_basis) == F)],
                                       colnames(tab_basis)
                                       
                      )
              )
  )
  
  # order columns and rows
  tab_basis <- tab_basis[order (rownames(tab_basis)),
                                order(colnames(tab_basis))]
  
  ; # return
  tab_basis
                
})

# checknames
colnames(table_data_basis[[1]]) == colnames(table_data_basis[[33]])
rownames(table_data_basis[[1]]) == rownames(table_data_basis[[33]])



# order cols
# arrange colnames by time
#table_data_basis <- table_data_basis[,match(bins$bin ,
#                                      colnames (table_data_basis))]
# there is no observation in 88% of the combinations between intervals and formations
# table(is.na(table_data_basis))/sum(table(is.na(table_data_basis)))
# summation of detections (each entry in the database is a detection)
#range(table_data_basis,na.rm=T)
# detections
# hist(as.numeric(data.matrix(table_data_basis))[is.na(as.numeric(data.matrix(table_data_basis)))!=T])

# since 1 to 1775 detections
range(as.numeric(data.matrix(table_data_basis[[1]]))[is.na(as.numeric(data.matrix(table_data_basis[[1]])))!=T])


# k = cynodontia_data[67]
# i = time_bins[26]
# 18

# do the table per interval and taxon

# table 
table_data <- lapply (cynodontia_data, function (k) # each species
  
                        
                     lapply (time_bins, function (i) { # each time bin
                      
                         
                         
                        tryCatch({
                                 
                              
                              
                         # basis data 
                         basis_data <- coll_occ_taxa_perm_cret %>% 
                           filter (unique_name %in% k) %>%
                           filter (bin %in% i)
  
                          # if no data is available for the k genus at time bin i, return an  matrix filled with NAs
                         if (nrow (basis_data) == 0) {
                           
                           table_data <- matrix (NA,
                                                nrow =length(cells),
                                                ncol =length(formations),
                                                dimnames = list (cells,
                                                                 formations))
                           
                         } else {
                         # pivot table
                           table_data <- cast (formula = cell_ID ~ formation.x,
                                           data = basis_data,
                                           value = "detection", 
                                           fun.aggregate = sum,
                                           drop=F,
                                           fill=NA)
                           
                           # problems with only one col and row
                           if (sum(colnames (table_data) != "cell_ID")== 1) {
                             name <- colnames (table_data)[2]
                             table_data<-data.frame (name = table_data[,2],
                                                     row.names = table_data$cell_ID)
                             colnames (table_data) <- name
                             
                           } else {
                            rownames (table_data)<- table_data$cell_ID 
                            table_data<- (table_data[,-1])
                           }
                         }
            
                         
                        # if just one column, set the name
                        # other wise it will be a dataframe
                      
                        #if (nrow (table_data) == 1 & ncol (table_data) ==1) {
                        #  name <- colnames (table_data)[2]
                        #  table_data<-data.frame (name = table_data[,2],
                        #                          row.names = table_data$cell_ID)
                        #  colnames (table_data) <- name
                        #  
                        #} else {
                        #
                        #  table_data
                        #
                        #}
                        
                        # bind missing formations (zeros; true absences)
                        input_formations <- matrix (NA, 
                                                    
                                                  nrow = nrow (table_data),
                                                  ncol = sum (formations %in% colnames(table_data) != T),
                                                  byrow=T,
                                                  
                                                  dimnames = list (
                                                      
                                                    rownames(table_data),
                                                    formations[which(formations %in% colnames (table_data) == F)]
                                                      
                                                       
                                                      )
                                          )
                        # bind 
                        table_data_formations <- cbind (table_data,
                                                        input_formations)
                        
                        # input cells (NAs, not assessed)
                        input_cells <- matrix (NA, 
                                                    
                                                   nrow = sum(cells %in% rownames (table_data_formations) != T),
                                                   ncol = ncol (table_data_formations),
                                                    
                                                    byrow=T,
                                                    
                                                    dimnames = list (
                                                      
                                                      cells [which(cells %in% rownames (table_data_formations) != T)],
                                                      colnames (table_data_formations)
                                                    )
                        )
                        
                        # bind 
                        table_data_formations_cells <- rbind (table_data_formations,
                                                              input_cells)
                        
                        
                        # order of formations
                        table_data_formations_cells <- table_data_formations_cells[order (rownames(table_data_formations_cells)),
                                                                                   order (colnames(table_data_formations_cells))]
                        
                        
                        ; # return
                        
                        table_data_formations_cells
                            }, # close tryCatch
                        
                        
                        error = function(e) return ("NULL"))
       


  }
)
)

table(table_data[[1]][[2]]==1)


# data to long format (each row one observation per time bin, taxon, formation, cell)

# check nrows and cols
# do.call(rbind, lapply (seq(1,length (table_data)), function (i)
#   range(unlist(lapply (table_data[[i]],nrow)))))
# lapply (table_data[[470]],nrow)

# matching with basis table
# if no taxa was detected in any time, set NA
# if at least one taxa was detected in a time, set zero

# i = 477
# k=24

table_data_long <- lapply (seq(1,length(table_data)), function (k)  # each species
  
              lapply (seq (1,length(table_data [[k]])), function (i) { # each time
  
                
                tryCatch({
                               
                     
                    # df frame with detections of each taxon and detections in the complete basis dataset 
                    test <- data.frame (det_taxon=as.numeric(unlist(table_data[[k]][[i]])),
                                        det_all_taxon=(as.numeric(unlist(table_data_basis[[i]]))))
                    
                    #test[which(test$det_taxon ==1),]
                    
                    # matching
                    test$test <- ifelse (is.na (test$det_all_taxon) &  is.na(test$det_taxon), # if both are NA
                            NA, # NA is the output
                            ifelse ((is.na (test$det_all_taxon) != T) &  is.na(test$det_taxon), # if at least one taxon was detected in any time and formation
                                    0, # zero is output (non-detection)
                                    test$det_taxon # otherwise the data will be filled with ones
                    ))
                  
                    # matrix
                    matrix_data <- matrix(as.numeric(test$test),
                                           nrow=nrow(table_data_basis[[1]]),
                                           ncol=ncol(table_data_basis[[1]]),
                                           byrow = F,
                                           dimnames = list(seq(1,length(rownames(table_data_basis[[1]]))),
                                                           seq(1,length(colnames(table_data_basis[[1]])))
                                                           )
                                          )
                  
                    # check if tables match
                    # table (is.na(matrix_data) == is.na(table_data_basis[[1]]))
                    
                    # melt to long format
                    table_data_long <- melt(matrix_data,as.is=T)
                    
                    # range(table((table_data_long$X1))) # check if it's ok
                    colnames(table_data_long) <- c("site", "form", "det")
                    
                    # removing NAs
                    table_data_long<- table_data_long[is.na(table_data_long$det)!=T,]
                    
                    # bind interval
                    table_data_long$int <- i
                    
                    # bind taxon 
                    table_data_long$taxon <- k
                    
                    # need to input zeros in the begging to help the model
                    #table_data_long <- rbind (data.frame (form = 500,
                    #                          int=seq (1,7),
                    #                          det=0),
                    #                          table_data_long)
                    
                    
                    ; # return
                    table_data_long }, # close tryCatch
                    
                    
                    error = function(e) return ("NULL"))
  
  })
)

# check formations
formations == colnames(table_data_basis[[1]])

#subcells[order(subcells)] == colnames(table_data_basis[[1]])

# bind lithology factor
#  to the table
table_data_long <- lapply (table_data_long, function (i) # for each spp
  
                      lapply (i , function (k) {# for each time bin

                        
                            k$lith <- formations [k$form]#sapply (strsplit (,"_"), "[[",2)
                          
                          
                          
                          ; # return
                          
                          k


                      }))


# melt at the first level (time bins)
table_data_long <- lapply (table_data_long, function (i)
  
              
                do.call(rbind , i)
                
)
    

# detection/non-detection data relative to the total N of possible observations
#mean((unlist(lapply(table_data_long,nrow))))/(nrow (table_data_basis)*ncol(table_data_basis))
# detection/non-detection data
#range(unlist(lapply(table_data_long,nrow)))

# detections of each genus
hist(unlist(
lapply (table_data_long, function (i)
  
  sum(i$det)
)
)
)


# find the number of detections per genus
#par(mar = rep(4,4))

# number of genus with 5 or more detections
#table(unlist(lapply (
  
#  table_data_long, function (i)
    
#    sum(i$det)
#)
#) >=2
#)

# select detection in at least one period and formation
sel_cynodontia <- unlist(lapply (
  
  table_data_long, function (i)
    
    sum(i$det)
)>=5
)
barplot(table(sel_cynodontia))



# table  for naive occupancy
table_naive_cynodontia <- table_data [which(sel_cynodontia==T)]
save (table_naive_cynodontia,  file = here ("processed_data", "table_naive_cynodontia_small_bins.RData"))


# filter
table_data_long <- table_data_long [which(sel_cynodontia==T)]

# melt
table_data_long_df <- do.call(rbind,table_data_long) # [1:2]

# match new names
df_match <- data.frame (names=unique(table_data_long_df$taxon),
            new_names=seq(1,length(unique(table_data_long_df$taxon))))
table_data_long_df$taxon_code <- df_match$new_names[match (table_data_long_df$taxon,df_match$names)]


naive<-cast (formula=taxon_code ~ int,
      data = table_data_long_df,
      value = "det",
      fun.aggregate = sum)
plot(rev(colSums(naive[,-1]>0)),type="b")


# edit lithologies
#table_data_long_df <-table_data_long_df %>%
#  
#  mutate (lith2 = recode (lith, 
#                         "mudstone"=1,
#                         #"NA"="NA",
#                         "siliciclastic"=2,
#                         "sandstone"=3,
#                         "claystone"=4,
#                         "siltstone"=5,
#                         "shale"=6,
#                         "conglomerate"=7,
#                         "marl"=8,
#                         "dolomite"=9,
#                         "mixed carbonate-siliciclastic"=10,
#                         "breccia"=11,
#                         "coal"=12,
#                         "framestone"=13,
#                         "carbonate"=14,
#                         "limestone"=15,
#                         "lime mudstone"=16,
#                         "tuff"=17))
#


# data checks


# data onto 0,1
table_data_long_df$det[which(table_data_long_df$det >0)]<-1

# observed data
obs_data <- (cast (data = table_data_long_df,
                   formula = int ~taxon_code,
                   value = "det",
                   fun.aggregate = max))
# time trend
data.frame (ntaxa=rev(rowSums(obs_data)),
            int = seq(1,33)) %>%
  ggplot (aes(x=int,y=ntaxa))+
  geom_point()+
  geom_line()+
  geom_smooth()+
  theme_bw()





# calculando o percentual de zeros nas réplicas dos sítios q tem pelo menos uma ocorrência

sites <- unique(table_data_long_df$site)[order(unique(table_data_long_df$site))]
spp<- unique(table_data_long_df$taxon)[order(unique(table_data_long_df$taxon))]


# calculate
perc_rep <- lapply (sites, function (i)
  
  lapply (spp, function (k)
    
    table (table_data_long_df [which(table_data_long_df$site == i
                                     &
                                       table_data_long_df$taxon == k)
                               
                               ,"det"])/sum(table (table_data_long_df [which(table_data_long_df$site == i
                                                                             &
                                                                               table_data_long_df$taxon == k)
                                                                       
                                                                       ,"det"]))
    
  )
  
)

# plot 
hist (unlist (
  lapply (perc_rep, function (i)
    
    do.call(rbind,i [which(unlist(lapply (i, length))>1)])[,1]
  )),
  xlab= "Proportion of zeros in sites with detections", main="")



# ---------------------------------------------------

# save to use in models

# save to use in models
save (coll_occ_taxa_perm_cret,
      cynodontia_data,
      sel_cynodontia,
      table_data_basis,
      table_data_long,
      table_data_long_df,
      formations,
      cynodontia_data, 
      time_bins, 
      cells,
      file = here ("processed_data","table_data_array_cynodontia_smaller_bins.RData"))


# save data to use in spatial analyses
save (coll_occ_taxa_perm_cret,
      grid_info,
      time_bins, 
      cells,
      file = here ("processed_data","table_space_smaller_bins.RData"))


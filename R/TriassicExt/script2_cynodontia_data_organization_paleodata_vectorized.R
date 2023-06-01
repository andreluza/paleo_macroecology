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
bins <- time_bins(interval = c("Permian", "Cretaceous"), rank = "stage")



# the max_ma.y and min_ma.y will be the original ages
#colnames(coll_occ_taxa)[which(colnames(coll_occ_taxa) == "max_ma.x")] <- "max_ma"
#colnames(coll_occ_taxa)[which(colnames(coll_occ_taxa) == "min_ma.x")] <- "max_ma"
#
## Get new numeric ages for named intervals using the interval key that is supplied with Palaeoverse
#coll_occ_taxa <- look_up(coll_occ_taxa, 
#                     early_interval = "early_interval.x",
#                     late_interval = "late_interval.x")
#
## Make sure that any values which could not be matched contain their original values
#coll_occ_taxa$interval_max_ma <- ifelse(is.na(coll_occ_taxa$interval_max_ma),
#                                        coll_occ_taxa$max_ma,
#                                        coll_occ_taxa$interval_max_ma)
#coll_occ_taxa$interval_min_ma <- ifelse(is.na(coll_occ_taxa$interval_min_ma),
#                                        coll_occ_taxa$min_ma,
#                                        coll_occ_taxa$interval_min_ma)
#coll_occ_taxa$interval_mid_ma <- (coll_occ_taxa$min_ma + coll_occ_taxa$max_ma)/2
coll_occ_taxa$interval_mid_ma <- (coll_occ_taxa$min_ma.x + coll_occ_taxa$max_ma.x)/2

# rename
#colnames(coll_occ_taxa)[which(colnames (coll_occ_taxa) %in% c("max_ma", "min_ma"))] <- c("old_max_ma", "old_min_ma")
#colnames(coll_occ_taxa)[which(colnames (coll_occ_taxa) %in% c("interval_max_ma", "interval_min_ma"))] <- c("max_ma", "min_ma")



# we're interested in permian to triassic data
# Remove occurrences that are younger than the time intervals we're interested in
coll_occ_taxa_perm_cret <- subset(coll_occ_taxa, min_ma.x > min(bins$min_ma))


# rename cols to feed the binning function
# the max_ma.y and min_ma.y will be the original ages
colnames(coll_occ_taxa_perm_cret)[which(colnames(coll_occ_taxa_perm_cret) == "max_ma.x")] <- "max_ma"
colnames(coll_occ_taxa_perm_cret)[which(colnames(coll_occ_taxa_perm_cret) == "min_ma.x")] <- "min_ma"


# Generate time bins
coll_occ_taxa_perm_cret <- bin_time(occdf = coll_occ_taxa_perm_cret,
                      bins = bins,
                      method = 'all')

# Make a table of potential number of bins
table(coll_occ_taxa_perm_cret$n_bins) # raw counts
table(coll_occ_taxa_perm_cret$n_bins) / nrow(coll_occ_taxa_perm_cret) # proportions


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
                              plot=T,return=T)

# object with the grid
grid_info <- coll_occ_taxa_perm_cret [2:4]
plot(grid_info$grid_base)
plot(grid_info$grid,add=T, col = rgb(0.1,0.1,0.5,alpha=0.2))

# occurrence data
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret$occdf


# cell coordinates 
cell_coordinates <- coll_occ_taxa_perm_cret %>% 
  group_by (cell_ID) %>%
  summarise (centroid_lat = mean (cell_centroid_lat),
             centroid_lng = mean(cell_centroid_lng),
             sd_lat = sd (cell_centroid_lat),
             sd_lng = sd (cell_centroid_lng))




# load spatial covariates 
load (here ("processed_data", "brick_rasters_elevation.RData"))
load (here ("processed_data", "brick_rasters_precipitation.RData"))
load (here ("processed_data", "brick_rasters_temp.RData"))


# load intervals
# find the interval of each raster
intervals <- openxlsx::read.xlsx (here ("data", "periods_intervals_ages.xlsx"))
# reverse the order to the past be the first period
rev_order <- rev(intervals$Interval)
intervals <- intervals [match (intervals$Interval, rev_order),]


# extract covariates for each cell

site_covs <- lapply (list(brick_rasters,
                        brick_rasters_paleoprec,
                        brick_rasters_paleotemp), function (layer) {
                          
                          # extract
                          
                          df_var <- (extract(
                              
                                  layer,
                                  st_as_sf(cell_coordinates[,c("centroid_lng","centroid_lat")],
                                        coords = c("centroid_lng", "centroid_lat"),
                                        crs = "+proj=longlat +datum=WGS84 +no_defs"),
                               
                               fun = "mean",
                               method = "simple",
                               cellnumbers=T,
                               df=T))
                          colnames(df_var) <- gsub ("index_", "",colnames(df_var)) # adjust colnames
                          rownames (df_var) <- cell_coordinates$cell_ID
                          ; #return
                          df_var
                          
                          }
                     )

            
names (site_covs) <- c("elevation", "prec", "temp")


# bind missing data

site_covs_input <- lapply (site_covs, function (i) { 
  
  # if not in the colnames, take the previous one
  
  sel_int <- intervals$Interval [-which(intervals$Interval %in% colnames (i)[-c(1:2)])]
  
  # for each missing interval made the input
  to_bind <- lapply (seq (1,length(sel_int)), function (int) {
    
    # interval to select
    sel_col <- intervals$Interval[which(intervals$Interval %in% sel_int[int])-1]
    to_bind <- i[which(colnames(i) %in% sel_col)]
    
    
    if (ncol (to_bind) > 0 ) {
      
      colnames(to_bind) <- sel_int[int]
      
      
    } else { # if stil not in the colnames, regret two time steps
      
      sel_col <- intervals$Interval[which(intervals$Interval %in% sel_int[int])-2]
      to_bind <- i[which(colnames(i) %in% sel_col)]
      colnames(to_bind) <- sel_int[int]
      
    }
    
    
    #colnames(to_bind) <- sel_int[int]
    # try also with NAs
    to_bind_NAs <- to_bind
    to_bind_NAs  <- ifelse (to_bind > -9999,NA,NA)
    output <- list (input_vals = to_bind,
                    input_NAs = to_bind_NAs)
    output
    })
  
  
})

# input elevation
site_covs[[1]] <- cbind (site_covs[[1]],
                         do.call(cbind, sapply (site_covs_input[[1]], "[[", "input_vals")))[,-c(1:2)]
site_covs[[1]] <- site_covs[[1]] [match (intervals$Interval, colnames(site_covs[[1]]))]

# input precipitation
site_covs[[2]] <- cbind (site_covs[[2]],
                         do.call(cbind, sapply (site_covs_input[[2]], "[[", "input_vals")))[,-c(1:2)]
site_covs[[2]] <- site_covs[[2]] [match (intervals$Interval, colnames(site_covs[[2]]))]

# input temp
site_covs[[3]] <- cbind (site_covs[[3]],
                         do.call(cbind, sapply (site_covs_input[[3]], "[[", "input_vals")))[,-c(1:2)]
site_covs[[3]] <- site_covs[[3]] [match (intervals$Interval, colnames(site_covs[[3]]))]

# -----------------------------------

#  observation covariates

# i= 1594 
# i = 1000
# i = 163


# extract altitude per occurrence point
paleo_env <- lapply (list(brick_rasters,
                          brick_rasters_paleoprec,
                          brick_rasters_paleotemp), function (layer)
      
      lapply (seq (1,nrow (coll_occ_taxa_perm_cret)), function (i) {
  
      # extract
      extracted_data <- (extract(layer,
                                coll_occ_taxa_perm_cret[i,c("paleolng2.x","paleolat2.x")],
                                fun = "mean",
                                method = "simple"))
      colnames(extracted_data) <- gsub ("index_", "",colnames(extracted_data)) # adjust colnames
      
      # select
      selected_bin_n <- coll_occ_taxa_perm_cret[i,"bin_assignment"]
      selected_bin <- bins[which(bins$bin %in% selected_bin_n),"interval_name"]
      output <- (extracted_data[,which(colnames(extracted_data) %in% selected_bin)])
      
      
      # if no output is reported
      if (length(output)==0){
        
        # gather the paloaltitude from the previous period in the same point
        selected_bin_n <- coll_occ_taxa_perm_cret[i,"bin_assignment"]
        selected_bin <- bins[which(bins$bin %in% (selected_bin_n-1)),"interval_name"]
        output <- (extracted_data[,which(colnames(extracted_data) %in% selected_bin)])
        
        
      } else { output }
      
      
      # if still the previous one misses data, regret two periods
      if (length(output)==0){
          
          # gather the paloaltitude from the previous period in the same point
          selected_bin_n <- coll_occ_taxa_perm_cret[i,"bin_assignment"]
          selected_bin <- bins[which(bins$bin %in% (selected_bin_n-2)),"interval_name"]
          output <- (extracted_data[,which(colnames(extracted_data) %in% selected_bin)])
          
          
        } else {
          
        
        output
        
        }
        
      # if still the previous one misses data, regret three periods
      if (length(output)==0){
        
        # gather the paloaltitude from the previous period in the same point
        selected_bin_n <- coll_occ_taxa_perm_cret[i,"bin_assignment"]
        selected_bin <- bins[which(bins$bin %in% (selected_bin_n-3)),"interval_name"]
        output <- (extracted_data[,which(colnames(extracted_data) %in% selected_bin)])
        
        
      } else {
        
        
        output
        
      }
        
  })
)
#paleoaltitude[sapply(paleoaltitude, function(x) length(x)==0L)] <- NULL


# all filled with data
which(unlist(lapply (paleo_env[[3]], function (i) length(i) == 0))==T)
which(unlist(lapply (paleo_env[[3]], function (i) is.na(i) == T))==T)


# melt and bind

coll_occ_taxa_perm_cret$elevation <- do.call(rbind,paleo_env[[1]])[,1]
coll_occ_taxa_perm_cret$precipitation <- do.call(rbind,paleo_env[[2]])[,1]
coll_occ_taxa_perm_cret$temperature <- do.call(rbind,paleo_env[[3]])[,1]

# mean(coll_occ_taxa_perm_cret$elevation[which(coll_occ_taxa_perm_cret$elevation>-1000)])
# sd(coll_occ_taxa_perm_cret$elevation[which(coll_occ_taxa_perm_cret$elevation>-1000)])

# some elevations are lower than zero
table(coll_occ_taxa_perm_cret$elevation<0)

# number of formations per cellid
hist(colSums(table (coll_occ_taxa_perm_cret$formation.x, coll_occ_taxa_perm_cret$cell_ID)>0))
mean(colSums(table (coll_occ_taxa_perm_cret$formation.x, coll_occ_taxa_perm_cret$cell_ID)>0))# mean

# Extract unique interval midpoints
midpoints <- sort(unique(coll_occ_taxa_perm_cret$bin_midpoint))

# Calculate the number of unique cells in each time bin
unique_cells <- unique(coll_occ_taxa_perm_cret[, c("bin_midpoint","cell_ID")])
spat.cov <- group_apply(unique_cells, group = "bin_midpoint", fun = nrow)



# taxonomic resolution 
# species level : Occurrences identified to a coarser taxonomic resolution than the desired level
# are retained if they belong to a clade which is not otherwise represented in the dataset 
coll_occ_taxa_perm_cret <- tax_unique(occdf = coll_occ_taxa_perm_cret, 
                      genus = "genus", 
                      family = "family",
                      order = "order", 
                      class = "class", 
                      name = "accepted_name",
                      resolution="species",
                      append=T)

# changes
length(unique(coll_occ_taxa_perm_cret$accepted_name))
length(unique(coll_occ_taxa_perm_cret$unique_name))

# some mistakes (only " sp.")
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$unique_name %in% " sp." ==F),]
# remove NAs
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret [is.na(coll_occ_taxa_perm_cret$unique_name)==F,]

# Get the names of unique genera per collection
unique_genera <- unique(coll_occ_taxa_perm_cret[, c("unique_name", "formation.x")]) 

# Calculate the number of unique genera per collection
coll_taxa <- group_apply(unique_genera, group = "formation.x", fun = nrow)

# Rename column names:
colnames(coll_taxa) <- c("n_taxa", "formation.x")

# Take the columns pertaining to collections and their ages in millions of years:
coll_info <- coll_occ_taxa_perm_cret[, c("formation.x",
                               "max_ma", 
                               "interval_mid_ma", 
                               "min_ma")]

# Remove duplicated collections based on the collection number (column 1)
coll_info <- coll_info[!duplicated(coll_info[1]),]

# Combine this dataframe with the dataframe from above
alpha_data <- merge(coll_info, coll_taxa, by = "formation.x")
# Take a look:
head(alpha_data)

# uses the cynodontia diversity data from above
# help here https://cran.r-project.org/web/packages/deeptime/vignettes/coord_geo.html
ggplot(alpha_data,aes(x = (interval_mid_ma),
                      y = n_taxa)) +
  
  geom_point() +
  
  geom_smooth(method = "loess", span = 0.8, alpha = 0.2) +
  ylab("Cynodontia species") +
  coord_geo(
    dat = list("periods", "eras"), 
    xlim = c(max(alpha_data$interval_mid_ma)+10, 
             min(alpha_data$interval_mid_ma)-10), 
    ylim = c(0, max(alpha_data$n_taxa)+10),
    pos = list("b", "b"), abbrv = list(TRUE, FALSE)
  ) +
  scale_x_reverse("Age (Ma)") +
  theme_classic() 


# exploring range size
# First, remove any occurrences without a genus
space_coll_occ_taxa_perm_cret <- subset(coll_occ_taxa_perm_cret, !is.na(unique_name)) # accepted_name


# Find temporal range of all genera
space_coll_occ_taxa_perm_cret <- group_apply(occdf = space_coll_occ_taxa_perm_cret,
                               group = c("bin_midpoint"),
                               fun = tax_range_space,
                               name = "unique_name",
                               lng = "paleolng2.x",
                               lat = "paleolat2.x",
                               method = "con")


# Have a look at the dataset
head(space_coll_occ_taxa_perm_cret)


# Find the average geographic range size for each time interval
space_coll_occ_taxa_perm_cret_mean <- group_apply(space_coll_occ_taxa_perm_cret, 
                                        "bin_midpoint",
                                        function(df) 
                                          mean(df$area))

colnames(space_coll_occ_taxa_perm_cret_mean) <- c("mean_area", "bin_midpoint")

# Create a plot of average range size through time
ggplot(space_coll_occ_taxa_perm_cret_mean,
       aes(x = as.numeric(bin_midpoint),
                      y = log(mean_area+1))) +
  
  geom_point() +
  
  geom_smooth(method = "loess", span = 0.8, alpha = 0.2) +
  ylab("Cynodontia species' range size") +
  coord_geo(
    dat = list("periods", "eras"), 
    xlim = c(max(as.numeric(space_coll_occ_taxa_perm_cret_mean$bin_midpoint))+10, 
             min(as.numeric (space_coll_occ_taxa_perm_cret_mean$bin_midpoint))-10), 
    #ylim = c(-100000, max(space_coll_occ_taxa_perm_cret_mean$mean_area)+100000),
    pos = list("b", "b"), abbrv = list(TRUE, FALSE)
  ) +
  scale_x_reverse("Age (Ma)") +
  theme_classic() 


# range size correlation with sample number 
cor.test(log(space_coll_occ_taxa_perm_cret_mean$mean_area+1),
         log(spat.cov$nrow), method = "spearman")


# data exploration according to feeding guidsl
# uses the coral occurrence data from above
coll_occ_taxa_perm_cret %>%
  filter(interval_mid_ma != "") %>%
  filter(diet != "") %>%
  mutate(diet=recode(diet, 
                    "durophage"="carnivore",
                    "insectivore" = "carnivore",
                    "piscivore" = "carnivore",
                    "carnivore, insectivore"="carnivore",
                    "carnivore, durophage"="carnivore",
                    "insectivore, herbivore" = "omnivore",
                    "omnivore, frugivore" = "omnivore",
                    "omnivore, carnivore" = "omnivore")) %>%
  group_by(diet, early_interval.x) %>%
  summarise(n = n_distinct(accepted_name),
            age = mean(interval_mid_ma)) %>%
  
  ggplot(aes(x = age, y = n)) +
        geom_line() +
        geom_point()+
        scale_x_reverse("Age (Ma)") +
        ylab("Cynodontia species") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c(260, 70), 
    ylim = c(0, 100),
    pos = list("b", "b"),
    size = list(4, 6),
    abbrv = list(TRUE, FALSE)
  ) +
        theme_classic() +
        facet_wrap(~diet, nrow = 3)





# ---------------------------------------------------------------------
# sites



# binned intervals per cell
hist (
  

  table(rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                       coll_occ_taxa_perm_cret$bin_assignment
       )>0)
       )
)


# formations per cell
hist (
  
  
  table(rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                       coll_occ_taxa_perm_cret$formation.x
  )>0)
  )
)


# formations per interval
hist (
  
  
  (colSums(
    table (coll_occ_taxa_perm_cret$formation.x,
           coll_occ_taxa_perm_cret$bin_assignment
  )>0)
  )
)




# --------------------------------------------------------------------



# create index as value in the table
coll_occ_taxa_perm_cret$detection <- 1



# match bins to have interval name
coll_occ_taxa_perm_cret<-cbind (coll_occ_taxa_perm_cret,
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


# --------------------------------------------------------------------

# select interesting cols (to downsize the dataset)

coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[,  c("unique_name", 
                                                        "genus",
                                                        "formation.x",
                                                        "lithology1.x",
                                                        "elevation",
                                                        "precipitation",
                                                        "temperature",
                                                        "bin_assignment",
                                                        "Interval",
                                                        "Period",
                                                        "min_ma",
                                                        "max_ma",
                                                        "cell_ID", 
                                                        "paleolng2.x",
                                                        "paleolat2.x",
                                                        "detection"
)]


# plot table
table_taxon_period <- cast (unique_name ~ bin_assignment,
                            data = coll_occ_taxa_perm_cret,
                            value = "detection",
                            fun.aggregate = sum)
rownames(table_taxon_period)<- table_taxon_period$unique_name; table_taxon_period<- table_taxon_period[,-1]
table_taxon_period[table_taxon_period>1] <- 1
# melt data to plot
melt_data <- melt (data.matrix(table_taxon_period))

#
melt_data %>%
  filter (value >0) %>%
  group_by(X1) %>%
  summarise (min_int = min (X2),
             max_int = max (X2)) %>%
  mutate (diff = max_int - min_int) %>%
  ggplot() +
  geom_linerange(mapping=aes(x = reorder(X1, diff), 
                             ymin = min_int, 
                             ymax = max_int), 
                 size = 1, alpha = 0.5,
                 position = position_dodge(width = 0.1)) +
  geom_point (aes(x = reorder(X1, diff), 
                  y = min_int),col = "orange",size=2)+
  geom_point (aes(x = reorder(X1, diff), 
                  y = max_int),col = "green3")+
  theme (axis.text.y =  element_text(size=1.5))+
  coord_flip()+
  xlab ("Taxon")+
  ylab ("Time bin")

ggsave (file = here ("output", "spp_range.pdf"))


# last check in time
tail(table(coll_occ_taxa_perm_cret$min_ma))
tail(table(coll_occ_taxa_perm_cret$max_ma))

ggplot(coll_occ_taxa_perm_cret)+
  geom_histogram(aes(x = max_ma)) + 
  geom_histogram(aes(x = min_ma),col="yellow", alpha =0.5) 

# --------------------------------------


# bind the interaction between formation and lithology
coll_occ_taxa_perm_cret$formation_lithology <-  paste (coll_occ_taxa_perm_cret$formation.x,
                                                       coll_occ_taxa_perm_cret$lithology1.x,
                                                       sep="_")




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
write_pyrate(PyRate_data, fname = "processed_data/pyrate_data", status = PyRate_data2$status,
             taxon = "unique_name", min_age = "min_ma", max_age = "max_ma")





# --------------------------------------




# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)



synapsida_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin_assignment)[order(unique(coll_occ_taxa_perm_cret$bin_assignment))]
cells <- unique(coll_occ_taxa_perm_cret$cell_ID)
formations <- unique(coll_occ_taxa_perm_cret$formation.x)[order(unique(coll_occ_taxa_perm_cret$formation.x))]
formations_lithologies <- unique(coll_occ_taxa_perm_cret$formation_lithology )[order( unique(coll_occ_taxa_perm_cret$formation_lithology ))]




# the basic table summarizing all detections across intervals and formations

table_data_basis <- lapply (time_bins, function (i) {
  
    
  tab_basis <- cast (formula = cell_ID ~ formation_lithology,
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$bin_assignment %in% i),],
                     value = "detection", 
                     fun.aggregate = sum,
                     drop=F,
                     fill=NA)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$cell_ID
  tab_basis <- tab_basis[,-1]
  
  # input formations not in this data
  tab_basis <- cbind (tab_basis, 
              matrix (NA, 
                      nrow = nrow (tab_basis),
                      ncol= sum(formations_lithologies %in% colnames (tab_basis) == F),
                      
                      dimnames = list (rownames(tab_basis),
                                       formations_lithologies [which(formations_lithologies  %in% colnames (tab_basis) == F)]
                                       
                      )
              )
  )
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


# k = synapsida_data[67]
# i = time_bins[26]
# 18

# do the table per interval and taxon

# table 
table_data <- lapply (synapsida_data, function (k) # each species
  
                        
                     lapply (time_bins, function (i) { # each time bin
                      
                         
                         
                        tryCatch({
                                 
                              
                              
                         # basis data 
                         basis_data <- coll_occ_taxa_perm_cret %>% 
                           filter (unique_name %in% k) %>%
                           filter (bin_assignment %in% i)
  
                          # if no data is available for the k genus at time bin i, return an  matrix filled with NAs
                         if (nrow (basis_data) == 0) {
                           
                           table_data <- matrix (NA,
                                                nrow =length(cells),
                                                ncol =length(formations_lithologies),
                                                dimnames = list (cells,
                                                                 formations_lithologies))
                           
                         } else {
                         # pivot table
                           table_data <- cast (formula = cell_ID ~ formation_lithology,
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
                                                  ncol = sum (formations_lithologies %in% colnames(table_data) != T),
                                                  byrow=T,
                                                  
                                                  dimnames = list (
                                                      
                                                    rownames(table_data),
                                                    formations_lithologies[which(formations_lithologies %in% colnames (table_data) == F)]
                                                      
                                                       
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
formations_lithologies == colnames(table_data_basis[[1]])

# bind lithology factor
#  to the table
table_data_long <- lapply (table_data_long, function (i) # for each spp
  
                      lapply (i , function (k) {# for each time bin

                        
                            k$lith <- sapply (strsplit (formations_lithologies [k$form],"_"), "[[",2)
                          
                          
                          
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
sel_synapsida <- unlist(lapply (
  
  table_data_long, function (i)
    
    sum(i$det)
)
) >=1
barplot(table(sel_synapsida))



# table  for naive occupancy
table_naive_synapsida <- table_data [which(sel_synapsida==T)]
save (table_naive_synapsida,  file = here ("processed_data", "table_naive_synapsida.RData"))


# filter
table_data_long <- table_data_long [which(sel_synapsida==T)]

# melt
table_data_long_df <- do.call(rbind,table_data_long) # [1:2]


# bind temperature to these data
# i = 987
# table_data_long_df[i,]
temperature_det <- lapply (seq (1,nrow (table_data_long_df)), function (i)


  mean(coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$cell_ID %in% rownames(table_data_basis[[1]]) [table_data_long_df [i,"site"]] 
                                  & 
                                   coll_occ_taxa_perm_cret$bin_assignment %in% time_bins [table_data_long_df [i,"int"]] 
                                 &
                                   coll_occ_taxa_perm_cret$formation_lithology %in% formations_lithologies [table_data_long_df [i,"form"]] 
                                   
                                   ),"temperature"])
)

# bind temperature
table_data_long_df <- cbind (table_data_long_df,
                             temp = unlist(temperature_det))


# edit lithologies
table_data_long_df <-table_data_long_df %>%
  
  mutate (lith2 = recode (lith, 
                         "mudstone"=1,
                         #"NA"="NA",
                         "siliciclastic"=2,
                         "sandstone"=3,
                         "claystone"=4,
                         "siltstone"=5,
                         "shale"=6,
                         "conglomerate"=7,
                         "marl"=8,
                         "dolomite"=9,
                         "mixed carbonate-siliciclastic"=10,
                         "breccia"=11,
                         "coal"=12,
                         "framestone"=13,
                         "carbonate"=14,
                         "limestone"=15,
                         "lime mudstone"=16,
                         "tuff"=17))



# ---------------------------------------------------

# save to use in models
save (site_covs,
      bins,
      file = here ("processed_data","site_covs.RData"))


# save to use in models
save (coll_occ_taxa_perm_cret,
      synapsida_data,
      sel_synapsida,
      table_data_basis,
      table_data_long,
      table_data_long_df,
      formations_lithologies,
      synapsida_data, 
      time_bins, 
      cells,
      file = here ("processed_data","table_data_array_synapsida.RData"))


# save data to use in spatial analyses
save (coll_occ_taxa_perm_cret,
      grid_info,
      time_bins, 
      cells,
      file = here ("processed_data","table_space.RData"))


# made animations to show climate and continents


# -----------------------------------------------
# produce animations to explore data

# paleomaps of formations - take a time ...
#maps_models <- getmap(ma = bins$mid_ma, 
#                      model = "GOLONKA")
#maps_models <- maps_models [order(bins$mid_ma)] # ordering
#save (maps_models, file = here ("processed_data", "paleomaps.RData")) 
load(here ("processed_data", "paleomaps.RData"))
names(maps_models) <- bins$interval_name

# plot points
lapply (bins$interval_name, function (i) {
  
  png (here ("processed_data", paste ("maps_", (i), ".png",sep ="")),
       width = 30,height = 15, res=300, 
       units = "cm")
  par (mar=rep (1,4))
  # triassic
  plot (maps_models[[which(names(maps_models) == (i))]],
        main = paste (i, " (",
                      bins[which( bins$interval_name == i),"mid_ma"],
                      " Ma)",sep=""),
        
        border = "gray",
        col = rgb (0.1,0.1,0.1,alpha=0.1))
  
  with (coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$Interval == i),], 
        
        points (paleolng2.x,
                paleolat2.x,
                pch=19,
                col = rgb(0.1,0.5,0,alpha=0.05))
  )
  
  
  dev.off()
})

#anime
list_img <- list.files(path = here ("processed_data"), full.names = T, pattern = "maps_")
# ordering
list_img <- list_img[match (names(ages),
                            gsub (".png", "",gsub ("D:/Pos_Doc_Paleonto_Macroecology/modeling/paleo_macroecology/processed_data/maps_","",list_img))
)]
##https://cran.r-project.org/web/packages/magick/vignettes/intro.html
a_image<-image_read(list_img)
animation <-  image_animate(a_image, fps = 1)
image_write(animation, here ("output","animation_data.gif"))




rm(list=ls())

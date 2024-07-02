
# ======================================================

# CMR design (species in the rows and intervals/time bin in the cols)
# covariates (file "site_covs.RData") will gathered from the code "script2_cynodontia_data_organization_paleodata_vectorized.R"

# organizing paleodata
# the metadata for all this
# https://taphonomy.doit.wisc.edu/data1.2/colls/single_doc.html
# geoplates
# https://www.gplates.org/docs/user-manual/Interacting_Features/#30-assign-plate-ids
# https://www.gplates.org/docs/user-manual/Reconstructions/#plate-ids


# =====================================================

# load packages
source("R/packages.R")
source("R/functions.R")

# --------------------------------------------------------------------

# load saved data
load(here ("processed_data", "PBDB_download.RData"))

# original size of the dataset
dim(pbdb_data_occ) # number of observations
length(unique(pbdb_data_occ$accepted_name))
length(unique(pbdb_data_collections$formation))
length(unique(pbdb_data_collections$collection_no))
max(pbdb_data_occ$max_ma)
min(pbdb_data_occ$min_ma)

# run a cleaning of coordinates
# help here : https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_PBDB_fossils_with_CoordinateCleaner.html

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

# result of data filtering
3904 - dim (coll_occ_taxa_perm_cret)[1]
3872 - dim (coll_occ_taxa_perm_cret)[1]

# --------------------------------------------------------------------


# adjust formation and lithology names
# remove quotes and other characters

coll_occ_taxa_perm_cret$formation <- gsub ("\\? ", "", coll_occ_taxa_perm_cret$formation)
coll_occ_taxa_perm_cret$formation <- gsub ("\"", "", coll_occ_taxa_perm_cret$formation)
coll_occ_taxa_perm_cret$formation <- gsub ("^[:alnum:]", "", coll_occ_taxa_perm_cret$formation)
coll_occ_taxa_perm_cret$formation <- gsub ("=", "", coll_occ_taxa_perm_cret$formation)
coll_occ_taxa_perm_cret$formation <-  (noquote(coll_occ_taxa_perm_cret$formation))
coll_occ_taxa_perm_cret$formation <- gsub("\"", '', coll_occ_taxa_perm_cret$formation)
coll_occ_taxa_perm_cret$formation <- gsub("'", '', coll_occ_taxa_perm_cret$formation)

# remove weird names
coll_occ_taxa_perm_cret<- coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$formation %in% c("", "=", "5","6") == F),]
unique(coll_occ_taxa_perm_cret$formation ) [order(unique(coll_occ_taxa_perm_cret$formation ))]

#as factor
coll_occ_taxa_perm_cret$formation <- as.factor (coll_occ_taxa_perm_cret$formation)

# recode
coll_occ_taxa_perm_cret <- data.frame (coll_occ_taxa_perm_cret) %>% 
  mutate (formation = recode (formation,
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




# --------------------------------------------------------

# add data from Agudotherium gassenae and Santagnathus mariensis
sp1 <- data.frame (matrix(NA,nrow=1,ncol = ncol(coll_occ_taxa_perm_cret),
                          dimnames = list(NA,
                                          colnames(coll_occ_taxa_perm_cret))))
#sp1$unique_name <- "Agudotherium sp."
sp1$accepted_name <- "Agudotherium gassenae"
sp1$identified_name <- "Agudotherium gassenae"
sp1$identified_rank <- "species"
sp1$accepted_rank <- "species"
sp1$formation <- "Santa Maria"
sp1$record_type.x <- "occ"
sp1$early_interval.x <- "Carnian"
sp1$early_interval.y <- "Carnian"
sp1$late_interval.x <- "Carnian"
sp1$late_interval.y <- "Carnian"
sp1$collection_name <- "Holotype CAPPA/UFSM 0262, a left lower jaw with teeth"
sp1$min_ma <- 227.000
sp1$interval_mid_ma <- 232.000
sp1$max_ma <- 237.000
sp1$class <- "Osteichthyes"
sp1$phylum <- "Chordata"
sp1$order <- "Therapsida"
sp1$genus <- "Agudotherium"
# Niemeyer, Agudo
sp1$paleolng2 <- coll_occ_taxa_perm_cret[grep ("Siriusgnathus",coll_occ_taxa_perm_cret$accepted_name),"paleolng2"][1]
sp1$paleolat2 <-  coll_occ_taxa_perm_cret[grep ("Siriusgnathus",coll_occ_taxa_perm_cret$accepted_name),"paleolat2"][1]
# set lat long
sp1$lng.x <- coll_occ_taxa_perm_cret[grep ("Siriusgnathus",coll_occ_taxa_perm_cret$accepted_name),"lng.x"][1]
sp1$lat.x <-  coll_occ_taxa_perm_cret[grep ("Siriusgnathus",coll_occ_taxa_perm_cret$accepted_name),"lat.x"][1]
# related publication:
# https://www.tandfonline.com/doi/full/10.1080/02724634.2020.1782415?fbclid=IwAR1zon2k-nfyXw3UnjS8vHqZpb_fr9mdrJQ9mnNIy0fyPlyC-9g3yjUfd7g
# rotate coordinates

# second species

# add data from Agudotherium gassenae and Santagnathus mariensis
sp2 <- data.frame (matrix(NA,nrow=1,ncol = ncol(coll_occ_taxa_perm_cret),
                          dimnames = list(NA,
                                          colnames(coll_occ_taxa_perm_cret))))
#sp2$unique_name <- "Agudotherium sp."
sp2$accepted_name <- "Santagnathus mariensis"
sp2$identified_name <- "Santagnathus mariensis"
sp2$identified_rank <- "species"
sp2$accepted_rank <- "species"
sp2$formation <- "Santa Maria"
sp2$record_type.x <- "occ"
sp2$early_interval.x <- "Carnian"
sp2$early_interval.y <- "Carnian"
sp2$late_interval.x <- "Carnian"
sp2$late_interval.y <- "Carnian"
sp2$collection_name <- "Holotype: Cranium with associated lower jaw, UFRGS-PV-1419-T"
sp2$min_ma <- 227.000
sp2$interval_mid_ma <- 232.000
sp2$max_ma <- 237.000
sp2$class <- "Osteichthyes"
sp2$phylum <- "Chordata"
sp2$order <- "Therapsida"
sp2$family <- "Traversodontidae"
sp2$genus <- "Santagnathus"
# coordinates : 29 41’43.21 S ; 53 47’45.19 W
sp2$lng.x <- -53.795886
sp2$lat.x <- -29.695336

#rotate
rotate_sp2 <- palaeorotate(
  sp2,
  lng = "lng.x",
  lat = "lat.x",
  age = "interval_mid_ma",
  model = "PALEOMAP",
  method = "point",
  uncertainty = TRUE,
  round = 3
)

# paste
sp2$paleolat2 <- rotate_sp2$p_lat
sp2$paleolng2 <- rotate_sp2$p_lng

# bind these occurrences to the dataset

coll_occ_taxa_perm_cret <- rbind (coll_occ_taxa_perm_cret,
                                  sp1,sp2)



# --------------------------------------------------------------------

# bin occurrences spatially into spatial bins


coll_occ_taxa_perm_cret <- bin_space(occdf = coll_occ_taxa_perm_cret %>%
                                       filter (is.na(paleolat2) != T),
                                     lng = 'paleolng2',
                                     lat = 'paleolat2',
                                     spacing = 2000,
                                     plot=T,return=T)


# object with the grid
grid_info <- coll_occ_taxa_perm_cret [2:4]

# change the class of the objects
grid_base <- (st_as_sf (grid_info$grid_base))
grid_data <- (st_as_sf (grid_info$grid))

# bind cell ID
grid_data <- cbind(grid_data,
           
           cell_ID = coll_occ_taxa_perm_cret$occdf$cell_ID )

# plot spatial data
plot(grid_info$grid_base,border="gray80")
plot(grid_info$grid,add=T, 
     border="gray80",
     col = rgb(0.1,0.1,0.5,alpha=0.1))
#plot(grid_info$sub_grid,add=T, col = rgb(0.5,0.1,0.1,alpha=0.2))

ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
        grid_data %>% 
          count (cell_ID),
        "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
      ), aes (fill=n),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis()


world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data = sf::st_transform(
    grid_data,
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ),aes (fill=cell_ID), alpha = 1) +
  
  geom_sf(data = sf::st_transform(
    grid_data[3635,],
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), alpha = 1,fill= "red") +

  geom_sf(data = world, alpha = 0.1) +
  scale_fill_viridis_d() +
  #ggtitle('H3 hexagons over County Ashe, NC', subtitle = 'Resolutions 6-10') +
  theme_minimal() +
  coord_sf() 


# save
save (grid_info,grid_base,grid_data, file = here ("processed_data", "grid_info.RData"))

# occurrence data
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret$occdf


# cell coordinates 
cell_coordinates <- coll_occ_taxa_perm_cret %>% 
  group_by (cell_ID) %>%
  summarise (centroid_lat = mean (cell_centroid_lat),
             centroid_lng = mean(cell_centroid_lng),
             sd_lat = sd (cell_centroid_lat),
             sd_lng = sd (cell_centroid_lng))

# A function to assign fossil occurrences to user-specified latitudinal bins.
bins_lat <- lat_bins(size = 20) # size = degrees
coll_occ_taxa_perm_cret <- bin_lat(coll_occ_taxa_perm_cret, 
                                   bins=bins_lat, 
                                   lat = "paleolat2",
                                   boundary = FALSE)


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



# save a table with PBDB taxonomy to complete
#PBDB_taxonomy <- tax_unique(occdf = coll_occ_taxa_perm_cret, 
#                            genus = "genus", 
#                            family = "family",
#                            order = "order", 
#                            class = "class", 
#                            name = "accepted_name",
#                            resolution="genus")
#
#  
## save
#write.csv (PBDB_taxonomy,
#           file = here ("processed_data", "taxonomy_PBDB.csv"))
#

# dim now : 3637

# remove others than sp.
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$unique_name != " sp."),]

3637 - nrow(coll_occ_taxa_perm_cret)

# bind the classification of genera made by hand
# load taxonomy
taxonomy <- read.xlsx (here ("processed_data", "table_taxonomy.xlsx"))

# match
coll_occ_taxa_perm_cret$clade <- taxonomy$clade3 [match(coll_occ_taxa_perm_cret$unique_name,
                                                        taxonomy$unique_name_PBDB_taxonomy)]
# remove synonyms
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[-which(coll_occ_taxa_perm_cret$clade %in%
                                                           c("REMOVER","Remover")),]

2763 - nrow(coll_occ_taxa_perm_cret)

# find the genus range size
# exploring range size (number of cells with detection)
range_area <- group_apply(occdf = coll_occ_taxa_perm_cret,
                                             group = c("interval_mid_ma"),
                                             fun = tax_range_space,
                                             name = "unique_name",
                                             lng = "paleolng2",
                                             lat = "paleolat2",
                                             method = "occ",
                                            spacing=50)


# aggregate per taxon(mean range area across time bins)
range_area_taxon <-group_apply(range_area, 
                         "taxon",
                         function(df) 
                          mean(df$n_cells))
colnames(range_area_taxon) <- c("range_area", "taxon")

# bind to the dataset
coll_occ_taxa_perm_cret$range_area <- range_area_taxon[match (coll_occ_taxa_perm_cret$unique_name,
                                                              range_area_taxon$taxon), "range_area"]


# Have a look at the dataset
# check 
coll_occ_taxa_perm_cret[grep ("Adelodelphys",coll_occ_taxa_perm_cret$unique_name),]


# Find the average geographic range size for each time interval
mean_range_area_bin <- group_apply(range_area, 
                                                  "interval_mid_ma",
                                                  function(df) 
                                                    mean(df$n_cells))

colnames(mean_range_area_bin) <- c("mean_area", "bin_midpoint")

# Create a plot of average range size through time
ggplot(mean_range_area_bin,
       aes(x = as.numeric(bin_midpoint),
           y = (mean_area))) +
  
  geom_point() +
  
  geom_smooth(method = "loess", span = 0.8, alpha = 0.2) +
  ylab("Cynodontia genus' range size\n(number of occupied cells, stage average)") +
  coord_geo(
    dat = list("periods", "eras"), 
    xlim = c(max(as.numeric(mean_range_area_bin$bin_midpoint))+10, 
             min(as.numeric (mean_range_area_bin$bin_midpoint))-10), 
    #ylim = c(-100000, max(space_coll_occ_taxa_perm_cret_mean$mean_area)+100000),
    pos = list("b", "b"), abbrv = list(TRUE, FALSE)
  ) +
  scale_x_reverse("Age (Ma)") +
  theme_classic() 

ggsave (file = here ("output","figures", "av_cells_stage.png"), height=20,width=11,units ="cm",dpi=600)


# plot of range size
range_area_taxon %>%
  ggplot() +
  geom_linerange(mapping=aes(x = reorder(taxon, range_area), 
                             ymin = 0, 
                             ymax = range_area), 
                 size = 1, alpha = 0.5,
                 position = position_dodge(width = 0.1)) +
  geom_point (aes(x =  reorder(taxon, range_area), 
                  y = range_area),col = "green3")+
  theme (axis.text.y =  element_text(size=1.5))+
  coord_flip()+
  xlab ("Taxon")

table(range_area_taxon$range_area ==1 )/sum(table(range_area_taxon$range_area ==1 ) ) 

ggsave (file = here ("output","figures", "genus_spatial_range.png"), height=20,width=11,units ="cm",dpi=600)


# -----------------------------------------------

# determine detection and intervals
# create index as value in the table
coll_occ_taxa_perm_cret$detection <- 1

# check preservation quality
# table(coll_occ_taxa_perm_cret$preservation_quality)
# filter(coll_occ_taxa_perm_cret, preservation_quality.x != "very poor")

# --------------------------------------------------------------------

# select interesting cols (to downsize the dataset)

coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[,  c("collection_no",
                                                        "unique_name", 
                                                        "genus",
                                                        "class",
                                                        "clade",
                                                        "formation",
                                                         "min_ma",
                                                        "interval_mid_ma",
                                                        "max_ma",
                                                        "cell_ID",
                                                        "cell_centroid_lat",
                                                        "lat_bin",
                                                        "lat_max",
                                                        "lat_mid",
                                                        "lat_min",
                                                        "lat.x",
                                                        "lng.x",
                                                        "paleolng2",
                                                        "paleolat2",
                                                        "range_area",
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

# paste period to the bins table
bins$Period <- coll_occ_taxa_perm_cret$Period[match (bins$interval_name,coll_occ_taxa_perm_cret$Interval)]


# plot table
table_taxon_period <- cast (unique_name ~ bin_assignment,
                            data = coll_occ_taxa_perm_cret,
                            value = "detection",
                            fun.aggregate = sum)
rownames(table_taxon_period)<- table_taxon_period$unique_name;
table_taxon_period<- table_taxon_period[,-1]
table_taxon_period[table_taxon_period>1] <- 1

# melt data to plot
melt_data <- melt (data.matrix(table_taxon_period))
melt_data$age <- bins$mid_ma [match (melt_data$X2, bins$bin) ]

# plot the stratigraphic range of each taxon
melt_data %>%
  filter (value >0) %>%
  group_by(X1) %>%
  summarise (min_int = min (age),
             max_int = max (age)) %>%
  mutate (diff = max_int - min_int) %>%
  arrange(-diff) %>%
  ggplot() +
  geom_linerange(mapping=aes(x = reorder(X1, diff), 
                             ymin = min_int, 
                             ymax = max_int), 
                 size = 1, alpha = 0.5,
                 position = position_dodge(width = 0.1)) +
  geom_point (aes(x = reorder(X1, diff), 
                  y = min_int),col = "red",size=2)+
  geom_point (aes(x = reorder(X1, diff), 
                  y = max_int),col = "orange")+
  theme (axis.text.y =  element_text(size=1.5))+
  coord_flip()+
  xlab ("Taxon")+
  ylab ("Time bin")+
  
  scale_y_reverse("Age (Ma)")


ggsave (file = here ("output","figures", "genus_range.png"), height=20,width=11,units ="cm",dpi=600)


# -----------------------------------------------
# produce animations to explore data


# paleomaps of formations - take a time ...
#maps_models <- getmap(ma = bins$mid_ma, 
#                      model = "PALEOMAP")
#maps_models <- maps_models [order(bins$mid_ma)] # ordering
#names(maps_models) <- rev(bins$interval_name)
#save (maps_models, file = here ("processed_data", "paleomaps.RData")) 

# load these maps
load(here ("processed_data", "paleomaps.RData"))


# get coordinates of each record
# coords to reconstruct
df_rec <- rbind (
  
          data.frame (Period = "Permian",
                      mid_ma = 253.021),
          
          data.frame (Period = "Triassic",
                      mid_ma = 204.900),
          
          data.frame (Period = "Jurassic",
                      mid_ma =  148.550),
          
          data.frame (Period = "Cretaceous",
                      mid_ma = 69.050)
)

# bind          
coll_occ_taxa_perm_cret$agr_rec <- df_rec [match (coll_occ_taxa_perm_cret$Period,df_rec$Period),"mid_ma"]


# rotate
coords_period <- lapply (rev(unique(coll_occ_taxa_perm_cret$Period)), function (i) {
  # rotate coords
    coords_records <- palaeorotate(
      coll_occ_taxa_perm_cret %>%
        filter (Period == i) ,
      lng = "lng.x",
      lat = "lat.x",
      age = "agr_rec",
      model = "PALEOMAP",
      method = "point",
      uncertainty = TRUE,
      round = 3
    )
    
    # sf object
    coords_records_sf <-st_as_sf(coords_records, coords = c("p_lng", "p_lat"),
                              crs = "+proj=longlat +datum=WGS84 +no_defs" )
    coords_records_sf

})

# set colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen")


# plot and save
pdf (here("output", "figures","maps_overlap.pdf"),width=15,height=10,onefile=T)

# plot

p_permian <- ggplot(data = st_transform( st_as_sf(maps_models[["Changhsingian"]]),
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +

  geom_sf(fill = "gray80", colour = 'gray20') +
  scale_y_continuous(breaks = seq(-90, 90, by = 10))+
  scale_x_continuous(breaks = seq(-180, 180, by = 50))+
  theme (panel.grid.major = element_line(colour = "gray", linewidth = 1),
         legend.position = "top",
         panel.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.5, linetype = "solid")) +
  # add records
  
  geom_sf(data = coords_period[[1]],
          #aes (colour = ifelse (coll_occ_taxa_perm_cret$paleolat2 >=0,1,2)),
          aes (colour = clade ),
          alpha=0.5,size=2) +
  
  scale_colour_manual(values = cols)+
  ggtitle ("Permian")

# triassic
p_triassic <- ggplot(data = st_transform( st_as_sf(maps_models[["Rhaetian"]]),
                                         "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  
  geom_sf(fill = "gray80", colour = 'gray20') +
  scale_y_continuous(breaks = seq(-90, 90, by = 10))+
  scale_x_continuous(breaks = seq(-180, 180, by = 50))+
  theme (panel.grid.major = element_line(colour = "gray", linewidth = 1),
         legend.position = "top",
         panel.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.5, linetype = "solid")) +
  
  # add records
  
  geom_sf(data = coords_period[[2]],
          #aes (colour = ifelse (coll_occ_taxa_perm_cret$paleolat2 >=0,1,2)),
          aes (colour = clade ),
          alpha=0.5,size=2) +
  
  scale_colour_manual(values = cols)+
  ggtitle ("Triassic")


# jurassic
p_jurassic <- ggplot(data = st_transform( st_as_sf(maps_models[["Tithonian"]]),
                                          "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  
  geom_sf(fill = "gray80", colour = 'gray20') +
  scale_y_continuous(breaks = seq(-90, 90, by = 10))+
  scale_x_continuous(breaks = seq(-180, 180, by = 50))+
  theme (panel.grid.major = element_line(colour = "gray", linewidth = 1),
         legend.position = "top",
         panel.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.5, linetype = "solid")) +
  
  # add records
  
  geom_sf(data = coords_period[[3]],
          #aes (colour = ifelse (coll_occ_taxa_perm_cret$paleolat2 >=0,1,2)),
          aes (colour = clade ),
          alpha=0.5,size=2) +
  
  scale_colour_manual(values = cols)+
  ggtitle ("Jurassic")


# jurassic
p_cretaceous <- ggplot(data = st_transform( st_as_sf(maps_models[["Maastrichtian"]]),
                                          "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  
  geom_sf(fill = "gray80", colour = 'gray20') +
  scale_y_continuous(breaks = seq(-90, 90, by = 10))+
  scale_x_continuous(breaks = seq(-180, 180, by = 50))+
  theme (panel.grid.major = element_line(colour = "gray", linewidth = 1),
         legend.position = "top",
         panel.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.5, linetype = "solid")) +
  
  # add records
  
  geom_sf(data = coords_period[[4]],
          #aes (colour = ifelse (coll_occ_taxa_perm_cret$paleolat2 >=0,1,2)),
          aes (colour = clade ),
          alpha=0.5,size=2) +
  
  scale_colour_manual(values = cols)+
  ggtitle ("Cretaceous")


grid.arrange(p_permian,
             p_triassic,
             p_jurassic,
             p_cretaceous,
             nrow=2)



dev.off()

# data for the legend
range(colSums(table (coll_occ_taxa_perm_cret$paleolng2,coll_occ_taxa_perm_cret$cell_ID)>0))

# animation per stage
# plot points
dir.create(here ("processed_data","animation"))
lapply (bins$interval_name[7:length(bins$interval_name)], function (i) {
  
  png (here ("processed_data","animation", paste ("maps_", (i), ".png",sep ="")),
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
  
  with (coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$Interval %in% i),], 
        
        points (paleolng2,
                paleolat2,
                pch=19,
                col = rgb(0.1,0.5,0,alpha=0.1))
  )
  
  
  dev.off()
})

#anime
list_img <- list.files(path = here ("processed_data","animation"), full.names = T, pattern = "maps_")
# ordering
list_img <- list_img[match (bins$interval_name[7:length(bins$interval_name)],
                            gsub (".png", "",gsub ("D:/Pos_Doc_Paleonto_Macroecology/modeling/paleo_macroecology/processed_data/animation/maps_","",list_img))
)]
##https://cran.r-project.org/web/packages/magick/vignettes/intro.html
a_image<-image_read(list_img)
animation <-  image_animate(a_image, fps = 1)
image_write(animation, here ("output","animation_data.gif"))


# ---------------------------------------------------------------------
# check


# adjust period order
coll_occ_taxa_perm_cret$Period <- factor(coll_occ_taxa_perm_cret$Period,
                                         levels = c("Permian",
                                                    "Triassic",
                                                    "Jurassic", 
                                                    "Cretaceous"))


png (here ("output", "figures", "descriptive_hists.png"),height=15,width=20,units ="cm",res=300)

par(mfrow=c(2,2),mar=c(4,4,4,4)) # (15,5,15,4)

# binned intervals per cell
hist (
  
  
  (rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                  coll_occ_taxa_perm_cret$bin
  )>0)
  ),
  main = "Stages per region (cell)",
  xlab = "Count"
)


# formations per cell
hist (
  
  
  (rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                  coll_occ_taxa_perm_cret$formation
  )>0)
  ),
  main= "Formations per region (cell)",
  xlab = "Count"
)


# formations per interval
hist (
  
  
  (colSums(
    table (coll_occ_taxa_perm_cret$formation,
           coll_occ_taxa_perm_cret$bin
    )>0)
  ),
  main = "Formations per stage",
  xlab = "Count"
)

# formations per lat bin
hist (
  
  
  (colSums(
    table (coll_occ_taxa_perm_cret$formation,
           coll_occ_taxa_perm_cret$lat_bin
    )>0)
  ),
  main = "Formations per latitudinal bin",
  xlab = "Count"
)



dev.off()


# latitude of records
# histogram
p1<-ggplot(data = coll_occ_taxa_perm_cret, 
       aes(x=paleolat2))+
  geom_histogram()+
  facet_grid(~Period)+
  theme_bw()+
  xlab ("Paleolatitude") +
  ylab ("Number of records")

# longitude of records
# histogram
p2<-ggplot(data = coll_occ_taxa_perm_cret, 
       aes(x=paleolng2))+
  geom_histogram()+
  facet_grid(~Period)+
  theme_bw()+
  xlab ("Paleolongitude")+
  ylab ("Number of records")

png (here ("output","figures", "hist_lat_period.png"),height=15,width=20,units ="cm",res=300)
grid.arrange(p1,p2,nrow=2,
             ncol=1)
dev.off()


# --------------------------------------
# final number of formatiosn adn collections
length(unique(coll_occ_taxa_perm_cret$formation))
length(unique(coll_occ_taxa_perm_cret$collection_no))


# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)

groups <- unique(coll_occ_taxa_perm_cret$clade)[order (unique(coll_occ_taxa_perm_cret$clade))]
cynodontia_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin_assignment)[order(unique(coll_occ_taxa_perm_cret$bin_assignment))]
cells <- unique(coll_occ_taxa_perm_cret$lat_bin)[order( unique(coll_occ_taxa_perm_cret$lat_bin))]
formations <- unique(coll_occ_taxa_perm_cret$formation)[order(unique(coll_occ_taxa_perm_cret$formation))]

# save
write.csv (data.frame(cynodontia_data),
           file = here ("processed_data", "table_taxonomy.csv"))


# number of detections per genera

tapply(coll_occ_taxa_perm_cret$detection,
       coll_occ_taxa_perm_cret$unique_name,
       sum) [order(tapply(coll_occ_taxa_perm_cret$detection,
                          coll_occ_taxa_perm_cret$unique_name,
                          sum)
       )]

# the basic table summarizing all detections across intervals and formations

#table_data_basis <- lapply (cells, function (i) {
#  
#  
#  tab_basis <- cast (formula = unique_name ~ bin_assignment,
#                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$lat_bin %in% i),],
#                     value = "detection", 
#                     fun.aggregate = sum,
#                     drop=F,
#                     fill=0)
#  
#  
#  # order names
#  rownames(tab_basis) <- tab_basis$unique_name
#
#    
#  # input formations/subcells not in this data
#  tab_basis <- cbind (tab_basis, 
#                      matrix (0, 
#                              nrow = nrow (tab_basis),
#                              ncol= length(time_bins [which(time_bins  %in% colnames (tab_basis) == F)]), # subcells or formations
#                              
#                              dimnames = list (rownames(tab_basis),
#                                               time_bins [which(time_bins  %in% colnames (tab_basis) == F)]
#                                               
#                              )
#                      )
#  )
#  tab_basis <- tab_basis[,-which(colnames(tab_basis) == "unique_name")]
#  
#  # input cells not in this data
#  tab_basis <- rbind (tab_basis, 
#                      matrix (0, 
#                              nrow = sum(cynodontia_data %in% rownames (tab_basis) == F),
#                              ncol= ncol (tab_basis),
#                              
#                              dimnames = list (cynodontia_data [which(cynodontia_data  %in% rownames (tab_basis) == F)],
#                                               colnames(tab_basis)
#                                               
#                              )
#                      )
#  )
#  
#  # order columns and rows
#  tab_basis <- tab_basis[order (rownames(tab_basis)),
#                         order(as.numeric(colnames(tab_basis)))]
#  
#  ; # return
#  tab_basis
#  
#})
#
## checknames
#colnames(table_data_basis[[1]]) == colnames(table_data_basis[[9]])
#rownames(table_data_basis[[1]]) == rownames(table_data_basis[[9]])
#
#rownames(table_data_basis[[1]]) == genus
#
## data to array
#array_genus_bin_site <- array (unlist (table_data_basis),
#                               dim = c(length(cynodontia_data),
#                                       length(time_bins),
#                                       length(table_data_basis)))
#
#array_genus_bin_site[array_genus_bin_site>0]<-1


# global scale array
# tapply option
array_genus_bin_site <- tapply (coll_occ_taxa_perm_cret$detection,
                                list(coll_occ_taxa_perm_cret$unique_name,
                                     coll_occ_taxa_perm_cret$bin_assignment,
                                      coll_occ_taxa_perm_cret$lat_bin,
                                     coll_occ_taxa_perm_cret$formation),
                        sum,na.rm=T)

# obtain the number of formations with detection per taxon
array_genus_bin_site[array_genus_bin_site>0]<-1
#range(array_genus_bin_site,na.rm=T)

# aggregate formations
array_genus_bin_site <- (apply (array_genus_bin_site,c(1,2,3), sum,na.rm=T))

# ---------------------------------------------

# number of formations per cell and interval
# the basic table summarizing all detections across intervals and formations

formations_per_site_interval<- (tapply (coll_occ_taxa_perm_cret$detection,
                                        list(coll_occ_taxa_perm_cret$bin_assignment,
                                             coll_occ_taxa_perm_cret$lat_bin,
                                             coll_occ_taxa_perm_cret$formation),
                                            sum,na.rm=T))
# aggregate formations
formations_per_site_interval<-(apply (formations_per_site_interval, c(1,2),sum,na.rm=T))

# check whether the number of formations with one genus is above the number of formations   
checkit<-lapply (seq(1,nrow(array_genus_bin_site)), function (i)
  table(array_genus_bin_site[i,,] > formations_per_site_interval)
)
table(unlist(checkit)) # good of all is false

# ------------------------------------------------------------------------------------

# global scale array
array_genus_bin <- tapply (coll_occ_taxa_perm_cret$detection,
        list(coll_occ_taxa_perm_cret$unique_name,
             coll_occ_taxa_perm_cret$bin_assignment,
             coll_occ_taxa_perm_cret$formation),
        sum,na.rm=T)

# obtain the number of formations with detection per taxon
array_genus_bin[array_genus_bin>0]<-1
#range(array_genus_bin_site,na.rm=T)

# aggregate formations
array_genus_bin <- (apply (array_genus_bin,c(1,2), sum,na.rm=T))

# clades
clades <- taxonomy$clade3 [match (rownames (array_genus_bin), taxonomy$genus)]
genus <- rownames (array_genus_bin)


# get the total number of detections per formation and interval

formations_per_interval <-  (tapply (coll_occ_taxa_perm_cret$detection,
        list(coll_occ_taxa_perm_cret$bin_assignment,
             coll_occ_taxa_perm_cret$formation),
        sum,na.rm=T))

# obtain the number of formations with detection per taxon
formations_per_interval[formations_per_interval>0]<-1
# aggregate
formations_per_interval <- apply (formations_per_interval, c(1),sum,na.rm=T)

# get the number of collections per interval
coll_per_interval <- tapply (coll_occ_taxa_perm_cret$detection,
                                   list(coll_occ_taxa_perm_cret$collection_no,
                                        coll_occ_taxa_perm_cret$bin_assignment),
                                   sum,na.rm=T)
coll_per_interval[coll_per_interval>0]<-1
coll_per_interval[is.na(coll_per_interval)]<-0
coll_per_interval <- data.frame (bin_assignment=colnames(coll_per_interval),
                                 coll_per_interval = apply (coll_per_interval,2,sum))
# save fig
png(here ("output", "figures","pull_of_the_recent.png"),height=20,width=22,units ="cm",res=300)
par(mfrow=c(2,1),mar=c(5,5,5,5))

# plot observed SR per site
plot(NA,ylim=c(range(apply (array_genus_bin_site>0,c(2,3),sum))[1],
               150),
     xlim=c(1,33),
     ylab = "Number of genus",
     xlab = "Intervals (0: Capitanian; 33: Maastrichtian)")

lapply (seq(1,dim(array_genus_bin_site)[3]), function (i)
  
  lines (apply (array_genus_bin_site>0,c(2,3),sum)[,i],
         col="gray")
  
)
# observed n of genus
lines (colSums(array_genus_bin>0), col="red", lwd=2)

# pull of the recent?
plot(NA,ylim=c(range(coll_per_interval$coll_per_interval)[1],
               range(coll_per_interval$coll_per_interval)[2]),
     xlim=c(1,33),
     ylab = "Number of formations (black)\nand collections (green)",
     xlab = "Intervals (0: Capitanian; 33: Maastrichtian)")
lines(formations_per_interval,type= "b",lwd=2)
lines(coll_per_interval$coll_per_interval,type= "b",lwd=2,col="green")

dev.off()

# correlation
cor(coll_per_interval$coll_per_interval,
    formations_per_interval)
plot(log(coll_per_interval$coll_per_interval),
    log(formations_per_interval))

# --------------------------------------------------

# try a vectorized version of the region-level dataset
vectorized_occ <-melt(array_genus_bin_site)
colnames (vectorized_occ) <- c("genus", "stage","site", "detection")
vectorized_occ$int <- paste (vectorized_occ$stage,vectorized_occ$site,sep="_")

# paste genus
vectorized_occ$genusID <- genus
vectorized_occ$cladeID <- clades
vectorized_occ$clade <- (as.numeric(as.factor(vectorized_occ$cladeID)))

# vectorized formations
vectorized_formations <- (melt (unname(formations_per_site_interval), as.is=T))
colnames (vectorized_formations) <- c("site","stage", "formations")
vectorized_formations$int <- paste (vectorized_formations$stage+6,vectorized_formations$site,sep="_")

# match
vectorized_occ$formations<- vectorized_formations [match (vectorized_occ$int,
       vectorized_formations$int
       ),"formations"]
# yes!
vectorized_occ[which(vectorized_occ$detection > vectorized_occ$formations),]

# paste period
vectorized_occ$periodID <- bins$Period [match (vectorized_occ$stage+6,bins$bin)]
vectorized_occ$period <- (as.numeric(factor(vectorized_occ$periodID,
                                            levels = c("Permian","Triassic","Jurassic","Cretaceous"))))


(tapply (vectorized_occ$detection,
        vectorized_occ$genusID,
        sum))[order((tapply (vectorized_occ$detection,
                             vectorized_occ$genusID,
                             sum))
        )]

# range size
vectorized_occ$range_size <- range_area_taxon [match (vectorized_occ$genusID,
                                                      range_area_taxon$taxon
),"range_area"]


# --------------------------------------------------
# vectorized version of the global-level dataset


# try a vectorized version of the region-level dataset
vectorized_occ_global <-melt(array_genus_bin)
colnames (vectorized_occ_global) <- c("genusID", "stage" , "detection")
# clade and genus code
vectorized_occ_global$clade <- (as.numeric(as.factor(vectorized_occ_global$cladeID)))
vectorized_occ_global$genus <- (as.numeric(as.factor(vectorized_occ_global$genusID)))

#  formations
vectorized_occ_global$formations<- formations_per_interval [match (vectorized_occ_global$stage,
                                                                   formations_per_interval$bin_assignment),
                                                            "formations_per_interval"]

# yes!
vectorized_occ_global[which(vectorized_occ_global$detection > vectorized_occ_global$formations),]

# paste period
vectorized_occ_global$periodID <- bins$Period [match (vectorized_occ_global$stage,bins$bin)]
vectorized_occ_global$period <- (as.numeric(factor(vectorized_occ_global$periodID,
                                                                levels = c("Permian","Triassic","Jurassic","Cretaceous"))))



(tapply (vectorized_occ_global$detection,
         vectorized_occ_global$genusID,
         sum))[order((tapply (vectorized_occ_global$detection,
                              vectorized_occ_global$genusID,
                              sum))
         )]

# range size
vectorized_occ_global$range_size <- range_area_taxon [match (vectorized_occ_global$genusID,
                                                            range_area_taxon$taxon
),"range_area"]



# save 
save (coll_occ_taxa_perm_cret,
      array_genus_bin, 
      clades,
      genus,
      array_genus_bin_site,
      formations_per_interval,
      coll_per_interval,
      formations_per_site_interval,
      cynodontia_data,
      time_bins ,
      cells ,
      bins_lat,
      formations,
      bins,
      range_area_taxon,
      # global
      vectorized_occ_global,
      # regional
      vectorized_occ,
      file= here("processed_data", "CMR_data.RData")
      )

rm(list=ls())




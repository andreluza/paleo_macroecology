

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

# =====================================================
# load data
# collections
#pbdb_data_collections <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_collections.csv"),
#                                   skip=18)
url_col <- "https://paleobiodb.org/data1.2/colls/list.csv?base_name=cynodontia&max_ma=259.5&min_ma=66"
pbdb_data_collections <- read.csv (url_col, sep=",")

# strata
#pbdb_data_strata <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_strata.csv"),
#                              skip=17)
url_strata <- "https://paleobiodb.org/data1.2/occs/strata.csv?base_name=cynodontia&max_ma=259.5&min_ma=66&pgm=gplates,scotese,seton"
pbdb_data_strata <-  read.csv (url_strata, sep=",")

#taxa
#pbdb_data_taxa <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_taxa.csv"),
#                            skip=15)
url_taxa <- "https://paleobiodb.org/data1.2/occs/taxa.csv?base_name=cynodontia&max_ma=259.5&min_ma=66"
pbdb_data_taxa <- read.csv(url_taxa,sep=",")

# occ
#pbdb_data_occ <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_occ.csv"),
#                           skip=18)
url_occ <- "https://paleobiodb.org/data1.2/occs/list.csv?base_name=cynodontia&max_ma=259.5&min_ma=66&pgm=gplates,scotese,seton&show=genus,paleoloc,class"
pbdb_data_occ <- read.csv(url_occ,sep=",")


# match collection and occurrence (lat long of collection)
pbdb_data_occ <- cbind (pbdb_data_occ,
                        pbdb_data_collections [match (pbdb_data_occ$collection_no,
                                                      pbdb_data_collections$collection_no),
                                               c("lat","lng")])

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


save (grid_info, file = here ("processed_data", "grid_info.RData"))

# occurrence data
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret$occdf


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


# find the genus range size
# exploring range size (number of cells with detection)
range_area <- group_apply(occdf = coll_occ_taxa_perm_cret,
                                             group = c("interval_mid_ma"),
                                             fun = tax_range_space,
                                             name = "unique_name",
                                             lng = "paleolng2",
                                             lat = "paleolat2",
                                             method = "occ",
                                            spacing=2000)

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
  ylab("Cynodontia genus' range size (number of occupied cells)") +
  coord_geo(
    dat = list("periods", "eras"), 
    xlim = c(max(as.numeric(mean_range_area_bin$bin_midpoint))+10, 
             min(as.numeric (mean_range_area_bin$bin_midpoint))-10), 
    #ylim = c(-100000, max(space_coll_occ_taxa_perm_cret_mean$mean_area)+100000),
    pos = list("b", "b"), abbrv = list(TRUE, FALSE)
  ) +
  scale_x_reverse("Age (Ma)") +
  theme_classic() 


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

ggsave (file = here ("output", "genus_spatial_range.pdf"))

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
                                                        "formation",
                                                         "min_ma",
                                                        "interval_mid_ma",
                                                        "max_ma",
                                                        "cell_ID",
                                                        "cell_centroid_lat",
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


ggsave (file = here ("output", "genus_range.pdf"))

# -----------------------------------------------
# produce animations to explore data

# paleomaps of formations - take a time ...
#maps_models <- getmap(ma = bins$mid_ma, 
#                      model = "GOLONKA")
#maps_models <- maps_models [order(bins$mid_ma)] # ordering
#save (maps_models, file = here ("processed_data", "paleomaps.RData")) 

load(here ("processed_data", "paleomaps.RData"))
names(maps_models) <- bins$interval_name

pdf (here("output", "maps_overlap.pdf"),width=15,height=10)
# test
plot(maps_models[[1]],
     border=rgb(0.1,0.1,0.1,alpha=0.03), 
     col=rgb(0.1,0.1,0.1,alpha=0.01)
)

# plot paleomaps
lapply (maps_models, plot,
        border=rgb(0.1,0.1,0.1,alpha=0.05), 
        col=rgb(0.1,0.1,0.1,alpha=0.03),add=T)
# 
plot(maps_models[[1]],
     border=rgb(0.1,0.5,0.7,alpha=1), 
     col=rgb(0.1,0.5,0.7,alpha=0.1),add=T
)
plot(maps_models[[39]],
     border=rgb(0.1,0.1,0.1,alpha=1), 
     col=rgb(0.1,0.1,0.1,alpha=0.1),add=T
)

# plot spatial data
plot(grid_info$grid_base,border="gray80",add=T)
plot(grid_info$grid,add=T, 
     border="gray80",
     col = rgb(0.1,0.1,0.5,alpha=0.1))
#plot(grid_info$sub_grid,add=T, col = rgb(0.5,0.1,0.1,alpha=0.2))

# map
points(
  coll_occ_taxa_perm_cret$paleolng2,
  coll_occ_taxa_perm_cret$paleolat2,
  col = rgb (0.8,0.1,0.1,alpha=0.1),pch=19)

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
                  coll_occ_taxa_perm_cret$formation
  )>0)
  ),
  main= "formations per cell"
)


# formations per interval
hist (
  
  
  (colSums(
    table (coll_occ_taxa_perm_cret$formation,
           coll_occ_taxa_perm_cret$bin
    )>0)
  ),
  main = "formations per interval"
)


# pull of the recent?
plot((colSums(
  table (coll_occ_taxa_perm_cret$formation,
         coll_occ_taxa_perm_cret$bin_assignment
  )>0)
),type="l",xlab= "intervals (older (0 - Capitanian) to the most recent (33 - Maastrichian))", ylab= "formations")



# --------------------------------------




# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)



cynodontia_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin_assignment)[order(unique(coll_occ_taxa_perm_cret$bin_assignment))]
cells <- unique(coll_occ_taxa_perm_cret$cell_ID)[order( unique(coll_occ_taxa_perm_cret$cell_ID))]
formations <- unique(coll_occ_taxa_perm_cret$formation)[order(unique(coll_occ_taxa_perm_cret$formation))]

# the basic table summarizing all detections across intervals and formations

table_data_basis <- lapply (cells, function (i) {
  
  
  tab_basis <- cast (formula = unique_name ~ bin_assignment,
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
  
  
  tab_basis <- cast (formula = bin_assignment ~ formation,
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$cell_ID %in% i),],
                     value = "detection", 
                     fun.aggregate = max,
                     drop=F,
                     fill=0)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$bin_assignment
  
  
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
  tab_basis <- tab_basis[,-which(colnames(tab_basis) == "bin_assignment")]
  
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
      bins,
      range_area_taxon,
      file= here("processed_data", "CMR_data.RData")
      )
rm(list=ls())

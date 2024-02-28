# ======================================================

# organizing covariates


# organizing paleodata
# the metadata for all this
# https://taphonomy.doit.wisc.edu/data1.2/colls/single_doc.html
# geoplates
# https://www.gplates.org/docs/user-manual/Interacting_Features/#30-assign-plate-ids
# https://www.gplates.org/docs/user-manual/Reconstructions/#plate-ids

# ======================================================

# load packages
source("R/packages.R")

# =====================================================
# load data
# load spatial covariates 
load (here ("processed_data", "brick_rasters_elevation.RData"))
load (here ("processed_data", "brick_rasters_precipitation.RData"))
load (here ("processed_data", "brick_rasters_temp.RData"))

# load fossil data
load(here ("processed_data","CMR_data.RData"))

# load intervals
# find the interval of each raster
intervals <- openxlsx::read.xlsx (here ("data", "periods_intervals_ages.xlsx"))
# reverse the order to the past be the first period
#rev_order <- rev(intervals$Interval)
#intervals <- intervals [match (intervals$Interval, rev_order),]

# bind missing stages
# apply to each raster
raster_imputation <- lapply (list(brick_rasters,
                                 brick_rasters_paleoprec,
                                 brick_rasters_paleotemp), function (layer) {
                                   
      stages_to_imput <- bins$interval_name [-which(bins$interval_name %in% names(layer))]
      # imputations 
      imputations <- lapply (stages_to_imput, function (stage) {
        each_stage_to_imput <- bins$interval_name[which(bins$interval_name %in% stage)-1]
        if (each_stage_to_imput %in% names (layer)) {
          t1<-layer[[each_stage_to_imput]]
          names (t1) <- stage
        
          } else {
            
            each_stage_to_imput <- bins$interval_name[which(bins$interval_name %in% stage)-2]
            t1<-layer[[each_stage_to_imput]]
            names (t1) <- stage
            
          
          }
        
        t1
      })
      
      # to brick
      imputed_brick <- brick (stack(layer, brick(imputations)))
      #names (imputed_brick [[match(bins$interval_name,names (imputed_brick))]]) == bins$interval_name
      imputed_brick <- (imputed_brick [[match(bins$interval_name,names (imputed_brick))]])
      imputed_brick

})

names (raster_imputation[[1]]) == names (raster_imputation[[2]])
names (raster_imputation[[1]]) == names (raster_imputation[[3]])
names (raster_imputation[[2]]) == names (raster_imputation[[3]])

# precipitation and temperature in land
elevation_raster_continuous <- raster_imputation[[1]]
elevation_raster <- raster_imputation[[1]]
precipitation_raster <- raster_imputation[[2]]
temperature_raster <- raster_imputation[[3]]

# create a mask of elevation
values (elevation_raster) <- ifelse (values (elevation_raster) >=0, 1,NA)
# weight temperature and precipitation
precipitation_raster <- precipitation_raster*elevation_raster
temperature_raster <- temperature_raster*elevation_raster
names(temperature_raster) <- names(raster_imputation[[3]])
names(precipitation_raster) <- names(raster_imputation[[2]])

# select stages we're studying
precipitation_raster <- precipitation_raster[[which(names(precipitation_raster) %in% 
                                                      bins$interval_name [7:length(bins$interval_name)])]]
temperature_raster <- temperature_raster[[which(names(temperature_raster) %in% 
                                                      bins$interval_name [7:length(bins$interval_name)])]]
elevation_raster <- elevation_raster[[which(names(elevation_raster) %in% 
                                             bins$interval_name [7:length(bins$interval_name)])]]
elevation_raster_continuous <- elevation_raster_continuous[[which(names(elevation_raster_continuous) %in% 
                                              bins$interval_name [7:length(bins$interval_name)])]]


# make a dataframe
table(is.na(values(temperature_raster[[1]])))
table(is.na(values(precipitation_raster[[1]])))

# area 
plot(apply (values(precipitation_raster),2, function (x) sum(x>0,na.rm=T))*prod(res(temperature_raster)),type="b")
plot(apply (values(precipitation_raster),2, function (x) sum(x>0,na.rm=T))*prod(res(temperature_raster)),type="b")

cbind(
  apply (values(precipitation_raster),2, function (x) sum(x>0,na.rm=T)) , 
  apply (values(temperature_raster),2, function (x) sum(x>0,na.rm=T)))

cor(  apply (values(precipitation_raster),2, function (x) sum(x>0,na.rm=T)) , 
      apply (values(temperature_raster),2, function (x) sum(x>0,na.rm=T)))


cor(  apply (values(precipitation_raster),2, function (x) sum(x>0,na.rm=T)) , 
      apply (values(elevation_raster),2, function (x) sum(x>0,na.rm=T)) )

cor(  apply (values(temperature_raster),2, function (x) sum(x>0,na.rm=T)) , 
      apply (values(elevation_raster),2, function (x) sum(x>0,na.rm=T)) )

lines(apply (values(elevation_raster),2, function (x) sum(x>0,na.rm=T)) ,type="b")



df_var <- lapply (seq(1,dim(temperature_raster)[3]), function (i){
  
  df_var <- data.frame ("temperature" = mean(values(temperature_raster[[i]]),na.rm=T),
                        "precipitation" = mean(values(precipitation_raster[[i]]),na.rm=T),
                        "temperature_sd" = sd(values(temperature_raster[[i]]),na.rm=T),
                        "precipitation_sd" = sd(values(precipitation_raster[[i]]),na.rm=T),
                        "stage" = names (precipitation_raster)[i],
                        "area" = sum (values(elevation_raster[[i]])>0,na.rm=T)*prod(res(elevation_raster[[i]])),
                        "coastalLine" = table(values(clamp(elevation_raster_continuous[[i]], 
                                                           
                                                           lower=-700, 
                                                           
                                                           upper=0, 
                                                           
                                                           useValues=F))
                                              <0)[2]
                        )
  
  }
  
)
time_covariates <- do.call(rbind,df_var)

# correlations
cor.test(time_covariates$temperature,time_covariates$precipitation)
cor.test(time_covariates$temperature,time_covariates$area)
cor.test(time_covariates$precipitation,time_covariates$area)
cor.test(time_covariates$coastalLine,time_covariates$area)


# -----------------------------------

#  observation covariates

# i= 1594 
# i = 1000
# i = 163


# extract altitude per occurrence point
paleo_env <- lapply (list(brick_rasters_paleoprec,
                          brick_rasters_paleotemp), function (layer)
      
      lapply (seq (1,nrow (coll_occ_taxa_perm_cret)), function (i) {
  
      # extract
      extracted_data <- raster::extract(layer,
                                 coll_occ_taxa_perm_cret[i,c("paleolng2","paleolat2")],
                                fun = "mean",
                                method = "simple")
      
      # select
      selected_bin_n <- coll_occ_taxa_perm_cret[i,"bin_assignment"]
      selected_bin <- bins[which(bins$bin %in% selected_bin_n),"interval_name"]
      output <- (extracted_data[,which(colnames(extracted_data) %in% selected_bin)])
      
      
      # if no output is reported
      if (length(output)==0){
        
        # gather the paleoaltitude from the previous period in the same point
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
which(unlist(lapply (paleo_env[[2]], function (i) length(i) == 0))==T)
which(unlist(lapply (paleo_env[[2]], function (i) is.na(i) == T))==T)


# melt and bind
precipitation_record <- do.call(rbind,paleo_env[[1]])[,1]
temperature_record <- do.call(rbind,paleo_env[[2]])[,1]

# observation covariates 
observation_covariates <- data.frame ("precipitation" = precipitation_record,
                                      "temperature" = temperature_record,
                                      "latitude" = coll_occ_taxa_perm_cret$paleolat2,
                                      "unique_name" = coll_occ_taxa_perm_cret$unique_name)

# aggregate
observation_covariates <- observation_covariates %>%
  group_by (unique_name) %>%
  summarise(precipitation = mean(precipitation),
            temperature = mean(temperature),
            latitude = mean(latitude))


# match the observation covariates with the vectorized dataset
vectorized_occ <-cbind (vectorized_occ,
       observation_covariates [match (vectorized_occ$genusID,
                               observation_covariates$unique_name),]
)
table(vectorized_occ$genusID == vectorized_occ$unique_name)

# time for preservation analysis
vectorized_occ$time <- bins$mid_ma[match (vectorized_occ$stage+6,bins$bin)]


# match the observation covariates with the vectorized dataset
vectorized_occ_global <-cbind (vectorized_occ_global,
                        observation_covariates [match (vectorized_occ_global$genusID,
                                                       observation_covariates$unique_name),]
)
table(vectorized_occ_global$genusID == vectorized_occ_global$unique_name)

# time for preservation analysis
vectorized_occ_global$time <- bins$mid_ma[match (vectorized_occ_global$stage,
                                                 bins$bin)]


# choose lat band
bins_lat <- data.frame (bin = 1,
                        max = 20,
                        mid = 0,
                        min = -20)

#  variables at site level (latitudinal bin)
#site_covs <-  lapply (seq(1,nrow(bins_lat[4:6,])), function (k)
site_covs <-  lapply (nrow(bins_lat), function (k)
  
                                do.call(rbind,
                                        lapply (seq(1,dim(precipitation_raster)[3]), function (i) {
                              
                                    # bounding box
                                    e <- extent(-180, 180,bins_lat[k,4],bins_lat[k,2])
                                    a.temp <- crop(temperature_raster, e)
                                    a.prec <- crop(precipitation_raster, e)
                                    a.elev <- crop(elevation_raster, e)
                                    a.elev.cont <- crop(brick_rasters, e)
                                    
                                    df_var <- data.frame ("temperature" = mean(values(a.temp[[i]]),na.rm=T),
                                                          "precipitation" = mean(values(a.prec[[i]]),na.rm=T),
                                                          "temperature_sd" = sd(values(a.temp[[i]]),na.rm=T),
                                                          "precipitation_sd" = sd(values(a.prec[[i]]),na.rm=T),
                                                          "area" =  sum (values(a.elev[[i]])>0,na.rm=T)*prod(res(a.elev[[i]])),
                                                          "coastalLine" = table(values(clamp(a.elev.cont[[i]], 
                                                                                             
                                                                                             lower=-700, 
                                                                                             
                                                                                             upper=0, 
                                                                                             
                                                                                             useValues=F))
                                                                                <0)[2],
                                                          "stage" = names (precipitation_raster)[i],
                                                          "lat_min" = bins_lat[k,4],
                                                          "lat_mid" = bins_lat[k,3],
                                                          "lat_max" = bins_lat[k,2])
                                    
                                    
                                   
                                    
                                    }
                                    )
                    )
                    
                    )

# extract covariates
region_temperature <- sapply (site_covs, "[[", "temperature")
region_temperature_sd <- sapply (site_covs, "[[", "temperature_sd")
region_precipitation <- sapply (site_covs, "[[", "precipitation")
region_precipitation_sd <- sapply (site_covs, "[[", "precipitation_sd")
region_latitude <- sapply (site_covs, "[[", "lat_mid")
region_area <- sapply (site_covs, "[[", "area")
region_coastalLine <- sapply (site_covs, "[[", "coastalLine")


# save  the updated dataset
save (coll_occ_taxa_perm_cret,
      observation_covariates,
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
      vectorized_occ,
      vectorized_occ_global,
      file= here("processed_data", "CMR_data_observation_cov.RData")
)


# save to use in models
save (time_covariates,
      region_precipitation_sd,
      region_precipitation,
      region_temperature,
      region_temperature_sd,
      region_latitude,
      region_area,
      region_coastalLine,
      file = here ("processed_data","site_covs.RData"))

# plots
par(mfrow=c(2,2))
plot(region_precipitation, type = "b", ylab = "Precipitation (mm/day)",xlab ="Stage\n(1=Capitanian, 33= Maastrichtian)",pch=19)
plot(region_temperature, type = "b", ylab = "Temperature (ºC)",xlab ="Stage",pch=19)
plot(region_area, type = "b", ylab = "Area (1x1º cells above sea level)",xlab ="Stage",pch=19)
plot(region_coastalLine, type = "b", ylab = "Area (1x1º cells from -700 to 0 m)",xlab ="Stage",pch=19)



rm(list=ls())


# -------------------------------------------------------






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
                      resolution="genus",
                      append=T)

# changes
length(unique(coll_occ_taxa_perm_cret$accepted_name))
length(unique(coll_occ_taxa_perm_cret$unique_name))
unique(coll_occ_taxa_perm_cret$unique_name)[order(unique(coll_occ_taxa_perm_cret$unique_name))]

# some mistakes (only " sp.")
coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$unique_name != " sp."),]

# remove NAs
# coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret [is.na(coll_occ_taxa_perm_cret$unique_name)==F,]

# Get the names of unique genera per collection
unique_genera <- unique(coll_occ_taxa_perm_cret[, c("unique_name", "formation")]) 

# Calculate the number of unique genera per collection
coll_taxa <- group_apply(unique_genera, group = "formation", fun = nrow)

# Rename column names:
colnames(coll_taxa) <- c("n_taxa", "formation.x")

# Take the columns pertaining to collections and their ages in millions of years:
coll_info <- coll_occ_taxa_perm_cret[, c("formation",
                               "max_ma", 
                               "interval_mid_ma", 
                               "min_ma")]

# Remove duplicated collections based on the collection number (column 1)
coll_info <- coll_info[!duplicated(coll_info[1]),]

# Combine this dataframe with the dataframe from above
alpha_data <- merge(coll_info, coll_taxa)#, by = "formation")
# Take a look:
head(alpha_data)

# uses the cynodontia diversity data from above
# help here https://cran.r-project.org/web/packages/deeptime/vignettes/coord_geo.html
ggplot(alpha_data,aes(x = (interval_mid_ma),
                      y = n_taxa)) +
  
  geom_point() +
  geom_line()+
  geom_smooth(method = "loess", span = 0.8, alpha = 0.2) +
  ylab("# Cynodontia genera") +
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
  ylab("Cynodontia genus' range size (log)") +
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
  summarise(n = n_distinct(unique_name),
            age = mean(interval_mid_ma)) %>%
  
  ggplot(aes(x = age, y = n)) +
        geom_line() +
        geom_point()+
        scale_x_reverse("Age (Ma)") +
        ylab("Cynodontia genera") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c(260, 70), 
    ylim = c(0, 50),
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
  

  (rowSums(table (coll_occ_taxa_perm_cret$cell_ID,
                       coll_occ_taxa_perm_cret$bin_assignment
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
           coll_occ_taxa_perm_cret$bin_assignment
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

# remove poor quality
#coll_occ_taxa_perm_cret<-coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$preservation_quality.x != "poor"),]


# change the order of bin assignmente (more recent first)
coll_occ_taxa_perm_cret$bin_assignment<-40-(coll_occ_taxa_perm_cret$bin_assignment)

# --------------------------------------------------------------------

# select interesting cols (to downsize the dataset)

coll_occ_taxa_perm_cret <- coll_occ_taxa_perm_cret[,  c("collection_no",
                                                        "unique_name", 
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
                                                        "cell_centroid_lat",
                                                        "paleolng2.x",
                                                        "paleolat2.x",
                                                        #"cell_ID_sub",
                                                        #"cell_centroid_lng_sub",
                                                        #"cell_centroid_lat_sub",
                                                        "detection"
)]


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

#
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
                  y = min_int),col = "orange",size=2)+
  geom_point (aes(x = reorder(X1, diff), 
                  y = max_int),col = "green3")+
  theme (axis.text.y =  element_text(size=1.5))+
  coord_flip()+
  xlab ("Taxon")+
  ylab ("Time bin")+
  
  scale_y_reverse("Age (Ma)")


ggsave (file = here ("output", "genus_range.pdf"))


# time trend
time_bins <- unique(coll_occ_taxa_perm_cret$bin_assignment)[order(unique(coll_occ_taxa_perm_cret$bin_assignment))]

data.frame (ntaxa=colSums(table_taxon_period),
            int = rev(bins$mid_ma[7:39])) %>%
#[match (time_bins,bins$bin)])   
  ggplot (aes(x=int,y=ntaxa))+
  ylab ("Number of genera")+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  scale_x_reverse("Age (Ma)") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 70,260), 
    ylim = c(0, 100),
    pos = list("b", "b"),
    size = list(3, 6),
    abbrv = list(TRUE, FALSE)
  ) 

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




# data to divDyn
divDyn_data <- coll_occ_taxa_perm_cret %>%
  select (unique_name, bin_assignment, 
          Interval,
          min_ma, max_ma, collection_no, paleolat2.x,paleolng2.x)

# write
write.csv(divDyn_data, file = here ("processed_data", "divDyn_data.csv"))


# --------------------------------------




# organize dataset to hierarchical models
# organize a table for each taxon in the data
# list of all variables (species, sites, intervals, formations, litologies)



cynodontia_data <- unique(coll_occ_taxa_perm_cret$unique_name)[order (unique(coll_occ_taxa_perm_cret$unique_name))]
time_bins <- unique(coll_occ_taxa_perm_cret$bin_assignment)[order(unique(coll_occ_taxa_perm_cret$bin_assignment))]
cells <- unique(coll_occ_taxa_perm_cret$cell_ID)[order( unique(coll_occ_taxa_perm_cret$cell_ID))]
formations <- unique(coll_occ_taxa_perm_cret$formation.x)[order(unique(coll_occ_taxa_perm_cret$formation.x))]
formations_lithologies <- unique(coll_occ_taxa_perm_cret$formation_lithology )[order( unique(coll_occ_taxa_perm_cret$formation_lithology ))]
#subcells <- unique(coll_occ_taxa_perm_cret$cell_ID_sub)[order(unique(coll_occ_taxa_perm_cret$cell_ID_sub))]


# the basic table summarizing all detections across intervals and formations

table_data_basis <- lapply (time_bins, function (i) {
  
    
  tab_basis <- cast (formula = cell_ID ~ formation.x,#cell_ID_sub,#
                     data = coll_occ_taxa_perm_cret[which(coll_occ_taxa_perm_cret$bin_assignment %in% i),],
                     value = "detection", 
                     fun.aggregate = sum,
                     drop=F,
                     fill=NA)
  
  
  # order names
  rownames(tab_basis) <- tab_basis$cell_ID
  tab_basis <- tab_basis[,-1]
  
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
                           filter (bin_assignment %in% i)
  
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
)>=3
)
barplot(table(sel_cynodontia))



# table  for naive occupancy
table_naive_cynodontia <- table_data [which(sel_cynodontia==T)]
save (table_naive_cynodontia,  file = here ("processed_data", "table_naive_cynodontia.RData"))


# filter
table_data_long <- table_data_long [which(sel_cynodontia==T)]

# melt
table_data_long_df <- do.call(rbind,table_data_long) # [1:2]

# match new names
df_match <- data.frame (names=unique(table_data_long_df$taxon),
            new_names=seq(1,length(unique(table_data_long_df$taxon))))
table_data_long_df$taxon_code <- df_match$new_names[match (table_data_long_df$taxon,df_match$names)]

# bind temperature to these data
# i = 987
# table_data_long_df[i,]
temperature_det <- lapply (seq (1,nrow (table_data_long_df)), function (i)


  mean(coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$cell_ID %in% rownames(table_data_basis[[1]]) [table_data_long_df [i,"site"]] 
                                  & 
                                   coll_occ_taxa_perm_cret$bin_assignment %in% time_bins [table_data_long_df [i,"int"]] 
                                 &
                                   coll_occ_taxa_perm_cret$formation.x %in% formations [table_data_long_df [i,"form"]] 
                                   
                                   ),"temperature"])
)

# paleolat
paleolat_det <- lapply (seq (1,nrow (table_data_long_df)), function (i)
  
  
  mean(coll_occ_taxa_perm_cret [which(coll_occ_taxa_perm_cret$cell_ID %in% rownames(table_data_basis[[1]]) [table_data_long_df [i,"site"]] 
                                      & 
                                        coll_occ_taxa_perm_cret$bin_assignment %in% time_bins [table_data_long_df [i,"int"]] 
                                      &
                                        coll_occ_taxa_perm_cret$formation.x %in% formations [table_data_long_df [i,"form"]] 
                                      
  ),"paleolat2.x"])
)


# bind temperature
table_data_long_df <- cbind (table_data_long_df,
                             temp = unlist(temperature_det),
                             paleolat=unlist(paleolat_det))

hist(unlist(paleolat_det))

naive<-cast (formula=taxon_code ~ int,
      data = table_data_long_df,
      value = "det",
      fun.aggregate = sum)
plot(rev(colSums(naive[,-1]>0)),type="l")


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
      file = here ("processed_data","table_data_array_cynodontia.RData"))


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

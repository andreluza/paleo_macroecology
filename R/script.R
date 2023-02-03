# =====================================================
# load packages
require(here);require(rgeos);require(rgdal);require(sp);require(raster); require(openxlsx)
# data processing
library(ggplot2); library(rasterVis); library(rgdal);library(viridis)



# ======================================================
#  Brazil data
tab_data_BR <- readOGR(here("data", "GeoSBG", "ocorrencia_fossilifera"),
         layer = "ocorrencia_fossilifera")
# (tab_data_BR@coords)


#matched_data_coordinates_mpas_institutions <- matched_data_coordinates_mpas_institutions[which(matched_data_coordinates_mpas_institutions$data == "function"),]
# load south america shapefile
# for mapping, cropping, etc
southAme<- readOGR(dsn= here("data","South_America"),encoding="latin1", 
                   layer="South_America")



# map resolution (cell size)
map_resolution<-1

# create a grid for mapping
# based on the extent of extracted data
grd_df <- expand.grid(x = seq(from = extent (southAme)[1]-8,
                              to = extent (southAme)[2]+8, 
                              by = map_resolution),
                      y = seq(from = extent (southAme)[3]-8,                                           
                              to = extent (southAme)[4]+8, 
                              by = map_resolution))  # expand points to grid

# Convert grd object to a matrix and then turn into a spatial
# points object
coordinates(grd_df) <- ~x + y

# Sp points into raster
grd_raster <- (raster(grd_df,resolution = map_resolution))
crs(grd_raster) <-crs(southAme)

# rasterize to count the number of points per cell
overlap_grid_dois <- rasterize((tab_data_BR@coords),
                                 grd_raster,
                                 fun="count")
#clamp
#overlap_grid_dois<-clamp(overlap_grid_dois, 
#                           lower=0,
#                           useValues=T)

# maps
# plot - number of sites per cell
plot1 <- gplot(overlap_grid_dois) +
  geom_tile(aes(x=x, y=y, fill=value), alpha=1) + 
  coord_fixed (xlim = c( -80, -25), 
               ylim = c(-40, 10), ratio = 1) +
  scale_fill_viridis(option="magma",direction=-1,begin=0,
                     breaks = seq (range(values(overlap_grid_dois),na.rm=T)[1],
                                   range(values(overlap_grid_dois),na.rm=T)[2],
                                   100),
                     limits=c(range(values(overlap_grid_dois),na.rm=T)[1]+1,
                              range(values(overlap_grid_dois),na.rm=T)[2]),
                     na.value=NA,
                     name="Count of coordinates") +
  ggtitle ("Fossil record,\nGeological Service of Brazil (2022)") + 
  theme_classic() +
  theme (legend.position = "top",
         legend.direction = "horizontal") + 
  xlab("Longitude")+
  ylab("Latitude")

# add south america map
plot1 <- plot1 + geom_polygon(data=southAme, 
                                    aes(x=long, y=lat, group=group),
                                    size = 0.1, fill="gray60", 
                                    colour="gray75",alpha=0.1) + 
  xlab("Longitude") + ylab("Latitude")


plot1



# ---------------------------------------------------
# paleodata

paleodata <- read.csv (here ("data", "PaleoData", "pbdb_data_collections.csv"),
                       skip=16)


# rasterize to count the number of points per cell
overlap_grid_dois <- rasterize(paleodata[,c("lng","lat")],
                               grd_raster,
                               fun="count")
#clamp
#overlap_grid_dois<-clamp(overlap_grid_dois, 
#                           lower=0,
#                           useValues=T)

# maps
# plot - number of sites per cell
plot2 <- gplot(overlap_grid_dois) +
  geom_tile(aes(x=x, y=y, fill=value), alpha=1) + 
  coord_fixed (xlim = c( -80, -25), 
               ylim = c(-40, 10), ratio = 1) +
  scale_fill_viridis(option="magma",direction=-1,begin=0,
                     breaks = seq (range(values(overlap_grid_dois),na.rm=T)[1],
                                   range(values(overlap_grid_dois),na.rm=T)[2],
                                   150),
                     limits=c(range(values(overlap_grid_dois),na.rm=T)[1]+1,
                              range(values(overlap_grid_dois),na.rm=T)[2]),
                     na.value=NA,
                     name="Count of coordinates") +
  ggtitle ("Fossil record,\nPaleobiology Database (2022)") + 
  theme_classic() +
  theme (legend.position = "top",
         legend.direction = "horizontal") + 
  xlab("Longitude")+
  ylab("Latitude")

# add south america map
plot2 <- plot2 + geom_polygon(data=southAme, 
                              aes(x=long, y=lat, group=group),
                              size = 0.1, fill="gray60", 
                              colour="gray75",alpha=0.1) + 
  xlab("Longitude") + ylab("Latitude")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank())


plot2

# arrange
gridExtra::grid.arrange(plot1,plot2,nrow=1)

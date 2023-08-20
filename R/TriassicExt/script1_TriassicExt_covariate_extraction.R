# using geographic information to gather cell covarites
# help webpasges: https://rpubs.com/boyerag/297592
# https://pjbartlein.github.io/REarthSysSci/netCDF.html#data-frame-to-array-conversionrectangular-to-raster


source(here ("R","packages.R"))


# ------------------------------------------------------------

# function do deal with netcdf data
# extract variables and return a raster 
function_deal_with_netcdf <- function (netcdf_data, name_long = "lon", name_lat = "lat",name="z") {
  
          # get lat, long, and the variable ("z" in this case)
          lon <- ncvar_get(netcdf_data, name_long)
          lat <- ncvar_get(netcdf_data, name_lat, verbose = F)
          t <- ncvar_get(netcdf_data, name)
          
          # get altitude
          dname <- name
          tmp_array <- ncvar_get(netcdf_data,dname)
          dlname <- ncatt_get(netcdf_data,dname,"long_name")
          dunits <- ncatt_get(netcdf_data,dname,"units")
          fillvalue <- ncatt_get(netcdf_data,dname,"_FillValue")
          # dim(tmp_array)
          # nc_close(nc_data)
          # replace netCDF fill values with NA's
          tmp_array[tmp_array==fillvalue$value] <- NA
          
          # quick map
          # require(RColorBrewer)
          # image(lon,lat,tmp_array, col=rev(brewer.pal(10,"RdBu")))
          # better map
          #require(lattice)
          #grid <- expand.grid(lon=lon, lat=lat)
          #cutpts <- seq (range(tmp_array)[1], range(tmp_array)[2],
          #               1500)
          #levelplot(tmp_array ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          #          col.regions=(rev(brewer.pal(10,"RdBu"))))
          
          # convert into raster 
          r <- raster(t(tmp_array), 
                      xmn=min(lon), xmx=max(lon), 
                      ymn=min(lat), ymx=max(lat), 
                      crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
          
          # We will need to transpose and flip to orient the data correctly. 
          # The best way to figure this out is through trial and error, but remember that most netCDF files record spatial data from the bottom left corner.
          r <- flip(r, direction='y')
          
          # return
          return (r)
}

# ------------------------------------------------------------

# paleoelevation


# list of files
list_ncdf_files <- list.files(here ("data", "paleoDEMlayers", 
                                    "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc",
                                    "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2"),
                              pattern = "nc$")[14:60] # kt boundary to early permian
# revert the order
list_ncdf_files <- rev(list_ncdf_files)

# load
nc_data <- lapply(list_ncdf_files, function (i) 
  
  nc_open(here ("data", "paleoDEMlayers", 
                "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc",
                "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2",
                i))
)
# print (nc_data[[10]])
# extract data          
raster_paleoalt <- lapply (nc_data, function_deal_with_netcdf, name="z")

# stack (or brick)
brick_rasters <- brick (raster_paleoalt)
# plot(function_paleoalt (nc_data[[47]]), col=rev(brewer.pal(10,"RdBu")))

# arrange a table of ages
ages <- lapply (list_ncdf_files, function (i){
  frst <- strsplit (i, "_")
  age <- as.numeric(gsub ("Ma.nc","", frst[[1]] [length(frst[[1]])]))
  data.frame (file = i,
              age=age)
})
ages<-do.call(rbind,ages)

# find the interval of each raster
intervals <- openxlsx::read.xlsx (here ("data", "periods_intervals_ages.xlsx"))
# reverse the order to the past be the first period
rev_order <- rev(intervals$Interval)
intervals <- intervals [match (intervals$Interval, rev_order),]

# extract inter val names
ages_intervals <- lapply (ages$age, function (i)

  intervals [which(intervals[,4] >= i & intervals[,5] <=   i)[1],  # there one reconstruction at a J-K boundary (choose J)
             c("Period", "Interval") ]

)
ages_intervals <- do.call(rbind, ages_intervals)
ages_intervals$age_paleoalt <- ages$age

# bind original intervals to check

cbind (intervals,
  int=ages_intervals[match (intervals$Interval,ages_intervals$Interval),"Interval"])

# some ages miss paleoalt
brick_rasters<-stackApply(brick_rasters, 
                          indices =  ages_intervals$Interval, 
                          fun = "mean", na.rm = T)
names(brick_rasters) <- gsub ("index_", "", names(brick_rasters))



# save 
save (brick_rasters, file = here ("processed_data", 
                                  "brick_rasters_elevation.RData"))



# -----------------------------------------------------------


# paleoprecipitation


# list of files
list_ncdf_files_prec <- list.files(here ("data", 
                                         "PaleoPrecipitation"),
                              pattern = "nc$")
# reverse the order
list_ncdf_files_prec <- rev(list_ncdf_files_prec)

# load
nc_data_prec <- lapply(list_ncdf_files_prec, function (i) 
  
  nc_open(here ("data", "PaleoPrecipitation",
                i))
)
names(nc_data_prec) <- list_ncdf_files_prec

# match to find age
metadata_prec <- read.xlsx (here ("data","PaleoPrecipitation", "list_ncdf_files_prec.xlsx"))
nc_data_prec <- nc_data_prec [match (metadata_prec$file_name,names(nc_data_prec))]

# print (nc_data_prec[[1]])

raster_paleoprec <- lapply (nc_data_prec, function_deal_with_netcdf, 
                            name_long = "longitude", name_lat= "latitude",
                            name="precip_mm_srf")

# plot
par (mfrow=c(2,3))
lapply (raster_paleoprec, plot)
plot(raster_paleoprec[[1]] - raster_paleoprec[[2]])


# stack (or brick)
brick_rasters_paleoprec <- brick (raster_paleoprec)
# plot(function_paleoalt (nc_data[[47]]), col=rev(brewer.pal(10,"RdBu")))


# find the interval of each raster
ages_intervals <- lapply ((metadata_prec$age), function (i)
  
  intervals [which(intervals[,4] > i & intervals[,5] <=   i), c("Period", "Interval") ]
  
)
ages_intervals <- do.call(rbind, ages_intervals)
ages_intervals$age_paleoprec <- (metadata_prec$age)

# bind original intervals to check

cbind (intervals,
       int=ages_intervals[match (intervals$Interval,ages_intervals$Interval),"Interval"])



# some ages miss paleoalt
brick_rasters_paleoprec<-stackApply(brick_rasters_paleoprec, 
                                    indices =  ages_intervals$Interval, 
                                    fun = "mean", na.rm = T)

# rotate to make it -180- 180 long
# https://stackoverflow.com/questions/70312970/change-extent-in-map-from-0-360-0-300-to-180-180-90-90
brick_rasters_paleoprec <- rotate(brick_rasters_paleoprec)


# resample to have the complete raster 
brick_rasters_paleoprec <- resample (brick_rasters_paleoprec, brick_rasters)
names(brick_rasters_paleoprec) <- gsub ("index_", "", names(brick_rasters_paleoprec))


# save 
save (brick_rasters_paleoprec, file = here ("processed_data", 
                                            "brick_rasters_precipitation.RData"))


# --------------------------------------------------------------

# paleotemperature (mean annual temperature)


# list of files
list_ncdf_files_temp <- list.files(here ("data", 
                                         "PaleoTemperature"),
                                   pattern = "nc$")
list_ncdf_files_temp <- rev(list_ncdf_files_temp)

# load
nc_data_temp <- lapply(list_ncdf_files_temp, function (i) 
  
  nc_open(here ("data", "PaleoTemperature",
                i))
)
names(nc_data_temp) <- list_ncdf_files_temp

# match to find age
metadata_temp <- read.xlsx (here ("data","PaleoTemperature", "list_ncdf_files_temp.xlsx"))
nc_data_temp <- nc_data_temp [match (metadata_temp$file_name,names(nc_data_temp))]

# print (nc_data_prec[[1]])

raster_paleotemp <- lapply (nc_data_temp, function_deal_with_netcdf, 
                            name_long = "lon", name_lat= "lat",
                            name="mat")

# plot
par (mfrow=c(2,3))
lapply (raster_paleotemp, plot)
plot(raster_paleotemp[[1]] - raster_paleotemp[[2]])


# stack (or brick)
brick_rasters_paleotemp <- brick (raster_paleotemp)
# plot(raster_paleotemp[[1]], col=rev(brewer.pal(10,"RdBu")))


# find the interval of each raster
ages_intervals <- lapply ((metadata_temp$age), function (i)
  
  intervals [which(intervals[,4] > i & intervals[,5] <=   i), c("Period", "Interval") ]
  
)
ages_intervals <- do.call(rbind, ages_intervals)
ages_intervals$age_paleoprec <- (metadata_temp$age)

# bind original intervals to check

cbind (intervals,
       int=ages_intervals[match (intervals$Interval,ages_intervals$Interval),"Interval"])



# some ages miss paleoalt
brick_rasters_paleotemp<-stackApply(brick_rasters_paleotemp, 
                                    indices =  ages_intervals$Interval, 
                                    fun = "mean", na.rm = T)
# rotate to make it -180- 180 long
# https://stackoverflow.com/questions/70312970/change-extent-in-map-from-0-360-0-300-to-180-180-90-90
brick_rasters_paleotemp <- rotate(brick_rasters_paleotemp)

# resample to have the complete raster 
brick_rasters_paleotemp <- resample (brick_rasters_paleotemp, brick_rasters)
names(brick_rasters_paleoprec) <- gsub ("index_", "", names(brick_rasters_paleoprec))

# save 
save (brick_rasters_paleotemp, file = here ("processed_data", 
                                            "brick_rasters_temp.RData"))
par (mfrow=c(1,1))

rm(list=ls())


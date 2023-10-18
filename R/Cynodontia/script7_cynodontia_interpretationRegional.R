# --------------------------------------------


#     Interpretation: Region-scale analysis


# ------------------------------------------------------------

# load packages
source("R/packages.R")
source("R/functions.R")

# load output and data
# LOAD BERNOULLI ESTIMATES
require(coda)
require(mcmcr)
require(scales)

# ------------------------------------------------------------ #

# load basic data

load(here ("processed_data", "grid_info.RData"))
load(here("processed_data", "CMR_data.RData"))
load(here("processed_data", "site_covs.RData"))
load(here ("processed_data", "paleomaps.RData"))


# load richness estimates
load(here("output", "interesting_params_occ.RData"))
load( here ("output", "interesting_params.RData") )


# build matrix
SR_matrix <- matrix( interesting_params_occ$SRexp$statistics[,"Mean"],
                     nrow=33,
                     ncol=44,
                     byrow=F)
colnames(SR_matrix) <- cells
rownames(SR_matrix) <- time_bins

         
# map maastritchian
p1 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[33,] [match (grid_data$cell_ID, names(SR_matrix[33,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Maastrichtian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[1]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)



# map cenomanian
p2 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[28,] [match (grid_data$cell_ID, names(SR_matrix[28,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Cenomanian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[5]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Kimmeridgian 
p3 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[20,] [match (grid_data$cell_ID, names(SR_matrix[20,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Kimmeridgian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[13]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Kimmeridgian 
p4 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[14,] [match (grid_data$cell_ID, names(SR_matrix[14,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Toarcian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[19]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Carnian  
p5 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[8,] [match (grid_data$cell_ID, names(SR_matrix[8,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Carnian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[25]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Olenekian
p6 <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrix[5,] [match (grid_data$cell_ID, names(SR_matrix[5,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,120,20),
                       limits=c(0,120))+
  ggtitle ("Olenekian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[28]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)


# uncertainty map
# build matrix
SR_matrixu <- matrix( interesting_params_occ$SRexp$statistics[,"SD"],
                     nrow=33,
                     ncol=44,
                     byrow=F)
colnames(SR_matrixu) <- cells
rownames(SR_matrixu) <- time_bins


# map maastritchian
p1u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[33,] [match (grid_data$cell_ID, names(SR_matrixu[33,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Maastrichtian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[1]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map cenomanian
p2u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[28,] [match (grid_data$cell_ID, names(SR_matrixu[28,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Cenomanian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[5]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Kimmeridgian 
p3u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[20,] [match (grid_data$cell_ID, names(SR_matrixu[20,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Kimmeridgian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[13]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Kimmeridgian 
p4u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[14,] [match (grid_data$cell_ID, names(SR_matrixu[14,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Toarcian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[19]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Carnian  
p5u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[8,] [match (grid_data$cell_ID, names(SR_matrixu[8,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Carnian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[25]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)

# map Olenekian
p6u <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           SR = SR_matrixu[5,] [match (grid_data$cell_ID, names(SR_matrixu[5,]))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=SR),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,70,15),
                       limits=c(0,70))+
  ggtitle ("Olenekian")+
  geom_sf(data = st_transform(st_as_sf(maps_models[[28]]),
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0"),alpha=0.1)
  

pdf(here ("output", "Maps_SR.pdf"),height=12,width=8)

grid.arrange(p1+theme(legend.position = "none"),p1u+theme(legend.position = "none"),
             p2+theme(legend.position = "none"),p2u+theme(legend.position = "none"),
             p3+theme(legend.position = "none"),p3u+theme(legend.position = "none"),
             p4+theme(legend.position = "none"),p4u+theme(legend.position = "none"),
             p5+theme(legend.position = "none"),p5u+theme(legend.position = "none"),
             p6,p6u,
             
             nrow=6,ncol=2)
dev.off()


# number of stages and probabilities
# build matrix
phi_matrix <- matrix( interesting_params_occ$phi$statistics[,"Mean"],
                     nrow=32,
                     ncol=44,
                     byrow=F)
gamma_matrix <- matrix( interesting_params_occ$gamma$statistics[,"Mean"],
                      nrow=32,
                      ncol=44,
                      byrow=F)
colnames(gamma_matrix) <- colnames(phi_matrix) <- cells
rownames(gamma_matrix) <- rownames(phi_matrix) <- time_bins[-1]

# transform matrix (discreticize)
disc_gamma <- colSums(gamma_matrix >= 0.02)
disc_phi <- colSums(phi_matrix >= 0.6)

# map maastritchian
p1Or <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           Nstages = disc_gamma [match (grid_data$cell_ID, names(disc_gamma))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=Nstages),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,32,6),
                       limits=c(0,32))+
  ggtitle ("Origination probability")

# map cenomanian
p2Pers <- ggplot(sf::st_transform(
  grid_base,
  "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
)) +
  geom_sf(fill = NA, colour = 'black') +
  geom_sf(data= sf::st_transform(
    cbind (grid_data,
           Nstages = disc_phi [match (grid_data$cell_ID, names(disc_phi))]),
    "+proj=moll +lon_0=0 +x_0=0 +y_0=0"
  ), aes (fill=Nstages),colour = 'black',alpha=0.5)+
  theme_void()+
  scale_fill_viridis_c(option = "magma",na.value= "white",
                       breaks=seq(0,32,6),
                       limits=c(0,32))+
  ggtitle ("Persistence probability")

pdf(file = here ("output", "figures", "stages_probabs.pdf"))
grid.arrange(p1Or,p2Pers,nrow=2)
dev.off()

# organize variables
# scale variables
elevation<-(scale (site_covs$elevation)) # scaled elevation
temperature <- (scale (site_covs$temp)) # scaled temp
precipitation <-  (scale (site_covs$prec))
lat <- ( (site_covs$paleolat))

# reverse the order of columns
elevation<-elevation[,rev(colnames(elevation))]
temperature<-temperature[,rev(colnames(temperature))]
precipitation<-precipitation[,rev(colnames(precipitation))]

# bins in data
elevation <- elevation [, which(colnames(elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(elevation)==bins [which(bins$bin %in% time_bins),"interval_name"]
temperature <- temperature [, which(colnames(temperature) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(temperature)==bins [which(bins$bin %in% time_bins),"interval_name"]
precipitation <- precipitation [, which(colnames(precipitation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
colnames(precipitation)==bins [which(bins$bin %in% time_bins),"interval_name"]

# subsetting of cells
precipitation<-(precipitation[which(rownames(precipitation) %in% cells),])
temperature<-(temperature[which(rownames(temperature) %in% cells),])
elevation<-(elevation[which(rownames(elevation) %in% cells),])
lat<-(lat[which(rownames(elevation) %in% cells)])

# check the colours
library(RColorBrewer)
library(scales)
# cols <- brewer_pal(pal = "RdBu")(5) # not valid in 1.1-2
cols <- brewer.pal(n = 8, name = "RdBu") 
cols
# [1] "#CA0020" "#F4A582" "#F7F7F7" "#92C5DE" "#0571B0"
# show_col(cols) # not valid in 1.1-2
display.brewer.pal(n = 8, name = "RdBu")

# build matrix
SR_matrix <- matrix( interesting_params_occ$SRexp$statistics[,"Mean"],
                     nrow=33,
                     ncol=44,
                     byrow=F)
colnames(SR_matrix) <- cells
rownames(SR_matrix) <- time_bins

# melt
SR_df <- cbind (melt (as.data.frame(SR_matrix)),
                age=bins$mid_ma[-c(1:6)])

SR_df$Latitude <- lat [match (SR_df$variable,cells)]


# plot
p3 <- ggplot (SR_df, aes (x=age,y=value, 
                          col=Latitude,group=variable)) + 
  geom_point() + 
  geom_line() + 
  scale_colour_gradientn(colours = cols, 
                         values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                         guide = "colorbar", 
                         limits=c(-51, 83))+
  theme_bw()+
  ylab ("Taxonomic diversity") + 
  xlab ("Stages")+
  theme_bw()+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.8,0.9),
         legend.direction = "horizontal")+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(0, 500),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 
ggsave (p3,filename=here ("output", "figures","SR_space.png"))


# variation in genus richness
# try to incorporate uncertainty
SR_matrix <- matrix(interesting_params_occ$SRexp$statistics[,"Mean"],
                    nrow=33,
                    ncol=44,
                    byrow=F)
SR_matrix_lwr <- matrix(interesting_params_occ$SRexp$quantiles[,"2.5%"],
                        nrow=33,
                        ncol=44,
                        byrow=F)
SR_matrix_upr <- matrix(interesting_params_occ$SRexp$quantiles[,"97.5%"],
                        nrow=33,
                        ncol=44,
                        byrow=F)
# set names
rownames(SR_matrix) <- rownames(SR_matrix_lwr)<-rownames(SR_matrix_upr) <- time_bins
colnames(SR_matrix) <- colnames(SR_matrix_lwr)<- colnames(SR_matrix_upr) <- cells

# melt each one and then bind
# average
SR_df <- cbind (melt (as.data.frame(SR_matrix)),
                data = "Average",
                age=bins$mid_ma[-c(1:6)])
SR_df$Latitude <- lat [match (SR_df$variable,cells)]
# lower interval
SR_df_lwd <- cbind (melt (as.data.frame(SR_matrix_lwr)),
                    data = "Lower.CI",
                    age=bins$mid_ma[-c(1:6)])
SR_df_lwd$Latitude <- lat [match (SR_df_lwd$variable,cells)]
# upper interval
SR_matrix_upr <- cbind (melt (as.data.frame(SR_matrix_upr)),
                        data = "Upper.CI",
                        age=bins$mid_ma[-c(1:6)])
SR_matrix_upr$Latitude <- lat [match (SR_matrix_upr$variable,cells)]

# plot with uncertainty
p6 <- ggplot (data.frame (SR_df,
                          Lower.CI=SR_df_lwd$value,
                          Upper.CI =SR_matrix_upr$value), 
              
              aes (x=Latitude,y=value, 
                   col=age),
              position =  position_jitter(width=2)) + 
  
  
  geom_errorbar(aes(x = Latitude, 
                    ymin = Lower.CI, 
                    ymax = Upper.CI,
                    col=age),
                position =  position_jitter(width=2),
                alpha=0.3,
                width=0.01,
                size=1)+
  geom_point(size=2,position =  position_jitter(width=2))+
  geom_smooth()+
  theme_bw()+
  ylab ("Estimated richness")+
  theme(axis.title = element_text(size=15)) + 
  scale_colour_viridis(option = "magma",direction = 1,end=0.8)

p6

ggsave (filename=here ("output", "figures","SR_lat.png"))

# interesting parameters
# temporal trends


# function to produce plots with uncertainty around a desired parameter
parameter = "SRexp"
function_plot_sites <- function (parameter, label_yaxis,ymax1=1,ymax2=1) {
  
  
  # calculate statistics
  s <- interesting_params_occ [which(names (interesting_params_occ) %in% parameter)]
  
  # variation in genus richness
  # try to incorporate uncertainty
  SR_matrix <- matrix(s[[1]]$statistics[,"Mean"],
                      nrow=33,
                      ncol=44,
                      byrow=F)
  SR_matrix_lwr <- matrix(s[[1]]$quantiles[,"2.5%"],
                          nrow=33,
                          ncol=44,
                          byrow=F)
  SR_matrix_upr <- matrix(s[[1]]$quantiles[,"97.5%"],
                          nrow=33,
                          ncol=44,
                          byrow=F)
  # set names
  rownames(SR_matrix) <- rownames(SR_matrix_lwr)<-rownames(SR_matrix_upr) <- time_bins
  colnames(SR_matrix) <- colnames(SR_matrix_lwr)<- colnames(SR_matrix_upr) <- cells
  
  # melt each one and then bind
  # average
  SR_df <- cbind (melt (as.data.frame(SR_matrix)),
                  data = "Average",
                  age=bins$mid_ma[-c(1:6)])
  SR_df$Latitude <- lat [match (SR_df$variable,cells)]
  # lower interval
  SR_df_lwd <- cbind (melt (as.data.frame(SR_matrix_lwr)),
                      data = "Lower.CI",
                      age=bins$mid_ma[-c(1:6)])
  SR_df_lwd$Latitude <- lat [match (SR_df_lwd$variable,cells)]
  # upper interval
  SR_matrix_upr <- cbind (melt (as.data.frame(SR_matrix_upr)),
                          data = "Upper.CI",
                          age=bins$mid_ma[-c(1:6)])
  SR_matrix_upr$Latitude <- lat [match (SR_matrix_upr$variable,cells)]
  
  # plot with uncertainty
  p6 <- ggplot (data.frame (SR_df,
                            Lower.CI=SR_df_lwd$value,
                            Upper.CI =SR_matrix_upr$value), 
                
                aes (x=age,y=(value), 
                     col=Latitude,group=variable),
                position =  position_jitterdodge()) + 
    
    geom_line(alpha=0.1)  + 
    
    geom_errorbar(aes(x = age, 
                      ymin = Lower.CI, 
                      ymax = Upper.CI,
                      col=Latitude),
                  position =  position_jitterdodge(),
                  alpha=0.5,
                  width=0.1,
                  size=1)+
    
    geom_point(size=2,position =  position_jitterdodge()) +
    
    scale_colour_gradientn(colours = cols, 
                           values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                           guide = "colorbar", 
                           limits=c(-51, 83))+
    theme_bw()+
    ylab (label_yaxis) + 
    xlab ("Stages")+
    theme_bw()+
    geom_line()+
    scale_x_reverse("Age (Ma)") +
    theme (legend.position = c(0.8,0.9),
           legend.direction = "horizontal")+
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,270), 
      ylim = c(0, ymax1),
      pos = list("b", "b"),
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) 
  
  
  # plot with uncertainty wihtout the permian AND begining of the triassic
  p7 <- ggplot (data.frame (SR_df,
                            Lower.CI=SR_df_lwd$value,
                            Upper.CI =SR_matrix_upr$value), 
                
                aes (x=age,y=(value), 
                     col=Latitude,group=variable),
                position =  position_jitterdodge()) + 
    geom_line(alpha=0.1)  + 
    
    geom_errorbar(aes(x = age, 
                      ymin = Lower.CI, 
                      ymax = Upper.CI,
                      col=Latitude),
                  position =  position_jitterdodge(),
                  alpha=0.5,
                  width=0.1,
                  size=1)+
    
    geom_point(size=2,position =  position_jitterdodge()) + 
    
    scale_colour_gradientn(colours = cols, 
                           values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                           guide = "colorbar", 
                           limits=c(-51, 83))+
    theme_bw()+
    ylab (label_yaxis) + 
    xlab ("Stages")+
    theme_bw()+
    geom_line()+
    scale_x_reverse("Age (Ma)") +
    theme (legend.position = c(0.8,0.9),
           legend.direction = "horizontal")+
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,240), 
      ylim = c(0, ymax2),
      pos = list("b", "b"),
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) 
  
  # arrange a grid
  
  
  pdf (here ("output", "figures", paste ("arrange_", parameter, ".pdf",sep="")),
       width=10,height=8)
  
  
  grid.arrange (p6,
                p7,
                ncol=1,
                nrow=2
  )
  
  dev.off()
  
}


# save origination prob
function_plot_sites(parameter = "SRexp", 
                    label_yaxis = "Taxonomic diversity",
                    ymax1=700,
                    ymax2=300)

# save origination prob
function_plot_sites(parameter = "gamma", 
                    label_yaxis = "Origination probability",
                    ymax1=0.7,ymax2=0.7)


function_plot_sites(parameter = "phi", 
                    label_yaxis = "Persistence probability",
                    ymax1=1,
                    ymax2=1)


# table of parameters


require(knitr)


rbind (
  
  data.frame (par =  "Origination",
              var = "Average (logit)",
              mean = interesting_params$intercept.gamma$statistics[1],
              lci = interesting_params$intercept.gamma$quantiles[1],
              uci = interesting_params$intercept.gamma$quantiles[2]),
  data.frame (par =  "Origination",
              var = "Elevation",
              mean = interesting_params$beta.gamma.elev$statistics[1],
              lci = interesting_params$beta.gamma.elev$quantiles[1],
              uci = interesting_params$beta.gamma.elev$quantiles[2]),
  
  data.frame (par =  "Origination",
              var = "Precipitation",
              mean = interesting_params$beta.gamma.prec$statistics[1],
              lci = interesting_params$beta.gamma.prec$quantiles[1],
              uci = interesting_params$beta.gamma.prec$quantiles[2]),
  data.frame (par =  "Origination",
              var = "Temperature",
              mean = interesting_params$beta.gamma.temp$statistics[1],
              lci = interesting_params$beta.gamma.temp$quantiles[1],
              uci = interesting_params$beta.gamma.temp$quantiles[2]),
  
  
  data.frame (par =  "Persistence",
              var = "Average (logit)",
              mean = interesting_params$intercept.phi$statistics[1],
              lci = interesting_params$intercept.phi$quantiles[1],
              uci = interesting_params$intercept.phi$quantiles[2]),
  data.frame (par =  "Persistence",
              var = "Elevation",
              mean = interesting_params$beta.phi.elev$statistics[1],
              lci = interesting_params$beta.phi.elev$quantiles[1],
              uci = interesting_params$beta.phi.elev$quantiles[2]),
  
  data.frame (par =  "Persistence",
              var = "Precipitation",
              mean = interesting_params$beta.phi.prec$statistics[1],
              lci = interesting_params$beta.phi.prec$quantiles[1],
              uci = interesting_params$beta.phi.prec$quantiles[2]),
  data.frame (par =  "Persistence",
              var = "Temperature",
              mean = interesting_params$beta.phi.temp$statistics[1],
              lci = interesting_params$beta.phi.temp$quantiles[1],
              uci = interesting_params$beta.phi.temp$quantiles[2]),
  
  
  data.frame (par =  "Detection",
              var = "Average (logit)",
              mean = interesting_params$intercept.p$statistics[1],
              lci = interesting_params$intercept.p$quantiles[1],
              uci = interesting_params$intercept.p$quantiles[2]),
  data.frame (par =  "Detection",
              var = "Time",
              mean = interesting_params$beta.p.time$statistics[1],
              lci = interesting_params$beta.p.time$quantiles[1],
              uci = interesting_params$beta.p.time$quantiles[2]),
  
  data.frame (par =  "Detection",
              var = "Range size",
              mean = interesting_params$beta.p.range$statistics[1],
              lci = interesting_params$beta.p.range$quantiles[1],
              uci = interesting_params$beta.p.range$quantiles[2]),
  data.frame (par =  "Detection",
              var = "Latitude",
              mean = interesting_params$beta.p.lat$statistics[1],
              lci = interesting_params$beta.p.lat$quantiles[1],
              uci = interesting_params$beta.p.lat$quantiles[2])
  
  
  
  
) %>%
  
  #kable(., format = "pipe", padding = 2,align="c") 
  
  
  # Show the between-S CI's in red, and the within-S CI's in black
  ggplot(aes(x=var, y=mean, group=par)) +
      geom_errorbar(width=.1, aes(ymin=lci, ymax=uci), colour="red") +
      geom_point(shape=21, size=3, fill="white") +
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Coefficient value")+
  xlab ("Model parameter")+
  
  coord_flip()+
  
  facet_wrap(~par)+
  theme_bw()

ggsave (here ("output", "figures","parameters.png"))

# fit statistics
# Rhat
data.frame(Rhat = unlist(do.call(rbind,rhat_params)[,"intercept.gamma"]),
           par = rownames(do.call(rbind,rhat_params))) %>% 
  ggplot(aes(x=Rhat, y=par)) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  xlab ("Rhat value")+
  ggtitle ("Convergence in the region-scale model")

ggsave (here ("output","figures", "RHat_region.png"))

# fit statistics
# Rhat

rhat_gamma <- data.frame(Rhat = melt(rhat_occ_params$gamma$gamma)) %>% 
  filter (is.na(Rhat.value) != T) %>%
  ggplot(aes(x=Rhat.value, y=paste (Rhat.X1,Rhat.X2,sep="."))) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  theme (axis.text.y = element_text(size=1))+
  xlab ("Rhat value")+
  ylab ("Region x stage combination") +
  ggtitle ("Convergence in the region-scale model (Origination probability)")

# phi
rhat_phi <- data.frame(Rhat = melt(rhat_occ_params$phi$phi)) %>% 
  filter (is.na(Rhat.value) != T) %>%
  ggplot(aes(x=Rhat.value, y=paste (Rhat.X1,Rhat.X2,sep="."))) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  theme (axis.text.y = element_text(size=1))+
  xlab ("Rhat value")+
  ylab ("Region x stage combination") +
  ggtitle ("Convergence in the region-scale model (Persistence probability)")



# SR
rhat_SR <- data.frame(Rhat = melt(rhat_occ_params$SRexp$SRexp)) %>% 
  filter (is.na(Rhat.value) != T) %>%
  ggplot(aes(x=Rhat.value, y=paste (Rhat.X1,Rhat.X2,sep="."))) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  theme (axis.text.y = element_text(size=1))+
  xlab ("Rhat value")+
  ylab ("Region x stage combination") +
  ggtitle ("Convergence in the region-scale model (Taxonomic diversity)")


png (here ("output", "figures", "convergence_occ.png"),height=30,width=15,units ="cm",res=300)
grid.arrange(rhat_gamma,
             rhat_phi,
             rhat_SR,
             nrow=3,ncol=1)
dev.off()


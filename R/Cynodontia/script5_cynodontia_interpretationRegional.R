# --------------------------------------------


#     Interpretation: Region-scale analysis


# ------------------------------------------------------------

# load packages
source("R/packages.R")
source("R/functions.R")

# load output and data
load (here ("output",
            "samples_cynodontia_CMR_sites.RData"))
load(here ("processed_data", "grid_info.RData"))
load(here("processed_data", "CMR_data.RData"))
load(here("processed_data", "site_covs.RData"))

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

# average covariate effects
sel_cols <- c("mean", "2.5%", "97.5%", "Rhat", "n.eff")

require(knitr)

data.frame (round (
  
    rbind (
  
    
      samples_cynodontia_CMR_sites$summary [grep ("intercept", rownames(samples_cynodontia_CMR_sites$summary)),sel_cols],
  
    
      samples_cynodontia_CMR_sites$summary [grep ("beta", rownames(samples_cynodontia_CMR_sites$summary)),sel_cols])
  

    ,2
)) %>%

kable(., format = "pipe", padding = 2,align="c") 


# posterior exceedance probability
# detection
table(samples_cynodontia_CMR_sites$sims.list$beta.p.lat>0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.p.form>0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)

# persistence
table(samples_cynodontia_CMR_sites$sims.list$beta.phi.elev<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.phi.prec<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.phi.temp<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.phi.lat<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)


# origination
table(samples_cynodontia_CMR_sites$sims.list$beta.gamma.elev<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.gamma.prec<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.gamma.temp<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)
table(samples_cynodontia_CMR_sites$sims.list$beta.gamma.lat<0)/length(samples_cynodontia_CMR_sites$sims.list$beta.p.form)

# fit statistics
# Rhat
data.frame(Rhat = samples_cynodontia_CMR_sites$summary[,"Rhat"]) %>% 
  ggplot (aes (x=Rhat)) +
  geom_density(fill="red",alpha=0.2)+
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  xlab ("Rhat value")


# effective sample size
data.frame(Rhat = samples_cynodontia_CMR_sites$summary[,"n.eff"]) %>% 
  ggplot (aes (x=Rhat/600)) +
  geom_density(fill="red",alpha=0.2)+
  geom_vline(xintercept=1) + 
  theme_bw() + 
  xlab ("Number of Effective Samples From MCMC")



# variation in gamma and beta over time, per site
gamma_matrix <- samples_cynodontia_CMR_sites$mean$gamma
colnames(gamma_matrix) <- cells
rownames(gamma_matrix) <- time_bins
# melt
gamma_df <- cbind (melt (as.data.frame(gamma_matrix)),
                   age=bins$mid_ma[-c(1:6)])
gamma_df$Latitude <- lat [match (gamma_df$variable,
                cells)]

# check the colours
library(RColorBrewer)
library(scales)
# cols <- brewer_pal(pal = "RdBu")(5) # not valid in 1.1-2
cols <- brewer.pal(n = 8, name = "RdBu") 
cols
# [1] "#CA0020" "#F4A582" "#F7F7F7" "#92C5DE" "#0571B0"
# show_col(cols) # not valid in 1.1-2
display.brewer.pal(n = 8, name = "RdBu")


# plot

p1 <- ggplot (gamma_df, aes (x=age,y=value, 
                       col=Latitude,group=variable)) + 
  geom_point() + 
  geom_line() + 
  scale_colour_gradientn(colours = cols, 
                         values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                         guide = "colorbar", 
                         limits=c(-51, 83))+
  theme_bw()+
  ylab ("Origination Probability") + 
  xlab ("Stages")+
  theme_bw()+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.8,0.9),
         legend.direction = "horizontal")+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(0, 0.035),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 

ggsave (p1,filename=here ("output","figures", "origination_space.png"))

# variation in extinction probability
# variation in gamma and beta over time, per site
phi_matrix <- samples_cynodontia_CMR_sites$mean$phi
colnames(phi_matrix) <- cells
rownames(phi_matrix) <- time_bins
# melt
phi_df <- cbind (melt (as.data.frame(phi_matrix)),
                   age=bins$mid_ma[-c(1:6)])

phi_df$Latitude <- lat [match (phi_df$variable,cells)]


# plot
p2<-ggplot (phi_df, aes (x=age,y=1-value, 
                       col=Latitude,group=variable)) + 
  geom_point() + 
  geom_line() + 
  scale_colour_gradientn(colours = cols, 
                         values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                         guide = "colorbar", 
                         limits=c(-51, 83))+
  theme_bw()+
  ylab ("Extinction Probability (1-phi)") + 
  xlab ("Stages")+
  theme_bw()+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.8,0.9),
         legend.direction = "horizontal")+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(0, 1),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 

ggsave (p2,filename=here ("output", "figures","extinction_space.png"))

# variation in genus richness
SR_matrix <- samples_cynodontia_CMR_sites$mean$SRexp
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
    ylim = c(0, 50),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 
ggsave (p3,filename=here ("output", "figures","SR_space.png"))


# plot
p4 <- ggplot (SR_df, aes (x=age,y=log(value), 
                          col=Latitude,group=variable)) + 
  geom_point() + 
  geom_line() + 
  scale_colour_gradientn(colours = cols, 
                         values = rescale(c(-50,-25, -10, 0, 10, 25,50, 75)),
                         guide = "colorbar", 
                         limits=c(-51, 83))+
  theme_bw()+
  ylab ("Taxonomic diversity (log)") + 
  xlab ("Stages")+
  theme_bw()+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.8,0.9),
         legend.direction = "horizontal")+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(0, 7.6),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 
ggsave (p4,filename=here ("output", "figures","SR_space_complete.png"))


# without the permian

# plot
p5 <- ggplot (SR_df, aes (x=age,y=(value), 
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
    xlim = c( 66,240), 
    ylim = c(0, 75),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 
ggsave (p5,filename=here ("output", "figures","SR_space_without_permian.png"))


# function to produce plots with uncertainty around a desired parameter

function_plot_sites <- function (parameter, label_yaxis,ymax1=1,ymax2=1) {

        # find the parameter
        parameter_to_search <-  which(names(samples_cynodontia_CMR_sites$sims.list) == parameter)
      
      # variation in genus richness
      # try to incorporate uncertainty
      SR_matrix <- apply (samples_cynodontia_CMR_sites$sims.list[[parameter_to_search]],
                          c(2,3),mean,na.rm=T)
      SR_matrix_lwr <- apply (samples_cynodontia_CMR_sites$sims.list[[parameter_to_search]],
                              c(2,3),quantile, 0.025,na.rm=T)
      SR_matrix_upr <- apply (samples_cynodontia_CMR_sites$sims.list[[parameter_to_search]],
                              c(2,3),quantile, 0.975,na.rm=T)
      # set names
      colnames(SR_matrix) <- colnames(SR_matrix_lwr)<-colnames(SR_matrix_upr) <- cells
      rownames(SR_matrix) <- rownames(SR_matrix_lwr)<- rownames(SR_matrix_upr) <- time_bins
      
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
                          col=Latitude,group=variable)) + 
        
        geom_line(alpha=0.2)  + 
        
        geom_errorbar(aes(x = age, 
                          ymin = Lower.CI, 
                          ymax = Upper.CI,
                          col=Latitude),
                          alpha=0.2,
                          width=0.1,
                      size=1)+
        
        geom_point(size=2) +
        
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
      
      
      # plot with uncertainty wihtout the permian AND being triassic
      p7 <- ggplot (data.frame (SR_df,
                                Lower.CI=SR_df_lwd$value,
                                Upper.CI =SR_matrix_upr$value), 
                    
                    aes (x=age,y=(value), 
                         col=Latitude,group=variable)) + 
        geom_line(alpha=0.2)  + 
        
        geom_errorbar(aes(x = age, 
                          ymin = Lower.CI, 
                          ymax = Upper.CI,
                          col=Latitude),
                      alpha=0.2,
                      width=0.1,
                      size=1)+
        
        geom_point(size=2) + 
        
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
function_plot_sites(parameter = "gamma", 
                    label_yaxis = "Origination probability",
                    ymax1=0.05,ymax2=0.05)

# save origination prob
function_plot_sites(parameter = "phi", 
                    label_yaxis = "Persistence probability",
                    ymax1=1,
                    ymax2=1)

# save origination prob
function_plot_sites(parameter = "SRexp", 
                    label_yaxis = "Taxonomic diversity",
                    ymax1=600,
                    ymax2=75)

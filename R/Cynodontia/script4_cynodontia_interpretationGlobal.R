

# --------------------------------------------


#     Interpretation: Global-scale analysis


# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")

# create a directory to receive the figures
dir.create (here ("output", "figures"))


# ------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "samples_paleo_cynodontia_covariates_CMR_stages_binomial.RData"))
binomial_output <- samples_paleo_cynodontia

# load model bernoulli output 
load (here ("output",
            "samples_paleo_cynodontia_covariates_CMR_stage.RData"))
bernoulli_output <- samples_paleo_cynodontia_stage

# load basic data
load (here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data.RData"))

# binning time intervals (from Paleoverse)
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)


# table of summary statistics of covariates

summary(colMeans(site_covs$elevation))
sd(colMeans(site_covs$elevation))
summary(colMeans(site_covs$temp))
sd(colMeans(site_covs$temp))
summary(colMeans(site_covs$prec))
sd(colMeans(site_covs$prec))
summary(site_covs$prec)


# bins in data
elevation <- site_covs$elevation [which(names(site_covs$elevation) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
temperature <- site_covs$temp [which(names(site_covs$temp) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]
precipitation <- site_covs$prec [which(names(site_covs$prec) %in% (bins [which(bins$bin %in% time_bins),"interval_name"]))]

# cells among the studied sites
elevation<-elevation[which(rownames (elevation) %in% cells),] 
temperature<-temperature[which(rownames (temperature) %in% cells),] 
precipitation<-precipitation[which(rownames (precipitation) %in% cells),] 

# scale variables
# for elevation, we will consider only cells above and at the sea level in a given stage
list_elevation <- as.vector(elevation)
mean_elevation_stage <- lapply (list_elevation, function (i)
  
  
  mean (i[which(i >=0)]) # mean elevation for cells above and at the sea level
  
)
#unlist
mean_elevation_stage <- unlist(mean_elevation_stage)
scaled_elevation <- scale (mean_elevation_stage)
# temperature, precipitation
scaled_temperature <- (scale (colMeans(temperature )))
scaled_precipitation <-  (scale (colMeans(site_covs$prec)))

# reverse the order of columns
scaled_elevation<-scaled_elevation[rev(rownames(scaled_elevation)),]
scaled_temperature<-scaled_temperature[rev(rownames(scaled_temperature)),]
scaled_precipitation<-scaled_precipitation[rev(rownames(scaled_precipitation)),]

# time ( preservation )
scaled_time_stage <- scale (bins[time_bins,"mid_ma"])

# genus spatial range
range_taxon_analysis <- (range_area_taxon[which(range_area_taxon$taxon %in% cynodontia_data),])
scaled_range <- scale(range_taxon_analysis$range_area)


# create a function to optimize plots
function_plots <- function (samples, samples_bin) {
        
  
  # plot covariates
       plot1<- rbind (
          data.frame (int = bins$mid_ma[7:39],
                      var= "Precipitation",
                      val = scaled_precipitation),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Temperature",
                      val = scaled_temperature),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Elevation",
                      val = scaled_elevation)
          
          
        ) %>%
          
          ggplot (aes(x=int,y=val, fill=var,group = var,col=var))+
          geom_point(position = position_jitter(width=0.1),size=3)+
          
          theme_bw()+
          geom_line()+
          scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_colour_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_x_reverse("Age (Ma)") +
          theme (legend.position = c(0.9,0.9))+
          ylab ("Standandized values of variables")+
          coord_geo(
            dat = list("stages", "periods"), 
            xlim = c( 66,270), 
            ylim = c(-5, 5),
            pos = list("b", "b"),
            size = list(2, 4),
            abbrv = list(TRUE, T)
          ) 
        
        
        
      # plot parameters
      dat <-   rbind (
          data.frame (int = bins$mid_ma[7:39],
                      var= "Origination probability",
                      model = "Bernoulli",
                      average = samples$mean$gamma,
                      lw=round(apply (samples$sims.list$gamma,2,quantile,0.025,na.rm=T),2),
                      up=round(apply (samples$sims.list$gamma,2,quantile,0.975,na.rm=T),2)),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Origination probability",
                      model = "Binomial",
                      average = samples_bin$mean$gamma,
                      lw=round(apply (samples_bin$sims.list$gamma,2,quantile,0.025,na.rm=T),2),
                      up=round(apply (samples_bin$sims.list$gamma,2,quantile,0.975,na.rm=T),2)),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Extinction probability",
                      model = "Bernoulli",
                      average = 1-samples$mean$phi,
                      lw=1-round(apply (samples$sims.list$phi,2,quantile,0.025,na.rm=T),2),
                      up=1-round(apply (samples$sims.list$phi,2,quantile,0.975,na.rm=T),2)),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Extinction probability",
                      model = "Binomial",
                      average = 1-samples_bin$mean$phi,
                      lw=1-round(apply (samples_bin$sims.list$phi,2,quantile,0.025,na.rm=T),2),
                      up=1-round(apply (samples_bin$sims.list$phi,2,quantile,0.975,na.rm=T),2))
          
        ) 
          
      # plot results
      # binmial model
      plot2 <- ggplot (data = dat %>% 
                    filter (model == "Binomial"),
                    
                    aes(x=int,y=average, fill=var,group = var,col=var))+
          geom_point(position = position_jitter(width=0.1),size=3)+
          #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_colour_viridis_d(option="magma",begin=0.1,end=0.5)+
            
          geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                        width=0.1,size=1,position = position_jitter(width=0.1))+
          
          geom_ribbon(aes(x=int, y=average, ymax=up, ymin=lw), 
                      alpha=0.2,fill="red") + 
            
            
            # bernoulli results
            # points - average
            geom_point(data = dat %>%
                         filter (model=="Bernoulli"),
                       position = position_jitter(width=0.1),size=3)+
            scale_colour_viridis_d(option="magma",begin=0.6,end=0.8)+
          
              # bars
            geom_errorbar(data = dat %>%
                            filter (model=="Bernoulli"),
                          aes(x = int, ymin = lw, ymax = up,col=var),
                          width=0.1,size=1,position = position_jitter(width=0.1))+
            
            # credible intervals
            geom_ribbon(data = dat %>%
                          filter (model=="Bernoulli"),
                        aes(x=int, y=average, ymax=up, ymin=lw), 
                        alpha=0.2, fill="cyan")+
            geom_point(position = position_jitter(width=0.1),size=3)+
            
            
            
          theme_bw()+
          geom_line()+
          scale_x_reverse("Age (Ma)") +
          theme (legend.position = c(0.9,0.9))+
          coord_geo(
            dat = list("stages", "periods"), 
            xlim = c( 66,270), 
            ylim = c(0, 0.8),
            pos = list("b", "b"),
            size = list(2, 4),
            abbrv = list(TRUE, T)
          ) 
        
       
      
      # plot of species richness
        
      dat1 <- rbind (
          data.frame (int = bins$mid_ma[7:39],
                      var= "Observed richness",
                      model="NA",
                      average = colSums(array_genus_bin),
                      lw=NA,
                      up=NA), 
          data.frame (int = bins$mid_ma[7:39],
                      var= "Expected richness",
                      model="Bernoulli",
                      average = samples$mean$SRexp,
                      lw=round(apply (samples$sims.list$SRexp,2,quantile,0.025),2),
                      up=round(apply (samples$sims.list$SRexp,2,quantile,0.975),2)),
          
          data.frame (int = bins$mid_ma[7:39],
                      var= "Expected richness",
                      model="Binomial",
                      average = samples_bin$mean$SRexp,
                      lw=round(apply (samples_bin$sims.list$SRexp,2,quantile,0.025),2),
                      up=round(apply (samples_bin$sims.list$SRexp,2,quantile,0.975),2))
          
          
          
          ) # %>%
          
          
      
      # plot the binomial results
      
     plot3 <-  ggplot (dat1 %>%
                filter (model=="Binomial"), aes(x=int,y=average, fill=var,group = var,col=var))+
        # points - average
        geom_point(position = position_jitter(width=0.1),size=3)+
       scale_colour_viridis_d(option="magma",begin=0.1,end=0.5)+
       
        # bars
        geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                        width=0.1,size=1,position = position_jitter(width=0.1))+
          # credible intervals
          geom_ribbon(data = dat1 %>%
                        filter (model=="Binomial"),
                        aes(x=int, y=average, ymax=up, ymin=lw), 
                      alpha=0.2, fill="red")+
        
        
        # bernoulli results
        # points - average
        geom_point(data = dat1 %>%
                     filter (model=="Bernoulli"),
                   position = position_jitter(width=0.1),size=3)+
       scale_colour_viridis_d(option="magma",begin=0.6,end=0.8)+
        # bars
        geom_errorbar(data = dat1 %>%
                        filter (model=="Bernoulli"),
                      aes(x = int, ymin = lw, ymax = up,col=var),
                      width=0.1,size=1,position = position_jitter(width=0.1))+
        
        # credible intervals
        geom_ribbon(data = dat1 %>%
                      filter (model=="Bernoulli"),
                    aes(x=int, y=average, ymax=up, ymin=lw), 
                    alpha=0.2, fill="cyan")+
        geom_point(position = position_jitter(width=0.1),size=3)+
       
       # include the observed data
       geom_point(data = dat1 %>%
                    filter (model=="NA"),
                  position = position_jitter(width=0.1),size=3, col ="black")+
       geom_line(data = dat1 %>%
                   filter (model=="NA"),col="black")+
       
        theme_bw()+
          geom_line()+
        theme (legend.position = c(0.1,0.9))+
          scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_colour_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_x_reverse("Age (Ma)") +
          coord_geo(
            dat = list("stages", "periods"), 
            xlim = c( 66,270), 
            ylim = c(0,max(samples_bin$mean$SRexp+30)),
            pos = list("b", "b"),
            size = list(2, 4),
            abbrv = list(TRUE, T)
          ) 
        
        # detection probability
        dat2 <- rbind (
          
                data.frame (int = bins$mid_ma[7:39],
                    var= "Detection probability in occupied sites",
                    model = "Bernoulli",
                    average = samples$mean$p,
                    lw=round(apply (samples$sims.list$p,2,quantile,0.025),2),
                    up=round(apply (samples$sims.list$p,2,quantile,0.975),2),
                    nForm = formations_per_interval),
                
                data.frame (int = bins$mid_ma[7:39],
                            var= "Detection probability in occupied sites",
                            model = "Binomial",
                            average = samples_bin$mean$p,
                            lw=round(apply (samples_bin$sims.list$p,2,quantile,0.025),2),
                            up=round(apply (samples_bin$sims.list$p,2,quantile,0.975),2),
                            nForm = formations_per_interval))
                
                
                
                
             # plot    
                
          
        plot4<-  ggplot (data = dat2 %>%
                           
                           # bernoulli model
                           filter (model == "Bernoulli"),
                           
                           aes(x=int,y=average, fill=var,group = var,col=var))+
          # pts
          geom_point(position = position_jitter(width=0.1),size=3)+
          # error
          geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                        width=0.1,size=1,position = position_jitter(width=0.1))+
          
          # credible intervals
          geom_ribbon(aes(x=int, y=average, ymax=up, ymin=lw), 
                      alpha=0.2,fill="cyan")+
          
          # binomial
          # pts
          geom_point(data = dat2 %>%
                       filter (model=="Binomial"),
                     position = position_jitter(width=0.1),size=3)+
          # error
          geom_errorbar(data = dat2 %>%
                          filter (model=="Binomial"),
                        aes(x = int, ymin = lw, ymax = up,col=var),
                        width=0.1,size=1,position = position_jitter(width=0.1))+
          
          # credible intervals
          geom_ribbon(data = dat2 %>%
                        filter (model=="Binomial"),
                      aes(x=int, y=average, ymax=up, ymin=lw), 
                      alpha=0.2, fill="red")+
          
          
          
          
          theme_bw()+
          geom_line()+theme (legend.position = "none")+
          scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_colour_viridis_d(option="magma",begin=0.3,end=0.7)+
          scale_x_reverse("Age (Ma)") +
          coord_geo(
            dat = list("stages", "periods"), 
            xlim = c( 66,270), 
            ylim = c(0,1),
            pos = list("b", "b"),
            size = list(2, 4),
            abbrv = list(TRUE, T)
          ) + 
          geom_line( aes(y=nForm/max(nForm))) +
          scale_y_continuous(
            
            # Features of the first axis
            name = "Detection probability in occupied sites",
            
            # Add a second axis and specify its features
            sec.axis = sec_axis(trans=~.*57, name="Number of geological formations per stage")
          )
        
        
        # covariates per site
        
        
        ret_panel <-gridExtra::grid.arrange (plot4,
                                             plot1,
                                             plot2,
                                             plot3)

        return (ret_panel)
        }

# save plot
pdf(here ("output","figures", "panel_params_binomial_bernoulli.pdf"),width=12,height=10)
function_plots(bernoulli_output, binomial_output)
dev.off()


# effect of covariates
# elevation

png(here ("output","figures", "panel_covariates.png"),width = 20, height = 15, units = "cm",res=300)

par(mfrow=c(2,3))

# bernoulli model
plot(seq(min(elevation), max(elevation),by=0.01),
     plogis(bernoulli_output$mean$intercept.gamma+
              bernoulli_output$mean$beta.gamma.elev*seq(min(elevation), max(elevation),by=0.01)
     ), type = "l", ylab = "Origination probability",xlab="",ylim = c(0,0.2))

legend("topright",legend = c("Binomial model",
                             "Bernoulli model"),
       lty=c(2,1),
       lwd=2,
       bty="n"
)

# add posterior distribution samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(elevation), max(elevation),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.gamma[i]+
                 bernoulli_output$sims.list$beta.gamma.elev[i]*
                 seq(min(elevation), max(elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01))
)

# add average
lines(seq(min(elevation), max(elevation),by=0.01),
      plogis(bernoulli_output$mean$intercept.gamma+
               bernoulli_output$mean$beta.gamma.elev*seq(min(elevation), max(elevation),by=0.01)
      ),lwd=2,col="black")


# binomial model results
# add posterior distribution samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(elevation), max(elevation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.elev[i]*
                 seq(min(elevation), max(elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01),
        lty=2)
)

# add average
lines(seq(min(elevation), max(elevation),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
               binomial_output$mean$beta.gamma.elev*seq(min(elevation), max(elevation),by=0.01)
      ),lwd=2,col="black",
      lty=2)




# -------------------------

# temperature 


plot(seq(min(temperature), max(temperature),by=0.01),
     plogis(bernoulli_output$mean$intercept.gamma+
              bernoulli_output$mean$beta.gamma.temp*seq(min(temperature), max(temperature),by=0.01)
     ), type = "l", ylab = "",xlab="",ylim = c(0,0.2))

# add samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(temperature), max(temperature),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.gamma[i]+
                 bernoulli_output$sims.list$beta.gamma.temp[i]*
                 seq(min(temperature), max(temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01))
)

# average
lines(seq(min(temperature), max(temperature),by=0.01),
      plogis(bernoulli_output$mean$intercept.gamma+
               bernoulli_output$mean$beta.gamma.temp*seq(min(temperature), max(temperature),by=0.01)
      ),lwd=2,col="black")


# binomial model

# add samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(temperature), max(temperature),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.temp[i]*
                 seq(min(temperature), max(temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01),
        lty=2)
)

# average
lines(seq(min(temperature), max(temperature),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
             binomial_output$mean$beta.gamma.temp*seq(min(temperature), max(temperature),by=0.01)
      ),lwd=2,col="black",
      lty=2)



# ----------------------------------------------

# precipitation


plot(seq(min(precipitation), max(precipitation),by=0.01),
     plogis(bernoulli_output$mean$intercept.gamma+
              bernoulli_output$mean$beta.gamma.prec*seq(min(precipitation), max(precipitation),by=0.01)
     ), type = "l", ylab = "",xlab="",ylim = c(0,0.2))


# posterior distribution samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(precipitation), max(precipitation),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.gamma[i]+
                 bernoulli_output$sims.list$beta.gamma.prec[i]*
                 seq(min(precipitation), max(precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01))
)

# average
lines(seq(min(precipitation), max(precipitation),by=0.01),
      plogis(bernoulli_output$mean$intercept.gamma+
               bernoulli_output$mean$beta.gamma.prec*seq(min(precipitation), max(precipitation),by=0.01)
      ),lwd=2,col="black")


# binomial model
# posterior distribution samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(precipitation), max(precipitation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.prec[i]*
                 seq(min(precipitation), max(precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01),lty=2)
)

# average
lines(seq(min(precipitation), max(precipitation),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
               binomial_output$mean$beta.gamma.prec*seq(min(precipitation), max(precipitation),by=0.01)
      ),lwd=2,col="black",
      lty=2)

#-------------------------------------------------------

# persistence probability


# elevation

plot(seq(min(elevation), max(elevation),by=0.01),
     plogis(bernoulli_output$mean$intercept.phi+
                bernoulli_output$mean$beta.phi.elev*seq(min(elevation), max(elevation),by=0.01)
     ), type = "l", ylab = "Persistence probability",xlab="Elevation (m)",ylim = c(0,1))



# samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(elevation), max(elevation),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.phi[i]+
                   bernoulli_output$sims.list$beta.phi.elev[i]*
                   seq(min(elevation), max(elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01))
)

# average
lines(seq(min(elevation), max(elevation),by=0.01),
      plogis(bernoulli_output$mean$intercept.phi+
                 bernoulli_output$mean$beta.phi.elev*seq(min(elevation), max(elevation),by=0.01)
      ),lwd=2,col="black")

# binomial model
# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(elevation), max(elevation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.elev[i]*
                   seq(min(elevation), max(elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01),lty=2)
)
# average
lines(seq(min(elevation), max(elevation),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
              binomial_output$mean$beta.phi.elev*seq(min(elevation), max(elevation),by=0.01)
      ),lwd=2,col="black",
      lty=2)


# -----------------------

# temperature 


plot(seq(min(temperature), max(temperature),by=0.01),
     plogis(bernoulli_output$mean$intercept.phi+
                bernoulli_output$mean$beta.phi.temp*seq(min(temperature), max(temperature),by=0.01)
     ), type = "l", ylab = "",xlab="Temperature (ºC)",ylim = c(0,1))

# samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(temperature), max(temperature),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.phi[i]+
                   bernoulli_output$sims.list$beta.phi.temp[i]*
                   seq(min(temperature), max(temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01))
)

# average
lines(seq(min(temperature), max(temperature),by=0.01),
      plogis(bernoulli_output$mean$intercept.phi+
                 bernoulli_output$mean$beta.phi.temp*seq(min(temperature), max(temperature),by=0.01)
      ),lwd=2,col="black")


# binomial
# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(temperature), max(temperature),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.temp[i]*
                   seq(min(temperature), max(temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01),lty=2)
)

# average
lines(seq(min(temperature), max(temperature),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
                 binomial_output$mean$beta.phi.temp*seq(min(temperature), max(temperature),by=0.01)
      ),lwd=2,col="black",lty=2)


# -------------------

# precipitation
plot(seq(min(precipitation), max(precipitation),by=0.01),
      plogis(bernoulli_output$mean$intercept.phi+
                 bernoulli_output$mean$beta.phi.prec*seq(min(precipitation), max(precipitation),by=0.01)
     ), type = "l", ylab = "",xlab="Precipitation (mm/day)",ylim = c(0,1))

# samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(precipitation), max(precipitation),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.phi[i]+
                   bernoulli_output$sims.list$beta.phi.prec[i]*
                   seq(min(precipitation), max(precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01))
)
# average
lines(seq(min(precipitation), max(precipitation),by=0.01),
      plogis(bernoulli_output$mean$intercept.phi+
                  bernoulli_output$mean$beta.phi.prec*seq(min(precipitation), max(precipitation),by=0.01)
      ),lwd=2,col="black")

# binomial model

# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(precipitation), max(precipitation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.prec[i]*
                   seq(min(precipitation), max(precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01),lty=2)
)
# average
lines(seq(min(precipitation), max(precipitation),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
                  binomial_output$mean$beta.phi.prec*seq(min(precipitation), max(precipitation),by=0.01)
      ),lwd=2,col="black",lty=2)





dev.off()



# plots of detection probability

# detection probability
scaled_formations <- scale(formations_per_interval)

# plot
par(mfrow=c(1,1))
plot(seq(min(scaled_formations), max(scaled_formations),by=0.01),
     plogis(
       bernoulli_output$mean$intercept.p+
              bernoulli_output$mean$beta.p.form*seq(min(scaled_formations), max(scaled_formations),by=0.01)
     ), type = "l", ylab = "Detection probability (p)",xlab="Number of formations",ylim = c(0,1))

# samples
lapply (seq(1,length(bernoulli_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(scaled_formations), max(scaled_formations),by=0.01),
        
        plogis(bernoulli_output$sims.list$intercept.p[i]+
                 bernoulli_output$sims.list$beta.p.form[i]*
                 seq(min(scaled_formations), max(scaled_formations),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01))
)

# average
lines(seq(min(scaled_formations), max(scaled_formations),by=0.01),
      plogis(bernoulli_output$mean$intercept.p+
               bernoulli_output$mean$beta.p.form*seq(min(scaled_formations), max(scaled_formations),by=0.01)
      ),lwd=2,col="black")


# origination probability
# posterior exceedance probability
sum(bernoulli_output$sims.list$beta.gamma.elev<0)/length(bernoulli_output$sims.list$beta.gamma.elev<0)
sum(bernoulli_output$sims.list$beta.gamma.temp >0)/length(bernoulli_output$sims.list$beta.gamma.elev<0)
sum(bernoulli_output$sims.list$beta.gamma.prec >0)/length(bernoulli_output$sims.list$beta.gamma.elev<0)

# posterior exceedance probability
sum(binomial_output$sims.list$beta.gamma.elev<0)/length(binomial_output$sims.list$beta.gamma.elev<0)
sum(binomial_output$sims.list$beta.gamma.temp >0)/length(binomial_output$sims.list$beta.gamma.elev<0)
sum(binomial_output$sims.list$beta.gamma.prec >0)/length(binomial_output$sims.list$beta.gamma.elev<0)

# persistence probability
# posterior exceedance probability
sum(bernoulli_output$sims.list$beta.phi.elev>0)/length(bernoulli_output$sims.list$beta.phi.elev<0)
sum(bernoulli_output$sims.list$beta.phi.temp <0)/length(bernoulli_output$sims.list$beta.phi.elev<0)
sum(bernoulli_output$sims.list$beta.phi.prec >0)/length(bernoulli_output$sims.list$beta.phi.elev<0)

# posterior exceedance probability
sum(binomial_output$sims.list$beta.phi.elev>0)/length(binomial_output$sims.list$beta.phi.elev<0)
sum(binomial_output$sims.list$beta.phi.temp <0)/length(binomial_output$sims.list$beta.phi.elev<0)
sum(binomial_output$sims.list$beta.phi.prec >0)/length(binomial_output$sims.list$beta.phi.elev<0)



par(mfrow=c(1,3),mar=c(7,4,15,1))

# change
# pch
plot(bernoulli_output$mean$propcH,type="b",ylim=c(-1,2))
lapply (seq(1,nrow(bernoulli_output$sims.list$propcH)), function (i)

  lines(bernoulli_output$sims.list$propcH[i,],col=rgb(0,1,1,alpha=0.01))
)
lines(bernoulli_output$mean$propcH,type="b")

# R0
plot(bernoulli_output$mean$R0,type="b",ylim=c(0,1))
lapply (seq(1,nrow(bernoulli_output$sims.list$R0)), function (i)
  
  lines(bernoulli_output$sims.list$R0[i,],col=rgb(0,1,1,alpha=0.01))
)
lines(bernoulli_output$mean$R0,type="b")

# RER
plot(bernoulli_output$mean$RER,type="b",ylim=c(0,40))
lapply (seq(1,nrow(bernoulli_output$sims.list$R0)), function (i)
  
  lines(bernoulli_output$sims.list$RER[i,],col=rgb(0,1,1,alpha=0.01))
)
lines(bernoulli_output$mean$RER,type="b")


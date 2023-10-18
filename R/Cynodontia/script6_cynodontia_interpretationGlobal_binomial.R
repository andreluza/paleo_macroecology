

# --------------------------------------------


#     Interpretation: Global-scale analysis


# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")



# ------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "interesting_params_global_occ.RData"))

# load basic data
load (here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))

# binning time intervals (from Paleoverse)
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)


# table of summary statistics of covariates

apply (time_covariates[,-5],2,mean,na.rm=T)
apply (time_covariates[,-5],2,sd,na.rm=T)


# create a function to optimize plots
# plot covariates
plot1<- ggplot (data = time_covariates %>%
                  cbind (group = 1,
                         time = bins$mid_ma[-c(1:6)]),
                
                aes(x=time,y=temperature, #fill=temperature,col=temperature,
                    group = group))+
  
   geom_point(position = position_jitter(width=0.1),size=3,col= "#F83E4B")+
   
   theme_bw()+
   geom_line()+
  geom_ribbon(aes(x=time, y=temperature, 
                  ymax=temperature+temperature_sd, 
                  ymin=temperature-temperature_sd), 
              alpha=0.2,fill="#F83E4B") + 
  
   scale_x_reverse("Age (Ma)") +
   theme (legend.position = "none")+
   coord_geo(
     dat = list("stages", "periods"), 
     xlim = c( 66,270), 
     ylim = c(-15, 45),
     pos = list("b", "b"),
     size = list(2, 4),
     abbrv = list(TRUE, T)
   ) + 
  geom_line( aes(y=precipitation*10),col = "black") +
  geom_point(aes(y=precipitation*10),size=3,col= "#00917C")+
  geom_ribbon(aes(x=time, y=precipitation*10, 
                  ymax=precipitation*10+precipitation_sd*10, 
                  ymin=precipitation*10-precipitation_sd*10), 
              alpha=0.2,fill="#00917C") + 
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Temperature (ºC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~./10, name="Precipitation (mm/day)")
    
  )

        
# plot parameters
dat <-   rbind (
    
    
    data.frame (int = bins$mid_ma[7:39][-1],
                Taxon = "Non-mammaliaform cynodonts",
                var= "Origination probability",
                model = "Binomial",
                average = matrix(interesting_params_occ$gamma$stat$statistics[,"Mean"], ncol=3,byrow=T)[,2],
                lw= matrix(interesting_params_occ$gamma$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,2],
                up=matrix(interesting_params_occ$gamma$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,2]),
    
    data.frame (int = bins$mid_ma[7:39][-1],
                Taxon = "Non-mammaliaform cynodonts",
                var= "Extinction probability",
                model = "Binomial",
                average = 1-matrix(interesting_params_occ$phi$stat$statistics[,"Mean"], ncol=3,byrow=T)[,2],
                lw=1-matrix(interesting_params_occ$phi$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,2],
                up=1-matrix(interesting_params_occ$phi$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,2]),
    
    data.frame (int = bins$mid_ma[7:39][-1],
                Taxon = "Non-mammalian Mammaliaformes",
                var= "Origination probability",
                model = "Binomial",
                average = matrix(interesting_params_occ$gamma$stat$statistics[,"Mean"], ncol=3,byrow=T)[,3],
                lw=matrix(interesting_params_occ$gamma$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,3],
                up=matrix(interesting_params_occ$gamma$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,3]),
    
    data.frame (int = bins$mid_ma[7:39][-1],
                Taxon = "Non-mammalian Mammaliaformes",
                var= "Extinction probability",
                model = "Binomial",
                average = 1-matrix(interesting_params_occ$phi$stat$statistics[,"Mean"], ncol=3,byrow=T)[,3],
                lw=1-matrix(interesting_params_occ$phi$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,3],
                up=1-matrix(interesting_params_occ$phi$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,3]),
    data.frame (int = bins$mid_ma[7:39][-1],
                
                Taxon = "Mammalia",
                var= "Origination probability",
                model = "Binomial",
                average = matrix(interesting_params_occ$gamma$stat$statistics[,"Mean"], ncol=3,byrow=T)[,1],
                lw=matrix(interesting_params_occ$gamma$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,1],
                up=matrix(interesting_params_occ$gamma$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,1]),
    
    data.frame (int = bins$mid_ma[7:39][-1],
                Taxon = "Mammalia",
                var= "Extinction probability",
                model = "Binomial",
                average = 1-matrix(interesting_params_occ$phi$stat$statistics[,"Mean"], ncol=3,byrow=T)[,1],
                lw=1-matrix(interesting_params_occ$phi$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,1],
                up=1-matrix(interesting_params_occ$phi$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,1])
    
  ) 
dat$Taxon <- factor (dat$Taxon,
                     levels = c("Non-mammaliaform cynodonts",
                                "Non-mammalian Mammaliaformes",
                                "Mammalia"))
# plot results
# binmial model
plot2 <- ggplot (data = dat ,
              
              aes(x=int,y=average, fill=var,group = var,col=var))+
    geom_point(position = position_jitter(width=0.1),size=3)+
    #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
    scale_colour_viridis_d(option="magma",begin=0.1,end=0.5)+
      
    geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                  width=0.1,size=1,position = position_jitter(width=0.1))+
    
    geom_ribbon(aes(x=int, y=average, ymax=up, ymin=lw), 
                alpha=0.2,fill="red") + 
      # other settings
    
  facet_wrap(~Taxon,scales = "fixed")+
theme_bw()+
  ylab ("Probability")+
    geom_line()+
    scale_x_reverse("Age (Ma)") +
    theme (legend.position = c(0.9,0.9))+
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,270), 
      ylim = c(0, 1),
      pos = list("b", "b"),
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) 
  

# plot of species richness
  
dat1 <- rbind (
    
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Non-mammaliaform cynodonts",
                var= "Observed richness",
                model="NA",
                average = colSums(array_genus_bin[which(clades %in% "Non-mammaliaform cynodonts"),-c(1,2)],na.rm=T),
                lw=NA,
                up=NA),
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Non-mammalian Mammaliaformes",
                var= "Observed richness",
                model="NA",
                average = colSums(array_genus_bin[which(clades %in% "Non-mammalian Mammaliaformes"),-c(1,2)],na.rm=T),
                lw=NA,
                up=NA),
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Mammalia",
                var= "Observed richness",
                model="NA",
                average = colSums(array_genus_bin[which(clades %in% "Mammalia"),-c(1,2)],na.rm=T),
                lw=NA,
                up=NA), 
    
    
    # expected values
    
    
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Non-mammaliaform cynodonts",
                var= "Expected richness",
                model="Binomial",
                average = matrix(interesting_params_occ$SRexp$stat$statistics[,"Mean"], ncol=3,byrow=T)[,2],
                lw=matrix(interesting_params_occ$SRexp$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,2],
                up=matrix(interesting_params_occ$SRexp$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,2]),
    
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Non-mammalian Mammaliaformes",
                var= "Expected richness",
                model="Binomial",
                average = matrix(interesting_params_occ$SRexp$stat$statistics[,"Mean"], ncol=3,byrow=T)[,3],
                lw=matrix(interesting_params_occ$SRexp$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,3],
                up=matrix(interesting_params_occ$SRexp$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,3]),
    data.frame (int = bins$mid_ma[7:39],
                Taxon = "Mammalia",
                var= "Expected richness",
                model="Binomial",
                average = matrix(interesting_params_occ$SRexp$stat$statistics[,"Mean"], ncol=3,byrow=T)[,1],
                lw=matrix(interesting_params_occ$SRexp$stat$quantiles[,"2.5%"], ncol=3,byrow=T)[,1],
                up=matrix(interesting_params_occ$SRexp$stat$quantiles[,"97.5%"], ncol=3,byrow=T)[,1])
    
    
    
    ) # %>%
    
dat1$Taxon <- factor (dat1$Taxon,
                     levels = c("Non-mammaliaform cynodonts",
                                "Non-mammalian Mammaliaformes",
                                "Mammalia"))

# plot the binomial results

plot3 <-  ggplot (dat1, 
                 aes(x=int,y=average, fill=var,group = var,col=var))+
 facet_wrap(~Taxon)+
  # points - average
  geom_point(position = position_jitter(width=0.1),size=3)+
 scale_colour_viridis_d(option="magma",begin=0.1,end=0.5)+
 
  # bars
  geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                  width=0.1,size=1,position = position_jitter(width=0.1))+
    # credible intervals
    geom_ribbon(data = dat1,
                  aes(x=int, y=average, ymax=up, ymin=lw), 
                alpha=0.2, fill="red")+
  
  
 
 # include the observed data
 geom_point(data = dat1 ,
            position = position_jitter(width=0.1),size=3, col ="black")+
 geom_line(data = dat1 ,col="black")+
 
  theme_bw()+
  geom_line()+
 
  theme (legend.position = c(0.15,0.9))+
    scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
    scale_colour_viridis_d(option="magma",begin=0.3,end=0.7)+
    scale_x_reverse("Age (Ma)") +
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,270), 
      ylim = c(0,550),
      pos = list("b", "b"),
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) +
 ylab ("Taxonomic diversity")
  
# detection probability
dat2 <- rbind (
    data.frame (int = bins$mid_ma[7:39],
                             Taxon = "Non-mammaliaform cynodonts",
                             var= "Detection probability in occupied sites",
                             model = "Binomial",
                             average = samples_bin[[2]]$mean$mean.p,
                             lw=round(apply (samples_bin[[2]]$sims.list$mean.p,2,quantile,0.025),2),
                             up=round(apply (samples_bin[[2]]$sims.list$mean.p,2,quantile,0.975),2),
                             nForm = formations_per_interval$formations_per_interval),
                 
                 data.frame (int = bins$mid_ma[7:39],
                             Taxon = "Non-mammalian Mammaliaformes",
                             var= "Detection probability in occupied sites",
                             model = "Binomial",
                             average = samples_bin[[3]]$mean$mean.p,
                             lw=round(apply (samples_bin[[3]]$sims.list$mean.p,2,quantile,0.025),2),
                             up=round(apply (samples_bin[[3]]$sims.list$mean.p,2,quantile,0.975),2),
                             nForm = formations_per_interval$formations_per_interval),
                 
                 data.frame (int = bins$mid_ma[7:39],
                      Taxon = "Mammalia",
                      var= "Detection probability in occupied sites",
                      model = "Binomial",
                      average = samples_bin[[1]]$mean$mean.p,
                      lw=round(apply (samples_bin[[1]]$sims.list$mean.p,2,quantile,0.025),2),
                      up=round(apply (samples_bin[[1]]$sims.list$mean.p,2,quantile,0.975),2),
                      nForm = formations_per_interval$formations_per_interval)
                 
                 
  )
                 
          
          
          
          
       # plot    
          
    
  plot4 <-  ggplot (data = dat2 ,
                     
                     aes(x=int,y=average, fill=Taxon,group = Taxon,col=Taxon))+
    # pts
    geom_point(position = position_jitter(width=0.1),size=3)+
    # error
    geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=var),
                  width=0.1,size=1,position = position_jitter(width=0.1))+
    
    # credible intervals
    geom_ribbon(aes(x=int, y=average, ymax=up, ymin=lw), 
                alpha=0.2,fill="red")+
    
    # other settings
    theme_bw()+
    geom_line()+theme (legend.position = c(0.2,0.9))+
    scale_fill_viridis_d(option="magma",begin=0.3,end=0.8)+
    scale_colour_viridis_d(option="magma",begin=0.3,end=0.8)+
    scale_x_reverse("Age (Ma)") +
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,270), 
      ylim = c(0,0.5),
      pos = list("b", "b"),
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) + 
    geom_line( aes(y=nForm/max(nForm)),col="black",size=1) +
    scale_y_continuous(
      
      # Features of the first axis
      name = "Detection probability in occupied sites",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans=~.*100, name="Number of geological formations per stage")
    )
  
  plot4
  
  # covariates per site
  
  
  ret_panel <-gridExtra::grid.arrange (plot1,
                                       plot4,
                                       plot2,
                                       plot3,
                                       nrow=4)


# save plot
pdf(here ("output","figures", "panel_params_binomialc.pdf"),width=15,height=20)
function_plots(samples_bin=binomial_output)
dev.off()


# coefficients

require(knitr)

table_coeff <- do.call(rbind, 
        lapply (seq (1,length(binomial_output)), function (i)
    
    rbind (
      
      data.frame (par =  "Origination",
                  taxon = unique(clades)[i],
                  var = "Average (logit)",
                  mean = binomial_output[[i]]$mean$intercept.gamma,
                  lci = quantile (binomial_output[[i]]$sims.list$intercept.gamma,c(0.025,0.975))[1],
                  uci =  quantile (binomial_output[[i]]$sims.list$intercept.gamma,c(0.025,0.975))[2]),
      
      data.frame (par =  "Origination",
                  taxon = unique(clades)[i],
                  var = "Precipitation",
                  mean = binomial_output[[i]]$mean$beta.gamma.prec,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.gamma.prec,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.gamma.prec,c(0.025,0.975))[2]),
      data.frame (par =  "Origination",
                  taxon = unique(clades)[i],
                  var = "Temperature",
                  mean = binomial_output[[i]]$mean$beta.gamma.temp,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.gamma.temp,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.gamma.temp,c(0.025,0.975))[2]),
      
      
      data.frame (par =  "Persistence",
                  taxon = unique(clades)[i],
                  var = "Average (logit)",
                  mean = binomial_output[[i]]$mean$intercept.phi,
                  lci = quantile (binomial_output[[i]]$sims.list$intercept.phi,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$intercept.phi,c(0.025,0.975))[2]),
      
      data.frame (par =  "Persistence",
                  taxon = unique(clades)[i],
                  var = "Precipitation",
                  mean = binomial_output[[i]]$mean$beta.phi.prec,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.phi.prec,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.phi.prec,c(0.025,0.975))[2]),
      
      data.frame (par =  "Persistence",
                  taxon = unique(clades)[i],
                  var = "Temperature",
                  mean = binomial_output[[i]]$mean$beta.phi.temp,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.phi.temp,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.phi.temp,c(0.025,0.975))[2]),
      
      
      data.frame (par =  "Detection",
                  taxon = unique(clades)[i],
                  var = "Average (logit)",
                  mean = binomial_output[[i]]$mean$intercept.p,
                  lci = quantile (binomial_output[[i]]$sims.list$intercept.p,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$intercept.p,c(0.025,0.975))[2]),
      
      data.frame (par =  "Detection",
                  taxon = unique(clades)[i],
                  var = "Time",
                  mean = binomial_output[[i]]$mean$beta.p.time,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.p.time,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.p.time,c(0.025,0.975))[2]),
      
      data.frame (par =  "Detection",
                  taxon = unique(clades)[i],
                  var = "Range size",
                  mean = binomial_output[[i]]$mean$beta.p.range,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.p.range,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.p.range,c(0.025,0.975))[2]),
      
      data.frame (par =  "Detection",
                  taxon = unique(clades)[i],
                  var = "Latitude",
                  mean = binomial_output[[i]]$mean$beta.p.lat,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.p.lat,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.p.lat,c(0.025,0.975))[2]),
      
      data.frame (par =  "Detection",
                  taxon = unique(clades)[i],
                  var = "Temperature",
                  mean = binomial_output[[i]]$mean$beta.p.temp,
                  lci = quantile (binomial_output[[i]]$sims.list$beta.p.temp,c(0.025,0.975))[1],
                  uci = quantile (binomial_output[[i]]$sims.list$beta.p.temp,c(0.025,0.975))[2])
      
    )     
  
  
)) 
  
#kable(., format = "pipe", padding = 2,align="c") 

# intercepts
data.frame (taxon = table_coeff [grep("logit",table_coeff$var),"taxon"],
            par = table_coeff [grep("logit",table_coeff$var),"par"],
            var = table_coeff [grep("logit",table_coeff$var),"var"],
            mean_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"mean"]),7),
          lci_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"lci"]),7),
          uci_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"uci"]),7))
  
  
# Show the between-S CI's in red, and the within-S CI's in black
ggplot(data= table_coeff,
       aes(x=var, y=mean, group=par)) +
facet_wrap(~taxon+par,scales="free")+
geom_errorbar(width=.1, aes(ymin=lci, ymax=uci), colour="red") +
geom_point(shape=21, size=3, fill="white") +
geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
ylab ("Coefficient value")+
xlab ("Model parameter")+

coord_flip()+

theme_bw()


# effect of covariates
# elevation

pdf(here ("output","figures", "panel_covariatesc.pdf"),width = 8, height = 8)

par(mfrow=c(3,3))

# bernoulli model
plot(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
     plogis(binomial_output$mean$intercept.gamma+
              binomial_output$mean$beta.gamma.elev*seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
     ), type = "l", ylab = "Origination probability",xlab="",ylim = c(0,1))


# binomial model results
# add posterior distribution samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.elev[i]*
                 seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01),
        lty=1)
)

# add average
lines(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
               binomial_output$mean$beta.gamma.elev*seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
      ),lwd=2,col="black",
      lty=1)


# -------------------------

# temperature 


plot(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
     plogis(binomial_output$mean$intercept.gamma+
              binomial_output$mean$beta.gamma.temp*seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
     ), type = "l", ylab = "",xlab="",ylim = c(0,0.2))

# binomial model
# add samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.temp[i]*
                 seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01),
        lty=1)
)

# average
lines(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
             binomial_output$mean$beta.gamma.temp*seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
      ),lwd=2,col="black",
      lty=1)



# ----------------------------------------------

# precipitation


plot(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
     plogis(binomial_output$mean$intercept.gamma+
              binomial_output$mean$beta.gamma.prec*seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
     ), type = "l", ylab = "",xlab="",ylim = c(0,0.2))


# binomial model
# posterior distribution samples
lapply (seq(1,length(binomial_output$sims.list$intercept.gamma)), function (i)
  
  lines(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.gamma[i]+
                 binomial_output$sims.list$beta.gamma.prec[i]*
                 seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01),lty=1)
)

# average
lines(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
      plogis(binomial_output$mean$intercept.gamma+
               binomial_output$mean$beta.gamma.prec*seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
      ),lwd=2,col="black",
      lty=1)

#-------------------------------------------------------

# persistence probability


# elevation

plot(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
     plogis(binomial_output$mean$intercept.phi+
              binomial_output$mean$beta.phi.elev*seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
     ), type = "l", ylab = "Persistence probability",xlab="Elevation (m)",ylim = c(0,1))

# binomial model
# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.elev[i]*
                   seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
        ), col = rgb(1,0,1,alpha=0.01),lty=1)
)
# average
lines(seq(min(scaled_elevation), max(scaled_elevation),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
              binomial_output$mean$beta.phi.elev*seq(min(scaled_elevation), max(scaled_elevation),by=0.01)
      ),lwd=2,col="black",
      lty=1)


# -----------------------

# temperature 


plot(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
     plogis(binomial_output$mean$intercept.phi+
              binomial_output$mean$beta.phi.temp*seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
     ), type = "l", ylab = "",xlab="Temperature (ºC)",ylim = c(0,1))

# binomial
# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.temp[i]*
                   seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
        ), col = rgb(1,0,0,alpha=0.01),lty=1)
)

# average
lines(seq(min(scaled_temperature), max(scaled_temperature),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
                 binomial_output$mean$beta.phi.temp*seq(min(scaled_temperature), max(scaled_temperature),by=0.01)
      ),lwd=2,col="black",lty=1)


# -------------------

# precipitation
plot(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
               binomial_output$mean$beta.phi.prec*seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
     ), type = "l", ylab = "",xlab="Precipitation (mm/day)",ylim = c(0,1))

# binomial model
# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.phi)), function (i)
  
  lines(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.phi[i]+
                   binomial_output$sims.list$beta.phi.prec[i]*
                   seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01),lty=1)
)
# average
lines(seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01),
      plogis(binomial_output$mean$intercept.phi+
                  binomial_output$mean$beta.phi.prec*seq(min(scaled_precipitation), max(scaled_precipitation),by=0.01)
      ),lwd=2,col="black",lty=1)



# ---------------------------------------------------------------

# plots of detection probability
# time
plot(seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01),
     plogis(
       binomial_output$mean$intercept.p+
         binomial_output$mean$beta.p.time*seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01)
     ), type = "l", ylab = "Detection probability (p)",xlab="Time",
     ylim = c(0,0.2))

# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.p)), function (i)
  
  lines(seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.p[i]+
                 binomial_output$sims.list$beta.p.time[i]*
                 seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01))
)

# average
lines(seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01),
      plogis(binomial_output$mean$intercept.p+
               binomial_output$mean$beta.p.time*seq(min(scaled_time_stage), max(scaled_time_stage),by=0.01)
      ),lwd=2,col="black")



# range size

plot(seq(min(scaled_range), max(scaled_range),by=0.01),
     plogis(
       binomial_output$mean$intercept.p+
         binomial_output$mean$beta.p.range*seq(min(scaled_range), max(scaled_range),by=0.01)
     ), type = "l", ylab = "Detection probability (p)",xlab="Range size",
     ylim = c(0,0.2))

# samples
lapply (seq(1,length(binomial_output$sims.list$intercept.p)), function (i)
  
  lines(seq(min(scaled_range), max(scaled_range),by=0.01),
        
        plogis(binomial_output$sims.list$intercept.p[i]+
                 binomial_output$sims.list$beta.p.range[i]*
                 seq(min(scaled_range), max(scaled_range),by=0.01)
        ), col = rgb(0,0,0.8,alpha=0.01))
)

# average
lines(seq(min(scaled_range), max(scaled_range),by=0.01),
      plogis(binomial_output$mean$intercept.p+
               binomial_output$mean$beta.p.range*seq(min(scaled_range), max(scaled_range),by=0.01)
      ),lwd=2,col="black")


dev.off()


# origination probability
# posterior exceedance probability
sum(binomial_output$sims.list$beta.gamma.elev>0)/length(binomial_output$sims.list$beta.gamma.elev)
sum(binomial_output$sims.list$beta.gamma.temp <0)/length(binomial_output$sims.list$beta.gamma.elev)
sum(binomial_output$sims.list$beta.gamma.prec <0)/length(binomial_output$sims.list$beta.gamma.elev)

# persistence probability
# posterior exceedance probability
sum(binomial_output$sims.list$beta.phi.elev<0)/length(binomial_output$sims.list$beta.phi.elev)
sum(binomial_output$sims.list$beta.phi.temp <0)/length(binomial_output$sims.list$beta.phi.elev)
sum(binomial_output$sims.list$beta.phi.prec >0)/length(binomial_output$sims.list$beta.phi.elev)

# detection
sum(binomial_output$sims.list$beta.p.time>0)/length(binomial_output$sims.list$beta.p.time)
sum(binomial_output$sims.list$beta.p.range>0)/length(binomial_output$sims.list$beta.p.range)

# summary statistics
plogis(binomial_output$mean$intercept.gamma)
plogis(quantile (binomial_output$sims.list$intercept.gamma,c(0.025,0.975)))
plogis(binomial_output$mean$intercept.phi)
plogis(quantile (binomial_output$sims.list$intercept.phi,c(0.025,0.975)))
# detection
plogis(binomial_output$mean$intercept.p)
plogis(quantile (binomial_output$sims.list$intercept.p,c(0.025,0.975)))
# n genus
binomial_output$mean$SRexp
quantile (binomial_output$sims.list$SRexp[,33],c(0.025,0.975))


# --------------------------------------------------------------

# change
# pch
plot(bernoulli_output$mean$propcH,type="b",ylim=c(-1,4))
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

# plot parameters
dat <-   rbind (
  data.frame (int = bins$mid_ma[7:39],
              Parameter= "Relative extinction rates (RER)",
              Average = binomial_output$mean$RER,
              lw=round(apply (binomial_output$sims.list$RER,2,quantile,0.025,na.rm=T),2),
              up=round(apply (binomial_output$sims.list$RER,2,quantile,0.975,na.rm=T),2)),
  
  data.frame (int = bins$mid_ma[7:39],
              Parameter= "Proportional change (PCh)",
              Average = binomial_output$mean$propcH,
              lw=round(apply (binomial_output$sims.list$propcH,2,quantile,0.025,na.rm=T),2),
              up=round(apply (binomial_output$sims.list$propcH,2,quantile,0.975,na.rm=T),2))
  
) 

# plot results
# binmial model
png (here ("output", "figures","change_indicators.png"),
     height=15,width=17,units ="cm",res=300)


ggplot (data = dat,
                 
                 aes(x=int,y=Average, fill=Parameter,group = Parameter,col=Parameter))+
  geom_point(position = position_jitter(width=0.1),size=3)+
  #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  scale_colour_viridis_d(option="magma",begin=0.1,end=0.5)+
  
  geom_errorbar(aes(x = int, ymin = lw, ymax = up,col=Parameter),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=int, y=Average, ymax=up, ymin=lw), 
              alpha=0.2,fill="red") + 
  # other settings
  theme_bw()+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.75,0.9))+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(-2,18),
    pos = list("b", "b"),
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) 

dev.off()

# --------------------------------------------------------------

# convergence (RHat!)

# fit statistics
# Rhat

sel_rows <- c(grep("intercept",rownames(samples_paleo_cynodontia$summary)),
  grep("beta",rownames(samples_paleo_cynodontia$summary)))

data.frame(Rhat = samples_paleo_cynodontia$summary [sel_rows, "Rhat"],
           par = rownames(samples_paleo_cynodontia$summary [sel_rows, ])) %>% 
  ggplot(aes(x=Rhat, y=par)) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Convergence in the global-scale model")

ggsave (here ("output", "figures","convergence_global_params.png"))

# dynamic params

sel_rows <- c(grep("intercept",rownames(samples_paleo_cynodontia$summary)),
              grep("beta",rownames(samples_paleo_cynodontia$summary)),
              grep("R0",rownames(samples_paleo_cynodontia$summary)),
              grep("mean",rownames(samples_paleo_cynodontia$summary)),
              grep("deviance",rownames(samples_paleo_cynodontia$summary)))


data.frame(Rhat = samples_paleo_cynodontia$summary [-sel_rows, "Rhat"],
           comb = rownames(samples_paleo_cynodontia$summary [-sel_rows, ]),
           par = substr(rownames(samples_paleo_cynodontia$summary [-sel_rows, ]),1,2)) %>% 
  ggplot(aes(x=Rhat, y=comb)) +
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  theme (axis.text = element_text(size=2))+
  xlab ("Rhat value")+
  facet_wrap(~par,nrow=3,scales = "free")+
  ggtitle ("Convergence in the global-scale model")

ggsave (here ("output", "figures","convergence_global.png"))


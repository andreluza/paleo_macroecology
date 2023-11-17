

# -----------------------------------------------------------------------


#     Interpretation: Global-scale analysis (SENSITIVITY ANALYSES - SUPP INFORMATION)


# -----------------------------------------------------------------------



# load packages
source("R/packages.R")
source("R/functions.R")


# define the standard of colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen")

# ------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "interesting_params_global_occ_no_cov.RData"))

# change
load (here ("output",
            "interesting_params_global_change_no_cov.RData"))

# load basic data
load(here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))

# binning time intervals (from Paleoverse)
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)


# table of summary statistics of covariates

apply (time_covariates[,-5],2,mean,na.rm=T)
apply (time_covariates[,-5],2,sd,na.rm=T)


# apply colors for the strip

bins <- cbind (bins,cols_strip = rep (c("yes", "no"),20)[-1])
cols_strip <- c("yes" = "gray", "no" = "white")


# plot of species richness

# find the period of first and last detection
function_stages <- function (x) {
  
  # the first stage with the taxon,    
  sel_cols <- if  (max(which(colSums (x)>0)) == 33) {
    
    a <- seq(
      ifelse (min(which(colSums (x)>0)) ==1,
              min(which(colSums (x)>0)),
              min(which(colSums (x)>0))-1),
      
      
      max(which(colSums (x)>0)),
      
      1)
    
  } else {
    
    a <- seq(ifelse (min(which(colSums (x)>0)) ==1,
                     min(which(colSums (x)>0)),
                     min(which(colSums (x)>0))-1),
             max(which(colSums (x)>0))+1,
             1)
  }
  
  return(a)
  
}

stages_cyn <- (function_stages (x = array_genus_bin[which(clades %in% "Non-mammaliaform cynodonts"),]))
stages_mammaf <- (function_stages (x = array_genus_bin[which(clades %in% "Non-mammalian Mammaliaformes"),]))
stages_mamm <- (function_stages (x = array_genus_bin[which(clades %in% "Mammalia"),]))


# plot parameters
dat <-   rbind (
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,][-1,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_occ_cynodonts$gamma$stat$statistics[,"Mean"],
              lw= interesting_params_occ_cynodonts$gamma$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_cynodonts$gamma$stat$quantiles[,"97.5%"]),
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,][-1,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-interesting_params_occ_cynodonts$phi$stat$statistics[,"Mean"],
              lw=1-interesting_params_occ_cynodonts$phi$stat$quantiles[,"2.5%"],
              up=1-interesting_params_occ_cynodonts$phi$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,][-1,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_occ_mammaliaformes$gamma$stat$statistics[,"Mean"],
              lw=interesting_params_occ_mammaliaformes$gamma$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_mammaliaformes$gamma$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,][-1,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-interesting_params_occ_mammaliaformes$phi$stat$statistics[,"Mean"],
              lw=1-interesting_params_occ_mammaliaformes$phi$stat$quantiles[,"2.5%"],
              up=1-interesting_params_occ_mammaliaformes$phi$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,][-1,],
              Taxon = "Mammalia",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_occ_mammalia$gamma$stat$statistics[,"Mean"],
              lw=interesting_params_occ_mammalia$gamma$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_mammalia$gamma$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,][-1,],
              Taxon = "Mammalia",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-interesting_params_occ_mammalia$phi$stat$statistics[,"Mean"],
              lw=1-interesting_params_occ_mammalia$phi$stat$quantiles[,"2.5%"],
              up=1-interesting_params_occ_mammalia$phi$stat$quantiles[,"97.5%"])
)


dat$Taxon <- factor (dat$Taxon,
                     levels = c("Non-mammaliaform cynodonts",
                                "Non-mammalian Mammaliaformes",
                                "Mammalia"))
# plot results
# binmial model
plot2 <- ggplot (data = dat ,
                 
                 aes(x=mid_ma,
                     y=average, 
                     fill=var,
                     group = var,
                     col=Taxon))+
  geom_point(aes (shape = var),
             position = position_jitter(width=0.1),size=3)+
  #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  #scale_colour_viridis_d(option="magma",begin=0.1,end=1)+
  #scale_fill_viridis_d(option="magma",begin=0.1,end=1)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  geom_errorbar(aes(x = mid_ma, ymin = lw, ymax = up,col=var),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=mid_ma, y=average, ymax=up, ymin=lw,fill=Taxon), 
              alpha=0.2) + 
  # other settings
  
  facet_wrap(~Taxon,scales = "fixed")+
  theme_bw()+
  ylab ("Probability")+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.8,0.9))+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(0, 1),
    pos = list("b", "b"),
    rot=90,
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) +
  
  geom_rect(aes(xmin = max_ma, 
                xmax = min_ma, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = cols_strip), 
            
            col=NA,
            alpha = 0.2)+
  scale_fill_manual(values = cols_strip)



plot2



# plot of species richness


dat1 <- rbind (
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammaliaform cynodonts"),],na.rm=T)[stages_cyn],
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammalian Mammaliaformes"),],na.rm=T)[stages_mammaf],
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,],
              Taxon = "Mammalia",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Mammalia"),],na.rm=T)[stages_mamm],
              lw=NA,
              up=NA), 
  
  
  # expected values
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_cynodonts$SRexp$stat$statistics[,"Mean"],
              lw=interesting_params_occ_cynodonts$SRexp$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_cynodonts$SRexp$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_mammaliaformes$SRexp$stat$statistics[,"Mean"],
              lw=interesting_params_occ_mammaliaformes$SRexp$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_mammaliaformes$SRexp$stat$quantiles[,"97.5%"]),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,],
              Taxon = "Mammalia",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_mammalia$SRexp$stat$statistics[,"Mean"],
              lw= interesting_params_occ_mammalia$SRexp$stat$quantiles[,"2.5%"],
              up= interesting_params_occ_mammalia$SRexp$stat$quantiles[,"97.5%"])
  
  
  
)# %>%

dat1$Taxon <- factor (dat1$Taxon,
                      levels = c("Non-mammaliaform cynodonts",
                                 "Non-mammalian Mammaliaformes",
                                 "Mammalia"))

# explore
dat1 %>%
  filter (var == "Expected richness" &
            Taxon == "Non-mammaliaform cynodonts" &
            lw >= 50)
dat1 %>%
  filter (var == "Expected richness" &
            Taxon == "Mammalia" &
            lw >= 100)

dat1 %>%
  filter (var == "Expected richness" &
            Taxon == "Non-mammaliaform cynodonts")

dat1 %>%
  filter (var == "Expected richness" &
            Taxon == "Non-mammalian Mammaliaformes")


dat1 %>%
  filter (var == "Expected richness" &
            Taxon == "Mammalia")


# plot the binomial results

plot3 <-  ggplot (dat1, 
                  aes(x=mid_ma,y=average, 
                      fill=var,
                      group = var,
                      col=Taxon))+
  facet_wrap(~Taxon)+
  # points - average
  geom_point(aes (shape = var),position = position_jitter(width=0.1),size=3)+
  #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  #scale_colour_viridis_d(option="magma",begin=0.1,end=1)+
  #scale_fill_viridis_d(option="magma",begin=0.1,end=1)+
  
  # bars
  geom_errorbar(aes(x = mid_ma, ymin = lw, ymax = up,col=Taxon),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  # credible intervals
  geom_ribbon(data = dat1,
              aes(x=mid_ma, y=average, ymax=up, ymin=lw,fill=Taxon), 
              alpha=0.2) +
  
  
  
  # include the observed data
  geom_point(data = dat1 , aes (shape=var),
             position = position_jitter(width=0.1),
             size=3)+
  geom_line(data = dat1)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  theme_bw()+
  geom_line()+
  
  theme (legend.position = c(0.15,0.9))+
  #scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  #scale_colour_viridis_d(option="magma",begin=0.3,end=0.7)+
  scale_x_reverse("Age (Ma)") +
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    #ylim = c(0,200),
    pos = list("b", "b"),
    rot=90,
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) +
  ylab ("Taxonomic diversity")+
  
  geom_rect(aes(xmin = max_ma, 
                xmax = min_ma, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = cols_strip),
            col=NA,
            alpha = 0.2)+
  scale_fill_manual(values = cols_strip)

plot3

# 
dat1 %>% filter  (Taxon == "Non-mammaliaform cynodonts" &
                    var == "Expected richness") %>%
  mutate (fact_SR = average)
# 
dat1 %>% filter  (Taxon == "Non-mammalian Mammaliaformes" &
                    var == "Expected richness") %>%
  mutate (fact_SR = average)
# 
dat1 %>% filter  (Taxon == "Mammalia" &
                    var == "Expected richness") %>%
  mutate (fact_SR = average)



# total N of genera

sum(apply(interesting_params_change_cynodonts$z$PP_mat$AA[[1]],1,max))
sum(apply(interesting_params_change_mammaliaformes$z$PP_mat$AA[[1]],1,max))
sum(apply(interesting_params_change_mammalia$z$PP_mat$AA[[1]],1,max))

# function  to calculate change in TD
TD_fc <- function (output) {
  
  change_TD <- lapply (seq(1,length(output$z$PP_mat)), function (k)
    
    lapply (output$z$PP_mat[[k]], function (i) {
      
      sum(apply(i,1,max))
      
    }))
  # melt
  change_TD<-unlist(change_TD)
  return(change_TD)
  
}

total_TD_cyn <- TD_fc (interesting_params_change_cynodonts)
total_TD_mammaliaforms <- TD_fc (interesting_params_change_mammaliaformes)
total_TD_mammals <- TD_fc (interesting_params_change_mammalia)

rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts",
              TD = total_TD_cyn),
  data.frame (Taxon = "Non-mammalian Mammaliaformes",
              TD = total_TD_mammaliaforms),
  data.frame (Taxon = "Mammalia",
              TD = total_TD_mammals)) %>%
  ggplot (aes(x=TD,group=Taxon,col= Taxon,fill=Taxon))+
  geom_density(alpha=0.5)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  xlim(c(50,500))+
  geom_vline(aes (xintercept=mean(total_TD_cyn)),col= "red")+
  geom_vline(aes (xintercept=mean(total_TD_mammaliaforms)),col= "blue")+
  geom_vline(aes (xintercept=mean(total_TD_mammals)),col= "green4")+
  theme(legend.position = "top")+
  theme_classic()



# total estimate
quantile(total_TD_cyn,c(0.025,0.5,0.975))
quantile(total_TD_mammaliaforms,c(0.025,0.5,0.975))
quantile(total_TD_mammals,c(0.025,0.5,0.975))

mean(total_TD_cyn)
mean(total_TD_mammaliaforms)
mean(total_TD_mammals)

# ==============================

# convergence dynamic params
dat_rhat <- rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              rhat = unlist(sapply (interesting_params_occ_cynodonts, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_occ_cynodonts, "[[", "rhat")))),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              rhat = unlist(sapply (interesting_params_occ_mammaliaformes, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_occ_mammaliaformes, "[[", "rhat")))),
  data.frame (Taxon = "Mammalia",
              rhat=unlist(sapply (interesting_params_occ_mammalia, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_occ_mammalia, "[[", "rhat")))))


# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))
# adjust names
dat_rhat <- rbind (dat_rhat [grep ("phi",dat_rhat$par),],
                   dat_rhat [grep ("gamma",dat_rhat$par),]) %>%
  filter (is.na(rhat) != T)
dat_rhat$par <- factor (dat_rhat$par,
                        levels = rev(unique(dat_rhat$par)))
# dynamic params
#plot
ggplot(data=  dat_rhat,
       aes(x=rhat, y=par)) +
  facet_wrap(~Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col= "red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Convergence of origination (gamma) and persistence (phi)\nin the global-scale model")

ggsave (here ("output", "figures","convergence_dyn_global_params.png"),width =10,height=10)



# ------------------------------------------------------


# function  to calculate change in TD
change_fc <- function (output) {
  
  change_TD <- lapply (seq(1,length(output$z$PP_mat)), function (k)
    
    do.call(rbind, lapply (output$z$PP_mat[[k]], function (i) {
      
      CH <- rep (NA,ncol (i))
      for (t in 2:ncol (i)) {
        
        CH[t] <- (sum(i[,t]) - sum(i[,t-1]))#/sum(i[,t-1])
        
        
      }
      
      CH
      
    })))
  # melt
  change_TD<-do.call(rbind,change_TD)
  return(change_TD)
  
}

change_TD_cyn <- change_fc (interesting_params_change_cynodonts)
change_TD_mammaliaforms <- change_fc (interesting_params_change_mammaliaformes)
change_TD_mammals <- change_fc (interesting_params_change_mammalia)

# bind data to plot
dat_change <- rbind (
  
  data.frame ('average' = (apply (change_TD_cyn,2,mean,na.rm=T)),
              'lci' = t(apply (change_TD_cyn,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (change_TD_cyn,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammaliaform cynodonts",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,]),
  data.frame ('average' = (apply (change_TD_mammaliaforms,2,mean,na.rm=T)),
              'lci' = t(apply (change_TD_mammaliaforms,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (change_TD_mammaliaforms,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammalian Mammaliaformes",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,])
  ,
  data.frame ('average' = (apply (change_TD_mammals,2,mean,na.rm=T)),
              'lci' = t(apply (change_TD_mammals,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (change_TD_mammals,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Mammalia",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,])
)


# order of groups
dat_change$Taxon <- factor (dat_change$Taxon,
                            levels = c("Non-mammaliaform cynodonts",
                                       "Non-mammalian Mammaliaformes",
                                       "Mammalia"))
# replacing inf by NA (all will be zero)
dat_change[is.infinite(dat_change$average),c("average","lci","uci")] <- NA

# plot results

plot_change <- ggplot (data = dat_change,
                       
                       aes(x=mid_ma,y=average, fill=Taxon,group = Taxon,col=Taxon))+
  geom_point(position = position_jitter(width=0.1),size=3)+
  #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  geom_errorbar(aes(x = mid_ma, ymin = lci, ymax = uci,col=Taxon),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=mid_ma, y=average, ymax=uci, ymin=lci,col=Taxon,fill=Taxon), 
              alpha=0.2) + 
  ylab ("Change in taxonomic diversity (ΔTDt)")+
  # other settings
  theme_bw()+
  facet_wrap(~Taxon,scales="free")+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = c(0.5,0.82))+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    #ylim = c(-30,85),
    pos = list("b", "b"),
    rot=90,
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) +
  
  geom_rect(aes(xmin = max_ma, 
                xmax = min_ma, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = cols_strip),
            col=NA,
            alpha = 0.2)+
  scale_fill_manual(values = cols_strip)

plot_change

# divesrification

# function  to calculate change in TD
div_fc <- function (output) {
  
  # gamma
  res_gamma<-output[which(names(output) == "gamma")]
  res_gamma<-rbind (
    do.call(rbind,res_gamma$gamma$PP_mat[[1]]),
    do.call(rbind,res_gamma$gamma$PP_mat[[2]]),
    do.call(rbind,res_gamma$gamma$PP_mat[[3]]))
  
  # phi
  res_phi<-output[which(names(output) == "phi")]
  res_phi<-rbind (
    do.call(rbind,res_phi$phi$PP_mat[[1]]),
    do.call(rbind,res_phi$phi$PP_mat[[2]]),
    do.call(rbind,res_phi$phi$PP_mat[[3]]))
  
  
  # calculate diversification
  div_rate <- lapply (seq(1,ncol (res_phi)), function (t) {
    
    RER <- (1-res_phi[,t])/res_gamma[,t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
    R0 <- res_gamma[,t]-(1-res_phi[,t]) ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
    output <- list(RER=RER,
                   R0=R0)
    output
    
  })
  
  res <- list (RER = sapply (div_rate, "[[", "RER"),
               R0 = sapply (div_rate, "[[", "R0"))
  
  return(res)
  
}

# run function
div_cyn <- div_fc (interesting_params_occ_cynodonts)
div_mammaliaforms <- div_fc (interesting_params_occ_mammaliaformes)
div_mammals <- div_fc (interesting_params_occ_mammalia)


# bind data to plot
dat_div <- rbind (
  
  data.frame ('average' = (apply (div_cyn$R0,2,mean,na.rm=T)),
              'lci' = t(apply (div_cyn$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_cyn$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammaliaform cynodonts",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,][-1,]),
  data.frame ('average' = (apply (div_mammaliaforms$R0,2,mean,na.rm=T)),
              'lci' = t(apply (div_mammaliaforms$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_mammaliaforms$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammalian Mammaliaformes",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,][-1,])
  ,
  data.frame ('average' = (apply (div_mammals$R0,2,mean,na.rm=T)),
              'lci' = t(apply (div_mammals$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_mammals$R0,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Mammalia",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,][-1,])
)

# order of groups
dat_div$Taxon <- factor (dat_div$Taxon,
                         levels = c("Non-mammaliaform cynodonts",
                                    "Non-mammalian Mammaliaformes",
                                    "Mammalia"))
# plot results
plot_div <- ggplot (data = dat_div,
                    
                    aes(x=mid_ma,y=average, fill=Taxon,group = Taxon,col=Taxon))+
  geom_point(position = position_jitter(width=0.1),size=3)+
  #  scale_fill_viridis_d(option="magma",begin=0.3,end=0.7)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  geom_errorbar(aes(x = mid_ma, ymin = lci, ymax = uci,col=Taxon),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=mid_ma, y=average, ymax=uci, ymin=lci,col=Taxon,fill=Taxon), 
              alpha=0.2) + 
  ylab ("Net diversification rate (NDRt)")+
  # other settings
  theme_bw()+
  facet_wrap(~Taxon,scales="free")+
  geom_line()+
  scale_x_reverse("Age (Ma)") +
  theme (legend.position = "none")+
  coord_geo(
    dat = list("stages", "periods"), 
    xlim = c( 66,270), 
    ylim = c(-1,1),
    pos = list("b", "b"),
    rot=90,
    size = list(2, 4),
    abbrv = list(TRUE, T)
  ) +
  
  geom_rect(aes(xmin = max_ma, 
                xmax = min_ma, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = cols_strip),
            col=NA,
            alpha = 0.2)+
  scale_fill_manual(values = cols_strip)

plot_div


# relationship between diversification and change
rel_div_change <-dat_div %>%
  #filter (Taxon == "Mammalia") %>%
  bind_cols(dat_change %>%
              
              #filter (Taxon == "Mammalia") %>%
              select (interval_name,average) %>%
              dplyr::rename("int_name"="interval_name",
                            "change" = "average") %>%
              filter (is.na(change)!=T) 
              
              
            )

# plot the relationship
  
ggplot (rel_div_change,aes (x=average, y=change, group=Taxon,col= Taxon))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols) +
  xlab ("Net diversification rate (NDRt)")+
  ylab ("Change in taxonomic diversity (ΔTDt)")
  



# insert a lag
plot(rel_div_change$average[-length(rel_div_change$average)],rel_div_change$change[-1])




# save plot
# params
pdf(here ("output","figures", "panel_diversity_SUpp_FIgS2-3.3.pdf"),width=15,height=17)
ret_panel <-gridExtra::grid.arrange (plot3,
                                     plot_change,
                                     plot_div,
                                     nrow=3)
dev.off()

rm(list=ls())

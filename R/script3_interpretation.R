
require(jagsUI)
require(ggplot2)
require(here)
require(reshape)

# load
load (here ("output","samples_paleo.RData"))
load (here ("output","table_data_array.RData"))

# naive occupancy
naive_occupancy <- lapply (unique(table_data_long_df$taxon), function (i) 
  
      reshape::cast(data = table_data_long_df[which(table_data_long_df$taxon == i),],
                      formula = form~int,
                      margins = "taxon",
                      value = "det",
                      fun.aggregate = max,
                      fill=NA))

apply (naive_occupancy[[1]], 2, max)

# observed number of genus per formation and interval 
par(mar=rep(5,4))
hist(as.numeric(as.matrix(test[,-1])),
     xlab="Number of genera",
     main="")


# list of intervals of interest
triassic <- c("Olenekian", "Anisian", "Ladinian", "Carnian", "Norian", "Rhaetian") # Triassic
# triassic %in% coll_occ_taxa$early_interval.x
jurassic <- c("Hettangian", "Sinemurian", "Pliensbachian", "Toarcian", "Aalenian","Bajocian", "Bathonian", "Callovian", "Oxfordian", "Kimmeridgian", "Tithonian")

# explore results
df_per <- data.frame (change = samples_paleo$mean$propcH[-1],
                      CI= t(apply (samples_paleo$sims.list$propcH[,-1],2, quantile,probs=c(0.05,0.95))),
                      interval=c(triassic,jurassic)[-1])
df_per$interval <- factor (df_per$interval,
                           levels = c(triassic,jurassic)[-1])
ggplot (df_per, aes (x=interval, 
                     y=change,
                     group= 1)) +
  theme_bw()+
  geom_line(size=2)+
  theme (axis.text.x = element_text(angle=90)) + 
  ylab ("Proportion of Change in Genus Number\nrelative to t-1") + 
  xlab ("Geological Interval") + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  geom_vline(xintercept=5.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.))+ 
  annotate (geom="text", x=3,y=0.3,label ="Triassic")+ 
  annotate (geom="text", x=12,y=0.3,label ="Jurassic")


# genus number
df_genus <- data.frame (ngenus = samples_paleo$mean$Ngen,
                        CI=t(apply (samples_paleo$sims.list$Ngen,2, quantile,probs=c(0.05,0.95))),
                        interval = c(triassic,jurassic))
df_genus$interval <- factor (df_genus$interval,
                             levels = c(triassic,jurassic))

ggplot (df_genus, aes (x=interval,y=ngenus,group= 1)) +
  theme_bw()+
  geom_line(size=2)+
  theme (axis.text.x = element_text(angle=90)) + 
  ylab ("Expected Number of Genera") + 
  xlab ("Geological Interval") + 
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.))+
  geom_vline(xintercept=6.5, 
             linetype="dashed", 
             color = "red", 
             size=1) + 
  geom_hline(yintercept=mean(samples_paleo$mean$Ngen), 
             linetype="dashed", 
             color = "red", 
             size=1)+
  annotate (geom="text", x=3,y=460,label ="Triassic")+ 
  annotate (geom="text", x=12,y=460,label ="Jurassic")

# extinction
df_pars <- data.frame (epslon=rowMeans(samples_paleo$mean$epslon[-1,],na.rm=T),
                       gamma=rowMeans(samples_paleo$mean$gamma[-1,],na.rm=T),
                       CI.gamma = t(apply (samples_paleo$sims.list$gamma[,-1,],2, quantile,probs=c(0.05,0.95))),
                       CI.epslon = t(apply (samples_paleo$sims.list$epslon[,-1,],2, quantile,probs=c(0.05,0.95))),
                       interval=c(triassic,jurassic)[-1])
df_pars$interval <- factor (df_pars$interval,
                            levels = c(triassic,jurassic)[-1])

ggplot(melt (df_pars, id.vars=c("interval","CI.gamma.5.","CI.gamma.95.","CI.epslon.5.","CI.epslon.95.")),
       aes (x=interval,y=value,
            color=variable,group=variable))+
  geom_line(size=2)+
  scale_colour_viridis_d(end=0.8)+
  theme_bw()+

  theme (axis.text.x = element_text(angle=90)) + 
  ylab ("Expected Number of Genus") + 
  xlab ("Geological Interval") + 
  #geom_linerange(aes(ymin=CI.gamma.5., ymax=CI.gamma.95.))+
  #geom_linerange(aes(ymin=CI.epslon.5., ymax=CI.epslon.95.))+
  geom_vline(xintercept=5.5, 
             linetype="dashed", 
             color = "red", 
             size=1) + 
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1) + 
  annotate (geom="text", x=3,y=0.52,label ="Triassic")+ 
  annotate (geom="text", x=12,y=0.52,label ="Jurassic")+
  annotate (geom="text", x=0.6,y=0.495,size=3,
            label =triassic[1], angle=90) 
  
  

# genus
dat_genus <- data.frame (genus = dimnames(table_data_array)[[3]],
                         psi.eq = samples_paleo$mean$psi.eq,
                         CI = t(apply (samples_paleo$sims.list$psi.eq,2, quantile,probs=c(0.05,0.95))))

ggplot (dat_genus, aes (x=reorder(genus, psi.eq),
                        y=psi.eq,group=1))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  theme(axis.text.x = element_text(angle=90,size=8))+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.1)


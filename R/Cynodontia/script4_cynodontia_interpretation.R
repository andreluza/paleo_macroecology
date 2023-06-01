# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")

# ------------------------------------------------------------
# load
load (here ("output","samples_paleo_genus.RData"))
load (here ("output","samples_paleo_family.RData"))
load (here ("output","table_data_array.RData"))
load (here ("output","table_data_array_family.RData"))
load (here ("output","table_naive.RData"))
load (here ("output","table_naive_family.RData"))

# ------------------------------------------------------------
# closest function
closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}

# ------------------------------------------------------------

# information about sampling effort
# number of time intervals per taxon
intervals_per_taxon<- lapply (table_naive, colSums,na.rm=T)
intervals_per_taxon <- (colSums(do.call(rbind, intervals_per_taxon)>0))

# family
intervals_per_taxon_family<- lapply (table_naive_family, colSums,na.rm=T)
intervals_per_taxon_family <- (colSums(do.call(rbind, intervals_per_taxon_family)>0))


# number of formations per interval

png (here ("output","effort.png"),res=300,unit="cm",height=15,width=18)
hist (colSums(table_data_basis>0,na.rm=T),
      xlab = "Number of formations",
      ylab = "Frequency (# of intervals)",
      main = "Number of formations per time interval")

# plot average
text (x=200,y=3,
      labels = paste ("Average+-SD=", 
                      round(mean(colSums(table_data_basis>0,na.rm=T)),2),
                      "+-",
                      round(sd(colSums(table_data_basis>0,na.rm=T)),2),
                      sep=" "))
# plot range
text (x=200,y=2.8,
      labels = paste ("Range=", 
                      range(colSums(table_data_basis>0,na.rm=T))[1],
                      ",",
                      range(colSums(table_data_basis>0,na.rm=T))[2],
                      sep=" "))


dev.off()

# ------------------------------------------------------------
# formations that cover more intervals
hist(rowSums(table_data_basis > 0,na.rm=T) [order(rowSums(table_data_basis > 0,na.rm=T))])

# number of formations per interval
png (here ("output","effort_formations.png"),res=300,unit="cm",height=15,width=18)
hist (rowSums(table_data_basis>0,na.rm=T),
      xlab = "Number of intervals",
      ylab = "Frequency (# of formations)",
      main = "Number of time intervals per formation")
# plot average
text (x=4,y=1200,
      labels = paste ("Average+-SD=", 
                      round(mean(rowSums(table_data_basis>0,na.rm=T)),2),
                      "+-",
                      round(sd(rowSums(table_data_basis>0,na.rm=T)),2),
                      sep=" "))
# plot range
text (x=4,y=1120,
      labels = paste ("Range=", 
                      range(rowSums(table_data_basis>0,na.rm=T))[1],
                      ",",
                      range(rowSums(table_data_basis>0,na.rm=T))[2],
                      sep=" "))


dev.off()

# ------------------------------------------------------------
# formations per taxon
taxon_per_formation <- lapply (table_naive, function (i) {
  
  rowSums (i>0,na.rm=T)
  
})
#  melt  
taxon_per_formation <- (do.call(cbind, taxon_per_formation))
# family
taxon_per_formation_family <- lapply (table_naive_family, function (i) {
  
  rowSums (i>0,na.rm=T)
  
})
#  melt  
taxon_per_formation_family <- (do.call(cbind, taxon_per_formation_family))


# number of formations per interval
png (here ("output","taxon_formations.png"),res=300,unit="cm",height=15,width=18)
par (mfrow=c(1,2))
# genera
hist (rowSums(taxon_per_formation>0,na.rm=T),
      xlab = "Number of genera",
      ylab = "Frequency (# of formations)",
      main = "Number of genera per formation")
# plot average
text (x=100,y=1200,
      labels = paste ("Average+-SD=\n", 
                      round(mean(rowSums(taxon_per_formation>0,na.rm=T)),2),
                      "+-",
                      round(sd(rowSums(taxon_per_formation>0,na.rm=T)),2),
                      sep=" "))
# plot range
text (x=100,y=1000,
      labels = paste ("Range=", 
                      range(rowSums(taxon_per_formation>0,na.rm=T))[1],
                      ",",
                      range(rowSums(taxon_per_formation>0,na.rm=T))[2],
                      sep=" "))

# family
hist (rowSums(taxon_per_formation_family>0,na.rm=T),
      xlab = "Number of families",
      ylab = "",
      main = "Number of families per formation")
# plot average
text (x=80,y=1200,
      labels = paste ("Average+-SD=\n", 
                      round(mean(rowSums(taxon_per_formation_family>0,na.rm=T)),2),
                      "+-",
                      round(sd(rowSums(taxon_per_formation_family>0,na.rm=T)),2),
                      sep=" "))
# plot range
text (x=80,y=1000,
      labels = paste ("Range=", 
                      range(rowSums(taxon_per_formation_family>0,na.rm=T))[1],
                      ",",
                      range(rowSums(taxon_per_formation_family>0,na.rm=T))[2],
                      sep=" "))



dev.off()

# ==================================================================

# list of intervals of interest
triassic <- c("Olenekian", "Anisian", "Ladinian", "Carnian", "Norian", "Rhaetian") # Triassic
# triassic %in% coll_occ_genera$early_interval.x
jurassic <- c("Hettangian", "Sinemurian", "Pliensbachian", "Toarcian", "Aalenian","Bajocian", "Bathonian", "Callovian", "Oxfordian", "Kimmeridgian", "Tithonian")


# test
(sum(samples_paleo_genus$mean$muZ[3,])-sum(samples_paleo_genus$mean$muZ[2,]))/sum(samples_paleo_genus$mean$muZ[2,])


# explore results
df_per <- data.frame (change_f = -1*samples_paleo_family$mean$propcH[-1],
                      CI_f= -1*t(apply (samples_paleo_family$sims.list$propcH[,-1],2, quantile,probs=c(0.05,0.95))),
                      change = -1*samples_paleo_genus$mean$propcH[-1],
                      CI= -1*t(apply (samples_paleo_genus$sims.list$propcH[,-1],2, quantile,probs=c(0.05,0.95))),
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
  annotate (geom="text", x=12,y=0.3,label ="Jurassic")+ 
  geom_line (aes (x=interval, 
                  y=change_f,
                  group= 1),colour="green3",size=1,
             )+
  geom_pointrange(aes(ymin=CI_f.5., ymax=CI_f.95.))

ggsave (here ("output", "prop_change.png"),
        dpi=300,
        units="cm",
        height=15,
        width=20)


# genus number
df_genus <- data.frame (ngenus = samples_paleo_genus$mean$Ngen,
                        ngenus_naive = intervals_per_taxon,
                        CI=t(apply (samples_paleo_genus$sims.list$Ngen,2, quantile,probs=c(0.05,0.95))),
                        nfam = samples_paleo_family$mean$Ngen,
                        nfam_naive = intervals_per_taxon_family,
                        CI_f=t(apply (samples_paleo_family$sims.list$Ngen,2, quantile,probs=c(0.05,0.95))),
                        interval = c(triassic,jurassic))
df_genus$interval <- factor (df_genus$interval,
                             levels = c(triassic,jurassic))

# plot
p1<-ggplot (df_genus, aes (x=interval,y=ngenus,group= 1)) +
  theme_bw()+
  geom_line(size=2)+
  theme (axis.text.x = element_text(angle=90)) + 
  ylab ("Expected Number of Genera and Families") + 
  xlab ("Geological Interval") + 
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.))+
  geom_vline(xintercept=6.5, 
             linetype="dashed", 
             color = "red", 
             size=1) + 
  geom_hline(yintercept=mean(samples_paleo_genus$mean$Ngen), 
             linetype="solid", 
             color = "black", 
             size=0.5)+
  annotate (geom="text", x=3,y=650,label ="Triassic")+ 
  annotate (geom="text", x=12,y=650,label ="Jurassic") + 
  geom_line(data = df_genus, 
            aes (x=interval,y=ngenus_naive,group= 1),
            alpha=0.35)+
  geom_hline(yintercept=mean(intervals_per_taxon), 
             linetype="solid", 
             color = "black", 
             size=0.5,
             alpha = 0.35)
# family
# estimated
p1 + geom_line(data=df_genus, aes (x=interval,y=nfam,group= 1),
               size=2,linetype="solid",colour = "green3")+
  theme (axis.text.x = element_text(angle=90)) + 
  geom_pointrange(aes(ymin=CI_f.5., ymax=CI_f.95.))+
  geom_hline(yintercept=mean(samples_paleo_family$mean$Ngen), 
             linetype="solid", 
             color = "green", 
             size=1)+
  # naive
  geom_line(data = df_genus, 
            aes (x=interval,y=nfam_naive,group= 1),
            alpha=0.7,colour="green3")+
  geom_hline(yintercept=mean(intervals_per_taxon_family), 
             linetype="solid", 
             color = "green3", 
             size=0.5,
             alpha = 0.35)
  

# save
ggsave (here ("output", "exp_Ngenera_fam.png"),
        dpi=300,
        units="cm",
        height=15,
        width=20)


# ------------------------------------------------------------
# detection probability
detection_df <- data.frame (genus = genus_data[sel_genus],
                            p = samples_paleo_genus$mean$p,
                            CI = t(apply (samples_paleo_genus$sims.list$p,2, quantile,probs=c(0.05,0.95))))

# family level
detection_df_fam <- data.frame (family = family_data [sel_family],
                            p_fam =  samples_paleo_family$mean$p,
                            CI_fam = t(apply (samples_paleo_family$sims.list$p,2, quantile,probs=c(0.05,0.95)))
)
                            



# number of formations per interval
png (here ("output","detection.png"),res=300,unit="cm",height=15,width=18)
par (mfrow=c(1,2))
hist (detection_df$p,
      xlab = "Detection probability (p)",
      ylab = "Number of genera",
      main = "Detection probability\nacross genera")
# plot average
text (x=0.1,y=300,
      labels = paste ("Average+-SD=\n", 
                      round (mean(detection_df$p),2),
                      "+-",
                      round (sd(detection_df$p),2),
                      sep=" "))
# plot range
text (x=0.1,y=260,
      labels = paste ("Range=", 
                      round (range(detection_df$p),2)[1],
                      ",",
                      round (range(detection_df$p),2)[2],
                      sep=" "))


# family
hist (detection_df_fam$p,
      xlab = "Detection probability (p)",
      ylab = "Number of families",
      main = "Detection probability\nacross families")
# plot average
text (x=0.25,y=450,
      labels = paste ("Average+-SD=\n", 
                      round (mean(detection_df_fam$p),2),
                      "+-",
                      round (sd(detection_df_fam$p),2),
                      sep=" "))
# plot range
text (x=0.25,y=400,
      labels = paste ("Range=", 
                      round (range(detection_df_fam$p),2)[1],
                      ",",
                      round (range(detection_df_fam$p),2)[2],
                      sep=" "))


dev.off()


# ------------------------------------------------------------

# extinction
df_pars <- data.frame (EPSLON=1-rowMeans(samples_paleo_genus$mean$phi[-1,],na.rm=T),
                       GAMMA=rowMeans(samples_paleo_genus$mean$gamma[-1,],na.rm=T),
                       #muZ=rowSums(samples_paleo_genus$mean$muZ[-1,]),
                       # relative extinction (epslon / gamma)
                       RER = rowMeans(samples_paleo_genus$mean$RER[-1,],na.rm=T),
                       CI.gamma = t(apply (samples_paleo_genus$sims.list$gamma[,-1,],2, quantile,probs=c(0.05,0.95))),
                       CI.epslon = 1-t(apply (samples_paleo_genus$sims.list$phi[,-1,],2, quantile,probs=c(0.05,0.95))),
                       CI.RER = 1-t(apply (samples_paleo_genus$sims.list$RER[,-1,],2, quantile,probs=c(0.05,0.95))),
                       interval=c(triassic,jurassic)[-1],
                       level = "Genus")
df_pars$interval <- factor (df_pars$interval,
                            levels = c(triassic,jurassic)[-1])
# family scale
# extinction
df_pars_fam <- data.frame (EPSLON=1-rowMeans(samples_paleo_family$mean$phi[-1,],na.rm=T),
                       GAMMA=rowMeans(samples_paleo_family$mean$gamma[-1,],na.rm=T),
                       #muZ=rowSums(samples_paleo$mean$muZ[-1,]),
                       # relative extinction (epslon / gamma)
                       RER = rowMeans(samples_paleo_family$mean$RER[-1,],na.rm=T),
                       CI.gamma = t(apply (samples_paleo_family$sims.list$gamma[,-1,],2, quantile,probs=c(0.05,0.95))),
                       CI.epslon = 1-t(apply (samples_paleo_family$sims.list$phi[,-1,],2, quantile,probs=c(0.05,0.95))),
                       CI.RER = 1-t(apply (samples_paleo_family$sims.list$RER[,-1,],2, quantile,probs=c(0.05,0.95))),
                       interval=c(triassic,jurassic)[-1],
                       level = "Family")
df_pars_fam$interval <- factor (df_pars_fam$interval,
                            levels = c(triassic,jurassic)[-1])

# bind
df_pars <- rbind (df_pars,
                  df_pars_fam)

# plot


ggplot(melt (df_pars, id.vars=c("level","interval","CI.gamma.5.","CI.gamma.95.","CI.epslon.5.","CI.epslon.95.","CI.RER.5.","CI.RER.95.")),
       aes (x=interval,y=value,group=1)) + 
  facet_wrap(~variable+level,scales = "free_y",nrow=3,ncol=2)+
  geom_line(size=2)+
  scale_colour_viridis_d(end=0.8)+
  theme_bw()+
  
  theme (axis.text.x = element_text(angle=90)) + 
  ylab ("Expected Number of Genus") + 
  xlab ("Geological Interval") + 
  #geom_linerange(aes(ymin=CI.gamma.5., ymax=CI.gamma.95.))+
  #geom_linerange(aes(ymin=CI.epslon.5., ymax=CI.epslon.95.))+
  #geom_linerange(aes(ymin=CI.RER.5., ymax=CI.RER.95.))+
  geom_vline(xintercept=5.5, 
             linetype="dashed", 
             color = "red", 
             size=1) 

ggsave (here ("output", "dyn_par.png"),
        dpi=300,
        units="cm",
        height=15,
        width=20)


# plots for each taxon
par(mfrow=c(1,1))
plot (NA,xaxt="n",xlim = c(1,16),ylim= c(0,1),
      ylab="Probability",
      xlab = "",
      main = "Genus level")

lapply (seq (1,ncol(samples_paleo_genus$mean$phi)), function (i) {

  # extinction
  lines (seq(1,16), 
         1-samples_paleo_genus$mean$phi[-1,i],
         col= rgb(0.10,0.5,0.1,alpha=0.1))
  # origination
  lines (seq(1,16), 
         samples_paleo_genus$mean$gamma[-1,i],
         col= rgb(0.10,0.1,0.5,alpha=0.1))
})

axis(side=1, at=seq (1,16), 
     labels=c(triassic,jurassic)[-1],las=2)  

legend ("topright",legend=c("Extinction", "Origination"),
        col = c("green3","blue"),
        lty=1,
        lwd=2,
        bty="n")


# the same plot for families
par(mfrow=c(1,1))
plot (NA,xaxt="n",xlim = c(1,16),ylim= c(0,1),
      ylab="Probability",
      xlab = "",
      main="Family level")

lapply (seq (1,ncol(samples_paleo_family$mean$phi)), function (i) {
  
  # extinction
  lines (seq(1,16), 
         1-samples_paleo_family$mean$phi[-1,i],
         col= rgb(0.10,0.5,0.1,alpha=0.1))
  # origination
  lines (seq(1,16), 
         samples_paleo_family$mean$gamma[-1,i],
         col= rgb(0.10,0.1,0.5,alpha=0.1))
})

axis(side=1, at=seq (1,16), 
     labels=c(triassic,jurassic)[-1],las=2)  

legend ("topright",legend=c("Extinction", "Origination"),
        col = c("green3","blue"),
        lty=1,
        lwd=2,
        bty="n")



# ------------------------------------------------------------
# posterior exceedance probability of psi.eq
# positive
pep_psi.eq.pos <- lapply (seq (1,ncol (samples_paleo_genus$sims.list$psi.eq )), function (i)
  
  sum(samples_paleo_genus$sims.list$psi.eq[,i]>0.5)/length(samples_paleo_genus$sims.list$psi.eq[,i])
  
)

# negative
pep_psi.eq.neg <- lapply (seq (1,ncol (samples_paleo_genus$sims.list$psi.eq )), function (i)
  
  sum(samples_paleo_genus$sims.list$psi.eq[,i]<0.5)/length(samples_paleo_genus$sims.list$psi.eq[,i])
  
)

# genus
dat_genus <- data.frame (genus = genus_data[sel_genus],
                         psi.eq = samples_paleo_genus$mean$psi.eq,
                         pep_pos = unlist(pep_psi.eq.pos),
                         pep_neg = unlist(pep_psi.eq.neg),
                         CI = t(apply (samples_paleo_genus$sims.list$psi.eq,2, quantile,probs=c(0.05,0.95))))

# match with taxonomic information
agg_genus <- coll_occ_taxa_jur_tri %>%
  group_by (genus) %>%
  summarise (env=getmode(taxon_environment),
             phylum=getmode(phylum))

# match
dat_genus <- cbind (dat_genus, 
                    env = agg_genus$env  [match (dat_genus$genus, agg_genus$genus)],
                    phylum =agg_genus$phylum  [match (dat_genus$genus, agg_genus$genus)]
)

# order
dat_genus <- dat_genus[order(dat_genus$psi.eq,decreasing=T),]

# factor to select taxa to plot
dat_genus$factor_plot <- NA
dat_genus$factor_plot[1:20] <- "Highest"
dat_genus$factor_plot[(nrow(dat_genus)-20):nrow(dat_genus)] <- "Lower"


# find the closest genera to maximum and minimum

# plot
ggplot (dat_genus %>% 
          filter (is.na(factor_plot) !=T) %>%
          filter (is.na(psi.eq) !=T) %>%
          filter (is.na(genus) !=T), 
        aes (x=reorder(genus, psi.eq),
                        y=psi.eq,group=1))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  theme(axis.text.x = element_text(angle=90,size=10))+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.5,col="gray20")



ggsave (here ("output", "psi_eq_genera.png"),
        dpi=300,
        units="cm",
        height=15,
        width=20)

# proportion
table(samples_paleo_genus$mean$psi.eq >=0.5)/sum(table(samples_paleo_genus$mean$psi.eq >=0.5))
table(samples_paleo_genus$mean$psi.eq ==0.5)/sum(table(samples_paleo_genus$mean$psi.eq ==0.5))

# considering posterior exceedance probability of 80%

# plot
ggplot (dat_genus %>% 
          filter (pep_pos > 0.8 | pep_neg > 0.8) %>%
          filter (phylum != "NO_PHYLUM_SPECIFIED") %>%
          #filter (CI.5. > 0.5 & CI.95. > 0.5) %>% # none
          filter (is.na(genus) !=T), 
        aes (x=reorder(genus, psi.eq),
             y=psi.eq,group=1))+
  geom_point(aes (colour=phylum))+
  #geom_point(aes (colour=env))+
  xlab ("Phylum")+
  ylab ("Equilibrium occupancy")+
  geom_line()+
  theme_bw()+
  #scale_colour_viridis_d()+
  facet_wrap(~phylum)+
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  theme(axis.text.x = element_blank(),
        legend.position = "none")#+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.5,col="gray20")

  # the same relationship for each environment
  # plot
  ggplot (dat_genus %>% 
            filter (pep_pos > 0.8 | pep_neg > 0.8) %>%
            filter (phylum != "NO_PHYLUM_SPECIFIED") %>%
            #filter (CI.5. > 0.5 & CI.95. > 0.5) %>% # none
            filter (is.na(genus) !=T), 
          aes (x=reorder(genus, psi.eq),
               y=psi.eq,group=1))+
    geom_point(aes (colour=env))+
    xlab ("Phylum")+
    ylab ("Equilibrium occupancy")+
    geom_line()+
    theme_bw()+
    #scale_colour_viridis_d()+
    facet_wrap(~env)+
    geom_hline(yintercept=0.5, 
               linetype="dashed", 
               color = "red", 
               size=1)+
    theme(axis.text.x = element_blank(),
          legend.position = "none")#+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.5,col="gray20")
  
  
ggsave (here ("output", "psi_eq_genera_PEP.png"),
        dpi=300,
        units="cm",
        height=15,
        width=30)

# ------------------------------------------------------------
# family level

# posterior exceedance probability of psi.eq
# positive
pep_psi.eq.pos <- lapply (seq (1,ncol (samples_paleo_family$sims.list$psi.eq )), function (i)
  
  sum(samples_paleo_family$sims.list$psi.eq[,i]>0.5)/length(samples_paleo_family$sims.list$psi.eq[,i])
  
)

# negative
pep_psi.eq.neg <- lapply (seq (1,ncol (samples_paleo_family$sims.list$psi.eq )), function (i)
  
  sum(samples_paleo_family$sims.list$psi.eq[,i]<0.5)/length(samples_paleo_family$sims.list$psi.eq[,i])
  
)

# fams
dat_fam <- data.frame (family = family_data[sel_family], 
                         psi.eq = samples_paleo_family$mean$psi.eq,
                         pep_pos = unlist(pep_psi.eq.pos),
                         pep_neg = unlist(pep_psi.eq.neg),
                         CI = t(apply (samples_paleo_family$sims.list$psi.eq,2, quantile,probs=c(0.05,0.95))))


# match with taxonomic information
dat_fam <- cbind (dat_fam, 
                    env = coll_occ_taxa_jur_tri$taxon_environment  [match (dat_fam$family, coll_occ_taxa_jur_tri$family)],
                    phylum =coll_occ_taxa_jur_tri$phylum  [match (dat_fam$family, coll_occ_taxa_jur_tri$family)]
)

# order
dat_fam <- dat_fam[order(dat_fam$psi.eq,decreasing=T),]

# factor to select taxa to plot
dat_fam$factor_plot <- NA
dat_fam$factor_plot[1:20] <- "Highest"
dat_fam$factor_plot[(nrow(dat_fam)-20):nrow(dat_fam)] <- "Lower"

# plot and save
png (here ("output","psi_eq.png"),res=300,unit="cm",height=15,width=18)
par (mfrow=c(1,2))
hist (samples_paleo_genus$mean$psi.eq,
      xlab = "Equilibrium occupancy",
      ylab= "Frequency (# of genera)",
      main="Genera")
abline (v=0.5,col = "black", lty = 2,lwd=2)
hist (samples_paleo_family$mean$psi.eq,
      xlab = "Equilibrium occupancy",
      ylab= "Frequency (# of families)",
      main="Family")
abline (v=0.5,col = "black", lty = 2,lwd=2)
dev.off()

# find the closest genera to maximum and minimum

# plot
ggplot (dat_fam %>% 
          filter (is.na(factor_plot) !=T) %>%
          filter (is.na(psi.eq) !=T) %>%
          filter (is.na(family) !=T), 
        aes (x=reorder(family, psi.eq),
             y=psi.eq,group=1))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  theme(axis.text.x = element_text(angle=90,size=10))+
  geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.5,col="gray20")



ggsave (here ("output", "psi_eq_family.png"),
        dpi=300,
        units="cm",
        height=15,
        width=20)

# proportion
table(samples_paleo_family$mean$psi.eq >=0.5)/sum(table(samples_paleo_family$mean$psi.eq >=0.5))
table(samples_paleo_family$mean$psi.eq ==0.5)/sum(table(samples_paleo_family$mean$psi.eq ==0.5))

# considering posterior exceedance probability of 80%

# plot
ggplot (dat_fam %>% 
          filter (pep_pos > 0.8 | pep_neg > 0.8) %>%
          filter (family != "NO_FAMILY_SPECIFIED") %>%
          filter (phylum != "NO_PHYLUM_SPECIFIED") %>%
          #filter (CI.5. > 0.5 & CI.95. > 0.5) %>% # none
          filter (is.na(family) !=T), 
        aes (x=reorder(family, psi.eq),
             y=psi.eq,group=1))+
  geom_point(aes (colour=phylum))+
  xlab ("Genera")+
  ylab ("Equilibrium occupancy")+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept=0.5, 
             linetype="dashed", 
             color = "red", 
             size=1)+
  theme(axis.text.x = element_text(angle=90,size=5))#+
geom_pointrange(aes(ymin=CI.5., ymax=CI.95.),alpha=0.5,col="gray20")


ggsave (here ("output", "psi_eq_family_PEP.png"),
        dpi=300,
        units="cm",
        height=15,
        width=30)


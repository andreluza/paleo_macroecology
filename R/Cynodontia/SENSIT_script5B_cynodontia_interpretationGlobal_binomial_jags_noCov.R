

# --------------------------------------------


#     Interpretation: Global-scale analysis (SENSITIVITY, model without covariates)


# ------------------------------------------------------------
# load packages
rm(list=ls())
source("R/packages.R")
source("R/functions.R")

# load output
load (here ("output","No_covariate_model",
            "CMR_global_binomial_no_covMammalia.RData"))

# mammalia output
mammalia_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","No_covariate_model",
            "CMR_global_binomial_no_covNon-mammalian Mammaliaformes.RData"))

# mammaliaformes output
mammaliaformes_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","No_covariate_model",
            "CMR_global_binomial_no_covNon-mammaliaform cynodonts.RData"))

# mammalia output
cynodonts_output <- samples_paleo_cynodontia_binomial


# define the standard of colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen")

# ------------------------------------------------------------

# load basic data
load(here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))

# bind region covariates on time covariates
time_covariates$trop_area <- region_area
time_covariates$trop_coast <- region_coastalLine

# show metacommunity size

dat_diversity<-rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              Diversity = cynodonts_output$sims.list$Ntotal),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              Diversity = mammaliaformes_output$sims.list$Ntotal),
  data.frame (Taxon = "Mammalia", 
              Diversity = mammalia_output$sims.list$Ntotal)
)
  

# data augmented
p_div <-   ggplot (dat_diversity,aes (x=Diversity,fill=Taxon))+
  
  geom_density(lwd = 1,
               linetype = 1,
               col="gray20",
               alpha = 0.5)+
geom_histogram(aes(y = ..density..),
                 colour = 1,
               alpha=0.5) +
  theme_classic()+
  scale_fill_manual(values = cols)+
  scale_x_break(c(450, 1000), expand=T,scales="fixed")+
  #scale_y_break(c(0.06, 0.1), expand=T,scales="fixed")+
  theme(legend.position = "none",legend.direction = "vertical") + 
  
  # median
  geom_vline (data =dat_diversity %>%
                filter (Taxon == "Non-mammaliaform cynodonts"),  
              aes (xintercept=median(Diversity)),
              color= cols[1],
              size=1,
              linetype="dashed") +
  geom_vline (data =dat_diversity %>%
                filter (Taxon == "Non-mammalian Mammaliaformes"),  
              aes (xintercept=median(Diversity)),
              color= cols[2],
              size=1,
              linetype="dashed") +
  geom_vline (data =dat_diversity %>%
                filter (Taxon == "Mammalia"),  
              aes (xintercept=median(Diversity)),
              color= cols[3],
              size=1,
              linetype="dashed") 
  
dat_diversity %>%
  group_by(Taxon) %>%
  summarise(median(Diversity))

dat_diversity %>%
  group_by(Taxon) %>%
  summarise(hdi(Diversity))


p_div

ggsave (here ("output","No_covariate_model", "density_Div.pdf"),
        p_div,
        width=8,height =3)


# omega
png(here ("output", "No_covariate_model", "omega.png"),units = "cm",width=20,height=10,res=300)
par (mfrow=c(1,3))
hist(cynodonts_output$sims.list$omega,main = "Cynodonts",xlab="")
hist(mammaliaformes_output$sims.list$omega,main = "Mammaliaformes",xlab="Omega")
hist(mammalia_output$sims.list$omega,main = "Mammals",xlab="")

dev.off()
par (mfrow=c(1,1))

quantile (mammalia_output$sims.list$omega,c(0.025,0.5, 0.975))

# psi
png(here ("output", "No_covariate_model", "psi.png"),units = "cm",width=20,height=10,res=300)
par (mfrow=c(1,3))
hist(cynodonts_output$mean$muZW,xlab="Psi[g,t]",ylab= "Cynodonts")
hist(mammaliaformes_output$mean$muZW,xlab="Psi[g,t]",ylab= "Mammaliaformes")
hist(mammalia_output$mean$muZW,xlab="Psi[g,t]", ylab="Mammals")
dev.off()
par (mfrow=c(1,1))

# --------------------------------------------------

# plot of averages of gamma, phi, and p

# show metacommunity size

dat_rates<-rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              Origination = apply (cynodonts_output$sims.list$gamma,1,median,na.rm=T),
              Extinction = apply (1-cynodonts_output$sims.list$phi,1,median,na.rm=T),
              Detection = apply (cynodonts_output$sims.list$p,1,median,na.rm=T)),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              Origination = apply (mammaliaformes_output$sims.list$gamma,1,median,na.rm=T),
              Extinction = apply (1-mammaliaformes_output$sims.list$phi,1,median,na.rm=T),
              Detection = apply (mammaliaformes_output$sims.list$p,1,median,na.rm=T)),
  data.frame (Taxon = "Mammalia", 
              Origination = apply (mammalia_output$sims.list$gamma,1,median,na.rm=T),
              Extinction = apply (1-mammalia_output$sims.list$phi,1,median,na.rm=T),
              Detection = apply (mammalia_output$sims.list$p,1,median,na.rm=T))
)

# melt
dat_rates <- melt(dat_rates)

# averages to plot
df2<- dat_rates %>%
  group_by(Taxon, variable) %>%
  summarise(val=median(value))


# data augmented
p_rates <-   ggplot (dat_rates,aes (x=value,
                                    fill=Taxon))+
  geom_density(lwd = 1,
               linetype = 1,
               col="gray20",
               alpha = 0.2)+
  geom_histogram(aes(y = ..density..),
                 colour = 1,
                 alpha=0.5) +
  theme_classic()+
  scale_fill_manual(values = cols)+
  #scale_x_break(c(420, 640), expand=T,scales="fixed")+
  #scale_y_break(c(0.06, 0.1), expand=T,scales="fixed")+
  theme(legend.position = "top",legend.direction = "vertical") + 
  
  # median
  geom_vline (data =df2,
              aes (xintercept=val,
                   col = Taxon),
              
              size=1,
              linetype="dashed") +
  scale_colour_manual(values = cols)+
  
  facet_wrap(~variable,scales = "free")


p_rates

ggsave (here ("output","No_covariate_model", "p_rates.pdf"),
        p_rates,
        width=10,height =5)


# --------------------------------------------------


# plot covariates


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


# -----------------------------------------------------------------

# dynamic rates

# plot parameters
dat <-   rbind (
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn[-length(stages_cyn)],],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Origination probability",
              model = "Binomial",
              average = apply (cynodonts_output$sims.list$gamma,2,median,na.rm=T),
              lw= apply (cynodonts_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=apply (cynodonts_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[2,]),
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn[-length(stages_cyn)],],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-apply (cynodonts_output$sims.list$phi,2,median),
              lw=1-apply (cynodonts_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=1-apply (cynodonts_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[2,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf[-length(stages_mammaf)],],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Origination probability",
              model = "Binomial",
              average = apply (mammaliaformes_output$sims.list$gamma,2,median,na.rm=T),
              lw= apply (mammaliaformes_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=apply (mammaliaformes_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[2,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf[-length(stages_mammaf)],],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-apply (mammaliaformes_output$sims.list$phi,2,median),
              lw=1-apply (mammaliaformes_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=1-apply (mammaliaformes_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[2,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm[-length(stages_mamm)],],
              Taxon = "Mammalia",
              var= "Origination probability",
              model = "Binomial",
              average = apply (mammalia_output$sims.list$gamma,2,median,na.rm=T),
              lw= apply (mammalia_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=apply (mammalia_output$sims.list$gamma,2,quantile, c(0.025,0.975),na.rm=T)[2,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm[-length(stages_mamm)],],
              Taxon = "Mammalia",
              var= "Extinction probability",
              model = "Binomial",
              average = 1-apply (mammalia_output$sims.list$phi,2,median),
              lw=1-apply (mammalia_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[1,],
              up=1-apply (mammalia_output$sims.list$phi,2,quantile, c(0.025,0.975),na.rm=T)[2,])
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
    
  facet_wrap(~var+Taxon,scales = "fixed",ncol=3)+
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


# plot of species richness -------------------------------------------

dat1 <- rbind (
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammaliaform cynodonts"),]>0,na.rm=T)[stages_cyn],
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammalian Mammaliaformes"),]>0,na.rm=T)[stages_mammaf],
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,],
              Taxon = "Mammalia",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Mammalia"),]>0,na.rm=T)[stages_mamm],
              lw=NA,
              up=NA), 
  
  
  # expected values
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Expected richness",
              model="Binomial",
              average = apply (cynodonts_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[2,],
              lw=apply (cynodonts_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[1,],
              up=apply (cynodonts_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[3,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Expected richness",
              model="Binomial",
              average = apply (mammaliaformes_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[2,],
              lw=apply (mammaliaformes_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[1,],
              up=apply (mammaliaformes_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[3,]),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,],
              Taxon = "Mammalia",
              var= "Expected richness",
              model="Binomial",
              average = apply (mammalia_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[2,],
              lw=apply (mammalia_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[1,],
              up=apply (mammalia_output$sims.list$SRexp,2,quantile, c(0.025,0.5,0.975),na.rm=T)[3,])
  
  
  
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

# comparison
a<-ggplot (dat1 %>%
          filter (var == "Expected richness"), 
        aes(x=mid_ma,
            y=average, 
            fill=Taxon,
            group = Taxon,
            col=Taxon))+
  geom_line(size=1.5)+
  # points - average
  geom_point(aes (shape = var),position = position_jitter(width=0.1),size=3)+
  geom_ribbon(data = dat1%>%
                filter (var == "Expected richness"),
              aes(x=mid_ma, y=average, ymax=up, ymin=lw,fill=Taxon), 
              alpha=0.2) +
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  theme_bw()+
  
  
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
  )+
  ylab ("Taxonomic diversity")+
  geom_rect(aes(xmin = max_ma, 
                xmax = min_ma, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = cols_strip),
            col=NA,
            alpha = 0.2)+
  scale_fill_manual(values = cols_strip)
  


ggsave (here ("output", "No_covariate_model", "comparison.pdf"),
        a,
        width=6,height=4)  
  

# -------------------------------------

# correlation between taxa regarding SR dynamics

cor_test <- lapply (seq(1,nrow(cynodonts_output$sims.list$SRexp)), function (i) {     
  
  # build the DF  
  dat_cor <- rbind (
      
      data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,],
                  Taxon = "Non-mammaliaform cynodonts",
                  SR = cynodonts_output$sims.list$SRexp[i,]),
      
      data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,],
                  Taxon = "Non-mammalian Mammaliaformes",
                  SR = mammaliaformes_output$sims.list$SRexp[i,]),
      
      data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm,],
                  Taxon = "Mammalia",
                  SR = mammalia_output$sims.list$SRexp[i,])
      )
      
  #dat_cor$SR[is.na(dat_cor$SR)] <- 0
      
  # organize table
  tab_df <- cast (data=dat_cor,
        formula = interval_name~Taxon,
        value = "SR",
        fill=0)
  # order
  tab_df <- tab_df[match(bins$interval_name[-c(1:6)],tab_df$interval_name),]
  
    # correlation test
    cor_test <- cor (
      tab_df[,-1],
      use="pairwise.complete.obs"
    )
    # melt
    cor_test <- melt(cor_test)
    cor_test


})

# melt list into a DF 
cor_test_DF<- sapply (cor_test, "[[",3)
cor_test_DF_summary <- apply (cor_test_DF, 1,quantile, c(0.025,  0.5, 0.975))
apply (cor_test_DF, 1,hdi)
apply (cor_test_DF, 1,median)

# set rownames
rownames(cor_test_DF) <- c("Cyn.Cyn", "Mammaf.cyn", "Mamm.cyn",
                           "Cyn.mammaf", "Mammaf.mammaf","Mamm.mammaf",
                           "Cyn.mamm", "Mammaf.mamm","Mamm.mamm")
# melt
plot_cor <- ggplot (melt(cor_test_DF,as.is=T) %>%
          
          filter (X1 %in% c("Cyn.mammaf",
                            "Cyn.mamm",
                            "Mammaf.mamm") )
          
          ,
  
              aes (x=value))+
  
  geom_density(lwd = 1,
               linetype = 1,
               col="gray20",
               alpha = 0.2)+
  
  geom_histogram(aes(y = ..density..),
                 colour = 1,
                 alpha=0.5) +
  
  theme_classic()+
  scale_fill_manual(values = cols)+
  
  
  facet_wrap(~X1,scales = "free",nrow=3)+
  #scale_x_break(c(420, 640), expand=T,scales="fixed")+
  #scale_y_break(c(0.06, 0.1), expand=T,scales="fixed")+
  theme(legend.position = "top",legend.direction = "vertical") + 
  
  # median
  geom_vline (xintercept = 0,
              size=1,
              linetype="dashed")

# arrange
a_plotcor <- grid.arrange(a,
             plot_cor,
             layout_matrix = rbind (c(1,1,2),
                                    c(1,1,2),
                                    c(1,1,2)))

ggsave (here ("output", "No_covariate_model", "comparison.pdf"),
        a_plotcor,
        width=12,height=6)  


# ==============================

# convergence dynamic params
dat_rhat <- rbind (
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn[-length(stages_cyn)],],
              Taxon = "Non-mammaliaform cynodonts", 
              rhat = melt(cynodonts_output$Rhat$gamma),
              par = "Origination"),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn[-length(stages_cyn)],],
              Taxon = "Non-mammaliaform cynodonts", 
              rhat = melt(cynodonts_output$Rhat$phi),
              par = "Persistence"),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf[-length(stages_mammaf)],],
              Taxon = "Non-mammalian Mammaliaformes", 
              rhat = melt(mammaliaformes_output$Rhat$gamma),
              par = "Origination"),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf[-length(stages_mammaf)],],
              Taxon = "Non-mammalian Mammaliaformes", 
              rhat = melt(mammaliaformes_output$Rhat$phi),
              par = "Persistence"),
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm[-length(stages_mamm)],],
              Taxon = "Mammalia",
              rhat = melt(mammalia_output$Rhat$gamma),
              par = "Origination"),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm[-length(stages_mamm)],],
              Taxon = "Mammalia",
              rhat = melt(mammalia_output$Rhat$gamma),
              par = "Persistence")
  
  
  
  )



# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))
# adjust factor
dat_rhat$par <- factor (dat_rhat$par,
                           levels = rev(unique(dat_rhat$par)))
dat_rhat$interval_name <- factor (dat_rhat$interval_name,
                        levels = bins$interval_name[-c(1:6)])


# dynamic params
#plot
ggplot(data=  dat_rhat%>%
         filter (is.na(rhat.value) != T),
       aes(x=rhat.value, y=(interval_name))) +
  facet_wrap(~par+Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col= "red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ylab("Stage")+
  ggtitle ("Convergence of dynamic parameters")

ggsave (here ("output", "No_covariate_model",
              "convergence_dyn_global_params.png"),width =10,height=10)


# rhat total richness

# order of groups
dat_rhatSR <- rbind (
  
  data.frame (Taxon = "Non-mammaliaform cynodonts",
                          Rhat = cynodonts_output$Rhat$Ntotal),
  data.frame (Taxon = "Non-mammalian Mammaliaformes",
              Rhat = mammaliaformes_output$Rhat$Ntotal),
  data.frame (Taxon = "Mammalia",
              Rhat = mammalia_output$Rhat$Ntotal)
)
  
  
dat_rhatSR$Taxon <- factor (dat_rhatSR$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))


#plot
ggplot(data=  dat_rhatSR,
       aes(x=Rhat, y=(Taxon))) +
  #facet_wrap(~par+Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col= "red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  
  ylab("Clade")+
  ggtitle ("Convergence of \ntotal richness estimate")

ggsave (here ("output", "No_covariate_model","convergence_total_richness.png"),width =5,height=3)

# ------------------------------------------------------


# function  to calculate change in TD
change_fc <- function (output) {

  change_TD <- lapply (seq(1,dim(output$sims.list$z)[1]), function (k){
    
            # empty vector to fill
              CH <- rep (NA,ncol (output$sims.list$z[k,,]))
             
            # calculate difference   
            for (t in 2:length(CH)) {
                  
                  CH[t] <- (sum(output$sims.list$z[k,,t]) - sum(output$sims.list$z[k,,t-1]))
                  
                  
                }
                CH
  }
  )
  
  
    # melt
  change_TD<-do.call(rbind,change_TD)
  return(change_TD)
    
}

change_TD_cyn <- change_fc (cynodonts_output)
change_TD_mammaliaforms <- change_fc (mammaliaformes_output)
change_TD_mammals <- change_fc (mammalia_output)

# bind data to plot
dat_change <- rbind (
  
          data.frame ('average' = t(apply (change_TD_cyn,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
                          'lci' = t(apply (change_TD_cyn,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
                          'uci' = t(apply (change_TD_cyn,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
                          "Taxon" = "Non-mammaliaform cynodonts",
                          bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn,]),
          data.frame ('average' =  t(apply (change_TD_mammaliaforms,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
                      'lci' = t(apply (change_TD_mammaliaforms,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
                      'uci' = t(apply (change_TD_mammaliaforms,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
                      "Taxon" = "Non-mammalian Mammaliaformes",
                      bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf,])
          ,
          data.frame ('average' = t(apply (change_TD_mammals,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
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


# diversification
# function  to calculate change in TD
div_fc <- function (output) {
  
  # gamma
  res_gamma<-output$sims.list[[which(names(output$sims.list) == "gamma")]]
  
  # phi
  res_phi<-output$sims.list[[which(names(output$sims.list) == "phi")]]
  
  # calculate diversification
  div_rate <- lapply (seq(1,nrow (res_phi)), function (i)
    
    do.call(rbind,
    
      lapply (seq(1,ncol (res_phi)), function (t) {
  
        RER <- (1-res_phi[i,t])/res_gamma[i,t] ## relative extinction rate (μ/λ; Rabosky 2018) of each time
        R0 <- res_gamma[i,t]-(1-res_phi[i,t]) ## net diversification rate (r= μ - λ; Rabosky 2018) of each time
        output <- cbind (RER=RER,
                        R0=R0)
        output
          
      })
    )
  )
  
  names(div_rate)
  
  res <- list (RER = do.call(rbind,lapply (div_rate, function (i) i[, "RER"])),
               R0 = do.call(rbind,lapply (div_rate, function (i) i[, "R0"])))
  
  return(res)
  
}

# run function
div_cyn <- div_fc (cynodonts_output)
div_mammaliaforms <- div_fc (mammaliaformes_output)
div_mammals <- div_fc (mammalia_output)

# bind data to plot
dat_div <- rbind (
  
  data.frame ('average' = t(apply (div_cyn$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
              'lci' = t(apply (div_cyn$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_cyn$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammaliaform cynodonts",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_cyn[-length(stages_cyn)],]),
  data.frame ('average' = t(apply (div_mammaliaforms$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
              'lci' = t(apply (div_mammaliaforms$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_mammaliaforms$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Non-mammalian Mammaliaformes",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mammaf[-length(stages_mammaf)],])
  ,
  data.frame ('average' = t(apply (div_mammals$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,2],
              'lci' = t(apply (div_mammals$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,1],
              'uci' = t(apply (div_mammals$RER,2,quantile,c(0.025,0.5,0.975),na.rm=T))[,3],
              "Taxon" = "Mammalia",
              bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][stages_mamm[-length(stages_mamm)],])
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
    #ylim = c(-1,.2), # R0
    ylim = c(0,1000),
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
                            "change" = "average") 
            
            
  )


# plot the relationship
ggplot (rel_div_change,aes (x=average, 
                            y=change, 
                            group=Taxon,
                            col= Taxon))+
  geom_point()+
  geom_smooth(method="lm",formula = y ~ x)+ # poly(x, 2)
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols) +
  xlab ("Net diversification rate (NDRt)")+
  ylab ("Change in taxonomic diversity (ΔTDt)")+
  theme_bw()+
  theme(legend.position="top",
        legend.direction = "vertical")

ggsave (here ("output", "No_covariate_model", "rel_div_change.png"))




# insert a lag
rel_div_change_lag <- rel_div_change
rel_div_change_lag<-rel_div_change_lag %>% 
  filter(Taxon == "Non-mammaliaform cynodonts") 
plot(rel_div_change_lag$average[-length(rel_div_change_lag$average)],
     rel_div_change_lag$change[-1])



ggplot (rel_div_change_lag,aes (x=average, 
                            y=change, 
                            group=Taxon,
                            col= Taxon))+
  geom_point()+
  geom_smooth(method="lm",formula = y ~ poly(x, 2))+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols) +
  xlab ("Net diversification rate (NDRt)")+
  ylab ("Change in taxonomic diversity (ΔTDt)")+
  theme_bw()+
  theme(legend.position="top",
        legend.direction = "vertical")


# save plot

# params
pdf(here ("output","No_covariate_model", "panel_diversity.pdf"),width=15,height=17)
ret_panel <-gridExtra::grid.arrange (plot3,
                                     plot_change,
                                     plot_div,
                                     nrow=3)
dev.off()

# end script
rm(list=ls())

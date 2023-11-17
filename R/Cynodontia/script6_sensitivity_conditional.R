

# -----------------------------------------------------------------------


#     Interpretation: Global-scale analysis (SENSITIVITY ANALYSES)


# -----------------------------------------------------------------------

# load packages
source("R/packages.R")
source("R/functions.R")


# define the standard of colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen")


# load basic data
load (here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))



# ------------------------------------------------------------

# load results
# load model binomial output without covariates


load (here ("output",
            "interesting_params_global_occ_conditional.RData"))

# interesting_params_global_occ_conditional
# interesting_params_global_occ_no_cov
# binning time intervals (from Paleoverse)
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)

# apply colors for the strip

bins <- cbind (bins,cols_strip = rep (c("yes", "no"),20)[-1])
cols_strip <- c("yes" = "gray", "no" = "white")




# plot of species richness

dat1 <- rbind (
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammaliaform cynodonts"),-c(1,2)],na.rm=T),
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Non-mammalian Mammaliaformes"),-c(1,2)],na.rm=T),
              lw=NA,
              up=NA),
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Mammalia",
              var= "Observed richness",
              model="NA",
              average = colSums(array_genus_bin[which(clades %in% "Mammalia"),-c(1,2)],na.rm=T),
              lw=NA,
              up=NA), 
  
  
  # expected values
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_cynodonts$SRexp$stat$statistics[,"Mean"],
              lw=interesting_params_occ_cynodonts$SRexp$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_cynodonts$SRexp$stat$quantiles[,"97.5%"]),
  
  data.frame (bins[-c(1:6, 33:39),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_mammaliaformes$SRexp$stat$statistics[,"Mean"],
              lw=interesting_params_occ_mammaliaformes$SRexp$stat$quantiles[,"2.5%"],
              up=interesting_params_occ_mammaliaformes$SRexp$stat$quantiles[,"97.5%"]),
  data.frame (bins[-c(1:13),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")],
              Taxon = "Mammalia",
              var= "Expected richness",
              model="Binomial",
              average = interesting_params_occ_mammalia$SRexp$stat$statistics[,"Mean"],
              lw= interesting_params_occ_mammalia$SRexp$stat$quantiles[,"2.5%"],
              up= interesting_params_occ_mammalia$SRexp$stat$quantiles[,"97.5%"])
  
  
  
) # %>%

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
  geom_errorbar(aes(x = mid_ma, ymin = lw, ymax = up,col=var),
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
    ylim = c(0,500),
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



# ------------------------------------------------------
# change


load (here ("output",
            "interesting_params_global_change_no_cov.RData"))


# plot parameters
dat_change <-   rbind (
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_change_cynodonts$CH$stat$statistics[,"Mean"],
              lw= interesting_params_change_cynodonts$CH$stat$quantiles[,"2.5%"],
              up=interesting_params_change_cynodonts$CH$stat$quantiles[,"97.5%"]),
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_change_mammaliaformes$CH$stat$statistics[,"Mean"],
              lw=interesting_params_change_mammaliaformes$CH$stat$quantiles[,"2.5%"],
              up=interesting_params_change_mammaliaformes$CH$stat$quantiles[,"97.5%"]),
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Mammalia",
              var= "Origination probability",
              model = "Binomial",
              average = interesting_params_change_mammalia$CH$stat$statistics[,"Mean"],
              lw=interesting_params_change_mammalia$CH$stat$quantiles[,"2.5%"],
              up=interesting_params_change_mammalia$CH$stat$quantiles[,"97.5%"])
  
  
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
  
  geom_errorbar(aes(x = mid_ma, ymin = lw, ymax = up,col=Taxon),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=mid_ma, y=average, ymax=up, ymin=lw,col=Taxon,fill=Taxon), 
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



# weighted diversification
dat_divW <-   rbind (
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Non-mammaliaform cynodonts",
              var= "Net diversification rate",
              model = "Binomial",
              average = 1-interesting_params_occ_cynodonts$phi$stat$statistics[,"Mean"]*
                ifelse (interesting_params_occ_cynodonts$SRexp$post_prob$PP_neg >= 0.55,
                        NA,1)[-1] - interesting_params_occ_cynodonts$gamma$stat$statistics[,"Mean"],
              lw= 1-interesting_params_occ_cynodonts$phi$stat$quantiles[,"2.5%"]*
                ifelse (interesting_params_occ_cynodonts$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_cynodonts$gamma$stat$quantiles[,"2.5%"],
              
              up= 1-interesting_params_occ_cynodonts$phi$stat$quantiles[,"97.5%"]*
                ifelse (interesting_params_occ_cynodonts$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_cynodonts$gamma$stat$quantiles[,"97.5%"]),
  
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Non-mammalian Mammaliaformes",
              var= "Net diversification rate",
              model = "Binomial",
              average = 1-interesting_params_occ_mammaliaformes$phi$stat$statistics[,"Mean"]*
                ifelse (interesting_params_occ_mammaliaformes$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammaliaformes$gamma$stat$statistics[,"Mean"],
              lw= 1-interesting_params_occ_mammaliaformes$phi$stat$quantiles[,"2.5%"]*
                ifelse (interesting_params_occ_mammaliaformes$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammaliaformes$gamma$stat$quantiles[,"2.5%"],
              up= 1-interesting_params_occ_mammaliaformes$phi$stat$quantiles[,"97.5%"]*
                ifelse (interesting_params_occ_mammaliaformes$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammaliaformes$gamma$stat$quantiles[,"97.5%"]),
  
  
  data.frame (bins[-c(1:6),c("interval_name","mid_ma","max_ma", "min_ma", "cols_strip")][-1,],
              Taxon = "Mammalia",
              var= "Net diversification rate",
              model = "Binomial",
              average = 1-interesting_params_occ_mammalia$phi$stat$statistics[,"Mean"]*
                ifelse (interesting_params_occ_mammalia$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammalia$gamma$stat$statistics[,"Mean"],
              lw=1-interesting_params_occ_mammalia$phi$stat$quantiles[,"2.5%"]*
                ifelse (interesting_params_occ_mammalia$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammalia$gamma$stat$quantiles[,"2.5%"],
              up= 1-interesting_params_occ_mammalia$phi$stat$quantiles[,"97.5%"]*
                ifelse (interesting_params_occ_mammalia$SRexp$post_prob$PP_neg > 0.55,
                        NA,1)[-1] - interesting_params_occ_mammalia$gamma$stat$quantiles[,"97.5%"]))


dat_divW$Taxon <- factor (dat_divW$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))

# replace average by zero
dat_divW$average[is.na(dat_divW$average)] <- 0

# plot results
# binmial model
plot_divW <- ggplot (data = dat_divW ,
                     
                     aes(x=mid_ma,
                         y=-1*average, 
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
  
  geom_errorbar(aes(x = mid_ma, ymin = -1*lw, ymax = -1*up,col=var),
                width=0.1,size=1,position = position_jitter(width=0.1))+
  
  geom_ribbon(aes(x=mid_ma, y=-1*average, ymax=-1*up, ymin=-1*lw,fill=Taxon), 
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
    #ylim = c(0, 1),
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

# http://127.0.0.1:14643/graphics/199aa5c1-e3ce-40e5-ae1f-6a9dbd09d73e.png
plot_divW


# save plot

# supporting info

# params
pdf(here ("output","figures", "panel_diversity_SUpp_FIgS2-3.3.pdf"),width=15,height=17)
ret_panel <-gridExtra::grid.arrange (plot3,
                                     plot_change,
                                     plot_divW,
                                     nrow=3)
dev.off()

rm(list=ls())

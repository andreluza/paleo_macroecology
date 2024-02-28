

# --------------------------------------------


#     Interpretation: Global-scale analysis

# SENSITIVYT ANALYSIS USING n=50 IMPUTED GENERA

dir.create(here ("output","figures_sens"))

# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")
require(ggbreak)

# load output
load (here ("output","global",
            "CMR_global_binomial_50spp_Mammalia.RData"))
# mammalia output
mammalia_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_50spp_Non-mammalian Mammaliaformes.RData"))

# mammaliaformes output
mammaliaformes_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_50spp_Non-mammaliaform cynodonts.RData"))

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
  scale_x_break(c(160, 410), expand=T,scales="fixed")+
  #scale_y_break(c(0.06, 0.1), expand=T,scales="fixed")+
  theme(legend.position = "top",legend.direction = "vertical") + 
  
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

ggsave (here ("output","figures_sens", "density_Div.png"),
        p_div,
        width=8,height =5)


# total N of genera in the z matrix


sum(apply(cynodonts_output$mean$z,1,max))
sum(apply(mammaliaformes_output$mean$z,1,max))
sum(apply(mammalia_output$mean$z,1,max))


# --------------------------------------------------

# plot of averages of gamma, phi, and p


# show metacommunity size

dat_rates<-rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              Origination = cynodonts_output$sims.list$gamma.u,
              Extinction = 1-cynodonts_output$sims.list$phi.u,
              Detection = cynodonts_output$sims.list$p.u),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              Origination = mammaliaformes_output$sims.list$gamma.u,
              Extinction = 1-mammaliaformes_output$sims.list$phi.u,
              Detection = mammaliaformes_output$sims.list$p.u),
  data.frame (Taxon = "Mammalia", 
              Origination = mammalia_output$sims.list$gamma.u,
              Extinction = 1-mammalia_output$sims.list$phi.u,
              Detection = mammalia_output$sims.list$p.u)
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

ggsave (here ("output","figures_sens", "p_rates.png"),
        p_rates,
        width=10,height =5)



# omega
png(here ("output", "figures_sens", "omega.png"),units = "cm",width=20,height=10,res=300)
par (mfrow=c(1,3))
hist(cynodonts_output$sims.list$omega,main = "Cynodonts",xlab="")
hist(mammaliaformes_output$sims.list$omega,main = "Mammaliaformes",xlab="Omega")
hist(mammalia_output$sims.list$omega,main = "Mammals",xlab="")

dev.off()



# --------------------------------------------


#     Interpretation: regression coefficients unifyed global and regional


# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")


# define the standard of colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen",
          "Cynodonts" = "red", "Mammaliaformes" = "blue", "Mammals" = "darkgreen")

# ------------------------------------------------------------



# load basic data
load (here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))




# ------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "interesting_params_global_coeff_matrix.RData"))

output_list <- list(interesting_params_cynodonts,
                    interesting_params_mammaliaformes,
                    interesting_params_mammalia)




# taxon list
tax_list <- c("Non-mammaliaform cynodonts",
              "Non-mammalian Mammaliaformes",
              "Mammalia")


require(knitr)
table_coeff <- lapply (seq(1,length(output_list)), function (i) {
  
  table_coeff<-rbind (
    
    data.frame (
      par =  "Origination",
      taxon = tax_list[i],
      var = "Intercept (logit)",
      mean = output_list[[i]]$intercept.gamma$stat$statistics["Mean"],
      lci = output_list[[i]]$intercept.gamma$stat$quantiles["2.5%"],
      uci =  output_list[[i]]$intercept.gamma$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Precipitation",
                mean = output_list[[i]]$beta.gamma.prec$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.prec$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.prec$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.gamma.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.temp$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Land area",
                mean = output_list[[i]]$beta.gamma.area$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.area$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.area$stat$quantiles["97.5%"]),
    
    
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Intercept (logit)",
                mean = output_list[[i]]$intercept.phi$stat$statistics["Mean"],
                lci = output_list[[i]]$intercept.phi$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$intercept.phi$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Precipitation",
                mean = output_list[[i]]$beta.phi.prec$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.prec$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.prec$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.phi.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.temp$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Land area",
                mean = output_list[[i]]$beta.phi.area$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.area$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.area$stat$quantiles["97.5%"]),
    
    
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Intercept (logit)",
                mean = output_list[[i]]$intercept.p$stat$statistics["Mean"],
                lci = output_list[[i]]$intercept.p$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$intercept.p$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Time",
                mean = output_list[[i]]$beta.p.time$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.time$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.time$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Range size",
                mean = output_list[[i]]$beta.p.range$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.range$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.range$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Latitude",
                mean = output_list[[i]]$beta.p.lat$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.lat$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.lat$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.p.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.temp$stat$quantiles["97.5%"])
  )     
  table_coeff
  
})

# melt
table_coeff<-do.call(rbind, table_coeff)

table_coeff$taxon <- factor (table_coeff$taxon ,
                             levels = c("Non-mammaliaform cynodonts",
                                        "Non-mammalian Mammaliaformes",
                                        "Mammalia"))

#kable(., format = "pipe", padding = 2,align="c") 

# intercepts
kable (data.frame (taxon = table_coeff [grep("logit",table_coeff$var),"taxon"],
                   par = table_coeff [grep("logit",table_coeff$var),"par"],
                   #var = table_coeff [grep("logit",table_coeff$var),"var"],
                   mean_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"mean"]),3),
                   lci_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"lci"]),3),
                   uci_prob = round (plogis(table_coeff [grep("logit",table_coeff$var),"uci"]),3)),
       format = "pipe", padding = 2,align="c")


# Show the between-S CI's in red, and the within-S CI's in black
ggplot(data= table_coeff,
       aes(x=var, y=mean, group=par)) +
  facet_wrap(~taxon+par,scales="free")+
  geom_errorbar(width=.1, aes(ymin=lci, ymax=uci,col=taxon)) +
  geom_point(aes (col=taxon,fill=taxon),alpha=0.5,shape=21, size=3) +
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Coefficient value")+
  xlab ("Model parameter")+
  
  coord_flip()+
  
  theme_bw()+
  theme(legend.position = "none")

# ggsave (here ("output", "figures", "coeffb.png"),width=9,height=7)
# ggsave (here ("output", "figures", "coeffb.pdf"),width=9,height=7)


# area
# origination
interesting_params_cynodonts$beta.gamma.area$post_prob
interesting_params_mammaliaformes$beta.gamma.area$post_prob
interesting_params_mammalia$beta.gamma.area$post_prob
# persistence
interesting_params_cynodonts$beta.phi.area$post_prob
interesting_params_mammaliaformes$beta.phi.area$post_prob
interesting_params_mammalia$beta.phi.area$post_prob


# temperature
# origination
interesting_params_cynodonts$beta.gamma.temp$post_prob
interesting_params_mammaliaformes$beta.gamma.temp$post_prob
interesting_params_mammalia$beta.gamma.temp$post_prob

# persistence
interesting_params_cynodonts$beta.phi.temp$post_prob
interesting_params_mammaliaformes$beta.phi.temp$post_prob
interesting_params_mammalia$beta.phi.temp$post_prob

# Precipitation
# origination
interesting_params_cynodonts$beta.gamma.prec$post_prob
interesting_params_mammaliaformes$beta.gamma.prec$post_prob
interesting_params_mammalia$beta.gamma.prec$post_prob
# persistence
interesting_params_cynodonts$beta.phi.prec$post_prob
interesting_params_mammaliaformes$beta.phi.prec$post_prob
interesting_params_mammalia$beta.phi.prec$post_prob

# detection
# time
interesting_params_cynodonts$beta.p.time$post_prob
interesting_params_mammaliaformes$beta.p.time$post_prob
interesting_params_mammalia$beta.p.time$post_prob

# temperature
interesting_params_cynodonts$beta.p.temp$post_prob
interesting_params_mammaliaformes$beta.p.temp$post_prob
interesting_params_mammalia$beta.p.temp$post_prob

# range
interesting_params_cynodonts$beta.p.range$post_prob
interesting_params_mammaliaformes$beta.p.range$post_prob
interesting_params_mammalia$beta.p.range$post_prob

# latitude
interesting_params_cynodonts$beta.p.lat$post_prob
interesting_params_mammaliaformes$beta.p.lat$post_prob
interesting_params_mammalia$beta.p.lat$post_prob



# ---------------------------------------------
# convergence (RHat!)

# fit statistics
# Rhat

dat_rhat <- rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              rhat = unlist(sapply (interesting_params_cynodonts, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_cynodonts, "[[", "rhat")))),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              rhat = unlist(sapply (interesting_params_mammaliaformes, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_mammaliaformes, "[[", "rhat")))),
  data.frame (Taxon = "Mammalia",
              rhat=unlist(sapply (interesting_params_mammalia, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_mammalia, "[[", "rhat"))))
)

# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))
# adjust names
dat_rhat$par <-  substr(dat_rhat$par,
                        1,
                        round (nchar (dat_rhat$par)/2))


# plot
ggplot(data= dat_rhat,
       aes(x=rhat, y=par)) +
  facet_wrap(~Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col= "red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Parameter convergence in the global-scale model")

ggsave (here ("output", "figures","convergence_global_params.png"),width =10)


# -------------------------------------------------------------

# regional scale

# load results
# load model binomial output 
load (here ("output",
            "interesting_params_regional_coeff_matrix.RData"))


# table of parameters


output_list <- list(interesting_params_cynodonts_reg,
                    interesting_params_mammaliaformes_reg,
                    interesting_params_mammalia_reg)

# taxon list
tax_list <- c("Non-mammaliaform cynodonts",
              "Non-mammalian Mammaliaformes",
              "Mammalia")


require(knitr)
table_coeff_reg <- lapply (seq(1,length(output_list)), function (i) {
  
  table_coeff<-rbind (
    
    data.frame (
      par =  "Origination",
      taxon = tax_list[i],
      var = "Intercept (logit)",
      mean = (output_list[[i]]$intercept.gamma$stat$statistics["Mean"]),
      lci = (output_list[[i]]$intercept.gamma$stat$quantiles["2.5%"]),
      uci =  (output_list[[i]]$intercept.gamma$stat$quantiles["97.5%"])),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Precipitation",
                mean = output_list[[i]]$beta.gamma.prec$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.prec$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.prec$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.gamma.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.temp$stat$quantiles["97.5%"]),
    
    
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Land area",
                mean = output_list[[i]]$beta.gamma.area$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.gamma.area$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.gamma.area$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Intercept (logit)",
                mean = (output_list[[i]]$intercept.phi$stat$statistics["Mean"]),
                lci = (output_list[[i]]$intercept.phi$stat$quantiles["2.5%"]),
                uci =  (output_list[[i]]$intercept.phi$stat$quantiles["97.5%"])),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Precipitation",
                mean = output_list[[i]]$beta.phi.prec$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.prec$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.prec$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.phi.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.temp$stat$quantiles["97.5%"]),
    
    
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Land area",
                mean = output_list[[i]]$beta.phi.area$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.phi.area$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.phi.area$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Intercept (logit)",
                mean = output_list[[i]]$intercept.p$stat$statistics["Mean"],
                lci = output_list[[i]]$intercept.p$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$intercept.p$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Time",
                mean = output_list[[i]]$beta.p.time$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.time$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.time$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Range size",
                mean = output_list[[i]]$beta.p.range$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.range$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.range$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Latitude",
                mean = output_list[[i]]$beta.p.lat$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.lat$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.lat$stat$quantiles["97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.p.temp$stat$statistics["Mean"],
                lci = output_list[[i]]$beta.p.temp$stat$quantiles["2.5%"],
                uci =  output_list[[i]]$beta.p.temp$stat$quantiles["97.5%"])
  )     
  table_coeff
  
})

# melt
table_coeff_reg<-do.call(rbind, table_coeff_reg)
# order
table_coeff_reg$taxon <- factor (table_coeff_reg$taxon ,
                                 levels = c("Non-mammaliaform cynodonts",
                                            "Non-mammalian Mammaliaformes",
                                            "Mammalia"))


# intercepts

# intercepts
kable (data.frame (taxon = table_coeff_reg [grep("logit",table_coeff_reg$var),"taxon"],
                   par = table_coeff_reg [grep("logit",table_coeff_reg$var),"par"],
                   #var = table_coeff [grep("logit",table_coeff$var),"var"],
                   mean_prob = round (plogis(table_coeff_reg [grep("logit",table_coeff_reg$var),"mean"]),3),
                   lci_prob = round (plogis(table_coeff_reg [grep("logit",table_coeff_reg$var),"lci"]),3),
                   uci_prob = round (plogis(table_coeff_reg [grep("logit",table_coeff_reg$var),"uci"]),3)),
       format = "pipe", padding = 2,align="c")

# Show the between-S CI's in red, and the within-S CI's in black
ggplot(data= table_coeff_reg,
       aes(x=var, y=mean, group=par)) +
  facet_wrap(~taxon+par,scales="free")+
  geom_errorbar(width=.1, aes(ymin=lci, ymax=uci,col=taxon)) +
  geom_point(aes (col=taxon,fill=taxon),alpha=0.5,shape=21, size=3) +
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Coefficient value")+
  xlab ("Model parameter")+
  
  coord_flip()+
  
  theme_bw()+
  theme(legend.position = "none")

# ggsave (here ("output", "figures", "coeff_regional.png"),width=9,height=7)
# ggsave (here ("output", "figures", "coeff_regional.pdf"),width=9,height=7)




# area
# origination
interesting_params_mammalia_reg$beta.gamma.area$post_prob
interesting_params_mammaliaformes_reg$beta.gamma.area$post_prob
interesting_params_cynodonts_reg$beta.gamma.area$post_prob
# persistence
interesting_params_mammalia_reg$beta.phi.area$post_prob
interesting_params_mammaliaformes_reg$beta.phi.area$post_prob
interesting_params_cynodonts_reg$beta.phi.area$post_prob


# temperature
# origination
interesting_params_mammalia_reg$beta.gamma.temp$post_prob
interesting_params_mammaliaformes_reg$beta.gamma.temp$post_prob
interesting_params_cynodonts_reg$beta.gamma.temp$post_prob
# persistence
interesting_params_mammalia_reg$beta.phi.temp$post_prob
interesting_params_mammaliaformes_reg$beta.phi.temp$post_prob
interesting_params_cynodonts_reg$beta.phi.temp$post_prob

# Precipitation
# origination
interesting_params_mammalia_reg$beta.gamma.prec$post_prob
interesting_params_mammaliaformes_reg$beta.gamma.prec$post_prob
interesting_params_cynodonts_reg$beta.gamma.prec$post_prob
# persistence
interesting_params_mammalia_reg$beta.phi.prec$post_prob
interesting_params_mammaliaformes_reg$beta.phi.prec$post_prob
interesting_params_cynodonts_reg$beta.phi.prec$post_prob


# detection
# time
interesting_params_mammalia_reg$beta.p.time$post_prob
interesting_params_mammaliaformes_reg$beta.p.time$post_prob
interesting_params_cynodonts_reg$beta.p.time$post_prob

# range
interesting_params_mammalia_reg$beta.p.range$post_prob
interesting_params_mammaliaformes_reg$beta.p.range$post_prob
interesting_params_cynodonts_reg$beta.p.range$post_prob

# latitude
interesting_params_mammalia_reg$beta.p.lat$post_prob
interesting_params_mammaliaformes_reg$beta.p.lat$post_prob
interesting_params_cynodonts_reg$beta.p.lat$post_prob


# temperature
interesting_params_mammalia_reg$beta.p.temp$post_prob
interesting_params_mammaliaformes_reg$beta.p.temp$post_prob
interesting_params_cynodonts_reg$beta.p.temp$post_prob

# rhat


# ---------------------------------------------
# convergence (RHat!)

# fit statistics
# Rhat

dat_rhat <- rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              rhat = unlist(sapply (interesting_params_cynodonts_reg, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_cynodonts_reg, "[[", "rhat")))),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              rhat = unlist(sapply (interesting_params_mammaliaformes_reg, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_mammaliaformes_reg, "[[", "rhat")))),
  data.frame (Taxon = "Mammalia",
              rhat=unlist(sapply (interesting_params_mammalia_reg, "[[", "rhat")),
              par = names(unlist(sapply (interesting_params_mammalia_reg, "[[", "rhat"))))
)

# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))
# adjust names
dat_rhat$par <-  substr(dat_rhat$par,
                        1,
                        round (nchar (dat_rhat$par)/2))

# plot
ggplot(data= dat_rhat,
       aes(x=rhat, y=par)) +
  facet_wrap(~Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col="red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Parameter convergence in the region-scale model")

ggsave (here ("output", "figures","convergence_region_params.png"))



# ----------------------------------------

# unified plot

unified_df <- rbind (table_coeff_reg %>%
                       cbind ("scale" = "region"),
                     table_coeff %>%
                       cbind ("scale" = "global"))
unified_df$par <- factor (unified_df$par,
                          levels = c("Origination",
                                     "Persistence",
                                     "Detection"))


# plot
ggplot(data= unified_df,
       aes(x=var, y=mean, group=scale,shape=scale)) +
  facet_wrap(~taxon+par,scales="free")+
  geom_errorbar(width=.1, aes(ymin=lci, ymax=uci,col=taxon), 
                position = position_dodge(width = 0.75)) +
  geom_point(aes (col=taxon,fill=taxon,shape=scale),alpha=0.5,size=3, position = position_dodge(width = 0.75)) +
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Coefficient value")+
  xlab ("Model parameter")+
  
  coord_flip()+
  
  theme_bw()+
  theme(legend.position = "top")

ggsave (here ("output", "figures", "coeff_unified.pdf"),width=9,height=7)



# intercepts
kable (data.frame (scale = unified_df [grep("logit",unified_df$var),"scale"], 
                   taxon = unified_df [grep("logit",unified_df$var),"taxon"],
                   par = unified_df [grep("logit",unified_df$var),"par"],
                   #var = table_coeff [grep("logit",table_coeff$var),"var"],
                   mean_prob = round (plogis(unified_df [grep("logit",unified_df$var),"mean"]),3),
                   lci_prob = round (plogis(unified_df [grep("logit",unified_df$var),"lci"]),3),
                   uci_prob = round (plogis(unified_df [grep("logit",unified_df$var),"uci"]),3)
                   ),
       format = "pipe", padding = 2,align="c")


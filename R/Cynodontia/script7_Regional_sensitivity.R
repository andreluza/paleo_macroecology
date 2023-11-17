

# --------------------------------------------


#     Interpretation: Region-scale analysis


# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")


# define the standard of colors

cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen",
          "Cynodonts" = "red", "Mammaliaformes" = "blue", "Mammals" = "darkgreen")

# ------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "interesting_params_regional_rdm_area_occ.RData"))

# load basic data
load (here ("processed_data", "site_covs.RData"))
load(here ("processed_data","CMR_data_observation_cov.RData"))


# ==============================

# convergence dynamic params
dat_rhat <- rbind (
  data.frame (Taxon = "Cynodonts", 
              rhat = melt(sapply (interesting_params_occ_cynodonts_reg, "[[", "rhat"))),
  data.frame (Taxon = "Mammaliaformes", 
              rhat = melt(sapply (interesting_params_occ_mammaliaformes_reg, "[[", "rhat"))),
  data.frame (Taxon = "Mammals",
              rhat=melt(sapply (interesting_params_occ_mammalia_reg, "[[", "rhat"))))


# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Cynodonts",
                                     "Mammaliaformes",
                                     "Mammals"))
# adjust names
dat_rhat <- rbind (dat_rhat [grep ("phi",dat_rhat$rhat.L1),],
                   dat_rhat [grep ("gamma",dat_rhat$rhat.L1),]) %>%
  filter (is.na(rhat.value) != T)

dat_rhat$par2 <-  sapply (strsplit( dat_rhat$rhat.L1, "[.]"), "[[",1)

# order
dat_rhat$par <- factor (dat_rhat$par,
                        levels = rev(unique(dat_rhat$par)))

# dynamic params
#plot
ggplot(data=  dat_rhat,
       aes(x=rhat.value,fill= Taxon)) +
  geom_histogram()+
  geom_vline(xintercept=1.1,col= "red") + 
  facet_wrap(~Taxon+par2+rhat.X2,ncol=9)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  
  
  theme_bw() + 
  
  
  xlab ("Rhat value")+
  ggtitle ("Convergence of origination (gamma) and persistence (phi)\nin the region-scale model")

ggsave (here ("output", "figures","convergence_dyn_region_params.png"),width =20,height=15)


# -----------------------------------------------------
# plot params


# check the colours
library(RColorBrewer)
library(scales)
# cols <- brewer_pal(pal = "RdBu")(5) # not valid in 1.1-2
cols <- brewer.pal(n = 9, name = "RdYlGn") 
cols[5] <- "black"


# [1] "#CA0020" "#F4A582" "#F7F7F7" "#92C5DE" "#0571B0"
# show_col(cols) # not valid in 1.1-2
display.brewer.pal(n = 9, name = "RdYlGn")



# binning time intervals (from Paleoverse)
bins <- time_bins(interval = c("Permian", "Cretaceous"), 
                  rank = "stage",
                  plot=T)


# apply colors for the strip

bins <- cbind (bins,cols_strip = rep (c("yes", "no"),20)[-1])
cols_strip <- c("yes" = "gray", "no" = "white")


# diversification
function_plot_div <- function (output,title,label_yaxis) {
  
  
        
      # variation in genus richness
      # try to incorporate uncertainty
      SR_matrix_div <- matrix(1-(output$phi$stat$statistics[,"Mean"]),
                             nrow=9,
                             ncol=32,
                             byrow=T) - matrix(output$gamma$stat$statistics[,"Mean"],
                          nrow=9,
                          ncol=32,
                          byrow=T) 
                        
      # lower bound
      SR_matrix_div_lwr <-  matrix(1 -output$phi$stat$quantiles[,"2.5%"],
                                nrow=9,
                                ncol=32,
                                byrow=T) -
                          matrix(output$gamma$stat$quantiles[,"2.5%"],
                              nrow=9,
                              ncol=32,
                              byrow=T)  
                          
      # upper bound
      SR_matrix_div_upr <- matrix(1-output$phi$stat$quantiles[,"97.5%"],
                                  nrow=9,
                                  ncol=32,
                                  byrow=T) - 
                      matrix(output$gamma$stat$quantiles[,"97.5%"],
                              nrow=9,
                              ncol=32,
                              byrow=T)
      # set names
      rownames(SR_matrix_div) <- rownames(SR_matrix_div_lwr)<-rownames(SR_matrix_div_upr) <- cells
      colnames(SR_matrix_div) <- colnames(SR_matrix_div_lwr)<- colnames(SR_matrix_div_upr) <- time_bins[-1]
      
      # melt each one and then bind
      # average
      div_df <- cbind (melt (SR_matrix_div,as.is=T),
                      data = "Average")
      div_df$Latitude <- bins_lat [match (div_df$X1,bins_lat$bin),"mid"]
      div_df$age <- bins [match (div_df$X2,bins$bin),"mid_ma"]
      
      # lower interval
      div_df_lwd <- cbind (melt (as.data.frame(SR_matrix_div_lwr)),
                          data = "Lower.CI")
      # upper interval
      div_df_upr <- cbind (melt (as.data.frame(SR_matrix_div_upr)),
                              data = "Upper.CI")
      
      # plot with uncertainty
      p6 <- ggplot (data.frame (div_df,
                                Lower.CI=div_df_lwd$value,
                                Upper.CI =div_df_upr$value)# %>%
                      #filter (Latitude == 0)
                      , 
                    
                    aes (x=age,y=-1*value, 
                         col=as.factor(Latitude)),
                    position =  position_jitter(width=2)) + 
        
        geom_point(size=2,position =  position_jitter(width=2))+
        geom_errorbar(aes(x = age, 
                          ymin = -1*Lower.CI, 
                          ymax = -1*Upper.CI,
                          col=as.factor(Latitude)),
                      position =  position_jitter(width=2),
                      alpha=0.3,
                      width=0.01,
                      size=1)+
        geom_line(aes(group=as.factor (Latitude)))+
        #geom_ribbon(aes(x=age, y=-1*value, ymax=-1*Upper.CI, 
        #                ymin=-1*Lower.CI), 
        #            alpha=0.2)+
      
        theme_bw()+
        ylab ("Estimated richness")+
        theme(axis.title = element_text(size=15)) + 
        #scale_colour_viridis(option = "magma",direction = 1,end=0.8)
        #scale_colour_gradientn(colours = cols, 
        #                      values = rescale(c(-80,-60,-40, -10, 0, 10, 40,60, 80)),
        #                      guide = "colorbar", 
        #                      limits=c(-80, 80))+
        scale_colour_manual(values = cols)+
        scale_x_reverse("Age (Ma)") +
        coord_geo(
          dat = list("stages", "periods"), 
          xlim = c( 66,270), 
          ylim = c(-1, 1),
          pos = list("b", "b"),
          rot=90,
          size = list(2, 4),
          abbrv = list(TRUE, T)
        ) 
      
      p6

      # summary per region
      # summary per region
      p2 <- ggplot (data.frame (div_df,
                                Lower.CI=div_df_lwd$value,
                                Upper.CI =div_df_upr$value), 
                    
                    aes (x=as.factor(Latitude),y=-1*value, 
                         col=Latitude,fill=Latitude),
                    alpha=0.5,
                    position =  position_jitter(width=2)) +
        geom_boxplot(alpha=0.7) + 
        xlab ("Latitude")+
        theme_bw()+
        theme(legend.position = "none")+
        ylab ("")+
        scale_colour_gradientn(colours = cols, 
                               values = rescale(c(-80,-60,-40, -20, 0, 20, 40,60, 80)),
                               guide = "colorbar", 
                               limits=c(-80, 80)) + 
        scale_fill_gradientn(colours = cols, 
                             values = rescale(c(-80,-60,-40, -20, 0, 20, 40,60, 80)),
                             guide = "colorbar", 
                             limits=c(-80, 80)) + 
        
        ggtitle(title)
      
      
      p2
      
      # arrange
      grid.arrange (p2,
                    p6,
                    ncol=5,
                    nrow=2,
                    layout_matrix = rbind (
                      
                      c(1,1,2,2,2),
                      c(1,1,2,2,2)
                      
                      
                    )
      )
}

#ggsave (filename=here ("output", "figures","diversification_regional.png"))

pdf (here ("output", "figures", "div_space_time_sensitivity.pdf"),onefile=T,width=15,height=11)
grid.arrange (
  
  function_plot_div(output = interesting_params_occ_cynodonts_reg,
                      title = "Non-mammaliaform cynodonts",
                      label_yaxis = "Net Diversification Rate")
  ,
  
  function_plot_div(output = interesting_params_occ_mammaliaformes_reg,
                      title = "Non-mammalian Mammaliaformes",
                      label_yaxis = "Net Diversification Rate")
  
  ,
  function_plot_div(output = interesting_params_occ_mammalia_reg,
                      title = "Mammalia",
                      label_yaxis = "Net Diversification Rate")
  , nrow=3
  
)
dev.off()




# interesting parameters
# temporal trends


# function to produce plots with uncertainty around a desired parameter
parameter = "SRexp"
function_plot_sites <- function (output,title,parameter, label_yaxis) {
  
  
  # calculate statistics
  s <- output [which(names (output) %in% parameter)]
  
  
  # set ncols (t)
  if (parameter %in% c("phi", "gamma","CH","CH1", "R0")) {
    
    ncols <- length(time_bins[-1])
    
  } else {ncols <- length(time_bins)}
    
  # variation in genus richness
  # try to incorporate uncertainty
  SR_matrix <- matrix(s[[1]]$stat$statistics[,"Mean"],
                      nrow=length(cells),
                      ncol=ncols,
                      byrow=T)
  SR_matrix_lwr <- matrix(s[[1]]$stat$quantiles[,"2.5%"],
                          nrow=length(cells),
                          ncol=ncols,
                          byrow=T)
  SR_matrix_upr <- matrix(s[[1]]$stat$quantiles[,"97.5%"],
                          nrow=length(cells),
                          ncol=ncols,
                          byrow=T)
  # set names
  rownames(SR_matrix) <- rownames(SR_matrix_lwr)<-rownames(SR_matrix_upr) <- cells
  
  if (parameter %in% c("phi", "gamma","CH","CH1", "R0")) {
  
    colnames(SR_matrix) <- colnames(SR_matrix_lwr)<- colnames(SR_matrix_upr) <- time_bins[-1]
  
  } else {
    
    colnames(SR_matrix) <- colnames(SR_matrix_lwr)<- colnames(SR_matrix_upr) <- time_bins
    
    
  }
  
  # melt each one and then bind
  # average
  SR_df <- cbind (melt (SR_matrix,as.is=T),
                  data = "Average")
  SR_df$Latitude <- bins_lat [match (SR_df$X1,bins_lat$bin),"mid"]
  SR_df$age <- bins [match (SR_df$X2,bins$bin),"mid_ma"]
  
  # lower interval
  SR_df_lwd <- cbind (melt (as.data.frame(SR_matrix_lwr)),
                      data = "Lower.CI")
  # upper interval
  SR_matrix_upr <- cbind (melt (as.data.frame(SR_matrix_upr)),
                          data = "Upper.CI")
  
  # plot with uncertainty
  p6 <- ggplot (data.frame (SR_df,
                            Lower.CI=SR_df_lwd$value,
                            Upper.CI =SR_matrix_upr$value), 
                
                aes (x=age,y=value, 
                     col=as.factor (Latitude)),
                position =  position_jitter(width=2)) + 
    
    
    geom_errorbar(aes(x = age, 
                      ymin = Lower.CI, 
                      ymax = Upper.CI,
                      col=as.factor (Latitude)),
                  position =  position_jitter(width=2),
                  alpha=0.3,
                  width=0.01,
                  size=1)+
    geom_point(size=2,position =  position_jitter(width=2))+
    geom_line(aes(group=as.factor (Latitude)))+
    theme_bw()+
    ylab (label_yaxis)+
    theme(axis.title = element_text(size=15)) + 
    #scale_colour_viridis(option = "magma",direction = 1,end=0.8)
    #scale_colour_gradientn(colours = cols, 
    #                       values = rescale(c(-80,-60,-40, -10, 0, 10, 40,60, 80)),
    #                       guide = "colorbar", 
    #                       limits=c(-80, 80))+
    scale_colour_manual(values = cols)+
    scale_x_reverse("Age (Ma)") +
    coord_geo(
      dat = list("stages", "periods"), 
      xlim = c( 66,270), 
      ylim = c(0, max(SR_matrix_upr$value)),
      pos = list("b", "b"),
      rot=90,
      size = list(2, 4),
      abbrv = list(TRUE, T)
    ) 
  
  p6 <- p6 + ggtitle(title)
  p6
  
  # summary per region
  p2 <- ggplot (data.frame (SR_df,
                      Lower.CI=SR_df_lwd$value,
                      Upper.CI =SR_matrix_upr$value), 
          
          aes (x=as.factor(Latitude),y=value, 
               col=Latitude,fill=Latitude),
          alpha=0.5,
          position =  position_jitter(width=2)) +
    geom_boxplot(alpha=0.7) + 
    xlab ("Latitude")+
    theme_bw()+
    theme(legend.position = "none")+
    ylab (label_yaxis)+
    scale_colour_gradientn(colours = cols, 
                           values = rescale(c(-80,-60,-40, -20, 0, 20, 40,60, 80)),
                           guide = "colorbar", 
                           limits=c(-80, 80)) + 
    scale_fill_gradientn(colours = cols, 
                           values = rescale(c(-80,-60,-40, -20, 0, 20, 40,60, 80)),
                           guide = "colorbar", 
                           limits=c(-80, 80)) + 
    
    ggtitle(title)
    
  
  p2
  
  #pdf (here ("output", "figures", paste ("arrange_", parameter, ".pdf",sep="")),
  #     width=10,height=8)
  
  
  grid.arrange (p2,
                p6,
                ncol=5,
                nrow=2,
                layout_matrix = rbind (
                  
                  c(1,1,2,2,2),
                  c(1,1,2,2,2)
                  
                  
                )
  )
  
  
  
}


# save origination prob

pdf (here ("output", "figures", "SR_space_time_sensitivity.pdf"),onefile=T,width=15,height=11)
grid.arrange (

  function_plot_sites(output = interesting_params_occ_cynodonts_reg,
                    title = "Non-mammaliaform cynodonts",
                    parameter = "SRexp", 
                    label_yaxis = "Taxonomic diversity")
  ,

  function_plot_sites(output = interesting_params_occ_mammaliaformes_reg,
                    title = "Non-mammalian Mammaliaformes",
                    parameter = "SRexp", 
                    label_yaxis = "Taxonomic diversity")

  ,
  function_plot_sites(output = interesting_params_occ_mammalia_reg,
                    title = "Mammalia",
                    parameter = "SRexp", 
                    label_yaxis = "Taxonomic diversity")
  , nrow=3

)
dev.off()


# -------------------------------------------------------------
# load results
# load model binomial output 
load (here ("output",
            "interesting_params_regional_rdm_area_coeff_matrix.RData"))


# table of parameters


output_list <- list(interesting_params_cynodonts_reg,
                    interesting_params_mammaliaformes_reg,
                    interesting_params_mammalia_reg)

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
      var = "Average (logit)",
      mean = (output_list[[i]]$intercept.gamma$stat$statistics[,"Mean"]),
      lci = (output_list[[i]]$intercept.gamma$stat$quantiles[,"2.5%"]),
      uci =  (output_list[[i]]$intercept.gamma$stat$quantiles[,"97.5%"])),
    
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
                lci = output_list[[i]]$beta.gamma.temp$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$beta.gamma.temp$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Origination",
                taxon = tax_list[i],
                var = "Latitude",
                mean = output_list[[i]]$beta.gamma.temp$stat$statistics[,"Mean"],
                lci = output_list[[i]]$beta.gamma.temp$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$beta.gamma.temp$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Average (logit)",
                mean = (output_list[[i]]$intercept.phi$stat$statistics[,"Mean"]),
                lci = (output_list[[i]]$intercept.phi$stat$quantiles[,"2.5%"]),
                uci =  (output_list[[i]]$intercept.phi$stat$quantiles[,"97.5%"])),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Precipitation",
                mean = output_list[[i]]$beta.phi.prec$stat$statistics[,"Mean"],
                lci = output_list[[i]]$beta.phi.prec$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$beta.phi.prec$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Temperature",
                mean = output_list[[i]]$beta.phi.temp$stat$statistics[,"Mean"],
                lci = output_list[[i]]$beta.phi.temp$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$beta.phi.temp$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Latitude",
                mean = output_list[[i]]$beta.phi.lat$stat$statistics[,"Mean"],
                lci = output_list[[i]]$beta.phi.lat$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$beta.phi.lat$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Detection",
                taxon = tax_list[i],
                var = "Average (logit)",
                mean = output_list[[i]]$intercept.p$stat$statistics[,"Mean"],
                lci = output_list[[i]]$intercept.p$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$intercept.p$stat$quantiles[,"97.5%"]),
    
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
# order
table_coeff$taxon <- factor (table_coeff$taxon ,
                             levels = c("Non-mammaliaform cynodonts",
                                        "Non-mammalian Mammaliaformes",
                                        "Mammalia"))


# table intercepts

table_intercepts <- lapply (seq(1,length(output_list)), function (i) {
  
  table_intercepts<-rbind (
    
    data.frame (
      par =  "Origination",
      taxon = tax_list[i],
      var = "Average (logit)",
      mean = output_list[[i]]$intercept.gamma$stat$statistics[,"Mean"],
      lci = output_list[[i]]$intercept.gamma$stat$quantiles[,"2.5%"],
      uci =  output_list[[i]]$intercept.gamma$stat$quantiles[,"97.5%"]),
    
    data.frame (par =  "Persistence",
                taxon = tax_list[i],
                var = "Average (logit)",
                mean = output_list[[i]]$intercept.phi$stat$statistics[,"Mean"],
                lci = output_list[[i]]$intercept.phi$stat$quantiles[,"2.5%"],
                uci =  output_list[[i]]$intercept.phi$stat$quantiles[,"97.5%"])
    )     
  table_intercepts
  
})


cols <- c("Non-mammaliaform cynodonts" = "red", "Non-mammalian Mammaliaformes" = "blue", "Mammalia" = "darkgreen",
          "Cynodonts" = "red", "Mammaliaformes" = "blue", "Mammals" = "darkgreen")

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

ggsave (here ("output", "figures", "coeff_regional.png"),width=9,height=7)
ggsave (here ("output", "figures", "coeff_regional.pdf"),width=9,height=7)


# random intercept
table_intercepts<- do.call(rbind, table_intercepts)
# order
table_intercepts$taxon <- factor (table_intercepts$taxon ,
                             levels = c("Non-mammaliaform cynodonts",
                                        "Non-mammalian Mammaliaformes",
                                        "Mammalia"))

# regions
table_intercepts$regions <- rep(bins_lat$mid,6)


# Show the between-S CI's in red, and the within-S CI's in black
ggplot(data= table_intercepts,
       aes(x=plogis(mean), 
           y=regions, 
           group=par)) +
  facet_wrap(~par+taxon,scales="free")+
  geom_errorbar(width=.1, aes(xmin=plogis(lci), 
                              xmax=plogis(uci),
                              col=taxon)) +
  geom_point(aes (col=taxon,fill=taxon),alpha=0.5,shape=21, size=3) +
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Latitude of the region (20-degree bin)")+
  xlab ("Probability")+
  
  coord_flip()+
  
  theme_bw()+
  theme(legend.position = "none")

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
                        nchar (dat_rhat$par)/2+1,
                        round (nchar (dat_rhat$par)))


# plot
ggplot(data= dat_rhat,
       aes(x=rhat, y=par)) +
  facet_wrap(~Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1) + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Parameter convergence in the region-scale model")

ggsave (here ("output", "figures","convergence_global_params.png"))


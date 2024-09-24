

# --------------------------------------------


#     Interpretation: regression coefficients


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


# load output
load (here ("output","global",
            "CMR_global_binomial1000sp_Mammalia.RData"))

# mammalia output
mammalia_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial1000sp_Non-mammalian Mammaliaformes.RData"))

# mammaliaformes output
mammaliaformes_output <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial1000sp_Non-mammaliaform cynodonts.RData"))

# cynodonts output
cynodonts_output <- samples_paleo_cynodontia_binomial


# extract summary statistics
output_list <- list(cynodonts_output$summary,
                    mammaliaformes_output$summary,
                    mammalia_output$summary)



# taxon list
tax_list <- c("Non-mammaliaform cynodonts",
              "Non-mammalian Mammaliaformes",
              "Mammalia")


# list coefs
coefs <- c("intercept.gamma",
  #"beta.gamma.prec",
  "beta.gamma.temp",
  "beta.gamma.area",
  "beta.gamma.area.t",
  #"beta.gamma.coast",
  #"beta.gamma.coast.t",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "beta.phi.area",
  "beta.phi.area.t",
  "beta.phi.coast",
  "beta.phi.coast.t",
  "intercept.p",
  "beta.p.time",
  "beta.p.range",
  "beta.p.lat",
  "beta.p.temp") 

require(knitr)
table_coeff <- lapply (seq(1,length(output_list)), function (i) 
  
  do.call(rbind, 
          
          lapply (coefs, function (k) {
  
          
          # extrat coeffs and CI
          table_coeff<-    data.frame (
              taxon = tax_list[i],
              var = k,
              mean = output_list[[i]] [which(rownames (output_list[[i]]) %in% k),"mean"],
              lci = output_list[[i]] [which(rownames (output_list[[i]]) %in% k),"2.5%"],
              uci =  output_list[[i]] [which(rownames (output_list[[i]]) %in% k),"97.5%"])
            
          
  
        }
    )
  )
)


# melt
table_coeff<-do.call(rbind, table_coeff)

# oragnize factor
table_coeff$taxon <- factor (table_coeff$taxon ,
                             levels = c("Non-mammaliaform cynodonts",
                                        "Non-mammalian Mammaliaformes",
                                        "Mammalia"))

#kable(., format = "pipe", padding = 2,align="c") 

# type of parameter
table_coeff$par <- c(rep("Origination",4),
  rep("Persistence",7),
  rep ("Detection",5))

# edit variable names
table_coeff<-table_coeff %>%
  mutate (var2 = recode (var, 
                         "intercept.gamma" = "Intercept",
                         #"beta.gamma.prec" = "Precipitation",
                         "beta.gamma.temp" = "Temperature",
                         "beta.gamma.area" = "Area",
                         "beta.gamma.area.t" = "Tr. area",
                         #"beta.gamma.coast" = "Fragmentation",
                         #"beta.gamma.coast.t" = "Tr. fragmentation",
                         "intercept.phi" = "Intercept",
                         "beta.phi.prec" = "Precipitation",
                         "beta.phi.temp" = "Temperature",
                         "beta.phi.area" = "Area",
                         "beta.phi.area.t" = "Tr. area",
                         "beta.phi.coast" = "Fragmentation",
                         "beta.phi.coast.t" = "Tr. fragmentation",
                         "intercept.p" = "Intercept",
                         "beta.p.time" = "Time", 
                         "beta.p.range" = "Range size",
                         "beta.p.lat" = "Latitude",         
                         "beta.p.temp" = "Temperature"
                  
                  ))

# oragnize factor
table_coeff$var2 <- factor (table_coeff$var2,
                             levels = rev(c("Intercept",
                                        "Area",
                                         "Fragmentation",
                                        "Tr. area",
                                        "Tr. fragmentation",
                                        "Precipitation",
                                        "Temperature",
                                        "Time",
                                        "Range size",
                                        "Latitude"
                             )))

# adjust paprams
table_coeff$par <- factor (table_coeff$par,
                            levels = c("Origination", "Persistence", "Detection"))

# intercepts
kable (data.frame (taxon = table_coeff [grep("intercept",table_coeff$var),"taxon"],
                   par = table_coeff [grep("intercept",table_coeff$var),"par"],
                   #var = table_coeff [grep("logit",table_coeff$var),"var"],
                   mean_prob = round (plogis(table_coeff [grep("intercept",table_coeff$var),"mean"]),3),
                   lci_prob = round (plogis(table_coeff [grep("intercept",table_coeff$var),"lci"]),3),
                   uci_prob = round (plogis(table_coeff [grep("intercept",table_coeff$var),"uci"]),3)),
       format = "pipe", padding = 2,align="c")


# Show the between-S CI's in red, and the within-S CI's in black
ggplot(data= table_coeff %>%
         filter (var2 != "Intercept"),
       aes(x=var2, y=mean, group=par)) +
  facet_wrap(~taxon+par,scales="free")+
  geom_errorbar(width=.1, aes(ymin=lci, ymax=uci,col=taxon)) +
  geom_point(aes (col=taxon,fill=taxon),alpha=0.5,shape=21, size=3) +
  #scale_colour_manual(values = cols)+
  #scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0, linewidth=1,alpha=0.3)+
  ylab ("Coefficient value")+
  xlab ("Model parameter")+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")


# ggsave (here ("output", "figures", "coeffb.png"),width=9,height=7)
ggsave (here ("output", "figures", "coefficients.pdf"),width=9,height=7)

#----------------------
# posterior exceedance probabilities

tab_pp <- lapply (coefs[-grep("intercept",coefs)], function (coefi) {
  
  # build df 
    df_pp <- data.frame (coef = coefi, 
          cyn = cynodonts_output$sims.list [which(names(cynodonts_output$sims.list) %in% coefi)],
          mammaf = mammaliaformes_output$sims.list [which(names(mammaliaformes_output$sims.list) %in% coefi)],
          mamm = mammalia_output$sims.list [which(names(mammalia_output$sims.list) %in% coefi)])
       
    # calcualte probabilities
    tab_pp <- rbind (
                  apply (df_pp[,-1],2,function(x) {(sum(x>0)/3000)}),
                  apply (df_pp[,-1],2,function(x) {(sum(x<0)/3000)}) 
    )
    tab_pp
    
  }
  
  
)

tab_pp



# ---------------------------------------------
# convergence (RHat!)

# fit statistics
# Rhat

dat_rhat <- rbind (
  data.frame (Taxon = "Non-mammaliaform cynodonts", 
              rhat = unlist(cynodonts_output$Rhat[coefs]),
              par = names(cynodonts_output$Rhat[coefs])),
  data.frame (Taxon = "Non-mammalian Mammaliaformes", 
              rhat = unlist(mammaliaformes_output$Rhat[coefs]),
              par = names(mammaliaformes_output$Rhat[coefs])),
  data.frame (Taxon = "Mammalia",
              rhat = unlist(mammalia_output$Rhat[coefs]),
              par = names(mammalia_output$Rhat[coefs]))
)

# order of groups
dat_rhat$Taxon <- factor (dat_rhat$Taxon,
                          levels = c("Non-mammaliaform cynodonts",
                                     "Non-mammalian Mammaliaformes",
                                     "Mammalia"))
# adjust names
dat_rhat$par <-  c("Intercept.gamma",
                   #"Beta.gamma.prec",
                   "Beta.gamma.temp",
                   "Beta.gamma.area",
                   "Beta.gamma.areaTr",
                  # "Beta.gamma.fragm",
                   #"Beta.gamma.fragmTr",
                   
                   "Intercept.phi",
                   "Beta.phi.prec",
                   "Beta.phi.temp",
                   "Beta.phi.area",
                   "Beta.phi.areaTr",
                   "Beta.phi.fragm",
                   "Beta.phi.fragmTr",
                   
                   "Intercept.p",
                   "Beta.p.time",
                   "Beta.p.range",
                   "Beta.p.latit",
                   "Beta.p.temp"
                   
                   
                   )
                   
                  

# plot
ggplot(data= dat_rhat %>%
         mutate (rhat = round (rhat,1))
         ,
       aes(x=rhat, y=par)) +
  facet_wrap(~Taxon)+
  geom_bar(stat="identity")+
  
  geom_vline(xintercept=1.1,col= "red") + 
  theme_bw() + 
  
  xlab ("Rhat value")+
  ggtitle ("Parameter convergence in the global-scale model")

ggsave (here ("output", "figures","convergence_global_coef_params.png"),width =10)


# ------------------------------------

# compare the magnitude of effect


magnitude_check <-lapply (tax_list, function (i)
  
  lapply (coefs[-grep("intercept", coefs)], function (x) {
        
        subs_table_coeff <- table_coeff %>%
          filter (taxon == i &
                    var == x)
        
        # check precipitation magnitude
        prec_gamma <- data.frame(high_gamma=(sum(abs(mammalia_output$sims.list$beta.gamma.prec) > abs(subs_table_coeff$mean))/3000),
                      var=x)
        prec_phi <- data.frame (  high_phi=sum(abs(mammalia_output$sims.list$beta.phi.prec) > abs(subs_table_coeff$mean))/3000,
                            var=x)
        
        # check area magnitude
        area_gamma <- data.frame (high_gamma=sum(abs(mammalia_output$sims.list$beta.gamma.area) > abs(subs_table_coeff$mean))/3000,
                      var=x)
        area_phi <- data.frame (high_phi=sum(abs(mammalia_output$sims.list$beta.phi.area) > abs(subs_table_coeff$mean))/3000,
                      var=x)
        
        # check temp magnitude
        temp_gamma <- data.frame (high_gamma=sum(abs(mammalia_output$sims.list$beta.gamma.temp) > abs(subs_table_coeff$mean))/3000,
                                  var=x)
        temp_phi <- data.frame (high_phi=sum(abs(mammalia_output$sims.list$beta.phi.temp) > abs(subs_table_coeff$mean))/3000,
                                var=x)
        
        # check area tropical magnitude
        area.t_gamma <- data.frame (high_gamma=sum(abs(mammalia_output$sims.list$beta.gamma.area.t) > abs(subs_table_coeff$mean))/3000,
                                  var=x)
        area.t_phi <- data.frame (high_phi=sum(abs(mammalia_output$sims.list$beta.phi.area.t) > abs(subs_table_coeff$mean))/3000,
                                var=x)
        
        # check fragmentation magnitude
        coast_gamma <- data.frame (high_gamma=sum(abs(mammalia_output$sims.list$beta.gamma.coast) > abs(subs_table_coeff$mean))/3000,
                                    var=x)
        coast_phi <- data.frame (high_phi=sum(abs(mammalia_output$sims.list$beta.phi.coast) > abs(subs_table_coeff$mean))/3000,
                                  var=x)
        
        # check tropical fragmentation magnitude
        coast.t_gamma <- data.frame (high_gamma=sum(abs(mammalia_output$sims.list$beta.gamma.coast.t) > abs(subs_table_coeff$mean))/3000,
                                   var=x)
        coast.t_phi <- data.frame (high_phi=sum(abs(mammalia_output$sims.list$beta.phi.coast.t) > abs(subs_table_coeff$mean))/3000,
                                 var=x)
        
        
        # detection
        # check tropical fragmentation magnitude
        time_det <- data.frame (high_P=sum(abs(mammalia_output$sims.list$beta.p.time) > abs(subs_table_coeff$mean))/3000,
                                   var=x)
        range_det <- data.frame (high_P=sum(abs(mammalia_output$sims.list$beta.p.range) > abs(subs_table_coeff$mean))/3000,
                              var=x)
        lat_det <- data.frame (high_P=sum(abs(mammalia_output$sims.list$beta.p.lat) > abs(subs_table_coeff$mean))/3000,
                               var=x)
        temp_det <- data.frame (high_P=sum(abs(mammalia_output$sims.list$beta.p.temp) > abs(subs_table_coeff$mean))/3000,
                             var=x)
        
        
        # bind res into a list
        res <- list(
          
                    prec_gamma=prec_gamma,
                    prec_phi=prec_phi,
                    area_gamma = area_gamma,
                    area_phi=area_phi,
                    area.t_gamma=area.t_gamma,
                    area.t_phi=area.t_phi,
                    temp_gamma=temp_gamma,
                    temp_phi=temp_phi,
                    coast_gamma=coast_gamma,
                    coast_phi=coast_phi,
                    coast.t_gamma = coast.t_gamma,
                    coast.t_phi=coast.t_phi,
                    time_det=time_det,
                    range_det=range_det,
                    lat_det=lat_det,
                    temp_det=temp_det                    
                    
                    )
        res
  }
  )
  )


# build plots for all clades

# and save
pdf (file = here ("output", "figures", "PP_corrplot.pdf"),onefile = T)

lapply (seq(1,length(tax_list)), function (i) {

          # melt
          test_out<-(melt(magnitude_check[[i]]))
          
          # var 1
          
          test_out$var1 <- table_coeff$var2 [match (test_out$var,table_coeff$var)]
          test_out<-test_out %>%
            mutate (var2 = recode (L2, 
                                   "prec_gamma" = "Precipitation",
                                   "temp_gamma" = "Temperature",
                                   "area_gamma" = "Area",
                                   "area.t_gamma" = "Tr. area",
                                   "coast_gamma" = "Fragmentation",
                                   "coast.t_gamma" = "Tr. fragmentation",
                                   "prec_phi" = "Precipitation",
                                   "temp_phi" = "Temperature",
                                   "area_phi" = "Area",
                                   "area.t_phi" = "Tr. area",
                                   "coast_phi" = "Fragmentation",
                                   "coast.t_phi" = "Tr. fragmentation",
                                   "time_det" = "Time", 
                                   "range_det" = "Range size",
                                   "lat_det" = "Latitude",         
                                   "temp_det" = "Temperature"
                                   
            ))
          
          
          # table to plot
          tab_for_corrplot <- (cast (formula = var2 ~ var1,
                value = "value",
                fun.aggregate = mean,
            data = test_out %>% 
              filter (grepl ("gamma", L2)) %>%
              filter (grepl ("gamma", var))
              
            ))
          
          rownames(tab_for_corrplot) <- tab_for_corrplot$var2; tab_for_corrplot <- tab_for_corrplot[,-1]
          # order
          tab_for_corrplot <- tab_for_corrplot[order(rownames(tab_for_corrplot),decreasing = T),
                                               order(colnames(tab_for_corrplot),decreasing = T)]
          # data to plot
          tab_for_corrplot <- data.matrix(tab_for_corrplot)
          
          # gamma
          require(corrplot)
          corrplot(tab_for_corrplot,
                   is.corr = FALSE,
                             method = "square",
                             diag=F,
                             col = COL2('RdYlBu', 8),
                             mar = rep(4,4),
                             title = paste ("Comparison of covariates, origination of\n",
                                            tax_list[i]))
          
          # ---------------
          # persistence
          
          # table to plot
          tab_for_corrplot <- (cast (formula = var2 ~ var1,
                                     value = "value",
                                     fun.aggregate = mean,
                                     data = test_out %>% 
                                       filter (grepl ("phi", L2)) %>%
                                       filter (grepl ("phi", var))
                                     
          ))
          
          rownames(tab_for_corrplot) <- tab_for_corrplot$var2; tab_for_corrplot <- tab_for_corrplot[,-1]
          # order
          tab_for_corrplot <- tab_for_corrplot[order(rownames(tab_for_corrplot),decreasing = T),
                                               order(colnames(tab_for_corrplot),decreasing = T)]
          # data to plot
          tab_for_corrplot <- data.matrix(tab_for_corrplot)
          
          # gamma
          require(corrplot)
          corrplot(tab_for_corrplot,
                   is.corr = FALSE,
                   method = "square",
                   diag=F,
                   col = COL2('RdYlBu', 8),
                   mar = rep(4,4),
                   title = paste ("Comparison of covariates, persistence of\n",
                                  tax_list[i]))
          
          # -------------------------------
          # detection
          # table to plot
          tab_for_corrplot <- (cast (formula = var2 ~ var1,
                                     value = "value",
                                     fun.aggregate = mean,
                                     data = test_out %>% 
                                       filter (grepl ("det", L2)) %>%
                                       filter (grepl ("\\b.p.\\b", var))
                                     
          ))
          
          rownames(tab_for_corrplot) <- tab_for_corrplot$var2; tab_for_corrplot <- tab_for_corrplot[,-1]
          # order
          tab_for_corrplot <- tab_for_corrplot[order(rownames(tab_for_corrplot),decreasing = T),
                                               order(colnames(tab_for_corrplot),decreasing = T)]
          # data to plot
          tab_for_corrplot <- data.matrix(tab_for_corrplot)
          
          # gamma
          require(corrplot)
          corrplot(tab_for_corrplot,
                   is.corr = FALSE,
                   method = "square",
                   diag=F,
                   col = COL2('RdYlBu', 8),
                   mar = rep(4,4),
                   title = paste ("Comparison of covariates, detection of\n",
                                  tax_list[i]))
          
})       

dev.off()

# end script
rm (list=ls())
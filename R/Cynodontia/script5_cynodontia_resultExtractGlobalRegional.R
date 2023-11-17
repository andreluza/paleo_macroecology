# -----------------------------------------------------------


#     Extract results Global- and Region-scale analysis


# ------------------------------------------------------------

# load packages
source("R/packages.R")
source("R/functions.R")

# load output and data
require(coda)
require(mcmcr)
require(scales)


# ------------------------------------------------------------ #


# find the parameter
test <- combineSaves(recoverSaves(
  
  paste ("output/region_rdm_area/CMR_regional_binomial_rdm",
         "Non-mammaliaform cynodonts",
         sep="_")), 
  
  params =  "intercept.gamma",thin = 10)

plot(test)
hist(unlist(test))
hist(plogis(unlist(test)))

# function to extract parameter 
ext_param <- function (parameter,group,path,nthin=1) {
  
  # find the parameter
  extract_param <- combineSaves(recoverSaves(
                                        
                                      paste (path,
                                             group,
                                             sep="_")), 
                                
                                params =  parameter,thin = nthin)
 
  # calculate statistics
  s <- summary(extract_param,quantiles = c(0.025, 0.975),na.rm=T)
  # function to extract RHat statistics
  est_rhat <- rhat(extract_param, by = "term",as_df=F)
  # function to extract parameter 
  # find the parameter
  
  PP<-apply (do.call(rbind,extract_param),2,function (i) {
    
    PP_pos <- sum(i >0) / length(i)
    PP_neg <- sum(i <=0) / length(i)
    PP <- data.frame (PP_neg=PP_neg,
                      PP_pos=PP_pos)
    PP
  })
  PP<-do.call(rbind,PP)
  
  # build a list
  output <- list (stat = s,
                  rhat = est_rhat,
                  post_prob = PP)
  
  return (output)
  
}


# extract posterior dist matrices of Z and muY

# function to extract parameter 
ext_mat <- function (parameter,group,path,nthin=1) {
  
  # find the parameter
  extract_param <- combineSaves(recoverSaves(
    
    paste (path,
           group,
           sep="_")), 
    
    params =  parameter,thin = nthin)
  
  # get matrices
  PP_mat <- lapply (extract_param, function (i) 
        lapply (seq(1,nrow (i)), function (k) {
        
          
        # obtain matrices
        # vector to matrix
    
        if (parameter %in% c("z", "muY")) {
          
          dims_mat <- gsub ("[^0-9,]" , "" , names(i[k,]),  perl = TRUE)
          nrow_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",1)))
          ncol_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",2)))
          
          
          # matrix
        
          mat_z <- matrix (i[k,],
               nrow = nrow_mat ,
               ncol = ncol_mat,
               byrow=F)
        
       
    
          } else {
            
            
            
            mat_z <- i[k,]
          
          
            
          }
          
          return (mat_z)
          
  }))
  
  # calculate statistics
  s <- summary(extract_param,quantiles = c(0.025, 0.975),na.rm=T)
  # function to extract RHat statistics
  est_rhat <- rhat(extract_param, by = "term",as_df=F)
  
  
  # build a list
  output <- list (stat = s,
                  rhat = est_rhat,
                  PP_mat = PP_mat)
  
  return (output)
  
}



# extract the interesting parameters (coefficients)
params <- c(
  "intercept.gamma",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "beta.gamma.area",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "beta.phi.area",
  "intercept.p",
  "beta.p.time",
  "beta.p.range",
  "beta.p.lat",
  "beta.p.temp"
)

# extract the interesting parameters (occurrence and dynamics)
params_occ <- c(
  "gamma",
  "phi",
  "SRexp",
  "p"
)

# extract the interesting parameters (change)
params_change <- c(
  "z",
  "muY"
)

# --------------------------------------------------------

# extract global scale

# run function for each clade
interesting_params_mammaliaformes <- lapply (params, ext_param, 
                                             group  = "Non-mammalian Mammaliaformes",
                                             path = "output/global/CMR_global_binomial")
interesting_params_cynodonts <- lapply (params, ext_param, 
                                        group  = "Non-mammaliaform cynodonts",
                                        path = "output/global/CMR_global_binomial")
interesting_params_mammalia <- lapply (params, ext_param, 
                                       group  = "Mammalia",
                                       path = "output/global/CMR_global_binomial")
# set names
names(interesting_params_mammaliaformes) <- names(interesting_params_cynodonts) <- names(interesting_params_mammalia) <- params

# save
save (params,
      interesting_params_mammaliaformes,
      interesting_params_cynodonts,
      interesting_params_mammalia,
      file = here ("output", "interesting_params_global_coeff_matrix.RData") )


# extract dynamic parameters
interesting_params_occ_mammaliaformes <- lapply (params_occ, ext_mat, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path = "output/global/CMR_global_binomial",
                                                 nthin = 50)
interesting_params_occ_cynodonts <- lapply (params_occ, ext_mat, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path = "output/global/CMR_global_binomial",
                                            nthin = 10)
interesting_params_occ_mammalia <- lapply (params_occ, ext_mat, group  = "Mammalia",
                                           path = "output/global/CMR_global_binomial",
                                           nthin = 10)
# set names
names(interesting_params_occ_mammaliaformes) <- names(interesting_params_occ_cynodonts) <- names(interesting_params_occ_mammalia) <- params_occ
save (params_occ,
      interesting_params_occ_mammaliaformes,
      interesting_params_occ_cynodonts,
      interesting_params_occ_mammalia,
      file = here ("output", "interesting_params_global_occ.RData") )



# extract change parameters
interesting_params_change_mammaliaformes <- lapply (params_change, ext_mat, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path = "output/global/CMR_global_binomial",
                                                 nthin = 50)
interesting_params_change_cynodonts <- lapply (params_change, ext_mat, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path = "output/global/CMR_global_binomial",
                                            nthin = 50)
interesting_params_change_mammalia <- lapply (params_change, ext_mat, group  = "Mammalia",
                                           path = "output/global/CMR_global_binomial",
                                           nthin = 50)
# set names
names(interesting_params_change_mammaliaformes) <- names(interesting_params_change_cynodonts) <- names(interesting_params_change_mammalia) <- params_change
save (params_occ,
      interesting_params_change_mammaliaformes,
      interesting_params_change_cynodonts,
      interesting_params_change_mammalia,
      file = here ("output", "interesting_params_global_change.RData") )


# -----------------------------------------------


# sensitivity analysis (no covariate)

# run function for each clade
# extract dynamic parameters
interesting_params_occ_mammaliaformes <- lapply (params_occ, ext_mat, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                                 nthin=10)
interesting_params_occ_cynodonts <- lapply (params_occ, ext_mat, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                            nthin=10)
interesting_params_occ_mammalia <- lapply (params_occ, ext_mat, group  = "Mammalia",
                                           path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                           nthin=10)
# set names
names(interesting_params_occ_mammaliaformes) <- names(interesting_params_occ_cynodonts) <- names(interesting_params_occ_mammalia) <- params_occ
save (params_occ,
      interesting_params_occ_mammaliaformes,
      interesting_params_occ_cynodonts,
      interesting_params_occ_mammalia,
      file = here ("output", "interesting_params_global_occ_no_cov.RData") )


# extract change parameters
interesting_params_change_mammaliaformes <- lapply (params_change, ext_mat, 
                                                    group  = "Non-mammalian Mammaliaformes",
                                                    path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                                    nthin=10)
interesting_params_change_cynodonts <- lapply (params_change, ext_mat, 
                                               group  = "Non-mammaliaform cynodonts",
                                               path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                               nthin=10)
interesting_params_change_mammalia <- lapply (params_change, ext_mat, group  = "Mammalia",
                                              path = "output/No_covariate_model/CMR_global_binomial_no_cov",
                                              nthin=10)
# set names
names(interesting_params_change_mammaliaformes) <- names(interesting_params_change_cynodonts) <- names(interesting_params_change_mammalia) <- params_change
save (params_occ,
      interesting_params_change_mammaliaformes,
      interesting_params_change_cynodonts,
      interesting_params_change_mammalia,
      file = here ("output", "interesting_params_global_change_no_cov.RData") )




# -----------------------------------------------------------------------


# region scale estimates

# redesign function to arrays

# extract posterior dist matrices of Z and muY

# function to extract parameter 
parameter = "z"
group = "Non-mammalian Mammaliaformes"
nthin= 100

ext_mat_reg <- function (parameter,group,path,nthin=1) {
  
  # find the parameter
  extract_param <- combineSaves(recoverSaves(
    
    paste (path,
           group,
           sep="_")), 
    
    params =  parameter,thin = nthin)
  
  # get matrices
  PP_mat <- lapply (extract_param, function (i) 
    lapply (seq(1,nrow (i)), function (k) {
      
      
      # obtain matrices
      # vector to matrix
      
      if (parameter %in% c("z", "muY")) {
        
        dims_mat <- gsub ("[^0-9,]" , "" , names(i[k,]),  perl = TRUE)
        nrow_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",1)))
        ncol_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",2)))
        nz_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",3)))
        
        
        # matrix
        
        mat_z <- array (i[k,],
                         dim = c (nrow_mat ,
                                  ncol_mat,
                                  nz_mat))
        
        
        
      } else {
        
        
        dims_mat <- gsub ("[^0-9,]" , "" , names(i[k,]),  perl = TRUE)
        nrow_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",1)))
        ncol_mat <- max(as.numeric(sapply (strsplit (dims_mat, ","), "[[",2)))
        
        
        # matrix
        
        mat_z <- matrix (i[k,],
                         nrow = nrow_mat ,
                         ncol = ncol_mat,
                         byrow=F)
        
        
        
      }
      
      return (mat_z)
      
    }))
  
  # calculate statistics
  s <- summary(extract_param,quantiles = c(0.025, 0.975),na.rm=T)
  # function to extract RHat statistics
  est_rhat <- rhat(extract_param, by = "term",as_df=F)
  
  
  # build a list
  output <- list (stat = s,
                  rhat = est_rhat,
                  PP_mat = PP_mat)
  
  return (output)
  
}



# extract the interesting parameters (occurrence and dynamics)
params_occ <- c(
  "SRexp",
  "p",
  "z",
  "muY"
)

# run function for each clade
interesting_params_mammaliaformes_reg <- lapply (params, ext_param, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path="output/region/CMR_regional_binomial")
interesting_params_cynodonts_reg <- lapply (params, ext_param, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path="output/region/CMR_regional_binomial")
interesting_params_mammalia_reg <- lapply (params, ext_param, group  = "Mammalia",
                                           path="output/region/CMR_regional_binomial")
# set names
names(interesting_params_mammaliaformes_reg) <- names(interesting_params_cynodonts_reg) <- names(interesting_params_mammalia_reg) <- params

# save
save (params,
      interesting_params_mammaliaformes_reg,
      interesting_params_cynodonts_reg,
      interesting_params_mammalia_reg,
      file = here ("output", "interesting_params_regional_coeff_matrix.RData") )


# extract dynamic parameters
interesting_params_occ_mammaliaformes_reg <- lapply (params_occ, ext_mat_reg, 
                                                     group  = "Non-mammalian Mammaliaformes",
                                                     path="output/region/CMR_regional_binomial",
                                                     nthin=50)
interesting_params_occ_cynodonts_reg <- lapply (params_occ, ext_mat_reg, group  = "Non-mammaliaform cynodonts",
                                                path="output/region/CMR_regional_binomial",
                                                nthin=50)
interesting_params_occ_mammalia_reg <- lapply (params_occ, ext_mat_reg, group  = "Mammalia",
                                               path="output/region/CMR_regional_binomial",
                                               nthin=50)
# set names
names(interesting_params_occ_mammaliaformes_reg) <- names(interesting_params_occ_cynodonts_reg) <- names(interesting_params_occ_mammalia_reg) <- params_occ
save (params_occ,
      interesting_params_occ_mammaliaformes_reg,
      interesting_params_occ_cynodonts_reg,
      interesting_params_occ_mammalia_reg,
      file = here ("output", "interesting_params_regional_occ.RData") )



# ----------------------------------------------



# regional scale (region_rdm_area)

# run function for each clade
interesting_params_mammaliaformes_reg <- lapply (params, ext_param, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path="output/region_rdm_area/CMR_regional_binomial_rdm")
interesting_params_cynodonts_reg <- lapply (params, ext_param, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path="output/region_rdm_area/CMR_regional_binomial_rdm")
interesting_params_mammalia_reg <- lapply (params, ext_param, group  = "Mammalia",
                                           path="output/region_rdm_area/CMR_regional_binomial_rdm")
# set names
names(interesting_params_mammaliaformes_reg) <- names(interesting_params_cynodonts_reg) <- names(interesting_params_mammalia_reg) <- params

# save
save (params,
      interesting_params_mammaliaformes_reg,
      interesting_params_cynodonts_reg,
      interesting_params_mammalia_reg,
      file = here ("output", "interesting_params_regional_rdm_area_coeff_matrix.RData") )




# extract dynamic parameters
interesting_params_occ_mammaliaformes_reg <- lapply (params_occ, ext_mat_reg, 
                                                     group  = "Non-mammalian Mammaliaformes",
                                                     path="output/region_rdm_area/CMR_regional_binomial_rdm",
                                               nthin=50)
interesting_params_occ_cynodonts_reg <- lapply (params_occ, ext_mat_reg, group  = "Non-mammaliaform cynodonts",
                                                path="output/region_rdm_area/CMR_regional_binomial_rdm",
                                                nthin=50)
interesting_params_occ_mammalia_reg <- lapply (params_occ, ext_mat_reg, group  = "Mammalia",
                                               path="output/region_rdm_area/CMR_regional_binomial_rdm",
                                               nthin=50)
# set names
names(interesting_params_occ_mammaliaformes_reg) <- names(interesting_params_occ_cynodonts_reg) <- names(interesting_params_occ_mammalia_reg) <- params_occ
save (params_occ,
      interesting_params_occ_mammaliaformes_reg,
      interesting_params_occ_cynodonts_reg,
      interesting_params_occ_mammalia_reg,
      file = here ("output", "interesting_params_regional_rdm_area_occ.RData") )





rm(list=ls())





# -----------------------------------------------

params_occ <- c(
  "gamma",
  "phi",
  #"SRexp",
  "p"#,
  #"FSS"
)

# sensitivity analysis (conditioning data to the first detection)

# run function for each clade
# extract dynamic parameters
interesting_params_occ_mammaliaformes <- lapply (params_occ, ext_param, 
                                                 group  = "Non-mammalian Mammaliaformes",
                                                 path = "output/res_det_conditional/CMR_global_binomial")
interesting_params_occ_cynodonts <- lapply (params_occ, ext_param, 
                                            group  = "Non-mammaliaform cynodonts",
                                            path = "output/res_det_conditional/CMR_global_binomial")
interesting_params_occ_mammalia <- lapply (params_occ, ext_param, group  = "Mammalia",
                                           path = "output/res_det_conditional/CMR_global_binomial")
# set names
names(interesting_params_occ_mammaliaformes) <- names(interesting_params_occ_cynodonts) <- names(interesting_params_occ_mammalia) <- params_occ
save (params_occ,
      interesting_params_occ_mammaliaformes,
      interesting_params_occ_cynodonts,
      interesting_params_occ_mammalia,
      file = here ("output","interesting_params_global_occ_conditional.RData") )



# ----------------------------------------------


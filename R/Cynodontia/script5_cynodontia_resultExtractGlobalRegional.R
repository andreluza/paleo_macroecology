# --------------------------------------------


#     Interpretation: Region-scale analysis


# ------------------------------------------------------------

# load packages
source("R/packages.R")
source("R/functions.R")

# load output and data
require(coda)
require(mcmcr)
require(scales)

# ------------------------------------------------------------ #

# load basic data
# extract the interesting parameters (coefficients)
params <- c(
  "intercept.gamma",
  "beta.gamma.prec",
  "beta.gamma.temp",
  "intercept.phi",
  "beta.phi.prec",
  "beta.phi.temp",
  "intercept.p",
  "beta.p.time",
  "beta.p.range",
  "beta.p.lat",
  "beta.p.temp"
)

# function to extract parameter 
ext_param <- function (parameter) {
  
  # find the parameter
  extract_param <- combineSaves(recoverSaves("output/CMR_global_binomial"), params =  parameter,thin = 1)
 
  # calculate statistics
  s <- summary(extract_param,quantiles = c(0.025, 0.975))
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


# extract
interesting_params <- lapply (params, ext_param)
names(interesting_params) <- params
save (params,interesting_params, file = here ("output", "interesting_params_global_coeff.RData") )



# extract the interesting parameters (occurrence and dynamics)
## long form
params_occ <- c(
  "gamma",
  "phi",
  "SRexp",
  "p",
  "FSS",
  "CH",
  "CH1",
  "psi1",
  "R0",
  "RER"
)


# extract
interesting_params_occ <- lapply (params_occ, ext_param)
names(interesting_params_occ) <- params_occ
save (params_occ,interesting_params_occ, file = here ("output", "interesting_params_global_occ.RData") )


rm(list=ls())

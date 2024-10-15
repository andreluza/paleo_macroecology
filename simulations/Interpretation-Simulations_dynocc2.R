# ----------------------------------------

# Load necessary libraries --------------------------------------
rm(list=ls())
library(rjags)
library(here)

# load output ---------------------------------
# results
filenames <- list.files(here ("simulations", "output2"), full.names=T, pattern = "sims_run")
results <- sapply(filenames, function(x) mget(load(x)), simplify = TRUE)

# simulated data
filenames <- list.files(here ("simulations", "output2"), full.names=T,pattern = "data")
simdata <- sapply(filenames, function(x) mget(load(x)), simplify = F)

# bind
df_res <- do.call(rbind,results)

# find variables I want
hat_gamma <- df_res [grep("gamma", rownames(df_res)),]
hat_gamma <- hat_gamma [-grep("intercept",rownames(hat_gamma)),] # remove intercept
hat_gamma <- hat_gamma [-grep("beta",rownames(hat_gamma)),] # remove coeffs

# plot
require(reshape)
require(dplyr)
require(ggplot2)

my_theme <- theme(legend.position = 'bottom', 
                  strip.text = element_text(size=12),
                  strip.text.y = element_text(color = 'black'),
                  strip.text.x = element_text(color = 'black'), 
                  #text = element_text(family="LM Roman 10"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle = 0, hjust = 1, size = 10), 
                  axis.text.y = element_text(size = 10),
                  axis.title = element_text(size=15))

# plot phi and gamma --------------------------------------------

p_gamma<-data.frame (hat_gamma = melt(matrix(hat_gamma[,"mean"],ncol=50,byrow=F)),
            true_gamma = melt(sapply (simdata, "[[","gamma" ))[,3]) %>%
  ggplot() +
  theme_light(base_size = 16) +
  geom_point(aes(x = true_gamma, y = hat_gamma.value,group=hat_gamma.X2),alpha=0.5,col="gray20",size=4)+
  geom_smooth(aes(x = true_gamma, y = hat_gamma.value,group=hat_gamma.X2), col = 'black', alpha = 1, se = F, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2,linewidth=1) +
  labs(y=bquote("Estimated origination probability "*(hat(gamma[t]))*""),
       x=bquote("True origination probability "*(gamma[t])*""))+
  ggtitle ("A")+my_theme

# phi ----------
hat_phi <- df_res [grep("phi", rownames(df_res)),]
hat_phi <- hat_phi [-grep("intercept",rownames(hat_phi)),] # remove intercept
hat_phi <- hat_phi [-grep("beta",rownames(hat_phi)),] # remove coeffs

# plot
p_phi<-data.frame (hat_phi = melt(matrix(hat_phi[,"mean"],ncol=50,byrow=F)),
                   true_phi = melt(sapply (simdata, "[[","phi" ))[,3]) %>%
  ggplot() +
  theme_light(base_size = 16) +
  geom_point(aes(x = true_phi, y = hat_phi.value,group=hat_phi.X2),alpha=0.5,col="gray20",size=4)+
  geom_smooth(aes(x = true_phi, y = hat_phi.value,group=hat_phi.X2), col = 'black', alpha = 1, se = F, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2,linewidth=1) +
  labs(y=bquote("Estimated persistence probability "*(hat(phi[t]))*""),
       x=bquote("True persistence probability "*(phi[t])*""))+
  ggtitle ("B")+my_theme



# plot p --------------

plot_p <- data.frame (p_est = melt(matrix(df_res [grepl("p\\[", ignore.case = T, rownames(df_res)),"mean"],nrow=30,byrow=F)),
                      lab_p=as.character(expression(paste(p))),
                      p_true = melt(sapply (simdata, "[[","p" ))[,3]) %>%
  ggplot() +
  theme_light(base_size = 16) +
  geom_point(aes(x = p_true, y = p_est.value,group=p_est.X2),alpha=0.5,col="gray20",size=4)+
  geom_smooth(aes(x = p_true, y = p_est.value,group=p_est.X2), col = 'black', alpha = 1, se = F, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2,linewidth=1) +
  labs(y=bquote("Estimated detection probability "*(hat(p[t]))*""),
       x=bquote("True detection probability "*(p[t])*""))+
  ggtitle ("C")+my_theme


# reegression coefs ---------------------------------

  coef_plot <- data.frame (intercept.gamma = (df_res [grepl("\\intercept.gamma\\b$", ignore.case = T, rownames(df_res)),"mean"]),
            intercept.phi  = (df_res [grepl("\\intercept.phi\\b$", ignore.case = T, rownames(df_res)),"mean"]),
            beta_phi2 = df_res [grepl("\\bbeta_phi2\\b$", ignore.case = T, rownames(df_res)),"mean"],
            beta_gamma1 = df_res [grepl("\\bbeta_gamma1\\b$", ignore.case = T, rownames(df_res)),"mean"],
            beta_gamma2 = df_res [grepl("\\bbeta_gamma2\\b$", ignore.case = T, rownames(df_res)),"mean"]) %>%
  melt () %>%
  ggplot (aes (y=value, x=variable)) +
  scale_x_discrete(labels=c(expression(paste(hat(beta), " "[gamma[0]])),
                            expression(paste(hat(beta), " "[phi[0]])),
                            expression(paste(hat(beta), " "[phi[2]])),
                            expression(paste(hat(beta), " "[gamma[1]])),
                            expression(paste(hat(beta), " "[gamma[2]]))
                            
                            
                            ))+
  labs(y="Coefficient value (logistical scale)", x="Parameter")+
  geom_violin() +
  geom_point (aes(x=1, y=qlogis(0.2)))+
  geom_point (aes(x=2, y=qlogis(0.3)))+
  geom_point (aes(x=3, y=-1))+
  geom_point (aes(x=4, y=0.5))+
  geom_point (aes(x=5, y=1))+
  ggtitle ("D")+my_theme


# coverage
# gamma
cbind (df_res [grepl("\\intercept.gamma\\b$", ignore.case = T, rownames(df_res)),c("2.5%", "97.5%")],
       qlogis(0.2) ) %>%
  data.frame () %>%
  mutate(is_in_range = qlogis(0.2) > X2.5. & qlogis(0.2) < X97.5.) %>% # Create TRUE/FALSE column
  summarise(total_TRUEs = sum(is_in_range)) # Sum number of TRUEs

# phi
cbind (df_res [grepl("\\intercept.phi\\b$", ignore.case = T, rownames(df_res)),c("2.5%", "97.5%")],
       qlogis(0.3)) %>%
  data.frame () %>%
  mutate(is_in_range = qlogis(0.3) > X2.5. & qlogis(0.3) < X97.5.) %>% # Create TRUE/FALSE column
  summarise(total_TRUEs = sum(is_in_range)) # Sum number of TRUEs


# arrange plot
require(gridExtra)
png(file = here("simulations", "figs", "sc2.png"),width = 20,height = 20,units = "cm",res=300)
grid.arrange (p_gamma, p_phi,
              plot_p,coef_plot,ncol=2)

dev.off()

# end

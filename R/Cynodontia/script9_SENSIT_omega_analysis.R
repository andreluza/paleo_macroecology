

# --------------------------------------------


#     Interpretation: Global-scale analysis

# SENSITIVYT ANALYSIS USING DIFFERENT NUMBER OF IMPUTED GENERA


# ------------------------------------------------------------
# load packages
source("R/packages.R")
source("R/functions.R")
require(ggbreak)

# load output

# imputation of 50 genera


load (here ("output","global",
            "CMR_global_binomial_50spp_Mammalia.RData"))
# mammalia output
mammalia_output50 <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_50spp_Non-mammalian Mammaliaformes.RData"))

# mammaliaformes output
mammaliaformes_output50 <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_50spp_Non-mammaliaform cynodonts.RData"))

# mammalia output
cynodonts_output50 <- samples_paleo_cynodontia_binomial



# imputation of 400 genera

load (here ("output","global",
            "CMR_global_binomial_Mammalia.RData"))
# mammalia output
mammalia_output400 <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_Non-mammalian Mammaliaformes.RData"))

# mammaliaformes output
mammaliaformes_output400 <- samples_paleo_cynodontia_binomial

# load output
load (here ("output","global",
            "CMR_global_binomial_Non-mammaliaform cynodonts.RData"))

# mammalia output
cynodonts_output400 <- samples_paleo_cynodontia_binomial


# imputation of 500 genera

load (here ("output","global",
            "CMR_global_binomial_500spp_Mammalia.RData"))
# mammalia output

mammalia_output500 <- samples_paleo_cynodontia_binomial

# imputation of 500 genera cynodonts

# load output
load (here ("output","global",
            "CMR_global_binomial_500spp_Non-mammaliaform cynodonts.RData"))

# mammalia output
cynodonts_output500 <- samples_paleo_cynodontia_binomial

# imputation of 800 genera

load (here ("output","global",
            "CMR_global_binomial_800spp_Mammalia.RData"))
# mammalia output
mammalia_output800 <- samples_paleo_cynodontia_binomial

# imputation of 1000 genera

load (here ("output","global",
            "CMR_global_binomial_1000spp_Mammalia.RData"))
# mammalia output
mammalia_output1000 <- samples_paleo_cynodontia_binomial

# omega
png(here ("output", "figures_sens", "omega.png"),units = "cm",width=20,height=25,res=300)

par (mfrow=c(4,3))
hist(cynodonts_output50$sims.list$omega,main = "Cynodonts (n=50)",xlab="")
hist(mammaliaformes_output50$sims.list$omega,main = "Mammaliaformes (n=50)",xlab="")
hist(mammalia_output50$sims.list$omega,main = "Mammals (n=50)",xlab="")
hist(cynodonts_output400$sims.list$omega,main = "Cynodonts (n=400)",xlab="")
hist(mammaliaformes_output400$sims.list$omega,main = "Mammaliaformes (n=400)",xlab="")
hist(mammalia_output400$sims.list$omega,main = "Mammals (n=400)",xlab="")
hist(cynodonts_output500$sims.list$omega,main = "Cynodonts (n=500)",xlab="")
hist(mammalia_output500$sims.list$omega,main = "Mammals (n=500)",xlab="Omega")
hist(mammalia_output800$sims.list$omega,main = "Mammals (n=800)",xlab="Omega")
hist(mammalia_output1000$sims.list$omega,main = "Mammals (n=1000)",xlab="Omega")

dev.off()



# total divesrity curves

# N total diversity
png(here ("output", "figures_sens", "tot_div.png"),units = "cm",width=20,height=25,res=300)

par (mfrow=c(4,3))
hist(cynodonts_output50$sims.list$Ntotal,main = "Cynodonts (n=50)",xlab="")
hist(mammaliaformes_output50$sims.list$Ntotal,main = "Mammaliaformes (n=50)",xlab="")
hist(mammalia_output50$sims.list$Ntotal,main = "Mammals (n=50)",xlab="")
hist(cynodonts_output400$sims.list$Ntotal,main = "Cynodonts (n=400)",xlab="")
hist(mammaliaformes_output400$sims.list$Ntotal,main = "Mammaliaformes (n=400)",xlab="")
hist(mammalia_output400$sims.list$Ntotal,main = "Mammals (n=400)",xlab="")
hist(cynodonts_output500$sims.list$Ntotal,main = "Cynodonts (n=500)",xlab="")
hist(mammalia_output500$sims.list$Ntotal,main = "Mammals (n=500)",xlab="Number of genera")
hist(mammalia_output800$sims.list$Ntotal,main = "Mammals (n=800)",xlab="Number of genera")
hist(mammalia_output1000$sims.list$Ntotal,main = "Mammals (n=1000)",xlab="Number of genera")

dev.off()

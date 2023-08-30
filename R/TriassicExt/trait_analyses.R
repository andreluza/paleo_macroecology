# Trait -based analyses 

# Why the niche of some taxa expand over time, while other shrink? Is this a matter of pure evolution, biotic interactions, or both? These questions 
# led much of research about ecological niche evolution over time, and give rise to the field of ecospace analysis...

# The tempo and mode of clade diversification are key topics in
# current macroevolutionary and biodiversity studies.

# organizing paleodata
# the metadata for all this
# https://taphonomy.doit.wisc.edu/data1.2/colls/single_doc.html
# geoplates
# https://www.gplates.org/docs/user-manual/Interacting_Features/#30-assign-plate-ids
# https://www.gplates.org/docs/user-manual/Reconstructions/#plate-ids

# ======================================================

# =====================================================

# load packages
source("R/packages.R")
source("R/functions.R")

# =====================================================


# create a folder for processed data
dir.create ("processed_data_TrJ")


# selected taxa
sel_taxa <- c("cynodontia","rhynchosauria", "pseudosuchia","dicynodontia","dinosauria")

# urls
# make a list of URLs
url_occ <- lapply (sel_taxa, function (i) 
  
  paste ("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",
         i,
         "&interval=Induan,Rhaetian&cc=AR,BR&pgm=gplates,scotese,seton&show=genus,subgenus,abund,coll,loc,paleoloc,env,acconly",
         sep="") # only Olenekian in the studied region
)

# measurements
url_measurement <- lapply (sel_taxa, function (i) 
  
  paste ("https://paleobiodb.org/data1.2/specs/measurements.csv?base_name=",
         i,"&interval=Induan,Maastrichtian&pgm=gplates,scotese,seton",
         sep="")
)


# specimens
url_spec <- lapply (sel_taxa, function (i) 
  
  paste ("https://paleobiodb.org/data1.2/specs/list.csv?base_name=",
         i,"&interval=Induan,Maastrichtian&pgm=gplates,scotese,seton&show=genus,env,acconly",
         sep="")
)

# load all at once
taxa_occurrence <- lapply (url_occ, function (i) 
  
  read.csv(file = i, sep=",")
  
)
taxa_occurrence <- do.call(rbind, taxa_occurrence)
taxa_occurrence <- taxa_occurrence%>%
  filter(genus != "")

# need to filter the occurrences
# here https://psmits.github.io/paleo_book/managing-and-processing-data-from-the-paleobiology-database.html

# load all at once
taxa_measurement <- lapply (url_measurement, function (i) 
  
  read.csv(file = i, sep=",")
  
)
taxa_measurement <- do.call(rbind, taxa_measurement)


taxa_measurement$measurement_type == unique(taxa_measurement$measurement_type)[6]

# load all at once
taxa_specimens <- lapply (seq(1,length(url_spec)), function (i) 
  
  read.csv(file = url_spec[[i]], sep=",") %>%
    cbind(taxon=sel_taxa[[i]])
  
)
taxa_specimens <- do.call(rbind, taxa_specimens)


# match speciments and measurements
taxa_specimens_measurements <- cbind (taxa_measurement,
                                      taxa_specimens [match (taxa_measurement$specimen_no,
                                                             taxa_specimens$specimen_no),])

# taxa_specimens_measurements
taxa_specimens_measurements_genus <- taxa_specimens_measurements %>%
  select(measurement_type,average,genus,taxon,environment) %>% # select cols
  mutate(average = as.numeric(average)) %>%
  group_by(genus, taxon, measurement_type) %>%
  summarise (val = mean(average,na.rm=T),
             env=getmode(environment)) %>%
  filter(genus != "") %>%
  spread(measurement_type, val) %>%
  mutate (env = as.numeric(as.factor(env)))
#%>%
  #mutate(env = recode(env, "0" = "All minus seafood",
  #                         "1" = "Seafood"
  #))
  


# aggregate at genus level
# Find temporal range of all genera
range_occ <- tax_range_time(
  taxa_occurrence,
  name = "genus",
  min_ma = "min_ma",
  max_ma = "max_ma",
  by = "FAD",
  plot = T
)

# binning time intervals
bins <- time_bins(interval = c("Permian", "Cretaceous"), rank = "stage")
bins$epoch <- NA
bins$epoch[which(bins$interval_name %in% "Changhsingian")] <- "Late Permian"
bins$epoch[which(bins$interval_name %in% c("Induan", "Olenekian"))] <- "Early Triassic"
bins$epoch[which(bins$interval_name %in% c("Anisian", "Ladinian"))] <- "Middle Triassic"
bins$epoch[which(bins$interval_name %in% c("Carnian", "Norian", "Rhaetian"))] <- "Late Triassic"

# Generate time bins
range_bins <- bin_time(occdf = taxa_occurrence,
                                    bins = bins,
                                    method = 'all')

# bins in data
bins<- bins[which(bins$bin %in% range_bins$bin_assignment),]
# match to have the epoch
range_bins <- cbind (range_bins, 
                     epoch = bins$epoch [match (range_bins$bin_assignment,
                                                bins$bin)]
)

# composition per bin
bin_composition <- lapply (unique(bins$epoch), function (i) 
  
   range_bins [which(range_bins$epoch %in% i),]
   
)

# -------------------------------------------------------
# ecospace analysis

# gawdis  
#install.packages("gawdis", lib = "C:/Program Files/R/R-4.1.2/library")
# help here: https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html

# standardize traits
require(vegan)
require(cluster)
require(ade4)
require(gawdis)


taxa_specimens_measurements_scaled<-data.frame (
  apply (taxa_specimens_measurements_genus[,-c(1:2)],2,scale)
  )

# remove rows with only NAs
range((rowSums (is.na (taxa_specimens_measurements_scaled)!= T)))
(colSums (is.na (taxa_specimens_measurements_scaled)!=T))

# distance matrix 
gower_matrix <- gawdis (taxa_specimens_measurements_scaled,
                        w.type="analytic") 

# principal coordinate analysis
# Building the functional space based on a PCOA 
pco<-dudi.pco(quasieuclid(gower_matrix), 
              scannf=F, 
              nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 

#barplot(pco$eig) # barplot of eigenvalues for each axis 
(Inertia2<-(pco$eig[1]+pco$eig[2]+pco$eig[3]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes

# estimate quality of f space
#quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
#                                       plot="quality_funct_space_I") # it will produce a plot (hosted in the root folder)



## only the frst axis
(Inertia.first <- (pco$eig[1]) /(sum(pco$eig)))
## only the frst axis
(Inertia.scnd <- (pco$eig[2]) /(sum(pco$eig)))
## only the frst axis
(Inertia.trd <- (pco$eig[3]) /(sum(pco$eig)))
Inertia.first+Inertia.scnd

# extracting scores
PCoA_scores <- pco$li
## bind taxon and species name
naxes_choose <- 3
PCoA_scores <- data.frame (PCoA_scores[,1:naxes_choose], 
                           Species = taxa_specimens_measurements_genus$genus,
                           Taxon=taxa_specimens_measurements_genus$taxon) 


## transforming SCORES to the long format
PCoA_scores_long <- melt(PCoA_scores, id.vars=c("Species", "Taxon"))

# and then to the wide format
require(reshape2)
PCoA_scores_wide <- dcast(data=PCoA_scores_long, 
                          formula=Species+Taxon ~ variable, 
                          value.var="value",
                          fun.aggregate = mean)

## extract the loadings (corretion of traits with each axis)
PCoA_loadings <- lapply (seq (1,naxes_choose), function (i) 
  
  cor ((taxa_specimens_measurements_scaled),PCoA_scores[,i],
       use = "pairwise.complete.obs")
  
)
PCoA_loadings <- do.call(cbind,PCoA_loadings)# melt    
PCoA_loadings<- round (PCoA_loadings,3)# round
colnames(PCoA_loadings)<- colnames(PCoA_scores[,1:naxes_choose])

# fortify (long format)    
# show the correlation of each trait with each axis
PCoA_loadings_wide <- fortify(data.frame (PCoA_loadings[,1:naxes_choose]))

# kernel density estimation for each period 
require(ks);require(hypervolume);require(vegan)

# load kernel density function cl
source("R/functions_kernel.R")

# apply to all periods
fator_de_correcao_arrow<-0.5 # fator pra corrigir o comprimento da seta

list_taxa <- list (unique(PCoA_scores$Taxon),
                   unique(PCoA_scores$Taxon)[1],
                   unique(PCoA_scores$Taxon)[2],
                   unique(PCoA_scores$Taxon)[3])

# obtain the trait space for each taxon
kde_taxon <- lapply (list_taxa, function (i) 
    
                lapply (bin_composition, function (k) {
  
  # subsetting 
  PCoA_scores_wide <- PCoA_scores[which(PCoA_scores$Taxon == i
                                        
                                        &
                                          
                                          PCoA_scores$Species %in% unique(k$genus)
                                          
                                        
                                        ),]
  
  # optimal bandwidth estimation
  hpi_mi_d1 <- Hpi(x = PCoA_scores[,c("A1","A2")])
  
  # kernel density estimation
  est_mi_d1 <- kde(x = PCoA_scores[,c("A1","A2")], 
                   H = hpi_mi_d1, 
                   compute.cont = TRUE)  
  
  # bandwidths for each point
  den_mi_d1 <- list(est_mi_d1$eval.points[[1]], est_mi_d1$eval.points[[2]], 
                    est_mi_d1$estimate)
  names(den_mi_d1) <- c("x", "y", "z")
  dimnames(den_mi_d1$z) <- list(den_mi_d1$x, den_mi_d1$y)
  dcc_mi_d1 <- melt(den_mi_d1$z)
  
  # 0.5 probability kernel
  # run kernel
  cl_50_mi_d1 <- cl(df = den_mi_d1, prob = 0.50)
  # 0.95 probability kernel
  cl_95_mi_d1 <- cl(df = den_mi_d1, prob = 0.95)
  # 0.99 probability kernel
  cl_99_mi_d1 <- cl(df = den_mi_d1, prob = 0.99)
  
  ## PCoA
  # colour palette
  col_pal <- colorRampPalette(c("red", "cyan", "white"))(100)
  
  # plot PCoA
  PCoA_plot_mi_d1 <- ggplot(dcc_mi_d1, aes(x = Var1, y = Var2)) +
    
    # coloured probabilty background
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colours = rev(col_pal)) +
    
    # points for species
    geom_point(data = PCoA_scores_wide, 
               aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
    
    # probability kernels
    geom_contour(aes(z = value), breaks = cl_50_mi_d1, colour = "grey30", size = 1) +
    geom_contour(aes(z = value), breaks = cl_95_mi_d1, colour = "grey60", size = 1) +
    geom_contour(aes(z = value), breaks = cl_99_mi_d1, colour = "grey70", size = 1) +
    coord_equal() +
    theme_classic()+
    
    # add arrows
    geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                        xend = A1, 
                                                                                        yend = A2), 
                 arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
    
    # add dashed arrows ends
    geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                                xend = -A1, 
                                                                                                yend = -A1), 
                 lty = 5, colour = "darkgrey") +
    # add arrow labels
    geom_text(data = PCoA_loadings_wide*fator_de_correcao_arrow, aes(x = A1, y = A2, 
                                                                     label = rownames(PCoA_loadings_wide)),
              size = 4, nudge_x = c(0, 0, 0, 0, -0.2), nudge_y = c(0.2, -0.2, -0.2, -0.2, 0.2)) +
    # axis labels - see comp_var
    #labs(x = paste ("Axis1 (",round(exp_axis[1],2),"%)",sep=""), 
    #     y = paste ("Axis2 (",round(exp_axis[2],2),"%)",sep="")) +
    labs (title=i)+
    # edit plot
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black",size =5),
          axis.title = element_text(colour = "black",size =7),
          legend.position = "none"#,
          #text = element_text(size = 20)
    )
  PCoA_plot_mi_d1
  
  ## adding the name of some species to check if this analysis makes sense
  ## affinity between species and axes (based on score value)
  range_of_values <- apply (PCoA_scores_wide [,c("A1", "A2")],2,range)
  
  # select
  species_to_project <- rbind (
    PCoA_scores_wide [which(PCoA_scores_wide [,"A1"] == range_of_values[1]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A1"] == range_of_values[2]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A2"] == range_of_values[3]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A2"] == range_of_values[4]),]
  )
  
  PCoA_plot_mi_d1 <- PCoA_plot_mi_d1 + 
    geom_text(data = species_to_project, aes(x = A1, y = A2, 
                                             label = Species),
              size = 3, nudge_x = c(0, 0, 0, 0, -0.2), nudge_y = c(0.2, -0.2, -0.2, -0.2, 0.2)) 
  
  
  # list of results
  res <- list (species_to_project = species_to_project,
               plot = PCoA_plot_mi_d1,
               density50 = cl_50_mi_d1,
               density95 = cl_95_mi_d1,
               density99 = cl_99_mi_d1,
               hpi_mi_d1 = hpi_mi_d1)
  
  ; # return
  # display plot
  res
})
)



# whole  trait space area 
all_pool <- data.frame (cbind (pco$li[,1:naxes_choose],
                   ext = F,
                   sp = (taxa_specimens_measurements_genus$genus),
                   taxon=taxa_specimens_measurements_genus$taxon))

# transform into numbers
all_pool<- all_pool%>%
  mutate (A1=as.numeric(A1),
          A2=as.numeric(A2),
          A3=as.numeric(A3))


# convex hull
a_pool <- data.frame (all_pool [chull(all_pool[,1:2], y = NULL),]) # its convex hull

# dinosauria
dinosauria <- all_pool [which (all_pool$taxon == "dinosauria"),]
chull_dinosauria <- dinosauria[chull(dinosauria[,1:2], y = NULL),] # its convex hull

# dicynodontia
dicynodontia <- all_pool [which (all_pool$taxon == "dicynodontia"),]
chull_dicynodontia <- dicynodontia[chull(dicynodontia[,1:2], y = NULL),] # its convex hull

# pseudosuchia
pseudosuchia <- all_pool [which (all_pool$taxon == "pseudosuchia"),]
chull_pseudosuchia <- pseudosuchia[chull(pseudosuchia[,1:2], y = NULL),] # its convex hull

# cynodontia
cynodontia <- all_pool [which (all_pool$taxon == "cynodontia"),]
chull_cynodontia <- cynodontia[chull(cynodontia[,1:2], y = NULL),] # its convex hull


# ========================================================
# trait spaces

# plots

# plot one particular period
require(gridExtra)

grid.arrange (kde_taxon[[1]][[1]]$plot + geom_polygon(data=a_pool, aes (A1,A2),
                                                 alpha=0.1,
                                                 fill="gray",
                                                 colour = "black",
                                                 size=1,
                                                 linetype = 2) ,
              kde_taxon[[2]][[1]]$plot + geom_polygon(data=chull_dinosauria, 
                                               aes (x=A1, y=A2),
                                               alpha=0.3,
                                               fill="green",
                                               colour = "green",
                                               size=1,
                                               linetype = 3),
              kde_taxon[[3]][[1]]$plot + geom_polygon(data=chull_cynodontia, 
                                                 aes (x=A1, y=A2),
                                                 alpha=0.3,
                                                 fill="yellow",
                                                 colour = "yellow",
                                                 size=1,
                                                 linetype = 3),
              kde_taxon[[4]][[1]]$plot + geom_polygon(data=chull_pseudosuchia, 
                                                 aes (x=A1, y=A2),
                                                 alpha=0.3,
                                                 fill="red",
                                                 colour = "red",
                                                 size=1,
                                                 linetype = 3),
              
              
              
              kde_taxon[[1]][[2]]$plot + geom_polygon(data=a_pool, aes (A1,A2),
                                                      alpha=0.1,
                                                      fill="gray",
                                                      colour = "black",
                                                      size=1,
                                                      linetype = 2) ,
              kde_taxon[[2]][[2]]$plot + geom_polygon(data=chull_dinosauria, 
                                                      aes (x=A1, y=A2),
                                                      alpha=0.3,
                                                      fill="green",
                                                      colour = "green",
                                                      size=1,
                                                      linetype = 3),
              kde_taxon[[3]][[2]]$plot + geom_polygon(data=chull_cynodontia, 
                                                      aes (x=A1, y=A2),
                                                      alpha=0.3,
                                                      fill="yellow",
                                                      colour = "yellow",
                                                      size=1,
                                                      linetype = 3),
              kde_taxon[[4]][[2]]$plot + geom_polygon(data=chull_pseudosuchia, 
                                                      aes (x=A1, y=A2),
                                                      alpha=0.3,
                                                      fill="red",
                                                      colour = "red",
                                                      size=1,
                                                      linetype = 3),
              
               
              nrow=4,ncol=2)


## plot A (complete space)
plotA <- ggplot(all_pool , 
                aes(x=A1, y=A2)) + 
  geom_point(size=2) + theme_bw()+
  geom_polygon(data=a_pool, aes (A1,A2),
               alpha=0.6,
               fill="gray",
               colour = "black",
               size=1,
               linetype = 2) + # complete space
  geom_polygon(data=chull_cynodontia, 
               aes (x=A1, y=A2),
               alpha=0.3,
               fill="yellow",
               colour = "yellow",
               size=1,
               linetype = 3) +
  geom_polygon(data=chull_dinosauria, 
               aes (x=A1, y=A2),
               alpha=0.3,
               fill="green",
               colour = "green",
               size=1,
               linetype = 3) + 
  geom_polygon(data=chull_pseudosuchia, 
               aes (x=A1, y=A2),
               alpha=0.3,
               fill="red",
               colour = "red",
               size=1,
               linetype = 3)



plotA + ggrepel::geom_text_repel(data = a_pool, 
                        aes (x=A1, y=A2, 
                             label=(sp)),
                        size=3)


kde_taxon[[1]]$plot 

  











head(taxa_measurement)
head(taxa_specimens)
head(taxa_occurrence)








# filter taxa
taxa_occ_filter <-   taxa_occurrence %>%
  janitor::clean_names() %>%           # standardizes names
  filter(accepted_rank %in% c('genus', 'species'))




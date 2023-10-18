
# =====================================================

# download data

# load packages
source("R/packages.R")
source("R/functions.R")

# =====================================================
# load data
# collections
#pbdb_data_collections <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_collections.csv"),
#                                   skip=18)
url_col <- "https://paleobiodb.org/data1.2/colls/list.csv?base_name=cynodontia&max_ma=259.5&min_ma=66"
pbdb_data_collections <- read.csv (url_col, sep=",")

# strata
#pbdb_data_strata <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_strata.csv"),
#                              skip=17)
url_strata <- "https://paleobiodb.org/data1.2/occs/strata.csv?base_name=cynodontia&max_ma=259.5&min_ma=66&pgm=gplates,scotese,seton"
pbdb_data_strata <-  read.csv (url_strata, sep=",")

#taxa
#pbdb_data_taxa <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_taxa.csv"),
#                            skip=15)
url_taxa <- "https://paleobiodb.org/data1.2/occs/taxa.csv?base_name=cynodontia&max_ma=259.5&min_ma=66"
pbdb_data_taxa <- read.csv(url_taxa,sep=",")

# occ
#pbdb_data_occ <- read.csv (here ("data", "PaleoData","cynodontia", "pbdb_data_occ.csv"),
#                           skip=18)
url_occ <- "https://paleobiodb.org/data1.2/occs/list.csv?base_name=cynodontia&max_ma=259.5&min_ma=66&pgm=gplates,scotese,seton&show=genus,paleoloc,class"
pbdb_data_occ <- read.csv(url_occ,sep=",")


# match collection and occurrence (lat long of collection)
pbdb_data_occ <- cbind (pbdb_data_occ,
                        pbdb_data_collections [match (pbdb_data_occ$collection_no,
                                                      pbdb_data_collections$collection_no),
                                               c("lat","lng")])


# save the original download
save (pbdb_data_collections,
      pbdb_data_strata,
      pbdb_data_taxa,
      pbdb_data_occ,
      file=here ("processed_data", "PBDB_download.RData"))

# save also the csvs
write.csv (pbdb_data_collections, file = here ("processed_data", "pbdb_data_collections.csv"))
write.csv (pbdb_data_strata, file = here ("processed_data", "pbdb_data_strata.csv"))
write.csv (pbdb_data_taxa, file = here ("processed_data", "pbdb_data_taxa.csv"))
write.csv (pbdb_data_occ, file = here ("processed_data", "pbdb_data_occ.csv"))

rm(list=ls())

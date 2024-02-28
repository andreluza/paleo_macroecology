# load packages
library(here);library(rgeos);library(rgdal);library(sp);library(raster); library(openxlsx)
# data processing
library(ggplot2);library(gridExtra); library(rasterVis); library(viridis); library(dplyr); library(magick); library(reshape)
library (purrr); library(janitor); library(knitr); library(kableExtra);library(tidyverse)

library(CoordinateCleaner)
# functions
library(palaeoverse)
library(deeptime)

# hierarchical models
library(jagsUI)
library(R2WinBUGS) 
#devtools::install_github("mikemeredith/saveJAGS")
require(saveJAGS)
require(bayestestR)

# load paleoDEM dataset
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(sp)
library(rgeos)
library(mapast)
library(sf)

# world map
require(rnaturalearth)
require(ggplot2)
require(ggrepel)
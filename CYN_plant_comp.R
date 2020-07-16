########### Cherkskiy North NDVI vs Plant Composition ###########

getwd()
setwd("/Users/elenaforbath/Downloads/loranty_lab/CYN_plant_comp_data")

## install packages 
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tiff")
install.packages("rtiff")
install.packages("raster") 
install.packages("sp")
install.packages("rgdal")

library(dplyr)
library(ggplot2)
library(tiff)
library(rtiff)
library(raster)
library(sp)
library(rgdal)


########## Extracted NDVI Values #########
##projecting rasters
FL016 <- raster("CYN_TR1_FL016M/RU_CYN_TR1_FL016B_index_ndvi.tif")
FL016b <- projectRaster(FL016, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                        method = "bilinear", 
                        alignOnly = FALSE)


FL020 <- raster("CYN_TR1_FL020M/NDVI.data.tif")
FL020b <- projectRaster(FL020, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                        method = "bilinear", 
                        alignOnly = FALSE)

##reading GPS coordinates
GPS <- na.omit(read.csv("CYN_plot_centers.csv"))
GPS2 <- subset(GPS, select= -c(plot, elevation))

## re order columns so longitude is first
GPS_order <- GPS2[,c("longitude", "latitude")]

## turn into spatial points
GPS_order2 <- SpatialPoints(GPS_order,
                            proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

## actually extract NDVI values from each flight 
ndvi_FL016 <- extract(FL016b, GPS_order2,
                      buffer = 0.25,
                      fun = mean,
                      df = TRUE, 
                      along = TRUE, 
                      sp = TRUE)
names(ndvi_FL016)[names(ndvi_FL016) == "RU_CYN_TR1_FL016B_index_ndvi"] <- "FL016_ndvi"
ndvi_FL016b <- as.data.frame(ndvi_FL016)


ndvi_FL020 <- extract(FL020b, 
                      GPS_order2,
                      buffer = 0.25,
                      fun = mean,
                      df = TRUE,
                      along = TRUE,
                      sp = TRUE)
names(ndvi_FL020)[names(ndvi_FL020) == "NDVI.data"] <- "FL020_ndvi"
ndvi_FL020b <- as.data.frame(ndvi_FL020)

## combine pre- and post-clipping data frames by lat and long
ndvi <- merge(ndvi_FL016b, ndvi_FL020b, by = c("longitude", "latitude")) ## just pre- and post-

ndvi_plots <- merge(ndvi, GPS, by = c("longitude", "latitude"))## add GPS points


## rearrage column orders for ease of reading
ndvi_plots2 <- ndvi_plots[,c("plot", "longitude", "latitude", "elevation", "FL016_ndvi", "FL020_ndvi")]


########## Read/Merge Treatments to Table ###########
treatments <- read.csv("treatments.csv")
treatments 

ndvi <- merge(ndvi_plots2, treatments, by = "plot")
ndvi

## plot by treatment (separate plots)
ndvi_CT <- subset(ndvi, treatment == "CT")
ndvi_GR <- subset(ndvi, treatment == "GR")
ndvi_SH <- subset(ndvi, treatment == "SH")
ndvi_GS <- subset(ndvi, treatment == "GS")


######### Read/Merge Percent Cover to Table ############
percent_cover <- na.omit(read.csv("percent_cover.csv"))
names(percent_cover)[names(percent_cover) == "Plot.ID"] <- "plot"

## aggregate data by functional group
percent_cover2 <- aggregate(percent_cover$percent.cover, 
                            by = list(percent_cover$plot, percent_cover$Functional.group), 
                            FUN = sum)
names(percent_cover2)[names(percent_cover2) == "Group.1"] <- "plot"
names(percent_cover2)[names(percent_cover2) == "Group.2"] <- "functional_group"
names(percent_cover2)[names(percent_cover2) == "x"] <- "percent_cover"


## subset by treatment 
pc_GR <- subset(percent_cover, Treatment == "GR")
pc_SH <- subset(percent_cover, Treatment == "SH")
pc_GS <- subset(percent_cover, Treatment == "G+S")

ndvi_GR$plot = gsub("P", "",ndvi_GR$plot)
ndvi_GS$plot = gsub("P", "",ndvi_GS$plot)
ndvi_SH$plot = gsub("P", "",ndvi_SH$plot)

## merge ndvi data with pc data
ndviGR <- merge(ndvi_GR, pc_GR, by = c("plot")) 
ndviSH <- merge(ndvi_SH, pc_SH, by = c("plot"))
ndviGS <- merge(ndvi_GS, pc_GS, by = c("plot"))

con <- subset(percent_cover, Functional.group == "CON")
evsh <- subset(percent_cover, Functional.group == "EVSH")
desh <- subset(percent_cover, Functional.group == "DESH")
gram <- subset(percent_cover, Functional.group == "GRAM")
forb <- subset(percent_cover, Functional.group == "FORB")
cwd <- subset(percent_cover, Functional.group == "CWD")
moss <- subset(percent_cover, Functional.group == "MOSS")
lichen <- subset(percent_cover, Functional.group == "LICH")
brg <- subset(percent_cover, Functional.group == "BRG")
litr <- subset(percent_cover, Functional.group == "LITR")
equ <- subset(percent_cover, Functional.group == "EQU")


### does change in NDVI correlate with veg type (ie does NDVI inc/dec with specific veg types)
pc_br <- merge(percent_cover, ndvi, by = c("plot")) ## merging percent cover with NDVI






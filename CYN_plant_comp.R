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
install.packages("tidyr")
install.packages("stringr")
install.packages("lme4")
install.packages("cowplot")
install.packages("lsmeans")
install.packages("DHARMa")

library(tiff)
library(rtiff)
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(stringr)
library(lme4)
library(cowplot)
library(lme4)
library(lsmeans)
library(DHARMa)



########## Extracted NDVI Values #########
##reading rasters (created in CYN_NDVI)
FL016 <- raster("CYN_TR1_FL016M/FL016.tif")  
writeRaster(FL016, "FL016.tif")
FL016

FL020 <- raster("CYN_TR1_FL020M/FL020.tif")
writeRaster(FL020, "FL020.tif")
FL020



##reading GPS coordinates
GPS <- na.omit(read.csv("CYN_plot_centers.csv"))
GPS$ID <- c(1:29) ## add ID number so it can be merged later with plot numbers

GPS2 <- subset(GPS, select= -c(plot, elevation))


## re order columns so longitude is first
GPS_order <- GPS2[,c("longitude", "latitude")]


## turn into spatial points
GPS_order2 <- SpatialPoints(GPS_order,
                            proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
CRS.new <- CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=56 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
GPS_final<- spTransform(GPS_order2, CRS.new)



## actually extract NDVI values from each flight 
## will not work with tidyr package; had to uninstall
ndvi_FL016 <- extract(FL016, GPS_final,
                      buffer = 0.25,
                      fun = mean,
                      df = TRUE,
                      along = TRUE, 
                      sp = TRUE)
names(ndvi_FL016)[names(ndvi_FL016) == "FL016"] <- "FL016_ndvi"
ndvi_FL016b <- as.data.frame(ndvi_FL016)


ndvi_FL020 <- extract(FL020, 
                      GPS_final,
                      buffer = 0.25,
                      fun = mean,
                      df = TRUE,
                      along = TRUE,
                      sp = TRUE)
names(ndvi_FL020)[names(ndvi_FL020) == "FL020"] <- "FL020_ndvi"
ndvi_FL020b <- as.data.frame(ndvi_FL020)

## combine pre- and post-clipping data frames by lat and long
ndvi <- merge(ndvi_FL016b, ndvi_FL020b, by = c("longitude", "latitude")) ## just pre- and post-
ndvi$ID <- c(1:29)

ndvi_plots <- merge(ndvi, GPS, by = c("ID"))## add GPS points
ndvi_plots2 <- subset(ndvi_plots, select= -c(longitude.y, latitude.y))

## rearrage column orders for ease of reading
ndvi_plots2 <- ndvi_plots2[,c("plot", "longitude.x", "latitude.x", "elevation", "FL016_ndvi", "FL020_ndvi")]

##### export GPS coordinates with plots #####
GPS_all <- ndvi_plots2[,c("plot", "longitude.x", "latitude.x", "elevation")]
write.csv(GPS_all, "GPS_all.csv")


########## Read/Merge Treatments to Table ###########
treatments <- read.csv("treatments.csv")
treatments 

ndvi <- merge(ndvi_plots2, treatments, by = "plot")
ndvi$plot = gsub("P", "", ndvi$plot)
ndvi


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


ndvi_pc <- merge(ndvi, percent_cover2, by = "plot")  ## data frame with ndvi and percent cover 



######### Read in Point Intercept Data #########
pt_int<- na.omit(read.csv("pt_intercept.csv"))

pt_int2<- merge(pt_int, ndvi, by = "plot")  ## data frame with ndvi and point intercept
pt_int2 <- pt_int2[ ,c("plot", "longitude.x", "latitude.x", "elevation",
                      "FL016_ndvi", "FL020_ndvi", "treatment", "functional_groups", 
                      "sum_hits", "percent_composition")]

write.csv(pt_int2, "plant_comp.csv", row.names = FALSE)

plant_comp <- subset(pt_int2, select= -c(longitude.x, latitude.x, elevation,
                                          FL020_ndvi, treatment, sum_hits)) 


## pivot dataframe (?) and create correlation ##
### correlation plots 
install.packages("ggcorrplot")
install.packages("GGally")
install.packages("tidyr")

library(tidyr)
library(ggcorrplot)
library(GGally)




veg = plant_comp %>% 
  dplyr::select(plot:percent_composition) %>% 
  pivot_wider(names_from = functional_groups, 
              values_from = percent_composition,
              values_fill = 0)

veg

names(veg)


veg_var = veg %>%
  dplyr::select(plot, FL016_ndvi, CON:FORB)
veg_var


corr <- round(cor(veg_var), 1)
corr

ggcorrplot(corr)




## subset percent cover by functional group
con <- subset(ndvi_pc, functional_group == "CON")
evsh <- subset(ndvi_pc, functional_group == "EVSH")
desh <- subset(ndvi_pc, functional_group == "DESH")
gram <- subset(ndvi_pc, functional_group == "GRAM")
forb <- subset(ndvi_pc, functional_group == "FORB")
cwd <- subset(ndvi_pc, functional_group == "CWD")
moss <- subset(ndvi_pc, functional_group == "MOSS")
lichen <- subset(ndvi_pc, functional_group == "LICH")
brg <- subset(ndvi_pc, functional_group == "BRG")
litr <- subset(ndvi_pc, functional_group == "LITR")
equ <- subset(ndvi_pc, functional_group == "EQU")



### box plot of NDVIs per treatment ###
boxplot(FL016_ndvi ~ treatment, 
        data = ndvi_pc, 
        main = "NDVI vs Treatment", 
        xlab = "Treatment", 
        ylab = "NDVI")


### box plot with just distribution of NDVI

Palette <- c("darkred","red","orange","yellow","lightgreen","green","blue",
                "lightblue","purple", "pink","brown","black")

Palette2 <-c("red", "orange", "yellow", "green")  

ggplot(data = pt_int2, aes(x ="" , y = FL016_ndvi, color = factor(""))) +
  geom_boxplot() +       
  geom_point() +
  scale_color_manual(values=Palette) +
  ylab("NDVI") +
  xlab("Cherskiy North")  +
  theme_bw() + theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(t = 0.3, r = 0.1, b = 0.3, l = 0.1), "cm")) +
  theme(axis.title.x = element_text(size = 11, hjust = 0.5, vjust = -0.1),
        axis.title.y = element_text(size = 11, hjust = 0.5, vjust = 1.1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA), title = "Transect ID"))



### boxplot of NDVi per treatment....looks a little funky
ggplot(data = pt_int2, aes(x = "", y = FL016_ndvi, color = factor(treatment))) +
  geom_boxplot() +       
  geom_point() +
  scale_color_manual(values=Palette2) +
  ylab("NDVI") +
  xlab("Cherskiy North")  +
  theme_bw() + theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(t = 0.3, r = 0.1, b = 0.3, l = 0.1), "cm")) +
  theme(axis.title.x = element_text(size = 11, hjust = 0.5, vjust = -0.1),
        axis.title.y = element_text(size = 11, hjust = 0.5, vjust = 1.1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA), title = "Transect ID"))




####### extracting all values from each plots (distribution of values) ######

## have to convert plot coordinates to shp file!

library(sf)
plot_locations <- st_as_sf(GPS_all, coords = c("longitude.x", "latitude.x"), crs = "+proj=aea +lat_1=50 +lat_2=70 +lat_0=56 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
st_crs(plot_locations)

ggplot() +
  geom_sf(data = plot_locations) +
  ggtitle("Map of Plot Locations")


st_write(plot_locations,
         "plot_locations.shp", driver = "ESRI Shapefile")

points<-shapefile("plot_locations.shp")


detach("package:tidyr", unload = TRUE) # extract wont ork with tidyr

ndvi_FL016_all <- extract(FL016, points,
                      buffer = 0.5, 
                      small = TRUE,
                      df = TRUE, 
                      factor = TRUE)
ndvi_FL016_all
ndvi_FL016_all <- as.data.frame(ndvi_FL016_all) 

all <- merge(ndvi_FL016_all, ndvi_plots, by = "ID")
all <- subset(all, select= -c(longitude.y, latitude.y, elevation, FL016_ndvi, FL020_ndvi))
all$plot = gsub("P", "", all$plot)


## for loop to create histogram of NDVI values for each plot 

sub <- split(all, all$plot) ## split data by plot (don't need this!!)
list2env(sub, envir= .GlobalEnv) ##separate into dataframes 


all_sub <- subset(all, select= -c(longitude.x, latitude.x, ID))
write.csv(all_sub, "ndvi_dist.csv", row.names = FALSE)

new <- read.csv("ndvi_dist_edit.csv") ## created new data frame bc I didn't know how to do it in r
names(new)


for (i in 1:length(new)) { # for every column in the "new" data frame
  x <- new[,i] # identifying columns (?)
  # Plot histogram of x
  jpeg(file = paste("dist", names((new)[i]), ".jpeg", sep = ""))
  hist(x,
       main = paste("NDVI distribution for", names((new)[i])), #paste name of column to the 
       xlab = "NDVI",
       xlim = c(0.2, 0.8),
       ylim = c(1, 60))
  dev.off()
}

## standard deviation for all plots?
FL016 = new  %>%
  summarise_all(sd, na.rm = TRUE)

boxplot(new, 
        main = "FL016 NDVI Distributions for All Plots",
        xlab = "Plot", 
        ylab = "NDVI", 
        las = 2)




.rs.restartR() ### too many graphs created; needed to restart r






###### REPEAT for FL020 ######
ndvi_FL020_all <- extract(FL020, points,
                          buffer = 0.5, 
                          small = TRUE,
                          df = TRUE, 
                          factor = TRUE)
ndvi_FL020_all
ndvi_FL020_all <- as.data.frame(ndvi_FL020_all) 

all <- merge(ndvi_FL020_all, ndvi_plots, by = "ID")
all <- subset(all, select= -c(longitude.y, latitude.y, elevation, FL016_ndvi, FL020_ndvi))
all$plot = gsub("P", "", all$plot)


## for loop to create histogram of NDVI values for each plot 

sub <- split(all, all$plot) ## split data by plot (don't need this??)
list2env(sub, envir= .GlobalEnv) ##separate into dataframes 


all_sub <- subset(all, select= -c(longitude.x, latitude.x, ID))
write.csv(all_sub, "ndvi_dist_post.csv", row.names = FALSE)

new2 <- read.csv("ndvi_dist_post_edit.csv") ## created new data frame bc I didn't know how to do it in r
names(new)


for (i in 1:length(new)) { # for every column in the "new" data frame
  y <- new[,i] # identifying columns (?)
  # Plot histogram of x
  jpeg(file = paste("FL020dist", names((new)[i]), ".jpeg", sep = ""))
  hist(y,
       main = paste("FL020 NDVI distribution for", names((new)[i])), #paste name of column to the 
       xlab = "NDVI",
       xlim = c(0.2, 0.8),
       ylim = c(1, 150))
  dev.off()
}

## standard deviation for all plots?
FL020 = new2  %>%
  summarise_all(sd, na.rm = TRUE)


FL016
FL020


boxplot(new2,
        main = "FL020 NDVI Distributions for All Plots",
        xlab = "Plot", 
        ylab = "NDVI", 
        las = 2)


##### Analyses #####
lm <- lm(FL016_ndvi ~ treatment, data = ndvi)
summary(lm)

lm <- lm(FL020_ndvi ~ treatment, data = ndvi)
summary(lm)

ndvi$ndvi_diff <- (ndvi$FL020_ndvi - ndvi$FL016_ndvi)
lm <- lm(ndvi_diff ~ treatment, data = ndvi)
summary(lm)






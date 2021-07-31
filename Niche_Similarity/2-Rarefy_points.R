#### MODIFIED FROM JERONYMO DALAPICOLLA, 2019 (https://github.com/jdalapicolla/) ########
#########################################################################################

#install packages
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("dismo" %in% rownames(installed.packages()) == FALSE){install.packages("dismo")
} else {print (paste0("'dismo' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("raster", "dismo")

##########################################################################################
################################# INPUT FILES ############################################
# 1- Enviromental layer 
# 2- Ocurrence points 

##########################################################################################
################################# REFERENCE LAYER ########################################
################################# Atlantic Forest ########################################

# select the directory with the bioclim layers to use as a reference

setwd ("~/AF_cut/")
layer_AF = stack(sapply(list.files(pattern='asc$', recursive = F), raster))

##################################### Pampas #############################################

setwd ("~/PA_cut/")
layer_PA = stack(sapply(list.files(pattern='asc$', recursive = F), raster))

##########################################################################################
################################ OCURRENCE DATA ##########################################

# select the directory with the occurence points
setwd(" ")

#### Occurrence data

points_raw <- read.table("occ_raw.txt", header=TRUE, sep ="\t")

# inspect the values of the file
head(points_raw)
# if you have more that long/lat columns remove all other columns
occ <- points_raw[,5:6]
head(occ)

# remove duplicated data based on latitude and longitude
dups <- duplicated(occ[c("Lat", "Long")])
occ_unique <- occ[!dups, ]
cat(nrow(occ) - nrow(occ_unique), "records are removed")

# make occ spatial
coordinates(occ_unique) <- ~Long + Lat
crs(occ_unique) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


################################# Atlantic Forest ########################################

#identify and selecting points within AF
ovr_AF = extract(layer_AF[[1]], occ_unique)
head(ovr_AF)
i_AF = which(ovr_AF != "-9999") #replace the value if the NoData value is different of -9999
i_AF #lines in the point_raw

#update the points raw and create a SpatialPoints for ocorrence points.
points_AF = points_raw[i_AF,]
occ_AF = points_AF[,5:6]
coordinates(occ_AF) = ~Long + Lat #name of the columns.
crs(occ_AF) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Eliminating sampling bias:
#create a buffer of 10Km around the points
buffer_AF = circles(occ_AF, d = 10000, lonlat=TRUE) #d is radius of each circle in meters
plot(buffer_AF)
class(buffer_AF)

#convert circles in polygon
buffer_AF = polygons(buffer_AF)
#rasterize the buffer as the layer
buffer_AF= rasterize(buffer_AF, layer_AF[[1]])

#select 1 point in each circle
sel_AF = gridSample(occ_AF, buffer_AF, n=1)

# Set working directory to save final results
setwd(" ")

## Save the results in a cvs file 
sel_AF = as.data.frame(sel_AF)

write.csv(sel_AF, "OCC_UNBIAS_AF.csv", row.names = FALSE)

##################################### Pampas #############################################

#identify and selecting points within PA
ovr_PA = extract(layer_PA[[1]], occ_unique)
head(ovr_PA)
i_PA = which(ovr_PA != "-9999") #replace the value if the NoData value is different of -9999
i_PA #lines in the point_raw

#update the points raw and create a SpatialPoints for ocorrence points.
points_PA = points_raw[i_PA,]
occ_PA = points_PA[,5:6]
coordinates(occ_PA) = ~Long + Lat #name of the columns.
crs(occ_PA) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Eliminating sampling bias:
#create a buffer of 10Km around the points
buffer_PA = circles(occ_PA, d = 10000, lonlat=TRUE) #d is radius of each circle in meters
plot(buffer_PA)
class(buffer_PA)

#convert circles in polygon
buffer_PA = polygons(buffer_PA)
#rasterize the buffer as the layer
buffer_PA= rasterize(buffer_PA, layer_PA[[1]])

#select 1 point in each circle
sel_PA = gridSample(occ_PA, buffer_PA, n=1)

# Set working directory to save final results
setwd(" ")

## Save the results in a cvs file 
sel_PA = as.data.frame(sel_PA)

write.csv(sel_PA, "OCC_UNBIAS_PA.csv", row.names = FALSE)

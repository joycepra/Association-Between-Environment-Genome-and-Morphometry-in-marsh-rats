#install packages
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

pacman::p_load("raster")

##########################################################################################
################################# INPUT FILES ############################################
# 1- Maxent Models 
# 2- Unbiased ocurrence points 

##########################################################################################
############################### OCCURRENCE POINTS ########################################
################################# Atlantic Forest ########################################

# select the directory with the unbiased occurence points
setwd(" ")

#unbiased points
occ_AF = read.csv("OCC_UNBIAS_AF.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ_AF = occ_AF[,2:3]
head(occ_AF)

##################################### Pampas #############################################

# select the directory with the unbiased occurence points
setwd(" ")

#unbiased points
occ_PA = read.csv("OCC_UNBIAS_PA.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ_PA = occ_PA[,2:3]
head(occ_PA)

##########################################################################################
################################# MAXENT MODELS ##########################################
################################# Atlantic Forest ########################################

# select the directory with the maxent results for Atlantic Forest
setwd(" ")

# Load the average model from the maxent results folder
layer_AF<- raster("model_avg.asc")
plot(layer_AF)

thresh_AF= 0.575 ### need to be changed acording to MAXENT results

## CREATE A BINARY MAP

#reclassify the model

pxReclass_AF = reclassify(layer_AF, rcl = c(0, thresh_AF, NA, thresh_AF, 1, 1))
plot(pxReclass_AF)

#save binary maps
dir.create("Binary_Maps")
setwd ("./Binary_Maps")

writeRaster(pxReclass_AF, "model_final_bin_AF.asc", driver='ascii', overwrite=TRUE, bylayer=T)


#Cut the models by threshold 

dir.create("Final_Models")
setwd ("./Final_Models")

cut_AF = mask(layer_AF, pxReclass_AF)
plot(cut_AF)

writeRaster(cut_AF, "Final_cut_AF.asc", driver='ascii', overwrite=TRUE, bylayer=T)

##################################### Pampas #############################################

# select the directory with the maxent results for Pampas
setwd(" ")

# Load the average model from the maxent results folder
layer_PA<- raster("model_avg.asc")
plot(layer_PA)

thresh_PA= 0.6286 ### need to be changed acording to MAXENT results

## CREATE A BINARY MAP

#reclassify the model

pxReclass_PA = reclassify(layer_PA, rcl = c(0, thresh_PA, NA, thresh_PA, 1, 1))
plot(pxReclass_PA)

#save binary maps
setwd ("./Binary_Maps")

writeRaster(pxReclass_PA, "model_final_bin_PA.asc", driver='ascii', overwrite=TRUE, bylayer=T)


#Cut the models by threshold 
setwd ("./Final_Models")

cut_PA = mask(layer_PA, pxReclass_PA)
plot(cut_PA)

writeRaster(cut_PA, "Final_cut_PA.asc", driver='ascii', overwrite=TRUE, bylayer=T)

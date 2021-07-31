#install packages
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("ENMeval" %in% rownames(installed.packages()) == FALSE){install.packages("ENMeval")
} else {print (paste0("'ENMeval' has already been installed in library"))}
if("rJava" %in% rownames(installed.packages()) == FALSE){install.packages("rJava")
} else {print (paste0("'rJava' has already been installed in library"))}
if("dismo" %in% rownames(installed.packages()) == FALSE){install.packages("dismo")
} else {print (paste0("'dismo' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("dplyr" %in% rownames(installed.packages()) == FALSE){install.packages("dplyr")
} else {print (paste0("'dplyr' has already been installed in library"))}

#load packages
pacman::p_load("ENMeval", "rJava", "raster", "dismo", "dplyr")


# maxent.jar need to be in the folder R/3.4/dismo/java/. To check it: 
jar = paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar) & require(rJava)){
  print("maxent.jar is in the correct place")
}else{
  print("maxent.jar file needs to be moved to dismo/java folder in R")
}

##########################################################################################
################################# INPUT FILES ############################################
# 1- Enviromental layer 
# 2- Unbiased ocurrence points 

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
####################### CHOOSING THE FEATURE CLASS AND BETA VALUES #######################
################################# Atlantic Forest ########################################

###create background points and save as a RDS object
bg_pnts_AF <- randomPoints(layer_AF, n = 1000)
colnames(bg_pnts_AF)[1:2] <- c("Long", "Lat")
saveRDS(bg_pnts_AF, file = "background_points_AF")
#bg_pnts_AF <- readRDS("background_points_AF")
head(bg_pnts_AF)

#testing the parameters

EVALme_AF <- ENMevaluate(occ_AF, layer_AF, bg.coords = bg_pnts_AF, 
                      RMvalues = seq(0.5, 3, 0.5), 
                      fc = c("L", "Q", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                      method = "block", clamp = TRUE, algorithm='maxent.jar',
                      numCores = 4)


# select the directory to save the results
setwd(" ")
dir.create("ENMeval") 
setwd ("./ENMeval")

#save the results
saveRDS(EVALme_AF, file = paste("Results_eval_","AF", ".rds", sep= ""))

#read the RDS file
#EVALme_AF = readRDS(paste("Results_eval_", ".rds", sep= ""))

#Best parameters
res_AF = arrange(EVALme_AF@results, EVALme_AF@results$AICc)
best_res_AF = res_AF[1,c(1:2, 16:17)]
best_res_AF

##################################### Pampas #############################################

###create background points and save as a RDS object
bg_pnts_PA <- randomPoints(layer_PA, n = 1000)
colnames(bg_pnts_PA)[1:2] <- c("Long", "Lat")
saveRDS(bg_pnts_PA, file = "background_points_PA")
#bg_pnts_PA <- readRDS("background_points_PA")
head(bg_pnts_PA)

#testing the parameters

EVALme_PA <- ENMevaluate(occ_PA, layer_PA, bg.coords = bg_pnts_PA, 
                      RMvalues = seq(0.5, 3, 0.5), 
                      fc = c("L", "Q", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                      method = "block", clamp = TRUE, algorithm='maxent.jar',
                      numCores = 4)


# select the directory to save the results

setwd ("./ENMeval")

#save the results
saveRDS(EVALme_PA, file = paste("Results_eval_","PA", ".rds", sep= ""))

#read the RDS file
#EVALme_PA = readRDS(paste("Results_eval_", biome, ".rds", sep= ""))

#Best parameters
res_PA = arrange(EVALme_PA@results, EVALme_PA@results$AICc)
best_res_PA = res_PA[1,c(1:2, 16:17)]
best_res_PA
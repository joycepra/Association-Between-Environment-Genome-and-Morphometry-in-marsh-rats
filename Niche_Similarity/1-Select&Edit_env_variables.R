#### MODIFIED FROM JERONYMO DALAPICOLLA, 2019 (https://github.com/jdalapicolla/) ########
#########################################################################################

#install packages
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("raster", "vegan")

##########################################################################################
################################# INPUT FILES ############################################
# 1- Environmental layers 
# 2- Mask for the study area 


##########################################################################################
##################################### MASK ###############################################


# Set working directory with the mask for the study area 
setwd(" ")

#load the mask for the study area
mask_AF<- shapefile("mask_AF.shp") #Atlantic forest
crs(mask_AF) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(mask_AF)

mask_PA<- shapefile("mask_PA.shp") #Pampas
crs(mask_PA) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(mask_PA)

##########################################################################################
################################ ENVIRONMENTAL LAYERS ####################################
################################# Atlantic Forest ########################################
###cut bioclimatic layers by mask 

# Set working directory 
setwd(" ")
dir.create("AF_cut")

#load bioclimatic layers
# Set working directory with the layers from the present
setwd(" ")
layers <- list.files(, pattern = '.tif')

outpath_AF <- "..../AF_cut"

outfiles_AF <- paste0(outpath_AF, layers)

for(i in 1:length(layers)) {
  e <- extent(mask_AF)
  r <-raster(layers[i])
  rc <- crop(r, e)
  rc <- extend(rc,mask_AF)
  rc<- mask(rc, mask_AF)
  rc <- writeRaster(rc, outfiles_AF[i], format = "ascii", overwrite=TRUE)
}

### re-scale variable when nedeed after the clip
# Set working directory where the cut variables are
setwd ("./AF_cut/")

layers <- list.files(, pattern='asc', full.names=TRUE )

r1<- raster(layers[12]) #choose variable to be re-scaled
crs(r1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
r2<- raster(layers[13]) #choose one variable in the correct scale
crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

r1resampled <- projectRaster(r1,r2,method = 'ngb')

writeRaster(r1resampled, "layer_resc.asc")

###################################### Pampas ############################################

###cut bioclimatic layers by mask 

# Set working directory 
setwd(" ")
dir.create("PA_cut")

#load bioclimatic layers
# Set working directory with the layers from the present
setwd(" ")
layers <- list.files(, pattern = '.tif')

outpath_PA <- "..../PA_cut"

outfiles_PA <- paste0(outpath_PA, layers)

for(i in 1:length(layers)) {
  e <- extent(mask_PA)
  r <-raster(layers[i])
  rc <- crop(r, e)
  rc <- extend(rc,mask_AF)
  rc<- mask(rc, mask_AF)
  rc <- writeRaster(rc, outfiles_PA[i], format = "ascii", overwrite=TRUE)
}

### re-scale variable when nedeed after the clip
# Set working directory where the cut variables are
setwd ("./PA_cut/")

layers <- list.files(, pattern='asc', full.names=TRUE )

r1<- raster(layers[12]) #choose variable to be re-scaled
crs(r1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
r2<- raster(layers[13]) #choose one variable in the correct scale
crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

r1resampled <- projectRaster(r1,r2,method = 'ngb')

writeRaster(r1resampled, "layer_resc.asc")

##########################################################################################
############################ SELECTING VARIABLES FOR ENMS ################################
####################################### PCA ##############################################
################################# Atlantic Forest ########################################

# Set working directory where the cut variables are
setwd ("~/AF_cut/")

# Calculating the Correlation among layers 
#load all layers
layers_AF = stack(sapply(list.files(pattern='asc$', recursive = F), raster))

# PCA of the layers

raw_values_AF = values(layers_AF)
raw_values_AF = na.omit(raw_values_AF) #remove the NAs
raw_values_AF = as.data.frame(raw_values_AF)
head(raw_values_AF)
summary(raw_values_AF) 

#standardize data 
values_trans_AF = decostand(raw_values_AF, method="standardize")
summary(values_trans_AF)

#PCA
pca_AF = prcomp (values_trans_AF)

#% PC
std_list_AF = pca_AF$sdev^2
std_list_pct_AF = std_list_AF / sum (std_list_AF) * 100
std_list_pct_sum_AF = round(std_list_pct_AF, 2)
write.csv(std_list_pct_sum_AF, file = "PC_variance_AF.csv")

contr_AF= as.data.frame(pca_AF$rotation)
write.csv(contr_AF, file = "variables_PCA_AF.csv")

# Look the number of PCs that explains 90% of the variation

n.pc_AF = length(std_list_pct_sum_AF)
sum_pca_AF = matrix(NA, n.pc_AF, 2)
for (i in 1:n.pc_AF) {sum_pca_AF[i,] = c(i, sum(std_list_pct_sum_AF[1:i]))}
colnames(sum_pca_AF) = c("PCs","% Variance")
sum_pca_AF
pc_AF = 5 ## CHANGE HERE DEPENDING OF THE RESULT
write.csv(sum_pca_AF, file = "PC_variance_sum_AF.csv", row.names = F)


## WHITIN THE PCs CHOSEN ABOVE, CHOOSE VARIABLES WHICH VALUES ARE ABOVE the fixed value of 
#0.32 (Based on Dormann et al., 2012 and Wang et al. 2016)


tab_AF = abs(contr_AF) 
lista_AF = list()
for (i in 1:pc_AF)
{
  linhas_AF=tab_AF[tab_AF[, i] > 0.32, ]
  lista_AF[[i]]=row.names(linhas_AF)
}

linhas.resultado_AF=unlist(lista_AF)
linhas.resultado_AF=unique(linhas.resultado_AF)
write.table(linhas.resultado_AF, file = "Variables_10%_PCA_AF.csv", sep = "\n", row.names = F, col.names = F)

###################################### Pampas ############################################

# Set working directory where the cut variables are
setwd ("~/PA_cut/")

# Calculating the Correlation among layers 
#load all layers
layers_PA = stack(sapply(list.files(pattern='asc$', recursive = F), raster))

# PCA of the layers

raw_values_PA = values(layers_PA)
raw_values_PA = na.omit(raw_values_PA) #remove the NAs
raw_values_PA = as.data.frame(raw_values_PA)
head(raw_values_PA)
summary(raw_values_PA) 

#standardize data 
values_trans_PA = decostand(raw_values_PA, method="standardize")
summary(values_trans_PA)

#PCA
pca_PA = prcomp (values_trans_PA)

#% PC
std_list_PA = pca_PA$sdev^2
std_list_pct_PA = std_list_PA / sum (std_list_PA) * 100
std_list_pct_sum_PA = round(std_list_pct_PA, 2)
write.csv(std_list_pct_sum_PA, file = "PC_variance_PA.csv")

contr_PA= as.data.frame(pca_PA$rotation)
write.csv(contr_PA, file = "variables_PCA_PA.csv")

# Look the number of PCs that explains 90% of the variation

n.pc_PA = length(std_list_pct_sum_PA)
sum_pca_PA = matrix(NA, n.pc_PA, 2)
for (i in 1:n.pc_PA) {sum_pca_PA[i,] = c(i, sum(std_list_pct_sum_PA[1:i]))}
colnames(sum_pca_PA) = c("PCs","% Variance")
sum_pca_PA
pc_PA = 5 ## change here depending os the result
write.csv(sum_pca_PA, file = "PC_variance_sum_PA.csv", row.names = F)


## WHITIN THE PCs CHOSEN ABOVE, CHOOSE VARIABLES WHICH VALUES ARE ABOVE the fixed value of 
#0.32 (Based on Dormann et al., 2012 and Wang et al. 2016)


tab_PA = abs(contr_PA) 
lista_PA = list()
for (i in 1:pc_PA)
{
  linhas_PA=tab_PA[tab_PA[, i] > 0.32, ]
  lista_PA[[i]]=row.names(linhas_PA)
}

linhas.resultado_PA=unlist(lista_PA)
linhas.resultado_PA=unique(linhas.resultado_PA)
write.table(linhas.resultado_PA, file = "Variables_10%_PCA_PA.csv", sep = "\n", row.names = F, col.names = F)

################################# CORRELATION ############################################
################################# Atlantic Forest ########################################

#Correlation analysis to select only the variables with less than 0.7 of correlation 

#load the variabled selected above (linhas.resultado_AF and linhas.resultado_PA)

setwd ("~/AF_cut/")
##Correlations among layers:
#load all layers
predictors_AF = stack(".asc",".asc", ...) #choose here according to the results (linhas.resultado_AF)


#calculete correlation:
correlation_AF = layerStats(predictors_AF, "pearson", na.rm = T)
correlation_AF = correlation_AF$`pearson correlation coefficient`
correlation_AF = round(correlation_AF,2)

write.csv(correlation_AF, file = "correlation_AF.csv")

###################################### Pampas ############################################

setwd ("~/PA_cut/")
##Correlations among layers:
#load all layers
predictors_PA = stack(".asc",".asc", ...) #choose here according to the results (linhas.resultado_PA)


#calculete correlation:
correlation_PA = layerStats(predictors_PA, "pearson", na.rm = T)
correlation_PA = correlation_PA$`pearson correlation coefficient`
correlation_PA = round(correlation_PA,2)

write.csv(correlation_PA, file = "correlation_PA.csv")


#install packages
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("RStoolbox" %in% rownames(installed.packages()) == FALSE){install.packages("RStoolbox")
} else {print (paste0("'RStoolbox' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("dismo" %in% rownames(installed.packages()) == FALSE){install.packages("dismo")
} else {print (paste0("'dismo' has already been installed in library"))}
if("hypervolume" %in% rownames(installed.packages()) == FALSE){install.packages("hypervolume")
} else {print (paste0("'hypervolume' has already been installed in library"))}
if("plot3D" %in% rownames(installed.packages()) == FALSE){install.packages("plot3D")
} else {print (paste0("'plot3D' has already been installed in library"))}
if("ecospat" %in% rownames(installed.packages()) == FALSE){install.packages("ecospat")
} else {print (paste0("'ecospat' has already been installed in library"))}

pacman::p_load("RStoolbox", "ggplot2", "raster", "dismo", "hypervolume", "plot3D", "ecospat")

##########################################################################################
################################# INPUT FILES ############################################
# 1- Enviromental layers selected and edited to build the ENM 
# 2- Binary maps from Maxent 
# 3- Unbiased occurence points

##########################################################################################
############################### ENVIRONMENTAL LAYER ######################################
##################################### Pampas #############################################

# select the directory with the bioclim layers selected to build the ENM
setwd(" ")

## Load the environmental layers
layers_PA = stack(sapply(list.files(pattern='asc$', recursive = F), raster)) 
crs(layers_PA) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers_PA

################################# Atlantic Forest ########################################

# select the directory with the bioclim layers selected to build the ENM
setwd(" ")

## Load the environmental layers
layers_AF = stack(sapply(list.files(pattern='asc$', recursive = F), raster)) 
crs(layers_AF) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers_AF

##########################################################################################
################################ ENVIRONMENTAL PCA #######################################

# select the directory for the hypervolume results
dir.create("Niche_Similarity")
setwd ("./Niche_Similarity")

##################################### Pampas #############################################



## Perform a PCA with the environmental data 
PCA_climate_WA_PA = rasterPCA(layers_PA, nSamples = NULL, nComp = 19, spca = TRUE, maskCheck = TRUE)

#plot the map with the 3 first PC
map_PA<- ggRGB(PCA_climate_WA_PA$map, 1, 2, 3, stretch="lin", q=0) + labs( x="Longitude", y = "Latitude") + theme_bw()

ggsave(filename="map_PA_dif.tiff", plot=map_PA, width=20, height=15, unit="cm", dpi=200)

#create a new object for the PCA map
map_PCA_WA_PA = PCA_climate_WA_PA$map
summary(PCA_climate_WA_PA$model) #choose the number of PCs that explain 90% of variance or more

#save the results
#save the PCA map in raster
writeRaster(PCA_climate_WA_PA$map, "PCA_WA_19PC_PA_RASTER.tif", driver='GTiff', bylayer=F)

#save the % of variance of each PC
std_dev_PA = PCA_climate_WA_PA$model$sdev
pr_var_PA = std_dev_PA^2
prop_varex_PA = pr_var_PA/sum(pr_var_PA)
prop_varex_PA = prop_varex_PA*100 # must be the same values than summary(PCA_climate_WA_PA$model)
write.csv(prop_varex_PA, file = "contribution_pc_eig_WA_PA.csv")

#save the variables contribution to the PC's
write.csv(PCA_climate_WA_PA$model$loadings, "contribuition_variables_WA_PA.csv")

################################# Atlantic Forest ########################################

## Perform a PCA with the environmental data 
PCA_climate_WA_AF = rasterPCA(layers_AF, nSamples = NULL, nComp = 19, spca = TRUE, maskCheck = TRUE)

#plot the map with the 3 first PC
map_AF<- ggRGB(PCA_climate_WA_AF$map, 1, 2, 3, stretch="lin", q=0) + labs( x="Longitude", y = "Latitude") + theme_bw()

ggsave(filename="map_AF_dif.tiff", plot=map_AF, width=20, height=15, unit="cm", dpi=200)

#create a new object for the PCA map
map_PCA_WA_AF = PCA_climate_WA_AF$map
summary(PCA_climate_WA_AF$model) #choose the number of PCs that explain 90% of variance or more

#save the results
#save the PCA map in raster
writeRaster(PCA_climate_WA_AF$map, "PCA_WA_19PC_AF_RASTER.tif", driver='GTiff', bylayer=F)

#save the % of variance of each PC
std_dev_AF = PCA_climate_WA_AF$model$sdev
pr_var_AF = std_dev_AF^2
prop_varex_AF = pr_var_AF/sum(pr_var_AF)
prop_varex_AF = prop_varex_AF*100 # must be the same values than summary(PCA_climate_WA_AF$model)
write.csv(prop_varex_AF, file = "contribution_pc_eig_WA_AF.csv")

#save the variables contribution to the PC's
write.csv(PCA_climate_WA_AF$model$loadings, "contribuition_variables_WA_AF.csv")

##########################################################################################
################################ MAXENT BINARY MAPS ######################################
##################################### Pampas #############################################

# select the directory with the maxent binary maps
setwd(" ")

# Load binary maps  
binary_PA = raster("model_final_bin_PA.asc")

#plot to verify
plot(binary_PA)

# Creating 1,000 random points inside the area of binary maps
randomPoints_PA = as.data.frame(randomPoints(binary_PA, 1000))
plot(randomPoints_PA, cex = 0.5)
colnames(randomPoints_PA) = c("Long", "Lat")
write.csv(randomPoints_PA, "randomPoints_PA.csv", row.names = F)

################################# Atlantic Forest ########################################

# select the directory with the maxent binary maps
setwd(" ")

# Load binary maps  
binary_AF = raster("model_final_bin_AF.asc")

#plot to verify
plot(binary_AF)

# Creating 1,000 random points inside the area of binary maps
randomPoints_AF = as.data.frame(randomPoints(binary_AF, 1000))
plot(randomPoints_AF, cex = 0.5)
colnames(randomPoints_AF) = c("Long", "Lat")
write.csv(randomPoints_AF, "randomPoints_AF.csv", row.names = F)

##########################################################################################
################################### HYPERVOLUME ##########################################
##################################### Pampas #############################################

#load random points saved if you need to:
#randomPoints_PA=read.csv("randomPoints_PA.csv")

# Creating the hypervolume with the number of pc equal or more than 90% of variation
#extracting values of variables for the Random Points  
values_PA = extract(map_PCA_WA_PA, randomPoints_PA)

#Use only the first 5 PC for analysis: ~0.945% of variation
values_PA=as.data.frame(values_PA[,1:5]) #change the number for the number of PC you need to use!

#hypervolume calculation 
hv_PA = hypervolume_gaussian(values_PA, name = "PAMPA", 
                              weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_PA))))/nrow(values_PA)),
                              kde.bandwidth = estimate_bandwidth(values_PA), 
                              sd.count = 3, 
                              quantile.requested = 0.95,
                              quantile.requested.type = "probability", 
                              chunk.size = 1000,
                              verbose = TRUE)

################################# Atlantic Forest ########################################

#load random points saved if you need to:
#randomPoints_AF=read.csv("randomPoints_AF.csv")

# Creating the hypervolume with the number of pc equal or more than 90% of variation
#extracting values of variables for the Random Points  
values_AF = extract(map_PCA_WA_AF, randomPoints_AF)

#Use only the first 5 PC for analysis: ~0.945% of variation
values_AF=as.data.frame(values_AF[,1:5]) #change the number for the number of PC you need to use!

#hypervolume calculation 
hv_AF = hypervolume_gaussian(values_AF, name = "ATLANTIC FOREST", 
                              weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_AF))))/nrow(values_AF)),
                              kde.bandwidth = estimate_bandwidth(values_AF), 
                              sd.count = 3, 
                              quantile.requested = 0.95,
                              quantile.requested.type = "probability", 
                              chunk.size = 1000,
                              verbose = TRUE)


############################# Hypervolume Comparisons ####################################

hv_wa = hypervolume_join(hv_PA, hv_AF)
centroid = get_centroid(hv_wa)
vol = get_volume(hv_wa)

write.csv(centroid, "centroid_WA.csv")
write.csv(vol, "volume_WA.csv")

##Comparing hypervolumes among different target goups

set_groups = hypervolume_set(hv_PA, hv_AF, check.memory=FALSE)
groups = hypervolume_overlap_statistics(set_groups)

varimp_FA = hypervolume_variable_importance(hv_PA, verbose = TRUE)
barplot(varimp_FA)
varimp_PA = hypervolume_variable_importance(hv_AF, verbose = TRUE)
barplot(varimp_PA)

# Save Figure

setEPS()
postscript("./Environmental_hypervolumes_WA.eps")
pdf("./Environmental_hypervolumes_WA.pdf", onefile = F ) 
plot.new()

plot(hv_wa,
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=F, show.density=F,show.data=F,
     names=NULL, show.legend=F, limits=NULL,
     show.contour=T, contour.lwd=2,
     contour.type="kde",
     contour.kde.level=0.01,
     show.centroid=TRUE, cex.centroid=2, pt.bg.centroid="black",
     colors= c("#77ab59", "#E7B800"), #choose the colors
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1.5,cex.axis=0.9,cex.names=2,cex.legend=2,
     num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE)

#create a legend box for your graph
legend =c("Atlantic Forest", "Pampa")
col = c("#77ab59", "#E7B800")
legend("bottomleft", legend = legend, cex = 1, bty="n",  col = "white", pt.bg = col, pch= 21, pt.cex=1.5)

dev.off()

# figure with the first three PCs

PCA_PA=as.data.frame(hv_PA@Data)
PCA_AF=as.data.frame(hv_AF@Data)

PCA_PA$Type <- "PA"
PCA_AF$Type <- "AF"

dat3D <- rbind(PCA_PA,PCA_AF)
dat3D$Type <- as.factor(dat3D$Type)

setEPS()
postscript("./Environmental_hypervolumes_WA_3D.eps")

pdf("./Environmental_hypervolumes_WA_3D.pdf", onefile = F ) #change the name of pdf
plot.new()

scatter3D(dat3D, x = dat3D$PC1, y = dat3D$PC2, z = dat3D$PC3, box=TRUE, pch=16, bty="b2", axes=TRUE, label=TRUE,
          xlab='PC1', ylab='PC2', zlab='PC3', colvar = as.integer(dat3D$Type), col = c("#77ab59", "#E7B800"), 
          theta = 60, phi = 0)

dev.off()

##########################################################################################
############################### OCCURRENCE POINTS ########################################
##################################### Pampas #############################################

# select the directory with the unbiased occurence points from Pampas
setwd(" ")

#unbiased points
occ_PA = read.csv("OCC_UNBIAS_PA.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ_PA = occ_PA[,2:3]
head(occ_PA)

extract_points_PA<-extract(map_PCA_WA_PA, occ_PA)
points_PA<- data.frame(occ_PA,extract_points_PA)



################################# Atlantic Forest ########################################

# select the directory with the unbiased occurence points from Atlantic Forest
setwd(" ")

#unbiased points
occ_AF = read.csv("OCC_UNBIAS_AF.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ_AF = occ_AF[,2:3]
head(occ_AF)

extract_points_AF<-extract(map_PCA_WA_AF, occ_AF)
points_AF<- data.frame(occ_AF,extract_points_AF)


##########################################################################################
##################################### PCA-ENV ############################################


extract_AF = na.omit(cbind(points_AF, rep(1, nrow(points_AF))))
extract_PA = na.omit(cbind(points_PA, rep(1, nrow(points_PA))))

colnames(extract_AF)[ncol(extract_AF)] = 'occ'
colnames(extract_PA)[ncol(extract_PA)] = 'occ'

extbgAF = na.omit(cbind(randomPoints_AF, extract(map_PCA_WA_AF, randomPoints_AF), rep(0, nrow(randomPoints_AF))))
extbgPA = na.omit(cbind(randomPoints_PA, extract(map_PCA_WA_PA, randomPoints_PA), rep(0, nrow(randomPoints_PA))))

colnames(extbgAF)[ncol(extbgAF)] = 'occ'
colnames(extbgPA)[ncol(extbgPA)] = 'occ'

colnames(extbgAF)[1] <- "Long"
colnames(extbgAF)[2] <- "Lat"
colnames(extbgPA)[1] <- "Long"
colnames(extbgPA)[2] <- "Lat"

dat_AF = rbind(extract_AF, extbgAF)
dat_PA = rbind(extract_PA, extbgPA)

#PCA-env

pca.env <- dudi.pca(
  rbind(dat_AF, dat_PA )[,3:7],
  scannf=FALSE,
  nf=2
)

#Variable contribution
pdf("var_contribution.pdf")
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
dev.off()

scores.globclim<-pca.env$li # PCA scores for the whole study area 

scores.AF <- suprow(pca.env,
                    extract_AF[which(extract_AF[,8]==1),3:7])$li # PCA scores for the Atlantic Forest biome

scores.PA <- suprow(pca.env,
                    extract_PA[which(extract_PA[,8]==1),3:7])$li # PCA scores for the Pampas biome

scores.clim_AF <- suprow(pca.env,dat_AF[,3:7])$li # PCA scores for the whole Atlantic Forest study area

scores.clim_PA <- suprow(pca.env,dat_PA[,3:7])$li # PCA scores for the whole Pampas study area

# gridding the Atlantic Forest niche

grid.clim_AF <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim_AF,
  sp = scores.AF,
  R = 100,
  th.sp = 0
)

# gridding the Pampas niche

grid.clim_PA <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim_PA,
  sp = scores.PA,
  R = 100,
  th.sp = 0
)

# Compute Schoener's D, index of niche overlap

D.overlap <- ecospat.niche.overlap (grid.clim_AF, grid.clim_PA, cor=T)$D
D.overlap

# Perform the Niche Equivalency Test 

eq.test<-ecospat.niche.equivalency.test(grid.clim_AF, grid.clim_PA, rep=2000, alternative = "greater", ncores = 4)
eq.test
eq.test$obs$D

pdf("equivalency.pdf")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

# Perform the Niche Similarity Test 

sim.test_AF_PA<-ecospat.niche.similarity.test(grid.clim_AF, grid.clim_PA, rep=2000, alternative = "greater", rand.type = 1, ncores = 4)
sim.test_AF_PA

pdf("Similarity_AF_PA.pdf")
ecospat.plot.overlap.test(sim.test_AF_PA, "D", "Similarity Atlantic Rainforest-> Pampa")
dev.off()

sim.test_PA_AF<-ecospat.niche.similarity.test(grid.clim_PA, grid.clim_AF, rep=2000, alternative = "greater", rand.type = 1, ncores = 4)
sim.test_PA_AF

pdf("Similarity_PA_FA.pdf")
ecospat.plot.overlap.test(sim.test_PA_AF, "D", "Similarity Pampa-> Atlantic Rainforest")
dev.off()

## Plot Niche

pdf("Niche_AF.pdf")
ecospat.plot.niche(grid.clim_AF, "Atlantic Rainforest", "PC1", "PC2", cor=FALSE)
dev.off()

pdf("Niche_PA.pdf")
ecospat.plot.niche(grid.clim_PA,"Pampas", "PC1", "PC2", cor=FALSE)
dev.off()



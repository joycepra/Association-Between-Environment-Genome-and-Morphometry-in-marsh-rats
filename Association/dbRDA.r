#################################################################################################################

#install packages
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("fossil" %in% rownames(installed.packages()) == FALSE){install.packages("fossil")
} else {print (paste0("'fossil' has already been installed in library"))}
if("ecodist" %in% rownames(installed.packages()) == FALSE){install.packages("ecodist")
} else {print (paste0("'ecodist' has already been installed in library"))}
if("RStoolbox" %in% rownames(installed.packages()) == FALSE){install.packages("RStoolbox")
} else {print (paste0("'RStoolbox' has already been installed in library"))}
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("vegan", "ecodist", "RStoolbox", "adegenet", "fossil", "ecodist", "raster")


#################################################################################################################
############################################# INPUT FILES #######################################################
# 1- structure file with the genomic data 
# 2- cvs file with the geographic coordinates of each individual
# 3- Skull and dental measurements
# 4- Environmental layers used in the ENMs 

#################################################################################################################
############################################ Geographic Distance ################################################
setwd("~/") #set working directory

coord_gen <- read.csv("coord_ind_gen.csv", header = T)
coord_gen <- coord_gen[,2:3] #select only the Long/Lat column

mDIST_gen <- earth.dist(coord_gen, dist = TRUE)
mDIST_gen <- as.dist(mDIST_gen)
mDIST_vec_gen <- as.vector(mDIST_gen)
hist(mDIST_vec_gen)

pcnm_gen<-pcnm(mDIST_gen, dist.ret = FALSE)

mDIST_gen2<-as.matrix(dist(mDIST_gen))
write.csv(mDIST_gen2, "geo_gen_dist.csv") # save the geographic distance of the genetic points

coord_mor <- read.csv("coord_ind_mor.csv", header = T)
coord_mor <- coord_mor[,2:3] #select only the Long/Lat column

mDIST_mor <- earth.dist(coord_mor, dist = TRUE)
mDIST_mor <- as.dist(mDIST_mor)
mDIST_vec_mor <- as.vector(mDIST_mor)
hist(mDIST_vec_mor)
mDIST_mor2<-as.matrix(dist(mDIST_mor))

pcnm_mor<-pcnm(mDIST_mor2, dist.ret = FALSE)

write.csv(mDIST_mor2, "geo_mor_dist.csv") # save the geographic distance of the morphometric points
write.csv(pcnm_mor$vectors, "geo_mor_pcnm.csv")

#################################################################################################################
############################################## Genomic Distance #################################################

setwd("~/") #set working directory with the genetic data

str <- read.structure('brasiliensis_final.str', n.ind = 20, n.loc = 24264, col.lab = 1, 
				col.pop = 2, onerowperind = FALSE, row.marknames = 0, NA.char= 0, 
				ask = FALSE)
str2 <- scaleGen (str, cent=TRUE,scale=TRUE,NA.method = "mean")

pca_gen <- prcomp(str2,center = TRUE,scale. =TRUE)
summary(pca_gen)

PCs_gen<- pca_gen$x[,1:11] ### ~70% of the variation

PCA_gen_dist<-ecodist::distance(PCs_gen, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_gen_dist2<-as.matrix(dist(PCA_gen_dist))

write.csv(PCA_gen_dist2, "PCA_gen_dist.csv") # save the PCA genetic distance

#################################################################################################################
########################################## Morphometric Distance ################################################

setwd(" ") # Set working directory with the measurements

data<- read.table("Morpho_H_brasiliensis.csv", header = T, sep= "," )

pca_mor<- prcomp(data[,36:56], center = TRUE, scale. = TRUE) #select the transformed data in log
summary(pca_mor)

PCs_mor<- pca_mor$x[,1:4] ### ~70% of the variation

PCA_mor_dist<-ecodist::distance(PCs_mor, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_mor_dist2<-as.matrix(dist(PCA_mor_dist))

write.csv(PCA_mor_dist2, "PCA_mor_dist.csv") # save the PCA genetic distance


#################################################################################################################
######################################### Environmental Variable ################################################

setwd("~/") #set working directory with the environmental variables for the entire area

predictors <- stack(sapply(list.files(pattern='asc$', recursive = F), raster))
names(predictors)

#Perform the PCA:
pca_env <- rasterPCA(predictors, nComp=1, norm=T, spca=F)
summary(pca_env$model) # PC1 81% of the total variance
pca_env$model$loadings

PC1_gen<-extract(pca_env$map$PC1, coord_gen)
PC1_mor<-extract(pca_env$map$PC1, coord_mor)

write.csv(PC1_gen, "PC1_env_gen.csv") # save the PCA env distance
write.csv(PC1_mor, "PC1_env_mor.csv") # save the PCA env distance

#################################################################################################################
########################################### Redundancy Analysis  ################################################

### Genetic

dbrda<-capscale(PCA_gen_dist2 ~ PC1_gen, dist = "man") 
dbrda
summary(dbrda)
sig<- anova(dbrda)
sig


dbrda2<-capscale(PCA_gen_dist2 ~ pcnm_gen$vectors, dist = "man") 
dbrda2
summary(dbrda2)
sig2<- anova(dbrda2)
sig2


dbrda3<-capscale(PCA_gen_dist2 ~ PC1_gen + Condition(pcnm_gen$vectors), dist = "man") 
dbrda3
summary(dbrda3)
sig3<- anova(dbrda3)
sig3

### Morphometric

dbrda4<-capscale(PCA_mor_dist2 ~ PC1_mor, dist = "man") 
summary(dbrda4)
sig4<- anova(dbrda4)
sig4

dbrda5<-capscale(PCA_mor_dist2 ~ pcnm_mor$vectors, dist = "man") 
dbrda5
summary(dbrda5)
sig5<- anova(dbrda5)
sig5

dbrda6<-capscale(PCA_mor_dist2 ~ PC1_mor + Condition(pcnm_mor$vectors), dist = "man") 
summary(dbrda6)
sig6<- anova(dbrda6)
sig6

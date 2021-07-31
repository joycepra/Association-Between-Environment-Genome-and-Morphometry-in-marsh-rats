#################################################################################################################

#install packages
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("AssocTests" %in% rownames(installed.packages()) == FALSE){install.packages("AssocTests")
} else {print (paste0("'AssocTests' has already been installed in library"))}
if("fossil" %in% rownames(installed.packages()) == FALSE){install.packages("fossil")
} else {print (paste0("'fossil' has already been installed in library"))}
if("ecodist" %in% rownames(installed.packages()) == FALSE){install.packages("ecodist")
} else {print (paste0("'ecodist' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("adegenet", "AssocTests", "fossil", "ecodist")

#################################################################################################################
############################################# INPUT FILES #######################################################
# 1- structure file with the genomic data 
# 2- cvs file with the geographic coordinates of each individual

#################################################################################################################
################################################# PCA ###########################################################

setwd("~/") #set working directory

str <- read.structure('brasiliensis_final.str', n.ind = 20, n.loc = 24264, col.lab = 1, 
				col.pop = 2, onerowperind = FALSE, row.marknames = 0, NA.char= 0, 
				ask = FALSE)
str2 <- scaleGen (str, cent=TRUE,scale=TRUE,NA.method = "mean")

pca <- prcomp(str2,center = TRUE,scale. =TRUE)
summary(pca)

pca$x

write.csv(pca$x, file = "scores_pca.csv") # save pca scores

###Plot PCA 1n2 #### 
tiff("pca1n2_moleculas.tiff", width=20, height=15, unit="cm", res=300)
quartz.options(height=10, width=12, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(5,5,1,1));
plot.window(xlim=c(-100,300), ylim=c(-100, 200));

points(pca$li[14:20,1],pca$li[14:20,2], col = 'slategray4', bg = "#77ab59", cex = 1.5, pch=21) # ATLANTIC FOREST
points(pca$li[1:13,1],pca$li[1:13,2], col = 'slategray4', bg = "#E7B800", cex = 1.5, pch=21) #PAMPA

axis(1, at=seq(-100, 300, by=130), cex.axis=1.15);
axis(2, at=seq(-100, 200, by=100), cex.axis=1.15, las=1);

mtext(side=1, text='PC1 (10.52%)',line=2.5, cex=1)
mtext(side=2, text='PC2 (8.79%)', line=2.8, cex=1)
legend<- c("Pampas", "Atlantic Forest")
col=c("#E7B800", "#77ab59")
legend("topright", legend = legend, cex=1, bty="n",  col = col, pch= 16)
dev.off()

#################################################################################################################
################################### TRACY-WIDOM TEST FOR EIGENVALUES ############################################

#Tracy CA and Widom H. (1994). Level spacing distributions and the bessel kernel. Commun Math Phys. 161 :289--309.
#Patterson N, Price AL and Reich D. (2006). Population structure and eigenanalysis. PLoS Genet. 2 :20.

eigenvalues<-pca$eig

eigenL<- length(eigenvalues)

#criticalpoint: a numeric value corresponding to the significance level. If the significance level is 0.05, 0.01, 
#0.005, or 0.001,the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. The default
# is 2.0234

tw<- tw(eigenvalues, eigenL, criticalpoint = 0.9793)

tw$SigntEigenL #the number of significant eigenvalues


#################################################################################################################
############################################# MANTEL TEST #######################################################

### Genomic distance #####

PCs<- pca$x[,1:11] ### ~70% of the variation

PCA_dist<-ecodist::distance(PCs, method = "mahalanobis", sprange=NULL, spweight=NULL) 

PCA_dist2<-as.matrix(dist(PCA_dist))

write.csv(PCA_dist2, "PCA_gen_dist.csv") # save the PCA genetic distance


### Geographic distance #####

coord <- read.csv("coord_ind_gen.csv", header = T)
coord <- coord[,2:3] #select only the Long/Lat column

mDIST <- earth.dist(coord, dist = TRUE)
mDIST <- as.dist(mDIST)
mDIST_vec <- as.vector(mDIST)
hist(mDIST_vec)

### Run the Mantel test #####

mantel <- mantel.rtest(PCA_dist, mDIST, 10000)

plot(mDIST, PCA_dist, pch = 16,  ylab = "PCA_dist", xlab = "mDIST", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist  ~ mDIST), col = "gray30", lwd = 3)

### Run the Mantel test for each biome ####

### Genomic distance for the Pampas ##### 

pca_pampa <- prcomp(str2[1:13,])
summary(pca_pampa)
PCs_pampa<- pca_pampa$x[,1:8] ### ~70% of the variation

PCA_dist_pampa<-ecodist::distance(PCs_pampa, method = "mahalanobis", sprange=NULL, spweight=NULL) 

### Geographic distance for the Pampas #####

coord <- read.csv("coord_ind.csv", header = T)
coord_pampa <- coord[1:13,2:3]
mDIST_pampa <- earth.dist(coord_pampa, dist = TRUE)
mDIST_pampa <- as.dist(mDIST_pampa )

### Run the Mantel test for the Pampas #####

mantel_pampa <- mantel.rtest(PCA_dist_pampa, mDIST_pampa, 10000)

plot(mDIST_pampa, PCA_dist_pampa,  pch = 16,  ylab = "PCA_dist", xlab = "mDIST", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist_pampa  ~ mDIST_pampa), col = "gray30", lwd = 3)


### Genomic distance for the Atlantic forest ##### 

pca_AF <- prcomp(str2[14:20,])
summary(pca_AF)
PCs_AF<- pca_AF$x[,1:4] ### ~70% of the variation

PCA_dist_FA<-ecodist::distance(PCs_FA, method = "mahalanobis", sprange=NULL, spweight=NULL) 

### Geographic distance for the Atlantic forest #####

coord_FA <- coord[14:20,2:3]
mDIST_FA <- earth.dist(coord_FA, dist = TRUE)
mDIST_FA <- as.dist(mDIST_FA )

### Run the Mantel test for the Atlantic forest #####

mantel_FA <- mantel.rtest(PCA_dist_FA, mDIST_FA, 10000)

plot( mDIST_FA, PCA_dist_FA,  pch = 16,  ylab = "PCA_dist", xlab = "mDIST", cex = 1.5, cex.lab=1.5, axes= T) 
abline(lm(PCA_dist_FA ~ mDIST_FA), col = "gray30", lwd = 3)






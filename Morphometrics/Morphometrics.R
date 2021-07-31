##############################MORPHOMETRIC STRUCTURE#####################################
#########################################################################################

#install packages
if("MASS" %in% rownames(installed.packages()) == FALSE){install.packages("MASS")
} else {print (paste0("'MASS' has already been installed in library"))}
if("QuantPsyc" %in% rownames(installed.packages()) == FALSE){install.packages("QuantPsyc")
} else {print (paste0("'QuantPsyc' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
                      
#load packages
pacman::p_load("MASS","QuantPsyc", "ggplot2")

##########################################################################################
################################# INPUT FILES ############################################
# 1- Skull and dental measurements
# 2- External measurements (body, length of tail, length of ears, length of hindfoot)

# Set working directory with the measurements
setwd(" ")

data<- read.table("Morpho_H_brasiliensis.csv", header = T, sep= "," )


##########################################################################################
################################## Normality #############################################

#Mardia’s Test (multivariate normality)

normality1<- mult.norm(data[,15:35])$mult.test #Skull and dental measurements

#Shapiro-Wilk’s test (univariate normality)
normality_LB<- shapiro.test(data[,11]) #body length
normality_LT<- shapiro.test(data[,12]) #tail length
normality_LH<- shapiro.test(data[,13]) #hindfoot length
normality_LE<- shapiro.test(data[,14]) #ears length

##########################################################################################
########################### Principal Component Analysis #################################

pca<- prcomp(data[,36:56], center = TRUE, scale. = TRUE) #select the transformed data in log

summary(pca)

pca_rot<-pca$rotation
write.table(pca_rot, "scores_pca_var.csv", sep=",") #save PCA variables scores

sc_pca<-data.frame(pca$x)
write.table(sc_pca, "scores_pca.csv", sep=",") #save PCA individuals scores
plot(pca$x)

# make the PCA plot
setEPS()
postscript("pca_morpho_data.eps")
#pdf("pca_morpho_data.pdf")
quartz.options(height=10, width=12, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(3,3,1,1));
plot.window(xlim=c(-8,10), ylim=c(-6,6));


points(pca$x[37:59,1],pca$x[37:59,2],  col = "#E7B800", cex = 1.5, pch=16) #PAMPAS
points(pca$x[1:36,1],pca$x[1:36,2],  col = "#77ab59", cex = 1.5, pch=16) #ATLANTIC FOREST

axis(1, at=seq(-8, 10, by=6), cex.axis=1);
axis(2, at=seq(-6, 6, by=4), cex.axis=1, las=1);

mtext(side=1, text='PC1(46.99%)',line=2.5, cex=1)
mtext(side=2, text='PC2(10.76%)', line=2.8, cex=1)
dev.off()

##########################################################################################
########################## Boxplot with Important Variables ##############################

# select biome column and the columns with the most important variables found in the PCA
box<- data.frame(data[,c(8,28,34,35)]) #LPB, LBB, CZL
fill <- c("#77ab59",  "#E7B800")

LPB<-ggplot(data, aes(x=data[,8], y=data[,28])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of Palatal Bridge") + 
  ylab("mm")+
  theme_classic()
LPB

ggsave("LPB.pdf",  plot =LPB, device = "pdf")

LBB<-ggplot(data, aes(x=data[,8], y=data[,34])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Lambdoidal Breadth") +
  ylab("mm")+
  theme_classic()
LBB

ggsave("LBB.pdf",  plot =LBB, device = "pdf")

CZL<-ggplot(data, aes(x=data[,8], y=data[,35])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of Condyle-Zygomatic") +
  ylab("mm")+
  theme_classic()
CZL

ggsave("CZL.pdf",  plot =CZL, device = "pdf")


##########################################################################################
############################## Discriminant Analysis #####################################

fit_disc <- lda(data[,8] ~ .,data= data[,36:56]) # select the biomes column and data in log

loadings_disc<-(fit_disc$scaling) #DFA variables scores

scores_disc <-predict(fit_disc)$x #DFA individuals scores

box_disc<- data.frame(data[,8],scores_disc)

fill <- c("#77ab59",  "#E7B800")

disc_plot<-ggplot(box_disc, aes(x=box_disc$data...8., y=LD1)) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("") +
  theme_classic()
disc_plot

ggsave("DFA.pdf",  plot =disc_plot, device = "pdf")


##########################################################################################
###################################### MANOVA ############################################

res.manova<- manova(cbind(LogLM,LogBM1,LogLIF,LogBIF,LogBIT,LogBP,LogLN,LogBN,LogLIB,LogBB,LogBZP,
                      LogDI,LogBOC,LogLPB,LogBOF,LogBR,LogBI,LogLI,LogBBU,LogLBB,LogCZL) ~ 
                  data[,8], data= data)
result_maov<-summary.aov(res.manova, test="Wilks")


##########################################################################################
################################## Student's t-test  #####################################

LB<-print(t.test(data[1:36,][[11]], data[37:61,][[11]]))

LT<-print(t.test(data[1:36,][[12]], data[37:61,][[12]]))

LH<-print(t.test(data[1:36,][[13]], data[37:61,][[13]]))


##########################################################################################
############################### Mann-Whitney  t-test  ####################################

LE<-wilcox.test(data[1:36,][[14]], data[37:61,][[14]]) 


##########################################################################################
########################## Boxplot with external Variables ###############################

LB_boxplot<-ggplot(data, aes(x=data[,8], y=data[,11])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of the Body") +
  ylab("mm")+
  theme_classic()
LB_boxplot

ggsave("LB_boxplot.pdf",  plot =LB_boxplot, device = "pdf")


LT_boxplot<-ggplot(data, aes(x=data[,8], y=data[,12])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of the Tail") +
  ylab("mm")+
  theme_classic()
LT_boxplot

ggsave("LT_boxplot.pdf",  plot =LT_boxplot, device = "pdf")

LH_boxplot<-ggplot(data, aes(x=data[,8], y=data[,13])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of the Hindfoot") +
  ylab("mm")+
  theme_classic()
LH_boxplot

ggsave("LH_boxplot.pdf",  plot =LH_boxplot, device = "pdf")

LE_boxplot<-ggplot(data, aes(x=data[,8], y=data[,14])) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  xlab("Length of the Ears") +
  ylab("mm")+
  theme_classic()
LE_boxplot

ggsave("LE_boxplot.pdf",  plot =LE_boxplot, device = "pdf")

##########################################################################################
############################## Descriptive Statistics ####################################

# Descriptive statistics including the sample size, mean, standard error, minimum and 
#maximum value for each skull and dental variable 

descript_stat<- data.frame(data[,c(8,15:35)]) # select the biomes column and raw data

AF<- descript_stat[1:36,]
PA<- descript_stat[37:61,]

### sample size
AF_n<- nrow(AF[1:36,])
PA_n<- nrow(PA[37:61,])

### means 
media_AF<- apply(AF[,-1], 2, mean)
media_PA<- apply(PA[,-1], 2, mean)

### standard error 
erro_padrao<- function(x) {
  2*(sd(x)/sqrt(length(x)))
}

AF_error<-apply(AF[,-1], 2, erro_padrao)
PA_error<-apply(PA[,-1], 2, erro_padrao)

#### maximum value  
max_AF<- apply(AF[,-1], 2, max)
max_PA<- apply(PA[,-1], 2, max)

#### minimum value
min_AF<- apply(AF[,-1], 2, min)
min_PA<- apply(PA[,-1], 2, min)

# Descriptive statistics including the sample size, mean, standard error, minimum and 
#maximum value for each external variable 

descript_stat_LB<- data.frame(na.omit(data[,c(8,11)])) # select the biomes column and raw external data
descript_stat_LT<- data.frame(na.omit(data[,c(8,12)])) # select the biomes column and raw external data
descript_stat_LH<- data.frame(na.omit(data[,c(8,13)])) # select the biomes column and raw external data
descript_stat_LE<- data.frame(na.omit(data[,c(8,14)])) # select the biomes column and raw external data

#### LB ######
AF_LB<- descript_stat_LB[1:9,]
PA_LB<- descript_stat_LB[10:19,]

### sample size
AF_n_LB<- nrow(AF_LB[1:9,])
PA_n_LB<- nrow(PA_LB[10:19,])

### means 
media_AF_LB<- mean(AF_LB[,-1])
media_PA_LB<- mean(PA_LB[,-1])

### standard error 

AF_error_LB<-erro_padrao(AF_LB[,-1])
PA_error_LB<-erro_padrao(PA_LB[,-1])

#### maximum value  
max_AF_LB<- max(AF_LB[,-1])
max_PA_LB<- max(PA_LB[,-1])

#### minimum value
min_AF_LB<- min(AF_LB[,-1])
min_PA_LB<- min(PA_LB[,-1])

#### LT ######
AF_LT<- descript_stat_LT[1:26,]
PA_LT<- descript_stat_LT[27:42,]

### sample size
AF_n_LT<- nrow(AF_LT[1:26,])
PA_n_LT<- nrow(PA_LT[27:42,])

### means 
media_AF_LT<- mean(AF_LT[,-1])
media_PA_LT<- mean(PA_LT[,-1])

### standard error 

AF_error_LT<-erro_padrao(AF_LT[,-1])
PA_error_LT<-erro_padrao(PA_LT[,-1])

#### maximum value  
max_AF_LT<- max(AF_LT[,-1])
max_PA_LT<- max(PA_LT[,-1])

#### minimum value
min_AF_LT<- min(AF_LT[,-1])
min_PA_LT<- min(PA_LT[,-1])

#### LH ######
AF_LH<- descript_stat_LH[1:24,]
PA_LH<- descript_stat_LH[25:40,]

### sample size
AF_n_LH<- nrow(AF_LH[1:26,])
PA_n_LH<- nrow(PA_LH[27:42,])

### means 
media_AF_LH<- mean(AF_LH[,-1])
media_PA_LH<- mean(PA_LH[,-1])

### standard error 

AF_error_LH<-erro_padrao(AF_LH[,-1])
PA_error_LH<-erro_padrao(PA_LH[,-1])

#### maximum value  
max_AF_LH<- max(AF_LH[,-1])
max_PA_LH<- max(PA_LH[,-1])

#### minimum value
min_AF_LH<- min(AF_LH[,-1])
min_PA_LH<- min(PA_LH[,-1])

#### LE ######
AF_LE<- descript_stat_LE[1:27,]
PA_LE<- descript_stat_LE[28:43,]

### sample size
AF_n_LE<- nrow(AF_LE[1:26,])
PA_n_LE<- nrow(PA_LE[27:42,])

### means 
media_AF_LE<- mean(AF_LE[,-1])
media_PA_LE<- mean(PA_LE[,-1])

### standard error 

AF_error_LE<-erro_padrao(AF_LE[,-1])
PA_error_LE<-erro_padrao(PA_LE[,-1])

#### maximum value  
max_AF_LE<- max(AF_LE[,-1])
max_PA_LE<- max(PA_LE[,-1])

#### minimum value
min_AF_LE<- min(AF_LE[,-1])
min_PA_LE<- min(PA_LE[,-1])
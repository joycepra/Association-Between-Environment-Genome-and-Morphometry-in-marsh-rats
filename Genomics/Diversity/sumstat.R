#########################################################################################

#install packages
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("BSDA" %in% rownames(installed.packages()) == FALSE){install.packages("BSDA")
} else {print (paste0("'BSDA' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("ggplot2", "BSDA")

##########################################################################################
################################# INPUT FILES ############################################
# 1- cvs file with genomic diversity summary statistics 

##########################################################################################
####################################### PLOT #############################################

setwd("~/") # set working directory

data <- read.csv("Sum_stat_20mis.csv", header=TRUE, sep=",")

plot<-ggplot(data, aes(x = data$Var, color = data$Pop)) +
  geom_errorbar(aes(ymax = data$Mean + data$StdErr, ymin = data$Mean - data$StdErr),
                position=position_dodge(1)) +
  scale_color_manual(values=c("#77ab59", "#E7B800"))+
  geom_point(aes(x=data$Var, y=data$Mean), size=4, shape=21, fill="white", position=position_dodge(1)) +
  theme_classic()

ggsave("sum_stat_20mis.pdf",  plot =plot, device = "pdf")


##########################################################################################
###################################### t TEST ############################################

### Obs_Het ### 

p1<-tsum.test(mean.x=0.11057,   s.x=0.144048603, n.x=13,
              mean.y=0.09704, s.y=0.147648231, n.y=7)
p1

### Exp_Het ### 

p2<-tsum.test(mean.x=0.12871,   s.x=0.02339, n.x=13,
              mean.y=0.13324, s.y=0.02644, n.y=7)
p2

### Pi ### 

p3<-tsum.test(mean.x=0.13444,   s.x=0.1597811, n.x=13,
              mean.y=0.14548, s.y=0.177651344, n.y=7)
p3

### Fis ### 


p4<-tsum.test(mean.x=0.07621,   s.x=0.250439613, n.x=13,
               mean.y=0.12911, s.y=0.325514977, n.y=7)
p4





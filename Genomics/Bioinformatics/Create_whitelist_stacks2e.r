#### MODIFIED FROM Andrea Thomaz (2015), Sarp Kaya (2017), Giorgia Auteri (2018) Cecilia Fiorini (2019)
#######################################################################################################

#install packages
if("plyr" %in% rownames(installed.packages()) == FALSE){install.packages("plyr")
} else {print (paste0("'plyr' has already been installed in library"))}
if("pegas" %in% rownames(installed.packages()) == FALSE){install.packages("pegas")
} else {print (paste0("'pegas' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}

#load packages
pacman::p_load("plyr", "pegas")

#########################################################################################################
############################################# INPUT FILES ###############################################
# 1- .vcf file from stacks output 

#########################################################################################################
#READ VCF
#setwd("~/") # set directory with the vcf file

data <- read.table("populations.snps.vcf", header = FALSE, sep = "\t")
head(data[1:10,1:10]) #look at the data

loci_num<-data$V1 #define the loci number
pos <- data$V2 #define the variable position
  
#creates dataframe with loci ID, the variable positions and the number of individuals in each loci
new_data <- data.frame(loci_ID = loci_num,
                       pos = pos,  
                       ind = rowSums(data[,10:length(data)] != "./."))


seq_len <- as.numeric(max(new_data$pos)) # check the sequence length

head(new_data)
min(new_data$pos)
max(new_data$pos)
table(new_data$pos)
length(unique(new_data$loci_ID)) #how many loci do I have
length(new_data$loci_ID) #how many snps do I have

plot.new()
par(mar = rep(2, 4))
#saving graph with frequency of variable sites along the loci
#pdf("./SNPdistr_pos140bp.pdf")
hist(new_data$pos, breaks = c(seq(1, seq_len , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites'); #for unmapped loci
hist(new_data$pos, xlim = c(-1,seq_len), breaks = c(seq(-1, seq_len, by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');

#play with first abline number to determine cutoff, move the lines around to visualize depending on the case
abline(1800, 0, col = "red") #helps to find where starts to increase toward the end, last positions have strong increase
abline(v = 139, col = "red") #helps to figure out where to cut off before increase in bad calls
abline(v = 5, col = "red") #helps to figure out where to cut off before increase in bad calls

# BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE *end*
to_del <- 18 #how many sites to delete in the end of the sequence 
seq_len_cut <- seq_len - to_del
# create a whitelist to exclude those (to_del) positions from the end and beginning (if the case)
whitelist <- subset(new_data, pos < seq_len_cut & pos > 5)

head(whitelist)

hist(whitelist$pos, xlim = c(0,seq_len_cut), breaks = c(seq(-1, seq_len_cut -1 , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');

table(whitelist$pos)
hist(whitelist$pos, breaks = c(seq(1, seq_len , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites after cut');

#calculating theta for all loci
var.sites <- count(whitelist, "loci_ID")
length(var.sites$loci_ID) 
max(var.sites$freq) #max variable sites in one loci; 
theta_calc <- merge(unique(whitelist[,-2]), var.sites, by = "loci_ID")
theta_calc$theta <- 0
head(theta_calc)
for (i in 1:length(theta_calc$theta)){
   theta_calc[i,4] <- theta.s(theta_calc$freq[i], theta_calc$ind[i])/seq_len
 }

#calculating the 95% quantile to exclude loci extremely variable
quant <- quantile(theta_calc$theta, probs=0.95) #set the value to be variability threshold
quant 
hist(theta_calc$theta)
abline(v = quant, col="red")

#what is the maximum number of mutations in a loci
max(theta_calc$freq) #max theta before 
x <- subset(theta_calc, theta < quant)
max(x$freq) #max theta after, make sure is realistic for the length of the sequence
#think about what mutation rate the spp might have

#saving whitelist for re-run populations in stacks (write blacklist and subtract it from whitelist)
blacklist <- subset(theta_calc, theta > quant)[,1]
#write.table(blacklist, file="blacklist.txt", sep = '\n', row.names = F, col.names = F)
#above blacklist would only be for highly variable loci

#removes the blacklist from the whitelist and write off whitelist
whitelist$blacklist <- match(whitelist$loci_ID, blacklist, nomatch = 0)
whitelist_final <- subset(whitelist, blacklist == 0)
blacklist_final <- subset(whitelist, blacklist > 0)

length(unique(whitelist_final$loci_ID)) # final number of unique loci, this is the number I need to get out with "write random loci"
length(whitelist_final$loci_ID) #final number of snps
write.table(whitelist_final[,1:2], file="whitelist.txt", sep = '\t', row.names = F, col.names = F)

length(unique(blacklist_final$loci_ID)) #final number of unique loci, this is the number I need to get out with "write random loci"
length(blacklist_final$loci_ID) #final number of snps
write.table(blacklist_final[,1:2], file="blacklist.txt", sep = '\t', row.names = F, col.names = F)


#########################################################################################################
############################################# INPUT FILES ###############################################
# 1- All Fastsimcoal Results 

#########################################################################################################

wd<- setwd("~/") # set directory with the Fastsimcoal results

files <- list.files(wd)
files

dir.create("results")


  for (i in 1:length(files)) {
    setwd(paste(wd,"/",i,"/", "bras_mig_assim", sep=""))
    f<-list.files() 
    for(iter in 1:length(f)){
      s.pattern <- f[iter];
      temp.p <- paste(i, '_', sep = '');
      rep.pattern <- paste(temp.p, s.pattern, sep = '');
      #new_name<-gsub(s.pattern, rep.pattern, f)
      file.rename(from= s.pattern, to = rep.pattern)
      f <- list.files()
      file.copy(from = f, 
       to = "~/results")
            }
  }


wd <- "~/results"

setwd(wd)

best <- list.files(wd, pattern = ".bestlhoods") ##Reading the bestlhoods files
best_DF <- data.frame(ANCSIZE = numeric(),NPOP=numeric(),TDIV=numeric(), MIG21=numeric(), #creating a dataframe
                      MIG12=numeric(),MaxEstLhood=numeric(),MaxObsLhood=numeric())

for (i in best){
  temp <- read.table(i, header = T)
  best_DF <- rbind(best_DF, temp)
}

names<- data.frame(best)
best2<- cbind(names, best_DF)

saveRDS(best2, "fastsimcoal_results")
write.csv(best2, "Fatsimcoal_res.csv")



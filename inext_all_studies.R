##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Lotte inext
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()
options(stringsAsFactors = F)

#install.packages(c("iNEXT", "dplyr", "parallel", "beepr", "xlsx"))

pac <- c("iNEXT", "ggplot2", "parallel", "beepr", "xlsx")
lapply(pac, require, character.only = T)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          read data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
studies <- list.files(path = "C://Users/ma22buky/Documents/Lotte_help/iNext/", pattern = ".csv")
dat_list <- c()

for(i in 1:length(studies))
{
  ## read the data into a list (only 1 way, could do everything 
  ## into one data frame with the study as seperate column in the beginning)
  ## beware! needed to change keuper cause I found a mistake in the dataset when changed i 
  ## forgott that i will change the seperator when I save it... 
  ## disable at your laptop the if else and the read after else
  if(studies[i] != "Keuper_2011a.csv") ## disable if you didn't change the studies
    dat_list[i] <- list(read.csv(paste0("C://Users/ma22buky/Documents/Lotte_help/iNext/", studies[i]), 
                                 sep = ";", header = T, dec = "."))
  
  
  else # disable at your pc if you didn't change the studies
    dat_list[i] <- list(read.csv(paste0("C://Users/ma22buky/Documents/Lotte_help/iNext/", studies[i]), 
                                 header = T, dec = "."))
  
  
  names(dat_list)[i] <- studies[i]
  
  ## make the columns into abundance data, where everytime there is an observation it gets 1
  dat_list[[i]][,15:ncol(dat_list[[i]])] <- ifelse(dat_list[[i]][,15:ncol(dat_list[[i]])] > 0 , 1, 0)
}



inext_extended <- function(b)
{
  set.seed(b)
  result_inext <- c()
  
  for(i in 1:length(dat_list))
  {
    
    ## get the number of plots per year per treatment for each study
    min_plot  <- unlist(lapply(dat_list, 
                               FUN = function(x) 
                                 aggregate(x$year, by = list(x$year, x$treatment, x$replic.id, x$site, 
                                                             x$treatment.direction, x$treatment.level, 
                                                             x$final.year, x$season), FUN = table)[,"x"]))
    
    ## get the minimum amount of plots over all studies for sampling purposis
    min_sample <- min(min_plot) ## change when keuper is clear
    
    ## subset data into different treatments of each study and the year (need year cause each year new plot...)
    dat_list[[i]]$end_year  <- ifelse(dat_list[[i]]$final.year == 1, dat_list[[i]]$year, 0)
    tr_yr  <- unique(dat_list[[i]][,c("replic.id", "site", "treatment", "treatment.direction",
                                      "treatment.level", "end_year", "season")])
    tr_yr <- subset(tr_yr, end_year == max(dat_list[[i]]$end_year))
    
    abundance_list <- c()
    inext_raw_dat <- c()
    
    for(j in 1:nrow(tr_yr))
    {
      ## subset data in the different groups per study
      abundance_list[j] <-  list(assign(paste(tr_yr[j,"replic.id"], tr_yr[j, "treatment"],
                                              max(tr_yr[,"end_year"]), 
                                              sep = "_"), subset(dat_list[[i]], 
                                                                 treatment == tr_yr[j, "treatment"] & final.year == 1 &
                                                                   site == tr_yr[j, "site"] &
                                                                   season == tr_yr[j, "season"] & 
                                                                   treatment.direction == tr_yr[j, "treatment.direction"] &
                                                                   treatment.level == tr_yr[j, "treatment.level"])))
      
      ## rename list entries so we can follow which is which
      names(abundance_list)[j] <- paste(tr_yr[j, "replic.id"], tr_yr[j, "treatment"], max(tr_yr[, "end_year"]), 
                                        sep = "_")
      
      rm(list = paste(tr_yr[j,"replic.id"], tr_yr[j, "treatment"], tr_yr[j, "end_year"], sep = "_"))
      
      ## sample the data after the minimum amount of groups we have in a treatment (3)
      inext_raw_dat[j] <- list(as.numeric(c(min_sample, 
                                            colSums(abundance_list[[j]][sample(nrow(abundance_list[[j]]),
                                                                               min_sample),
                                                                        15:(ncol(abundance_list[[j]])-1)]))))
      
      ## rename list entries so we can follow which is which
      names(inext_raw_dat)[j] <-  paste(tr_yr[j, "replic.id"], tr_yr[j, "treatment"], tr_yr[j, "end_year"], 
                                        sep = "_")
    }
    ## use iNext function and store in a list
    result_inext[i] <- list(iNEXT(inext_raw_dat, datatype = "incidence_freq", size = 1:3, se = F))
    
    ## rename list entries so we can follow which is which
    names(result_inext)[i] <- names(dat_list)[i]
  }
  
  return(result_inext)
}


## set the number in makeCluster to your wishing (never use all cores 
## if you work on a server or you plan on doing somehting else on your pc),
## if you have 8 cores this does the following, it creates a cluster with 7 sockets: 8 - 1
## (basically it opens R 7 times... )
cl <- makeCluster(detectCores()- 1)

##need to export everything into the other sockets... since R is newly opened on those
clusterExport(cl, c("inext_extended", "dat_list"))
clusterEvalQ(cl, c(library("iNEXT"), inext_extended, dat_list))

sta <- proc.time()
inext_out <- parLapply(cl, 1:10, fun = inext_extended)
proc.time() - sta

## close the other R processes
stopCluster(cl)

beepr::beep(sound = 3)

ggiNEXT(inext_out[[2]][[1]])


delist <- c()
study_vec <- c()
study_dat <- c()

for(i in 1:length(inext_out))
  for(j in 1:length(inext_out[[i]]))
  {
    study_vec[j] <- length(inext_out[[i]][[j]]$iNextEst)
    study_dat[j] <- list(data.frame(study_treatment = names(inext_out[[i]][[j]]$iNextEst), rep = i,
                                    study = rep(names(inext_out[[i]])[j], study_vec[j])))
    
    out <- do.call("rbind", study_dat)
    for(h in 1:nrow(out))
    {
      out[h,"obs_spe"]<- inext_out[[i]][[paste(out[h, "study"])]]$iNextEst[[paste(out[h, "study_treatment"])]][3, "qD"]
    }
    
    delist[i] <- list(out)
  }

final_out <- do.call("rbind", delist)

final <- aggregate(final_out$obs_spe, by = list(final_out$study_treatment), FUN = mean)
final[,"sd"] <- aggregate(final_out$obs_spe, by = list(final_out$study_treatment), FUN = sd)[,"x"]

ggplot(final, aes(x = Group.1, y = x)) + geom_point() + theme_classic() + 
  geom_errorbar(aes(ymin = x - sd, ymax = x + sd)) + theme(axis.text.x = element_text(angle = 90))
## mean über n3 und sd 

## write the resulting table
write.xlsx(final, "C:/Users/ma22buky/Documents/Lotte_help/iNext/lotte_out_all.xlsx")
write.xlsx(final, "C:/Users/ma22buky/Documents/Lotte_help/iNext/lotte_out_mean_sd.xlsx")

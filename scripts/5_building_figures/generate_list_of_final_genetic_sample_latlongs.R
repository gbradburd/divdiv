#idea: generate list of genetic sample points used in final analyses aka made it thru popgen
#with lat/longs and tax info

#load libraries -------
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



#set path to where popgen outputs with all sample names in them live
wmpath = "data/popgen/input_and_working/ALL_r80_gendiv_data/"
#set path to where og genetic sample renaming keys live - that match sample0 to SRR IDs, so we can link lat/longs and popgen
sampkeypath = "data/all_samplenamekeys"



#get all lat/longs -------------
datasets <- list.files(path = "data/abiotic/input_and_working/lat_long_tables_per_dataset",
                       pattern = "lat_long_table", full.names = TRUE)
sites <- data.frame("link"=NA, "run_acc_sra"=NA, "lat"=NA, "long"=NA, "coordinateUncertaintyInMeters"=NA, "origin"=NA)
for ( d in datasets ) {
  out <- read.delim(d)
  sites <- rbind(sites,out)
}
sites <- sites %>% filter(is.na(link)==F) %>% dplyr::rename("lon" = "long")



#filter to just those that made it thru popgen ----------------
wmlist <- list.files(path = wmpath, 
                     pattern="initPars", 
                     full.names = TRUE)

sites.kept <- data.frame(sampid_assigned_for_bioinf = NA, run_name = NA, run_acc_sra = NA)
for (i in 1:length(wmlist)) {
  
  wmFile <- wmlist[i]
  dataset = wmFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("WMfitwishart-","",.) %>% gsub("_initPars.Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  
  load(wmFile)
  if (dim(initPars$parHom)[1] != dim(initPars$parHom)[2]) {"error - matrix dims are not the same"}
  df.i = rownames(initPars$parHom) %>% as.data.frame() %>% dplyr::rename("sampid_assigned_for_bioinf"=".") %>% 
    mutate(run_name = run_name)
  
  sampkeyFile <- paste0(sampkeypath,"/samplenamekey-",run_name,".txt")
  sampkey <- read.delim(sampkeyFile)
  df.i <- merge(df.i, 
                sampkey %>% dplyr::select(sampid_assigned_for_bioinf, run_acc_sra), 
                by = "sampid_assigned_for_bioinf", all.x = TRUE)
  sites.kept <- rbind(sites.kept, df.i)
  
}
sites.kept <- sites.kept %>% filter(is.na(sampid_assigned_for_bioinf)==F)

#are run_acc_sra's all unique in two dfs about to use/filter by?
#should return 0 if so
nrow(sites) - sites %>% dplyr::select(run_acc_sra) %>% distinct() %>% nrow()
nrow(sites.kept) - sites.kept %>% dplyr::select(run_acc_sra) %>% distinct() %>% nrow()

sites <- sites %>% filter(run_acc_sra %in% sites.kept$run_acc_sra)



#and filter by datasets in master list of datasets we're using as final check ------------
using <- read.csv("data/master_df.csv") %>% mutate(link = gsub("bioprj_","",run_name))
using <- using %>% dplyr::select(species,link)
sites <- sites %>% filter(link %in% using$link)
#tack on species names
sites <- merge(sites, using, by = "link", all.x = TRUE)



#save ----------
write.csv(sites, "data/final_genetic_latlongs.csv", row.names = FALSE)





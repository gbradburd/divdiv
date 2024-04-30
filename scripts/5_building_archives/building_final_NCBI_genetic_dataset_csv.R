#idea: build final archive for NCBI genetic datasets that we used


#load libraries
library(dplyr)
library(tidyr)
library(googledrive)
library(googlesheets4)


rm(list = ls())
gc()


setwd("/Users/rachel/divdiv")




#get final datasets we used
using <- read.csv("data/master_df.csv") %>% dplyr::select(run_name)

#get master working (trt/bookkeeping) spreadsheet from Drive
drive  <- googledrive::shared_drive_find(pattern = "^divdiv$")
drive  <- googledrive::drive_ls(path = drive , pattern = "working_datasheets", recursive = FALSE)
drive  <- googledrive::drive_ls(path = drive , pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
drive$name
drive <- googlesheets4::range_read(drive , sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#fix one taxonomic thing so records will match up
drive <- drive %>% filter(organism_biosamp != "Exaiptasia diaphana")
drive$organism_biosamp[drive$organism_biosamp == "Exaiptasia pallida"] = "Exaiptasia diaphana"


#keep info from google drive we want for final datasets --------------
drive.slim <- drive %>% dplyr::select(run_name,project_acc_bioprj,common_name,
                                      organism_biosamp,has_latlong_in_NCBI,
                                      link_to_published_citation,datathon_project_index)

df <- drive.slim %>% filter(run_name %in% using$run_name) %>% distinct()




#check if where lat/long came from matches up with datathon id -------------------
df %>% mutate(check = ifelse(is.na(datathon_project_index)==F,"in datathon","not in datath")) %>% 
  group_by(has_latlong_in_NCBI, check) %>% summarise(n=n())
  #should not be any that are y for has_latlong_in_NCBI and in datathon for check col.

#get lat/long origin from lat/long files themselves, to be sure we get right info
filelist <- list.files(path = "data/abiotic/input_and_working/lat_long_tables_per_dataset/", pattern = "*.txt", full.names = TRUE)
genetic.georefs <- data.frame("link" = NA, "georef.origin"=NA)
for (loop.iter in 1:length(filelist)) {
  
  print(paste0("starting iter ",loop.iter," of ",length(filelist)))
  file <- filelist[loop.iter]
  
  tmp <- read.delim(file)
  tmp <- tmp %>% dplyr::select(link,origin) %>% distinct() %>% rename("georef.origin" = "origin")
  
  genetic.georefs <- rbind(genetic.georefs, tmp)
  
}
genetic.georefs <- genetic.georefs %>% filter(is.na(link)==F)
#find any datasets with mixed origin
genetic.georefs %>% group_by(link) %>%
  mutate(n=1:n()) %>% filter(n > 1)
#should return 0 if no mixed origin datasets
genetic.georefs %>% distinct(link) %>% nrow() - genetic.georefs %>% nrow()
genetic.georefs <- genetic.georefs %>% mutate(run_name = paste0("bioprj_",link))

#merge on to master
df <- merge(df, genetic.georefs, by = "run_name", all.x = T)

#one more check 
df %>% group_by(has_latlong_in_NCBI, georef.origin) %>% summarise(n=n())
#just keep one of these cols
df <- df %>% dplyr::select(-has_latlong_in_NCBI) %>% rename("latlong_source" = "georef.origin")




#get genetic samples we kept/used to very end in end -------------
filelist <- list.files(path = "data/popgen/input_and_working/ALL_r80_popgen_data/", pattern = "*.Robj", full.names = TRUE)

popgen.samps <- data.frame("link" = NA, "samps.in.r80popgen"=NA)
for (loop.iter in 1:length(filelist)) {
  
  print(paste0("starting iter ",loop.iter," of ",length(filelist)))
  file <- filelist[loop.iter]
  
  dataset = file %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    strsplit(., "\\.") %>% as.data.frame() %>% .[4,]
  link = dataset %>% gsub("_popgenstats","",.) %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_") %>% gsub("bioprj_","",.)
  
  load(file)
  tmp <- popgenstats$pwp
  tmp <- colnames(tmp) %>% as.data.frame() %>% rename("samps.in.r80popgen"=".")
  tmp <- tmp %>% mutate(link = link) %>% dplyr::select(link,samps.in.r80popgen)
  
  popgen.samps <- rbind(popgen.samps,tmp)
  
}
popgen.samps <- popgen.samps %>% filter(is.na(link)==F) %>% 
  mutate(matchid = paste0(link,".",samps.in.r80popgen))

#match bioinf genetic samp names up to to NCBI IDs
filelist <- list.files(path = "data/all_samplenamekeys", pattern = "*.txt", full.names = TRUE)

sampkeys <- data.frame("link"=NA, "run_acc_sra"=NA, "sampid_assigned_for_bioinf"=NA, "pop"=NA, 
                   "identifiers_biosamp"=NA, "infraspecies_biosamp"=NA, "read_type_sra"=NA, 
                   "total_number_reads_sra"=NA, "total_reads_written_to_final_fastq"=NA)

for (loop.iter in 1:length(filelist)) {
  
  print(paste0("starting iter ",loop.iter," of ",length(filelist)))
  file <- filelist[loop.iter]
  tmp <- read.delim(file)
  sampkeys <- rbind(sampkeys,tmp)
  
}
sampkeys <- sampkeys %>% filter(is.na(link)==F)

sampkeys.slim <- sampkeys %>% dplyr::select(link, run_acc_sra, sampid_assigned_for_bioinf, identifiers_biosamp) %>% 
  mutate(matchid = paste0(link,".",sampid_assigned_for_bioinf))

#filter sampnamekeys to just samples we used/kept thru popgen
sampkeys.slim <- sampkeys.slim %>% filter(matchid %in% popgen.samps$matchid)




#concatenate SRR IDs into one row per bioprj ---------------
srrs <- sampkeys.slim %>% group_by(link) %>%
  summarize(all_SRA_sequence_ids_used = paste(unique(run_acc_sra), collapse = "|")) %>% 
  mutate(run_name = paste0("bioprj_",link)) %>% dplyr::select(-link)

#merge 
df <- merge(df, srrs, by = "run_name", all.x = T)




#concatenate SRA IDs into one row per bioprj ---------------
#check we split correctly - if so should all be BioSample:
sampkeys.slim %>% separate(., identifiers_biosamp, into = c("x","y"), sep = ";") %>% 
  separate(., x, into = c("y","biosamp"), sep = " ") %>% 
  group_by(y) %>% summarise(n=n())
#check all biosamp IDs start with SAM - if so should all be "good"
sampkeys.slim %>% separate(., identifiers_biosamp, into = c("x","y"), sep = ";") %>% 
  separate(., x, into = c("y","biosamp"), sep = " ") %>%
  mutate(check = ifelse(grepl("SAM*",biosamp)==TRUE, "good", "weird")) %>%
  group_by(check) %>% summarise(n=n())
#build data
biosamps <- sampkeys.slim %>% separate(., identifiers_biosamp, into = c("x","y"), sep = ";") %>% 
  separate(., x, into = c("y","biosamp"), sep = " ") %>%
  group_by(link) %>% 
  summarize(all_BioSample_ids_used = paste(unique(biosamp), collapse = "|")) %>% 
  mutate(run_name = paste0("bioprj_",link)) %>% dplyr::select(-link)

#merge 
df <- merge(df, biosamps, by = "run_name", all.x = T)



# !!! START HERE !!!

#count how many biosamp ids there are and how many sra ids per datasets above
#compare all above match, they should
#keep final cols we want and do any renaming, save! done!

read.csv("data/methodological/methodological_predictors-wide.csv")


#then add two additional pieces to methods final df
#single end vs paired end
#library strategy (e.g., WGS vs RAD)
#rerun methods code, save methods df
#rebuild master df (shouldn't change)
#commit/push everything









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
using <- read.csv("data/master_df.csv")
using <- using %>% dplyr::select(link) %>% rename("run_name" = "link")

#get master working (trt/bookkeeping) spreadsheet from Drive
drive  <- googledrive::shared_drive_find(pattern = "^divdiv$")
drive  <- googledrive::drive_ls(path = drive , pattern = "working_datasheets", recursive = FALSE)
drive  <- googledrive::drive_ls(path = drive , pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
drive$name
drive <- googlesheets4::range_read(drive , sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#fix one taxonomic thing so records will match up
drive <- drive %>% filter(organism_biosamp != "Exaiptasia diaphana")
drive$organism_biosamp[drive$organism_biosamp == "Exaiptasia pallida"] = "Exaiptasia diaphana"


#keep info from google drive we want for final datasets
drive.slim <- drive %>% dplyr::select(run_name,project_acc_bioprj,common_name,
                                      organism_biosamp,has_latlong_in_NCBI,
                                      link_to_published_citation,datathon_project_index)

df <- drive.slim %>% filter(run_name %in% using$run_name) %>% distinct()

#check if where lat/long came from matches up with datathon id
df %>% mutate(check = ifelse(is.na(datathon_project_index)==F,"in datathon","not in datath")) %>% 
  group_by(has_latlong_in_NCBI, check) %>% summarise(n=n())

# !!! START HERE !!!

# use lat/long scripts/data to figure out where lat/long data came from for each dataset



#get whether or not we removed adapters or not (kept track of / entered manually into an Excel)
adpts <- read.csv("data/methodological/input_and_working/master_bookkeeping_sheet-preStacks.csv")
adpts.slim <- adpts %>% dplyr::rename("name_of_adapter_removed"="adapter_names_to_remove")





#get genetic samples we kept/used to very end in end
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
popgen.samps <- popgen.samps %>% filter(is.na(link)==F)

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

sampkeys.slim <- sampkeys %>% dplyr::select(link, run_acc_sra, sampid_assigned_for_bioinf, identifiers_biosamp)

#concatenate SRR IDs into one row per bioprj
popgen.samps %>% group_by(link) %>%
  summarize(SRR_all = paste(unique(samps.in.r80popgen), collapse = "|")) %>% 
  View()





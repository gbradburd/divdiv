#idea: write out lat/long table (with SRA sample IDs) for each dataset

#load libraries ------
library(dplyr)          #data handling
library(tidyr)          #data handling
#library(googledrive)    #google sheets
#library(googlesheets4)  #google sheets
#library(readr)         #to write out txt files

rm(list=ls())
gc()

#define a few variables ----------
#define directory to write lat/long tables to 
outdir = "/Users/rachel/divdiv/data/abiotic/input_and_working/lat_long_tables_per_dataset/"
#make dir 
dir.create(outdir)

#define master list file of individual fq.gz samples to toss (for specific reasons we have decided)
samplestotoss = read.csv("/Users/rachel/divdiv/scripts/master_keys/master_tossed_samples_file.csv", header = T)

#define file with lat/long from NCBI
ncbi.input = "/Users/rachel/divdiv/data/NCBI_metadata/11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"

#define file with lat/long from datathon
datathon.input = "/Users/rachel/divdiv/data/abiotic/input_and_working/datathon_metadata_all.csv"



#get list of datasets we want to write keys for
#get master spreadsheet from Drive
using <- googledrive::shared_drive_find(pattern = "^LSA-divdiv$")
using <- googledrive::drive_ls(path = using, pattern = "working_datasheets", recursive = FALSE)
using <- googledrive::drive_ls(path = using, pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
using$name
using <- googlesheets4::range_read(using, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#keep cols we want and only datasets that are "in"
using <- using %>% dplyr::select(organism_biosamp,run_name,link,keepinrunning_YN) %>%
  filter(grepl("yes|Yes",keepinrunning_YN))

#get lat/long full dataframes ---------
#from NCBI, only those with lat/long
df.1 <- read.csv(ncbi.input, header = T) %>% 
  filter(is.na(lat_long_biosamp)==F) %>% filter(lat_long_biosamp != "unknown")
#reformat or remove some weird lat/long entries (fine bc they aren't marine/ones we are using)
df.1 <- df.1 %>% filter(!grepl("Campanula-americana",link)) %>% filter(!grepl("Drosophila-melanogaster",link)) %>% 
  filter(lat_long_biosamp != "SP526068") %>% filter(!grepl("Sceloporus-clarkii",link))

df.1.slim <- df.1 %>% dplyr::select(link,project_acc_bioprj,biosample_acc_sra,run_acc_sra,lat_long_biosamp)

#convert lat/long from N/S/E/W to +/- (so they will plot nice)
df.1.slim <- df.1.slim %>% mutate(lat_long_working = lat_long_biosamp) %>% 
  separate(.,lat_long_working,into = c("lat","lat_dir","long","long_dir"), sep = " ", remove = F)
#check that there are only directions in lat/long_dir cols
df.1.slim %>% group_by(lat_dir) %>% summarise(n=n())
df.1.slim %>% group_by(long_dir) %>% summarise(n=n())
#make values numeric
df.1.slim <- df.1.slim %>% mutate(lat = as.numeric(lat), long = as.numeric(long))
#north and east are + ; south and west are - 
df.1.slim <- df.1.slim %>% mutate(lat = ifelse(lat_dir == "N", lat, lat*-1)) %>% 
  mutate(long = ifelse(long_dir == "E", long, long*-1)) %>% dplyr::select(-lat_dir,-long_dir,-lat_long_biosamp,-lat_long_working)
#add on blank coordinateUncertaintyInMeters col (so it will match datathon records)
df.1.slim <- df.1.slim %>% mutate(coordinateUncertaintyInMeters = "")
#add on origin id
df.1.slim <- df.1.slim %>% mutate(origin = "ncbi")



#get lat/long from datathon ------------
df.2 <- read.csv(datathon.input)
df.2.slim <- df.2 %>% dplyr::select(dataset_id,project_acc_bioprj,biosample_acc_sra,run_acc_sra,decimalLatitude,decimalLongitude,coordinateUncertaintyInMeters) %>% 
  rename("link" = "dataset_id") %>% rename("lat" = "decimalLatitude") %>% rename("long" = "decimalLongitude")
#add on origin id
df.2.slim <- df.2.slim %>% mutate(origin = "datathon")



#bind together --------
df <- rbind(df.1.slim, df.2.slim)

#remove samples we don't want (from master list of individual samples to toss for X reason)
df <- df %>% filter(!run_acc_sra %in% samplestotoss$run_acc_sra)

#there are some entries that are in both NCBI list and made it into datathon 
#(bc as a whole dataset might have had some lat/long and been missing others, so full dataset made it into datathon)
#if there are duplicates, just keep datathon version (that someone checked by hand)
#generate list of run_acc_sra IDs in both datathon and ncbi lists
duplist <- df %>% group_by(run_acc_sra,origin) %>% summarise(n=n()) %>% group_by(run_acc_sra) %>% summarise(n=n()) %>% filter(n > 1)
#remove these from NCBI list, so we're only keeping records from datathon
df.1.slim.nodatathon <- df.1.slim %>% filter(!run_acc_sra %in% duplist$run_acc_sra)
#rebind
df <- rbind(df.1.slim.nodatathon, df.2.slim)
#remove samples we don't want (from master list of individual samples to toss for X reason)
df <- df %>% filter(!run_acc_sra %in% samplestotoss$run_acc_sra)
#check no run_acc_sra IDS remain that have both ncbi and datathon origin
#should return 0 if so
df %>% group_by(run_acc_sra,origin) %>% summarise(n=n()) %>% group_by(run_acc_sra) %>% summarise(n=n()) %>% filter(n > 1) %>% nrow()

#important - fix data mistakes / add some on manually that we looked up by hand
df$long[df$run_acc_sra == "SRR4291672"] = 10.70
df$lat[df$run_acc_sra == "SRR9694602"] = 9.203925336683351
df$long[df$run_acc_sra == "SRR9694602"] = -81.71610023083699  

  


#write out lat/long tables ------------
#select columns we want in lat/long files
key <- df %>% dplyr::select(link,run_acc_sra,lat,long,coordinateUncertaintyInMeters,origin) %>% distinct()

#filter df to just include datasets we want to write keys for
key <- key %>% filter(link %in% using$link)

#write out lat/long tables, per dataset aka bioprj/sp combo
key %>% group_by(link) %>%
  do(readr::write_delim(., paste(outdir, "lat_long_table-bioprj_", unique(.$link), ".txt", sep = ""),
                 delim = "\t"))



# optional
# generate list of datasets to process for marmap/ecogregion code --------
write.table(using %>% dplyr::select(run_name) %>% mutate(x = "X"), 
            "list-marmap.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


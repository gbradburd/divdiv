#idea: generate files to use to download and save fastq sequences for datasets of interest (via hpcc)

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)          #used to write out text files

rm(list = ls())
gc()


# SET UP ENVIRO -------------------------------
#define folder to write all the lists to use to keep track of which projects and samples we are downloading
outputdir = "/Users/rachel/divdiv/data/bioinformatics/lists_to_download"
#define folder that holds all the master_keys (to use to specify which datasets to run X step of the pipeline on)
masterkeysdir = "/Users/rachel/divdiv/scripts/master_keys"

#define name for master list of all datasets that we are processing/preparing here/this time around (and will subsequently get processed as a set on the hpcc)
outlist_of_lists = "list_of_datasets-XXX.txt"

#make sub directories where output will be stored
dir.create(outputdir)
dir.create(masterkeysdir)

#define path to master bioinf list of datasets we already have downloaded and read it in
df.have <- read.csv("/Users/rachel/divdiv/data/methodological/input_and_working/master_bookkeeping_sheet-preStacks.csv") %>% 
  mutate(link = gsub("bioprj_","",run_name))

#define path to NCBI full metadata
ncbi.metadata <- "/Users/rachel/divdiv/data/NCBI_metadata/11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"

#define path to master list of samples to toss
df.samples.to.toss <- "/Users/rachel/divdiv/scripts/master_keys/master_tossed_samples_file.csv"



# ********************************************************************************************************
# ********* should not need to edit below this line *********

# READ IN LIST OF DATASETS TO DOWNLOAD -------------------------------

# Load data we want to use (usually at link - bioprj_species - level)

#get list of datasets to download (at bioprj/species level)
#get path to spreadsheet in GDrive
df <- googledrive::shared_drive_find(pattern = "^divdiv$")
df <- googledrive::drive_ls(path = df, pattern = "working_datasheets", recursive = TRUE)
df <- googledrive::drive_ls(path = df, pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = TRUE)
df$name
#read spreadsheet into df
df <- googlesheets4::range_read(df, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#keep the good/vetted datasets
df <- df %>% dplyr::select(organism_biosamp,run_name,link,keepinrunning_YN) %>%
  filter(grepl("yes|Yes",keepinrunning_YN))

df %>% group_by(keepinrunning_YN) %>% summarise(n=n())

#find datasets in working master Drive divdiv list that we don't have yet
df <- df %>% filter(!link %in% df.have$link) 

#take out these datasets for now that according to NCBI metadata have both single and paired end reads
df <- df %>% filter(!link %in% c("PRJNA286089_Sphyrna-tiburo","PRJNA511386_Ostrea-lurida",
                                 "PRJNA560239_Sebastes-diaconus", "PRJNA453553_Panulirus-argus"))
#notes  - PRJNA286089_Sphyrna-tiburo download paired and single end and just keep forward reads from everyone (authors used paired to build "ref genome" then aligned single to it)
#       - PRJNA511386_Ostrea-lurida download paired and single and just keep forward reads from everyone (that's what authors did)
#       - PRJNA560239_Sebastes-diaconus - remove samples: SRR9968853 (Sto17_055), SRR9968791 (Sil17_023), and SRR9968796 (Sea17_090) which are poor quality and never analyzed by authors/not in the paper, all reads are paired even tho 2 say single (typo)
#       - PRJNA453553_Panulirus-argus - only 10 indivs at 1 location and that location was poorly resolved

#round 2 from Eric - 10-15-2020
df <- data.frame(link = c("PRJDB7819_Fibramia-amboinensis","PRJDB7819_Lutjanus-fulvus","PRJDB7819_Zenarchopterus-dunckeri",
                          "PRJNA280898_Lagenorhynchus-acutus","PRJNA280898_Lagenorhynchus-albirostris","PRJNA311981_Engraulis-encrasicolus",
                          "PRJNA340157_Callinectes-sapidus","PRJNA356786_Lateolabrax-maculatus","PRJNA359404_Sebastiscus-marmoratus",
                          "PRJNA392526_Sebastiscus-marmoratus","PRJNA437462_Fundulus-grandis","PRJNA451040_Sebastes-paucispinis",
                          "PRJNA451040_Sebastes-pinniger","PRJNA451040_Sebastes-ruberrimus","PRJNA453151_Gasterosteus-aculeatus",
                          "PRJNA477377_Pungitius-pungitius","PRJNA508986_Ctenolabrus-rupestris","PRJNA563695_Amphiprion-clarkii",
                          "PRJNA564121_Galaxias-maculatus"))

#stragglers - 10-8-2021
df <- data.frame(link = c("PRJNA210803_Pygoscelis-adeliae", 
                                     "PRJNA295681_Anguilla-rostrata", "PRJNA311981_Engraulis-encrasicolus", 
                                     "PRJNA369717_Nautilus-pompilius", "PRJNA382467_Platichthys-flesus", 
                                     "PRJNA450328_Halichoerus-grypus-atlantica", "PRJNA523574_Eudyptes-chrysolophus", 
                                     "PRJNA544154_Pygoscelis-adeliae", "PRJNA560239_Sebastes-diaconus", 
                                     "PRJNA593014_Acropora-millepora", "PRJNA632874_Anthopleura-elegantissima", 
                                     "PRJNA659918_Phocoena-sinus"))


#TEST DATA - RT Icap. PhD GBS data
df <- data.frame(link = c("PRJNA524160_Impatiens-capensis"))

#Keegan paired end project
df <- data.frame(link = c("PRJNA371817_Ostrea-lurida"))

#stragglers
df <- data.frame(link = c("PRJNA210803_Pygoscelis-adeliae", 
                          "PRJNA369717_Nautilus-pompilius", 
                          "PRJNA593014_Acropora-millepora", 
                          "PRJNA659918_Phocoena-sinus"))

# ************************************************************************************

# READ IN FULL NCBI METADATA -------------------------------

#bring in all original seq. info, at level of sequence and/or indiv sample
df.all <- read.csv(ncbi.metadata) %>% dplyr::select(-X)
#important - one dataset accidentally listed as single but it is paired (confirmed with the authors directly)
df.all$read_type_sra[df.all$link == "PRJNA437462_Fundulus-grandis"] = "PAIRED"
df.all$read_type_sra[df.all$link == "PRJNA560239_Sebastes-diaconus"] = "SINGLE"
df.all$read_type_sra[df.all$link == "PRJNA369717_Nautilus-pompilius"] = "SINGLE" #this one author said "they think" they only uploaded the forward but said reads were paired
#get subset of cols with ID # type info and species names
df.all.slim <- df.all %>% dplyr::select(link,project_acc_bioprj,organism_biosamp,run_acc_sra,experiment_acc_sra,read_type_sra,lat_long_biosamp,identifiers_biosamp,infraspecies_biosamp,total_number_reads_sra) 


#get list of samples in above bioprjs/species combos that we want
dim(df.all.slim) #full list of orig. data at indiv/seq. level
dim(df) #list to keep at bioprj/species level
df.using <- df.all.slim %>% filter(link %in% df$link)
dim(df.using) #indiv/seq. records that match bioprj/species we want to use
#check if any rows in lat_long col. have NO numbers - meaning prob not lat/long coords
#df.using %>% filter(grepl("^[^0-9]+$",lat_long_biosamp)) 
  #if 0, all lat/long entries have at least 1 number in them - and therefore are prob actually lat/long data

# OPTIONAL ************
#manually take out samples we have decided to toss for specific reasons
#read in master list of samples to toss
df.toss <- read.csv(df.samples.to.toss, header = T)
head(df.toss)
dim(df.toss)
dim(df.using)
df.using <- df.using %>% filter(!run_acc_sra %in% df.toss$run_acc_sra)
dim(df.using)
#******************************************

#check if there are multiple types of reads (single and paired) within any datasets
df.using %>% group_by(link) %>% summarise(n=n()) %>% nrow() #number of datasets
df.using %>% group_by(link,read_type_sra) %>% summarise(n=n()) %>% nrow() #number of dataset/readtype combos
  #above 2 values should be equal if each dataset contains only paired or only single end reads
  #mix/matched datasets will mess up (bioinformatics) code later

#view which datasets have mixed types of reads if needed
df.using %>% group_by(link,read_type_sra) %>% summarise(n=n()) %>% spread(.,read_type_sra,n) %>% 
  mutate(check=ifelse(is.na(PAIRED)==F&is.na(SINGLE)==F,"mixed_types_of_reads","one_read_type")) %>%
  arrange(check)

#for now - filter out two with mixed paired-end and single-end data
df.using <- df.using %>% filter(!link %in% c("PRJNA659918_Phocoena-sinus"))


# GENERATE LISTS OF RUN ACCESSION IDS (text file per dataset with all SRA numbers to download) -----------------------------------------------------

# NOTE - we want to do actual fastq download from command line using sratoolkit
#(haven't figured out an elegant/clean way to download fastq's via R)
#we want to use run_acc as the accession ID that we pass to sratoolkit
#this is the NCBI accession ID number assigned directly to the/a sequence file itself

#from NCBI:
#There are four hierarchical levels of SRA entities and their accessions:
#STUDY with accessions in the form of SRP#, ERP#, or DRP#
#SAMPLE with accessions in the form of SRS#, ERS#, or DRS#
#EXPERIMENT with accessions in the form of SRX#, ERX#, or DRX#
#RUN with accessions in the form of SRR#, ERR#, or DRR#


#do we only have unique experiment acces. numbers (sequence ID level)?
acclist <- df.using
unique(acclist) %>% as.data.frame() %>% dim()
unique(acclist$experiment_acc_sra) %>% as.data.frame() %>% dim()
unique(acclist$run_acc_sra) %>% as.data.frame() %>% dim()
  #if row numbers here are equal, yes



#write out text files of RUN ACCESSION IDs per dataset (bioprj/sp) to pass through download bash script -------------------
#grouped by "link" col. (aka bioprj/species unqiue combos) per text file

#get pared down version of dataframe
acclist.slim <- acclist %>% dplyr::select(run_acc_sra,link)

#write out a .csv for each bioprj/species combo (can't figure out how to not include grouping variable in .csv and we want just list of SRA values so doing this in two steps) 
acclist.slim %>% group_by(link) %>%
  do(readr::write_csv(., paste(outputdir,"/bioprj_", unique(.$link), ".csv", sep = "")))

#get list of all .csv files that we just wrote out
list_of_csvs = list.files(path = paste(outputdir,"/", sep = ""), pattern="bioprj_(.*).csv")
list_of_csvs

#now write out just the SRA numbers (not whole df that is in .csv) to text file
#and deleted the .csvs
for( file in list_of_csvs ) {
  temp <- read.csv(paste(outputdir,"/", file, sep = ""))
  file_name <- unique(temp$link)
  write.table(temp$run_acc_sra, 
              paste(outputdir,"/bioprj_", file_name, ".txt", sep = ""),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  unlink(paste(outputdir,"/", file, sep = ""))
} 



# WRITE LIST OF LISTS (list of datasets (bioprj_species IDs) to process) -----------------------------------------------------------------------------------------------
list_of_datasets <- list_of_csvs %>% gsub(pattern = ".csv", replacement = ".txt") #if you want to download just the datasets you are working with now and not all datasets currently in /lists_to_download
list_of_datasets
#list_of_datasets = list.files(path = "marine/sequence_level_data/lists_to_download/", pattern="bioprj*")
#list_of_datasets
#add whether project is single or paired end (so bioinf pipelines know what to do with it)
list_of_datasets <- merge(list_of_datasets %>% as.data.frame() %>% rename(temp = ".") %>% 
                            mutate(link = gsub(".txt","",temp)) %>% mutate(link = gsub("bioprj_","",link)), 
                          df.using %>% dplyr::select(link,read_type_sra) %>% unique(),
                          by = "link", all.x = T) %>% dplyr::select(-link)
list_of_datasets

write.table(list_of_datasets, 
            paste(masterkeysdir, "/", outlist_of_lists, sep = ""),
            row.names = FALSE, col.names = FALSE, quote = FALSE)



# WRITE AND STORE LIST OF ALL EXPECTED READ COUNTS ------------------------------------------------------------------------------
#make list of expected read counts (what NCBI metadata says read counts should be) for each/all sequence data file, 
#so we can check that downloading worked correctly

#easiest to just have one master list of all read counts for all NCBI reads we might ever care about
#(instead of generating a new one for each dataset)

acclist.readcnts <- df.all %>% select(run_acc_sra, total_number_reads_sra)

# NOTE - file size does not seem to match between downloads and what NCBI says it should be, but read count does
#pull out fastq file size into clean col. 
#acclist.readcnts <- acclist.readcnts %>% separate(.,expxml,into = c("x","file_size_temp") , sep = "total_size", remove = FALSE) %>% 
#  dplyr::select(-x) %>% separate(.,file_size_temp,into = c("x","file_size") , sep = "\"", remove = TRUE) %>%
#  dplyr::select(-x)

write.table(acclist.readcnts,paste(masterkeysdir, "/", "allexpectedreadcounts-6-30-2021.txt", sep = ""), row.names = FALSE, quote = FALSE)



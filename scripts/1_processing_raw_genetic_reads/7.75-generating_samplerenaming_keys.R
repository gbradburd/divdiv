#idea: generate sample renaming scripts, so we can rename fq.gz sample files from NCBI ID names to sample0, sample1, etc. for Stacks

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)    #used to write out text files


rm(list = ls())
gc()


# SET UP ENVIRO -------------------------------
#define folder to write all the sample name key files to
outputdir = "/Users/rachel/divdiv/data/all_samplenamekeys"

#define folder to write all (dummy) popmaps for Stacks to
popmapdir = "/Users/rachel/divdiv/data/bioinformatics/popmaps"

#define path to files with read counts at end of read filtering/cleaning pipeline
finalreadcountdir = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results"

#define master list file of individual fq.gz samples to toss (for specific reasons we have decided)
samplestotoss = "/Users/rachel/divdiv/scripts/master_keys/master_tossed_samples_file.csv"

#define path to NCBI full metadata
ncbi.metadata <- "/Users/rachel/divdiv/data/NCBI_metadata/11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"

#define path to folder that holds text files of SRR IDs downloaded per dataset from NCBI
lists_to_download <- "/Users/rachel/divdiv/data/bioinformatics/lists_to_download"

#read in list of datasets to make rename keys for
list_of_datasets <- read.delim("/Users/rachel/divdiv/scripts/master_keys/list-finalroundof77.txt", sep = " ", header = F) %>% 
  rename("run_name" = "V1", "read_type" = "V2") %>% 
  separate(., run_name, into = c("garbage","link"), sep = "bioprj_", remove = F) %>% 
  dplyr::select(link, everything()) %>% 
  mutate(link = gsub(".txt","",link)) %>%
  dplyr::select(-garbage)



# ********************************************************************************************************
# ********* should not normally need to edit below this line *********
# ************************************************************************************



# READ IN FULL NCBI METADATA -------------------------------

#bring in all original seq. info, at level of sequence and/or indiv sample
df.all <- read.csv(ncbi.metadata) %>% dplyr::select(-X)
#important - these datasets accidentally listed as wrong read type but it is paired (confirmed with the authors directly)
df.all$read_type_sra[df.all$link == "PRJNA437462_Fundulus-grandis"] = "PAIRED"
df.all$read_type_sra[df.all$link == "PRJNA560239_Sebastes-diaconus"] = "SINGLE"
#get subset of cols with ID # type info and species names
df.all.slim <- df.all %>% dplyr::select(link,project_acc_bioprj,organism_biosamp,run_acc_sra,experiment_acc_sra,read_type_sra,lat_long_biosamp,identifiers_biosamp,infraspecies_biosamp,total_number_reads_sra) 

#read in lists of all SRA IDs downloaded
df.sras <- data.frame(V1 = NA)
files <- list.files(path = lists_to_download, full.names = TRUE)
for ( f in files ) {
  out <- read.delim(f, header = FALSE)
  df.sras <- rbind(df.sras, out)
}
df.sras <- df.sras %>% filter(is.na(V1)==F) %>% dplyr::rename("run_acc_sra"="V1")
rm(files)



# GET METADATA FOR SAMPLES WE ARE WORKING WITH (in above bioprj/species combos) --------------

dim(df.all.slim) #full list of orig. data at indiv/seq. level
dim(list_of_datasets) #list to keep at bioprj/species level
df.using <- df.all.slim %>% filter(link %in% list_of_datasets$link)
dim(df.using) #indiv/seq. records that match bioprj/species we want to use

#keep just SRAs actually in our lists to download (in case master df has more records that we don't want)
dim(df.sras)
df.using <- df.using %>% filter(run_acc_sra %in% df.sras$run_acc_sra)
dim(df.using)

#manually take out samples we have decided to toss for specific reasons
df.toss <- read.csv(samplestotoss, header = T)
head(df.toss)
dim(df.toss)
dim(df.using)
df.using <- df.using %>% filter(!run_acc_sra %in% df.toss$run_acc_sra)
dim(df.using)



# GET FINAL READ COUNTS FOR SAMPLES GOING INTO STACKS --------------

#read in data
filenames <- list.files(finalreadcountdir, pattern="summary_of_readcounts_posttrim", full.names = T)
filenames

for( file in filenames ) {
  
  #read in file
  x <- read.csv(file)
  
  #write out dataframe of results
  name <- x$run_name %>% unique()
  assign(paste0(name, sep=""), x) #give df a unique name for each vcf file
  
}

rm(list = c("name","file","x","filenames"))

#bind all dataframes together
df.posttrim <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE))) %>% 
  dplyr::select(run_name,sample_name,total_reads_written_to_final_fastq)
#reformat data
rownames(df.posttrim) <- NULL
df.posttrim <- df.posttrim %>% separate(., run_name, into = c("garbage","link"), sep = "bioprj_", remove = F) %>% 
  dplyr::select(link, everything()) %>% dplyr::select(-garbage) %>%
  mutate(sample_name = gsub("noadaptersremoved_","",sample_name)) %>% 
  mutate(sample_name = gsub("_1.fastq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub("_1.fq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub(".fq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub(".fastq.gz","",sample_name)) %>%
  mutate(sample_name = gsub("readyforreadlengthcheck_","",sample_name)) %>% 
  rename("run_acc_sra" = "sample_name")
head(df.posttrim)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))

#merge final pre-Stacks read counts onto rest of NCBI metadata
dim(df.using)
dim(df.posttrim)
df.using <- merge(df.using, df.posttrim %>% dplyr::select(run_acc_sra,total_reads_written_to_final_fastq), by = "run_acc_sra", all.x = T)
dim(df.using)
#check that no read counts are missing
df.using %>% filter(is.na(total_number_reads_sra)==T | is.na(total_reads_written_to_final_fastq)==T) %>% nrow()




# WRITE SAMPLE RENAMING KEYS (for each dataset aka bioprj/species combo) ------------------

#make file/key to use to rename samples later (to sample0, sample1, etc. so our code works)
#  NOTE - not obvious right now where to consistently get sample names and/or pop names from across all of our datasets
#   also not sure if we need to keep track of which pop samples are from or not
#   for now, just writing "dummy" in popmap so that we don't have to worry about getting pop info for each dataset and all samples will be treated as if they are from 1 pop in Stacks
# NOTE - this is also where we will probably want to pass a list of any individual samples that we don't want to use in Stacks (e.g. bc they have low read counts)

#assign sample0-sampleN within each bioprj/species combo
samps <- df.using %>% group_by(link) %>% mutate(sampid_assigned_for_bioinf =paste0("sample", 0:(n()-1), sep = "")) %>% mutate(pop = "dummy")

#select columns we want in key
samps <- samps %>% dplyr::select(link,run_acc_sra,sampid_assigned_for_bioinf,pop,
                                         identifiers_biosamp,infraspecies_biosamp,
                                         read_type_sra,total_number_reads_sra,total_reads_written_to_final_fastq)

#write out renaming files, per bioprj/sp combo
samps %>% group_by(link) %>%
  do(write_delim(., paste(outputdir, "/sample_name_key-bioprj_", unique(.$link), ".txt", sep = ""),
                 delim = "\t"))





# WRITE (DUMMY) POPMAPS FOR DATASETS (for each dataset aka bioprj/species combo) --------------------------

popmaps <- samps %>% select(link,sampid_assigned_for_bioinf,pop)

#write out a .csv for each bioprj/species combo (can't figure out how to not include grouping variable in .csv and we want just list of SRA values so doing this in two steps)
popmaps %>% group_by(link) %>%
  do(write_csv(., paste(popmapdir, "/popmap_bioprj_", unique(.$link), ".csv", sep = "")))

#get list of all .csv files that we just wrote out
list_of_csvs = list.files(path = popmapdir, pattern="popmap_bioprj_(.*).csv")
list_of_csvs

#now write out just sampleX and pop columns (not whole df that is in .csv) to final text file for Stacks
#and deleted the .csvs
for( file in list_of_csvs ) {
  temp <- read.csv(paste(popmapdir, file, sep = ""))
  file_name <- unique(temp$link)
  temp %>% select(sampid_assigned_for_bioinf,pop) %>% 
    write.table(., 
                paste(popmapdir, "/popmap_bioprj_", file_name, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  unlink(paste(popmapdir, file, sep = ""))
}



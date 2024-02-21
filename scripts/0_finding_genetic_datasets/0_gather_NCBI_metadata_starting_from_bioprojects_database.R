#idea - search the NCBI databases through R in reproducible way (as opposed to thru web interface)
# to get genetic sequence data for DivDiv project

# NOTES !!!
# bioprj PRJNA420112 didn't come up in our intial list in Nov. 2019 but we can't figure out why
# bioprjs PRJNA392526 and PRJNA359404 both did come up and it seems all 3 of these bioprjs
# are on same species and published in same 1 paper

#- NCBI pipeline, lost about 120 bioprj records that we had in Nov. 2019 and then didn't get in Nov. 2020
#- realized in investigating that if project contains Data Type = RNA and DNA we don't get it OR exome and DNA, etc.
#- phenotype genotype category (usually QTL, domestic/model org. e.g. red vs brown cows, samples not collected randomly, but sometimes people seem to pick this if they have phenotypic data somewhere)
#- target loci category (usually what it sounds like, but sometimes people seem to call GBS targeted)
#- not sure if we are getting projects where "Data Type" = "Genome sequencing and assembly"
#- also lost a few individual sequences - but checked and these sequences were in fact deleted from NCBI
#- only 3 or 4 of ~120 bioprjs we had in Nov. 2019 and lost in Nov. 2020 are ones we might actually want
#- for merging Nov. 2019 and Nov. 2020 runs of this code, 
# decided to drop all "lost" records and keep "already had" and "gained" and add:
# "date last downloaded" (Nov 2020) and "date originally added to our records" (Nov 2020 or Nov 2019 depending on row)

#load libraries ------
library(rentrez)  #searching NCBI database
library(dplyr)    #data handling
library(tidyr)    #data handling
library(stringr)  #counting number of characters in a string
library(ggplot2)  #graphing

#define directory to save intermediate files to
workdir = "/Users/rachel/divdiv/data/NCBI_metadata/input_and_working"

#define directory to save final file to
outdir = "/Users/rachel/divdiv/data/NCBI_metadata"

#get a list of all the NCBI databases
entrez_dbs()

#get brief description of one of the databases
entrez_db_summary("bioproject")

#get list of search terms allowed in a database
entrez_db_searchable("bioproject")


#NOTE - entrez_X() functions only allow you to send a max of 3 requests per second
#requesting and using an API key bumps it up to 10 requests per second
#RToczydlowski's current key is: api_key = "3137dd8892dc5b2e9be454fe5e42eab93309"

#perform a search for BioProjects ----------
# DOWNLOAD BIOPROJECT DATA FROM NCBI ----------
#retmax is the number of ID numbers that the search will store in a vector
#max within one search is 99,999 (but can loop thru with repeat if more than this many records)
#printing return object will tell you total number of matches tho (regardless of retmax value)
#retmax doesn't matter/isn't used if using webhistory feature instead - for large requests
#searching through bioproject database bc it has best list of IDs to link to other databases
#then can look up bioproject ID in SRA database to get all sequence data associated with each project
#SRA database also returns a bunch of IDs including biosample ID
#then we can finally use biosample ID from SRA to search in Biosample database to get biosamp attributes including lat/long

rm(list=ls())
gc()

#set up search terms
query <- paste("\"scope multiisolate\"[Filter] OR \"scope multispecies\"[Filter] OR \"scope other\"[Filter] NOT \"org human\"[Filter] NOT \"org archaea\"[Filter] NOT \"org viruses\"[Filter] NOT \"org bacteria\"[Filter] NOT \"metagenome\"[Filter] NOT \"targeted locus loci\"[Filter] NOT \"clone ends\"[Filter] NOT \"metagenomic assembly\"[Filter] NOT \"transcriptome gene expression\"[Filter] NOT \"proteome\"[Filter] NOT \"epigenomics\"[Filter] NOT \"exome\"[Filter] NOT \"bioproject protein\"[Filter] NOT \"bioproject gds\"[Filter] NOT \"phenotype genotype\"[Filter] NOT \"metagenome\"[All Fields] NOT \"map\"[Filter]")
#view search terms to make sure they look right
cat(query)
#run search
search_results <- entrez_search(db="bioproject", term=query, use_history = TRUE) # 361 seems to be the max number of IDs you can pass thru entrez_summary() at once, at least for our search in BioSamples
  #note - biosample has lat/long searchable attribute, bioproject and SRA do not, but trying to link biosamps to biorpjs and SRA starting from biosamps is a mess
#view summary of search results
search_results
#for large searches, we need to save the list of IDs to the NCBI web server as a web history object (use_history = TRUE above)
#view info for the web history, which is the IDs from our search saved to the NCBI web server
search_results$web_history



# extract attribute data from BioProject search ---------
#this is syntax for getting sample info without using webhistory
#id variable refers to id element of the list returned by entrez_search(), where length of list will be = to smaller of retmax and number of hits
#search_summary <- entrez_summary(db="bioproject", id=search_results$ids) 

#get the total number of hits returned by search and stored in webhistory
n_records = search_results$count
n_records

#set number of samples to grab at a time
batch_size = 200
batch_size

paste("Requesting to run ", round(n_records/batch_size, 0), " downloads", sep = "")


#seq_start is a variable that takes each number in the sequence that seq() generates and tells entrez_search which index number of 
#the web history image to start on with the next download request
#NOTE! - first sample in web history image is referenced/indexed by 0, not 1

for( seq_start in seq(from = 0, to = n_records-1, by = batch_size)) {
  
  #fetch a batch of BioSamples
  records <- entrez_summary(db="biosample", web_history=search_results$web_history, retmax=batch_size, retstart=seq_start)
  
  #extract a dataframe of the specified elements/attributes present in each BioSample
  df <- as.data.frame(t(extract_from_esummary(records, elements = c("uid", "taxid", "project_id", "project_acc", "project_type", 
                                                                    "project_data_type", "sort_by_projecttype", "sort_by_datatype", 
                                                                    "sort_by_organism", "project_subtype", "project_target_scope", 
                                                                    "project_target_material", "project_target_capture", "project_methodtype", 
                                                                    "project_method", "registration_date", 
                                                                    "project_name", "project_title", "project_description", "keyword", 
                                                                    "relevance_agricultural", "relevance_medical", "relevance_industrial", 
                                                                    "relevance_environmental", "relevance_evolution", "relevance_model", 
                                                                    "relevance_other", "organism_name", "organism_strain", "organism_label", 
                                                                    "sequencing_status", "supergroup"))))
  #make it into prettier df
  df <- data.frame(t(matrix(unlist(df), nrow=length(df), byrow=T)))
  
  names(df) <- c("uid", "taxid", "project_id", "project_acc", "project_type", 
                 "project_data_type", "sort_by_projecttype", "sort_by_datatype", 
                 "sort_by_organism", "project_subtype", "project_target_scope", 
                 "project_target_material", "project_target_capture", "project_methodtype", 
                 "project_method", "registration_date", 
                 "project_name", "project_title", "project_description", "keyword", 
                 "relevance_agricultural", "relevance_medical", "relevance_industrial", 
                 "relevance_environmental", "relevance_evolution", "relevance_model", 
                 "relevance_other", "organism_name", "organism_strain", "organism_label", 
                 "sequencing_status", "supergroup")
  
  assign(paste("batch.number", seq_start, sep = "."),as.data.frame(df)) 
  print(paste("BioSample", seq_start, "through", seq_start+(batch_size - 1), "downloaded of", n_records,"hits total", sep = " "))
}



#bind all of the metadata for BioProjects together
full.data <- do.call(rbind, mget(grep("batch.number",names(.GlobalEnv),value=TRUE)))
write.csv(full.data, paste0(workdir, "/1-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOPROJECTS.csv"))

#Note - tested that samples that run over are just not grabbed from NCBI, e.g. if there are 315 hits and sequence to grab samples
#in the loop runs to 320, you still only get 315 unique samples (aka download is correct - grabs all the samples but no more)

#clean up data briefly - if something seems weird and you need to check which iteration/s weirdness happened at
#full.data <- full.data %>% mutate(run = rownames(.))
#full.data <- full.data %>% separate(.,run, into = c("a","b","batch.no","iteration.no"), sep = '\\.', remove = F)
#full.data$batch.no <- as.numeric(full.data$batch.no)
#full.data$iteration.no <- as.numeric(full.data$iteration.no)
str(full.data)
n_records

#check that all UIDs are unique
dim(full.data)
unique(full.data) %>% as.data.frame() %>% nrow()
unique(full.data$project_acc) %>% as.data.frame() %>% nrow()
unique(full.data$uid) %>% as.data.frame() %>% nrow()
  #above 3 lines should return the same value
  #both uid and PRJ numbers are totally unique - yay!


# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CLEAN UP NCBI BioProject DOWNLOAD --------------
rm(list=ls())
gc()

#bring in full dataset that we downloaded from NCBI (not human, bacterial, viral, metagenome)
#we have not applied any filtering to this dataset (this is what/how it came from NCBI)
df.bioprj <- read.csv(paste0(workdir, "/1-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOPROJECTS.csv"))

#filter bioprojects more based on cols of data we got back - especially project title col.
df.bioprj %>% group_by(project_data_type) %>% dplyr::summarise(n=n())
df.bioprj.f <- df.bioprj %>% filter(!project_data_type %in% c("RefSeq Transcriptome or Gene expression", "Targeted loci cultured", 
                                                                "RefSeq Genome", "RefSeq Genome sequencing and assembly", "RefSeq Other",
                                                                "RefSeq Targeted Locus (Loci)")) %>% 
  filter(!project_target_material %in% c("Phenotype", "Proteome", "Purified Chromosome", "Reagent", "Transcriptome", "Targeted loci environmental")) %>% 
  filter(!project_target_capture %in% c("Clone Ends", "Exome")) %>% 
  filter(!grepl("metagenome",project_title)) %>% filter(!grepl("Metagenome",project_title)) %>% 
  filter(!grepl("Metagenomic",project_title)) %>% filter(!grepl("metagenomic",project_title)) %>% 
  filter(!grepl("transcriptome",project_title)) %>% filter(!grepl("Transcriptome",project_title)) %>% 
  filter(!grepl("transcriptomic",project_title)) %>% filter(!grepl("RNA",project_title)) %>% 
  filter(!grepl("16S",project_title)) %>% filter(!grepl("16s",project_title)) %>% 
  filter(!grepl("18S",project_title)) %>% filter(!grepl("18s",project_title)) %>% filter(!grepl("CRISPR",project_title)) %>% 
  filter(!grepl("gene expression",project_title)) %>% filter(!grepl("Gene Expression",project_title)) %>% filter(!grepl("Gene expression",project_title)) %>% 
  filter(!grepl("metabarcoding",project_title)) %>% filter(!grepl("viral",project_title)) %>% filter(!grepl("microbiome",project_title)) %>% 
  filter(!grepl("Microbiome",project_title)) %>% filter(!grepl("bacteria",project_title)) %>% filter(!grepl("Bacteria",project_title)) %>% 
  filter(!grepl("microbial diversity",project_title)) %>% filter(!grepl("methylation",project_title)) %>% filter(!grepl("microbial community",project_title)) %>% 
  filter(!grepl("microbial",project_title)) %>% filter(!grepl("microbiota",project_title)) %>% filter(!grepl("Plasmodium falciparum",project_title)) %>% 
  filter(!grepl("metagenome",organism_name))


#remove cols that don't have anything in them (or basically nothing)
df.bioprj.f <- df.bioprj.f %>% dplyr::select(-contains("sort_by")) %>% dplyr::select(-project_subtype, -project_method, -keyword)

#uid and project_id look like the same thing - check if they are identical/redundant cols.
df.bioprj.f %>% mutate(check = ifelse(uid == project_id, "exact match", "not a match")) %>% group_by(check) %>% dplyr::summarise(n=n())
  #if all rows have value of exact match, uid and project_id are totally redundant
#only keep uid col. (remove project_id)
df.bioprj.f <- df.bioprj.f %>% dplyr::select(-project_id)

#add bioprj suffix to all col names so we can keep track of where data is coming from
colnames(df.bioprj.f) <- paste(colnames(df.bioprj.f), "bioprj", sep = "_")


write.csv(df.bioprj.f, paste0(workdir, "/2-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPROJECTS"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD ASSOCIATED SRA DATA FROM NCBI --------------
rm(list=ls())
gc()

#get bioproject data
df <- read.csv(paste0(workdir, "/2-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPROJECTS")) %>% dplyr::select(-X)

#check that every row has a (unique) bioproject number
#two below lines should return same number of rows if so
dim(df)
unique(df$project_acc_bioprj) %>% as.data.frame() %>% dim()
df$project_acc_bioprj <- as.character(df$project_acc_bioprj)
  #they do

#see what search terms we can use in SRA database
entrez_db_searchable("sra")

# pull sequence attribute metadata from SRA database ---------
#figure out what specific sequence data are available using the SRA database - search for bioproject ID in SRA database
#this is so we can get a better sense for the type of sequence data using the metadata attached to SRA sequences
#(especially to filter out RNA-seq data)

n_bioprjs = nrow(df)
n_bioprjs  #row of dataframe to end loop on
start_row = 1 #row of dataframe to start loop on
PRJcolno <- which( colnames(df)=="project_acc_bioprj" ) #get col. number of PRJ BioProject number to use below in search

#build dummy dataframe to store output of loop in
df.building.bioprjs <- data.frame("uid" = NA ,"expxml" = NA ,"runs" = NA ,"extlinks" = NA , "createdate" = NA, "PRJ" = NA, "iteration" = NA)


starttime <- date()
for(i in start_row:n_bioprjs) {
  
  PRJ <- df[i, PRJcolno]
  query <- paste(PRJ,"[GPRJ]",sep = "")
  cat(query)
  search_results <- entrez_search(db="sra", term=query, use_history = TRUE)
  
      #get the total number of hits returned by search and stored in webhistory
      n_sras = search_results$count
  
      #set number of samples to grab at a time
      batch_size = 50  
  
      print(paste("Requesting to run ", round(n_sras/batch_size, 1), " downloads", sep = ""))
      
      #build dummy dataframe to store output of SRA loop in for each bioproject
      df.building.sras <- data.frame("uid" = NA ,"expxml" = NA ,"runs" = NA ,"extlinks" = NA , "createdate" = NA, "PRJ" = NA, "iteration" = NA)
      
      #***************
      for( seq_start in seq(from = 0, to = abs(n_sras-1), by = batch_size )) {
        tryCatch({ 
          #fetch a batch of SRAs (sequence data attributes) for each BioProject
          records <- entrez_summary(db="sra", web_history=search_results$web_history, retmax=batch_size, retstart=seq_start)
        
          #extract a dataframe of the specified attributes present for each SRA record
          temp.sras <- as.data.frame(t(extract_from_esummary(records, elements = c("uid","expxml","runs","extlinks","createdate"))))
          temp.sras <- data.frame(t(matrix(unlist(temp.sras), nrow=length(temp.sras), byrow=T)))
          names(temp.sras) <- c("uid","expxml","runs","extlinks","createdate")
          temp.sras <- temp.sras %>% mutate(PRJ = query) %>% mutate(iteration = i)
          df.building.sras <- rbind(df.building.sras,temp.sras)
        }, error=function(e){cat("No SRAs found for BioProject",PRJ,":",conditionMessage(e), "\n")})
        Sys.sleep(0.5)
      }
      
  df.building.bioprjs <- rbind(df.building.bioprjs,df.building.sras)
  print(paste("SRAs for BioProject", i, "downloaded of", n_bioprjs,"BioProjects total", sep = " "))
  Sys.sleep(0.12)
  
}
endtime <- date()  

write.csv(df.building.bioprjs, paste0(workdir, "/3-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-SRAresultsforBioProjects"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CHECK IF WE DOWNLOADED ALL SRA METADATA --------------

# code often breaks because the NCBI server seems to be glitchy when asking for SRA metadata
# we can run the SRA loop to grab just the number of SRA IDs per each bioprj - and not actually do the next step of downloading the metadata
# then we can compare this list to the number of SRA IDs that we actually download metadata for to make sure we don't have any errors
# and to download the SRAs again from bioprjs where expected and observed SRA counts per bioprj don't match
# note - we redownload all SRAs in a bioprj, so if some SRAs worked and some failed in first round we will have duplicates when 
# we merge final SRA download files together (so need to remove duplicate rows when we do this)

rm(list=ls())
gc()

# get expected SRA counts per bioproject ID ---------
#get bioproject data
df <- read.csv(paste0(workdir, "/2-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPROJECTS")) %>% dplyr::select(-X)

n_bioprjs = nrow(df)
n_bioprjs  #row of dataframe to end loop on
start_row = 1 #row of dataframe to start loop on
PRJcolno <- which( colnames(df)=="project_acc_bioprj" ) #get col. number of PRJ BioProject number to use below in search

#build dummy dataframe to store output of loop in
df.checkingwegotallSRAs <- data.frame("PRJ" = NA, "iteration" = NA, "n_sra_hits" = NA)

starttime <- date()
for(i in start_row:n_bioprjs) {
  
  PRJ <- df[i, PRJcolno]
  query <- paste(PRJ,"[GPRJ]",sep = "")
  cat(query)
  search_results <- entrez_search(db="sra", term=query, use_history = TRUE)
  
  #get the total number of hits returned by search and stored in webhistory
  n_sras = search_results$count
  
  #write to dataframe and keeping appending results onto it
  df.temp <- data.frame("PRJ" = NA, "iteration" = NA, "n_sra_hits" = NA)
  df.temp <- df.temp %>% mutate(PRJ = query) %>% mutate(iteration = i) %>% mutate(n_sra_hits = n_sras)
  df.checkingwegotallSRAs <- rbind(df.checkingwegotallSRAs,df.temp)
  print(paste("There are", n_sras, "SRA hits for bioprj", PRJ, " ;  completed bioprj", i, "of", n_bioprjs,"total", sep = " "))
  Sys.sleep(0.12)
  
}
endtime <- date()  

df.checkingwegotallSRAs.clean <- df.checkingwegotallSRAs %>% filter(is.na(iteration)==F)

write.csv(df.checkingwegotallSRAs.clean, paste0(workdir,"/expected_sra_counts.csv"))


# compare the expected SRA counts to the actual SRA counts we just downloaded -----------
rm(list=ls())
gc()

#read in expect sra counts
df.exp <- read.csv(paste0(workdir,"/expected_sra_counts.csv"), header = T) %>% dplyr::select(-X)

#read in actual sra counts from our search/downloading
df.obs <- read.csv(paste0(workdir, "/3-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-SRAresultsforBioProjects"), header = T) %>% dplyr::select(-X) %>% 
  filter(is.na(iteration)==F)
df.obs <- df.obs %>% group_by(PRJ) %>% summarise(n_sras_grabbed=n())

#merge, find bio prjs where expected and observed number of SRA IDs don't match
df <- merge(df.exp, df.obs, by = "PRJ", all = T)
df$n_sras_grabbed[is.na(df$n_sras_grabbed)==T] = 0  

df <- df %>% mutate(check = ifelse(n_sra_hits == n_sras_grabbed, "yay they match", "something went wrong"))

df <- df %>% filter(check == "something went wrong")

#these bioprjs were returning more SRAs than expected (and looking at NCBI online) said they should 
#not sure what is going on, removing them (likely not studies we would use anyways = lab strains)
df.obs <- df.obs %>% filter(PRJ != "PRJNA574869[GPRJ]") %>% filter(PRJ != "PRJNA497274[GPRJ]") %>% filter(PRJ != "PRJNA251548[GPRJ]")

#write out final complete sra results (still in messy/original download form)
write.csv(df.obs, paste0(workdir, "/3-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-SRAresultsforBioProjects"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CLEAN UP SRA RESULTS SOME AND EXTRACT METADATA -------------------
rm(list=ls())
gc()

#read in full but messy sra download
df.sra <- read.csv(paste0(workdir, "/3-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-SRAresultsforBioProjects"), header = T) %>% dplyr::select(-X)

#remove any remaining dummy NA rows (that the loop put in and we didn't already get rid of)
dim(df.sra)
df.sra <- df.sra %>% filter(is.na(PRJ)==F)
dim(df.sra)

#pull out library strategy into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","library_strat") , sep = "LIBRARY_STRATEGY", remove = FALSE) %>% 
  dplyr::select(-x) 
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = "&lt;/", replacement = "", df.sra$library_strat))
df.sra  <- df.sra  %>% mutate(library_strat = gsub(pattern = "&gt;", replacement = "", df.sra$library_strat))
df.sra  <- df.sra %>% mutate(library_strat = gsub(pattern = "xmlns=\"\"", replacement = "", df.sra$library_strat))
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = ">", replacement = "", df.sra$library_strat))
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = "<", replacement = "", df.sra$library_strat))
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = "/", replacement = "", df.sra$library_strat))
#do some factor level consolidation
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = "other", replacement = "OTHER", df.sra$library_strat))
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = " ", replacement = "", df.sra$library_strat))
df.sra <- df.sra %>% mutate(library_strat = gsub(pattern = "xmlns:xsi=\"http:www.w3.org2001XMLSchema-instance\"", replacement = "", df.sra$library_strat))
df.sra %>% group_by(library_strat) %>% dplyr::summarise(n=n()) %>% as.data.frame()

#pull out organism name into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","scientific_name.working") , sep = "ScientificName=\"", remove = FALSE) %>% 
  dplyr::select(-x) 
df.sra <- df.sra %>% separate(.,scientific_name.working,into = c("scientific_name","x") , sep = "\"", remove = FALSE) %>% 
  dplyr::select(-x,-scientific_name.working) 
df.sra %>% group_by(scientific_name) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out taxID # into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","taxID.working") , sep = "Organism taxid=\"", remove = FALSE) %>% 
  dplyr::select(-x) %>% separate(.,taxID.working,into = c("taxid","x") , sep = "\"", remove = TRUE) %>% dplyr::select(-x)
df.sra %>% group_by(taxid) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out study accession # into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","study_acc.working") , sep = "Study acc=\"", remove = FALSE) %>% 
  dplyr::select(-x) %>% separate(.,study_acc.working,into = c("study_acc","x") , sep = "\"", remove = TRUE) %>% dplyr::select(-x)
df.sra %>% group_by(study_acc) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out sample accession # into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","sample_acc.working") , sep = "Sample acc=\"", remove = FALSE) %>% 
  dplyr::select(-x) %>% separate(.,sample_acc.working,into = c("sample_acc","x") , sep = "\"", remove = TRUE) %>% dplyr::select(-x)
df.sra %>% group_by(sample_acc) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out experiment accession # into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","expt_acc.working") , sep = "Experiment acc=\"", remove = FALSE) %>% 
  dplyr::select(-x) %>% separate(.,expt_acc.working,into = c("experiment_acc","x") , sep = "\"", remove = TRUE) %>% dplyr::select(-x)
df.sra %>% group_by(experiment_acc) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out BioSample into clean col.
df.sra <- df.sra %>% separate(.,expxml,into = c("x","biosample_acc.working","y") , sep = "Biosample>", remove = FALSE) %>% 
  dplyr::select(-x,-y) %>% separate(.,biosample_acc.working,into = c("biosample_acc","y") , sep = "<", remove = TRUE) %>% 
  dplyr::select(-y)
df.sra %>% group_by(biosample_acc) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#pull out library selection category - so we can maybe better filter "other" sequencing category
df.sra <- df.sra %>% separate(.,expxml,into = c("x","library_selection.working") , sep = "<LIBRARY_SELECTION>", remove = FALSE) %>% 
  dplyr::select(-x) %>% separate(.,library_selection.working,into = c("library_selection","y") , sep = "<", remove = TRUE) %>% 
  dplyr::select(-y)
df.sra %>% group_by(library_selection) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)
#do some factor level consolidation
df.sra <- df.sra %>% mutate(library_selection = gsub(pattern = "other", replacement = "OTHER", df.sra$library_selection))
df.sra <- df.sra %>% mutate(library_selection = gsub(pattern = " ", replacement = "", df.sra$library_selection))
df.sra %>% group_by(library_selection) %>% summarise(n=n()) %>% as.data.frame(.) %>% head(25)

#clean up PJR col
df.sra <- df.sra %>% mutate(PRJ = gsub(pattern = "\\[GPRJ\\]", replacement = "", df.sra$PRJ))

#add sra suffix to all col names so we can keep track of where data is coming from
colnames(df.sra) <- paste(colnames(df.sra), "sra", sep = "_")

#write out results of SRA (to link back to biosample data for filtering later)
write.csv(df.sra, paste0(workdir,"/4-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-SRAresultsforBioProjects.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# FILTER SRA DOWNLOAD BY SEQ TYPE -------------------
rm(list=ls())
gc()

#bring SRA data in
df.sra <- read.csv(paste0(workdir,"/4-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-SRAresultsforBioProjects.csv"))
df <- df.sra #so we don't have to re-read in huge df that takes time if we screw up below and need to go backwards

#now we want to filter to include only SRA entries with types of sequence data that we want
df %>% group_by(library_strat_sra) %>% dplyr::summarise(n=n()) %>% as.data.frame()

#filter
dim(df)
df.f <- df %>% filter(library_strat_sra %in% c("RAD-Seq","WCS", "WGA", "WGS", "OTHER"))
dim(df.f)
df.f %>% group_by(library_strat_sra) %>% dplyr::summarise(n=n()) %>% as.data.frame()

#write out file of SRAs (sequence data entries) that have some sort of (public) sequence data that we might want
write.csv(df.f, paste0(workdir, "/5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-SRAresultsWITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# JOIN BioProjects AND SRA DATA AND FILTER OUT SOME SPECIES WE DEF DON'T WANT -------------------
#many BioProjects don't have sequence data type listed, or any (public) sequence data associated, and/or don't have organism name filled
#so, we looked up BioProject number in SRA database (which you can do via rentrez searches) as most of the data in above line IS in SRA database
#now we want to join this additional metadata from SRA back to BioProjects 
rm(list=ls())
gc()

#bring bioproject data in
df.bioprj <- read.csv(paste0(workdir, "/2-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPROJECTS.csv")) %>% dplyr::select(-X)

#bring SRA data in
df.sra <- read.csv(paste0(workdir, "/5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-SRAresultsWITHseqdataofinterest.csv")) %>% dplyr::select(-X, -X.1)

#at what level/s are 2 above df's unique?
dim(df.bioprj)
unique(df.bioprj) %>% dim()
unique(df.bioprj$project_acc_bioprj) %>% as.data.frame() %>% dim()

dim(df.sra)
unique(df.sra) %>% dim() #unique
unique(df.sra$uid_sra) %>% as.data.frame() %>% dim() #unique
unique(df.sra$expxml_sra) %>% as.data.frame() %>% dim() #unique
unique(df.sra$biosample_acc_sra) %>% as.data.frame() %>% dim()
unique(df.sra$experiment_acc_sra) %>% as.data.frame() %>% dim() #unique
unique(df.sra$sample_acc_sra) %>% as.data.frame() %>% dim()
unique(df.sra$study_acc_sra) %>% as.data.frame() %>% dim()
unique(df.sra$scientific_name_sra) %>% as.data.frame() %>% dim()
unique(df.sra$library_strat_sra) %>% as.data.frame() %>% dim()
unique(df.sra$runs_sra) %>% as.data.frame() %>% dim()
unique(df.sra$PRJ_sra) %>% as.data.frame() %>% dim()


#merge
df.bp.sra <- merge(df.bioprj, df.sra, by.x = "project_acc_bioprj", by.y = "PRJ_sra", all.x = T)
dim(df.bp.sra)

#filter out bioprojects that don't have any (public) sequence data associated with them
df.bp.sra.f <- df.bp.sra %>% filter(is.na(uid_sra)==F)
dim(df.bp.sra.f)

#filter out metagenomes that have still snuck through
df.bp.sra.f <- df.bp.sra.f %>% filter(!grepl("metagenome",scientific_name_sra)) %>% filter(!grepl("Metagenome",scientific_name_sra))

#filter out synthetic/artifical sequences, unidentified organisms, metagenomes still not kicked out
#we looked up these taxids by hand in NCBI online to know they are weird things
df.bp.sra.f <- df.bp.sra.f %>% filter(taxid_sra != "2293429") %>% filter(taxid_sra != "1440148") %>% filter(taxid_sra != "32630") %>% filter(taxid_sra != "1427524") %>% 
  filter(taxid_sra != "81077") %>% filter(taxid_sra != "32644") %>% filter(taxid_sra != "496923")

#filter out bacteria that snuck through 
#we looked up these taxids by hand in NCBI online to know they are weird things
df.bp.sra.f <- df.bp.sra.f %>% filter(taxid_sra != "77133") %>% filter(taxid_sra != "669196") %>% filter(taxid_sra != "2282120") %>% 
  filter(taxid_sra != "668369") %>% filter(taxid_sra != "727") %>% filter(taxid_sra != "2590021") %>% 
  filter(taxid_sra != "28448") %>% filter(taxid_sra != "272943") %>% filter(taxid_sra != "2021969") %>% 
  filter(taxid_sra != "527802") %>% filter(taxid_sra != "477819") %>% filter(taxid_sra != "562") %>% 
  filter(taxid_sra != "227") %>% filter(taxid_sra != "83333") %>% filter(taxid_sra != "316407") %>% 
  filter(taxid_sra != "703612") %>% filter(taxid_sra != "224308")

#filter out some environmental samples that snuck through
#we looked up these taxids by hand in NCBI online to know they are weird things
df.bp.sra.f <- df.bp.sra.f %>% filter(taxid_sra != "1229511") %>% filter(taxid_sra != "1126252") %>% filter(taxid_sra != "1126251") 

#and add a unique ID for each row in bioproject/SRA database so we can use it for joining biosamp data from loop on in next chunk
df.bp.sra.f <- df.bp.sra.f %>% mutate(our_unique_ID = paste("ouruniquesampID", 100, 1:nrow(.), sep = ""))

#write out file of SRAs (sequence data entries) that have some sort of (public) sequence data that we might want
write.csv(df.bp.sra.f, paste0(workdir, "/6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BioProjectANDSRArecords-WITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD BIOSAMPS (for each entry returned from SRA) - SO WE CAN HAVE LAT/LONG INFO AND SPECIES -------
# Note - beginning of iterative chunk for downloading BioSamples -------------
#idea - went through this chunk of code multiple times until all that was left in bs_get dataframe are BioSample IDs that don't return any hits
#     - read in our full working list of records (df), read in BioSamples we have downloaded so far, compare, send list of missing BioSample IDs back through NCBI, repeat 
#according to convo with NCBI team, the organism name returned by biosamp database is best species name to go by for a sample (most curated/reliable)

rm(list=ls())
gc()


#get bioprj/sra records
df.want <- read.csv(paste0(workdir, "/6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BioProjectANDSRArecords-WITHseqdataofinterest.csv"), header = T, stringsAsFactors = FALSE) 

#get list of biosamps we have downloaded so far (so as to not duplicate timely step)
df.have1 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-25-2020a.csv"), header = T) %>% 
  dplyr::select(-X)

df.have2 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020a.csv"), header = T) %>% 
  dplyr::select(-X)

df.have3 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020b.csv"), header = T) %>% 
  dplyr::select(-X)

df.have4 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020c.csv"), header = T) %>% 
  dplyr::select(-X)

str(df.have1)
str(df.have2)
str(df.have3)
str(df.have4)

df.have <- rbind(df.have1,df.have2,df.have3,df.have4)

#get full list of unique BioSamples that we want
df.want <- unique(df.want$biosample_acc_sra) %>% as.data.frame() %>% rename(biosample_acc_sra = ".") %>% mutate(x = "full list we want")

#get list of unique BioSamples that we already have downloaded
df.have <- unique(df.have$biosamp_id_passed_to_loop) %>% as.data.frame() %>% rename(biosample_acc_sra = ".") %>% mutate(y = "we have this biosamp record")

#make list of which biosamps we still need to get
df.get <- merge(df.want, df.have, by = "biosample_acc_sra", all.x = T, all.y = T)
df.get <- df.get %>% filter(x == "full list we want") %>% filter(is.na(y) == T)

df <- df.get %>% dplyr::select(biosample_acc_sra)


#view list of search terms allowed in a database
entrez_db_searchable("biosample")


#build dummy dataframe to store output of loop in
df.building <- data.frame("uid" = NA,"title" = NA,"accession" = NA,"date" = NA,"publicationdate" = NA,
                          "modificationdate" = NA,"organization" = NA,"taxonomy" = NA,
                          "organism" = NA,"sourcesample" = NA,"sampledata" = NA,
                          "identifiers" = NA,"infraspecies" = NA,"package" = NA,
                          "biosamp_id_passed_to_loop" = NA, "iteration" = NA)


n_records = nrow(df)
n_records  #row of dataframe to end loop on
start_row = 1 #row of dataframe to start loop on
BioSampcolno <- which( colnames(df)=="biosample_acc_sra" ) #get col. number of biosample ID col. to use below

#look up each row (unique biosample ID) in BioSample database (using biosamp ID that we got from SRA database)
starttime <- date()
for(i in start_row:n_records) {
  tryCatch({  
    biosamp <- df[i, BioSampcolno]
    biosamp <- str_replace(biosamp, pattern = " ", replacement = "")
    query <- biosamp
    print(query)
    search_results <- entrez_search(db="biosample", term=query, use_history = FALSE, api_key = "3137dd8892dc5b2e9be454fe5e42eab93309")
    
    #fetch a batch of BioSamples
    records <- entrez_summary(db="biosample", id=search_results$ids, api_key = "3137dd8892dc5b2e9be454fe5e42eab93309")
    
    #extract a dataframe of the specified attributes present in each BioSample
    temp <- as.data.frame(t(extract_from_esummary(records, elements = c("uid","title","accession","date","publicationdate",
                                                                        "modificationdate","organization","taxonomy",
                                                                        "organism","sourcesample","sampledata",
                                                                        "identifiers","infraspecies","package"))))
    temp <- temp %>% mutate(biosamp_id_passed_to_loop = query) %>% mutate(iteration = i)
    temp <- data.frame(t(matrix(unlist(temp), nrow=length(temp), byrow=T)))
    names(temp) <- c("uid","title","accession","date","publicationdate",
                     "modificationdate","organization","taxonomy",
                     "organism","sourcesample","sampledata",
                     "identifiers","infraspecies","package",
                     "biosamp_id_passed_to_loop","iteration")
    df.building <- rbind(df.building,temp)
    Sys.sleep(time=0.12)
    print(paste("BioSample ID", i, "looked up and grabbed of", n_records,"records total", sep = " "))
  }, error=function(e){cat("ERROR with BioSample",i,":",conditionMessage(e), "\n")})
}
endtime <- date()
#tryCatch() prints an error and then allows loop to continue - e.g. when BioSample ID doesn't return any search results  

write.csv(df.building, "working_files_to_build_NCBI_records/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020c.csv") 

#Note - end of iterative chunk for downloading BioSample info -------------



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# JOIN AND SAVE ALL BioSample DOWNLOAD CSVs -----------------

rm(list=ls())
gc()

#get all the lists of biosamps we have downloaded
df1 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-25-2020a.csv"), header = T) %>% 
  dplyr::select(-X)

df2 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020a.csv"), header = T) %>% 
  dplyr::select(-X)

df3 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020b.csv"), header = T) %>% 
  dplyr::select(-X)

df4 <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-10-26-2020c.csv"), header = T) %>% 
  dplyr::select(-X)

str(df1)
str(df2)
str(df3)
str(df4)

df <- rbind(df1,df2,df3,df4)

#save a full version of all biosample downloads in one file
write.csv(df, paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-ALL.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CLEAN UP NCBI BioSample DOWNLOAD -----------------

#clean up BioSample results some and extract metadata
rm(list=ls())
gc()
df <- read.csv(paste0(workdir, "/7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOSAMPLEresultsfor_SRA-WITHseqdataofinterest-ALL.csv"), header = T) %>% dplyr::select(-X)

df <- df %>% dplyr::select(-iteration)
dim(df)
df <- df %>% filter(is.na(uid)==F)
dim(df)
df <- unique(df)
dim(df)


#extract lat/long into new clean col.
#figure out how many times "display_name=\"latitude and longitude\">" occurs in sampledata col. (should be once) - so our code to extract it will work
ncolumns <- df %>% dplyr::select(uid,title,organism,sampledata,package,identifiers) %>% 
  mutate(n_latlong = stringr::str_count(df$sampledata, "display_name=\"latitude and longitude\">"))
ncolumns %>% ggplot(aes(x=n_latlong)) + geom_histogram(binwidth = 1) + theme_bw() + theme_classic()
ncolumns %>% filter(n_latlong > 1) %>% dim() #should be 0
rm(ncolumns)

#extract lat/long into new clean col.
df <- df %>% separate(., sampledata, into = c("garbage","latlong.working"), sep = "display_name=\"latitude and longitude\">", remove = FALSE) %>% 
  select(-garbage) %>% separate(.,latlong.working, into = c("lat_long","garbage"), sep = "<", remove = TRUE) %>% 
  dplyr::select(-garbage)

df %>% filter(grepl("^[^0-9]+$",lat_long)) %>% group_by(lat_long) %>% summarise(n=n()) %>% View()

#consolidate "missing" labels in lat/long col.
df <- df %>% mutate(lat_long = gsub(pattern = "missing", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Missing", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "MISSING", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "n/a", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "N/A", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "na", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "NA", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "None", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "not applicable", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Not applicable", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Not Applicable", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "NOT APPLICABLE", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "not available", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Not available", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "not collected", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Not collected", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "not determined", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "not recorded", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "Unknown", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "UNKNOWN", replacement = "unknown", df$lat_long))
df <- df %>% mutate(lat_long = gsub(pattern = "unspecified", replacement = "unknown", df$lat_long))

#view a df quick to visually check that all lat/long entries that don't start with a number are coded as "unknown" now
df %>% filter(grepl("^[^0-9]+$",lat_long)) %>% group_by(lat_long) %>% summarise(n=n()) %>% View()
#see if all lat/long have directions in them
df %>% filter(grepl("^[^W,E,S,N]+$",lat_long)) %>% group_by(lat_long) %>% summarise(n=n()) %>% View()

#add bioprj suffix to all col names so we can keep track of where data is coming from
colnames(df) <- paste(colnames(df), "biosamp", sep = "_")
df$accession_biosamp <- as.character(df$accession_biosamp)
df$biosamp_id_passed_to_loop_biosamp <- as.character(df$biosamp_id_passed_to_loop_biosamp)

#check if biosample acc that we pulled down always matches biosamp ID that we passed to BioSample database in our search
df %>% mutate(check = ifelse(accession_biosamp == biosamp_id_passed_to_loop_biosamp, "biosamps match", "biosamps don't match")) %>% filter(check == "biosamps don't match") %>% View()
  #hmm, about 1100 instances where we searched one biosamp ID and a different biosamp ID was returned - will worry about this later if still an issue - majority are virus, humans, metagenomes, mouse

#write out clean biosample records
write.csv(df, paste0(workdir, "/8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOSAMPLEresultsforSRA-WITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD ASSOCIATED LINEAGE IDS FROM NCBI TAXONOMY DATABASE--------------
rm(list=ls())
gc()

#this section looks up lineage info by NCBI taxID number and returns full lineage info plus NCBI taxID for each level (e.g. genus)

#get list of NCBI BioProjects/SRA that have some sort of (public) sequence data that we might want
df <- read.csv(paste0(workdir, "/8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOSAMPLEresultsforSRA-WITHseqdataofinterest.csv"), header = T, stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)

#extract out list of unique taxid's - using taxID from biosamp because NCBI staff told us organism names in biosamp are the most currated/checked by NCBI staff
df <- unique(df$taxonomy_biosamp) %>% as.data.frame() %>% rename(taxonomy_biosamp = ".")


#loop through and get metadata for each species from the NCBI taxonomy database
#doing this mainly so we have full lineage (kingdom,phylum,class,etc. for each species/SRA record we have)

n_taxids = nrow(df)
n_taxids  #row of dataframe to end loop on
start_row = 1 #row of dataframe to start loop on
taxidcolno <- which( colnames(df)=="taxonomy_biosamp" ) #get col. number of PRJ BioProject number to use below in search

#build dummy dataframe to store output of loop in
df.building <- data.frame("taxid"=NA, "scientific_name"=NA, "myrank"=NA, "division"=NA, "full_lineage"=NA, 
                          "class"=NA, "family"=NA, "genus"=NA, "kingdom"=NA, "order"=NA, "phylum"=NA, "class_taxid"=NA, 
                          "family_taxid"=NA, "genus_taxid"=NA, "kingdom_taxid"=NA, "order_taxid"=NA, 
                          "phylum_taxid"=NA, "iteration"=NA, "taxid_passed_to_loop"=NA)

starttime <- date()
for(i in start_row:n_taxids) {
  tryCatch({  
    #get taxid  
    taxid <- df[i, taxidcolno]
    query <- paste(taxid,"[UID]",sep = "")
    cat(query)
    search_results <- entrez_search(db="taxonomy", term=query, use_history = FALSE, api_key = "3137dd8892dc5b2e9be454fe5e42eab93309")
    
    #retrieve data
    record <- entrez_fetch(db="taxonomy", id=search_results$ids, rettype="xml", parsed=T, api_key = "3137dd8892dc5b2e9be454fe5e42eab93309")
    record <- XML::xmlToList(record)
    #tidy up data
    temp <- data.frame(placeholder = NA) %>% mutate(taxid = record$Taxon$TaxId, scientific_name = record$Taxon$ScientificName, 
                                                    myrank = record$Taxon$Rank, division = record$Taxon$Division, full_lineage = record$Taxon$Lineage) %>% dplyr::select(-placeholder)
    temp.lineage <- record$Taxon$LineageEx
    temp.lineage <- data.frame(matrix(unlist(temp.lineage), nrow=length(temp.lineage), byrow=T)) %>% filter(X3 %in% c("kingdom","phylum","class","order",
                                                                                                                      "family","genus"))
    temp.lineage1 <- temp.lineage %>% dplyr::select(-X1) %>% spread(.,X3,X2)
    temp.lineage2 <- temp.lineage %>% dplyr::select(-X2) %>% spread(.,X3,X1)
    names(temp.lineage2) <- paste(names(temp.lineage2),"taxid", sep = "_")
    
    #bind all data for this taxid together
    temp.lineage <- cbind(temp, temp.lineage1, temp.lineage2)
    temp <- temp.lineage %>% mutate(iteration = i, taxid_passed_to_loop = gsub(pattern = "\\[UID\\]", replacement = "", query))
    
    #bind onto running dataframe (fill missing cols. with NAs)
    df.building <- plyr::rbind.fill(df.building,temp)
    
    Sys.sleep(time=0.12)
    print(paste("Taxonomic ID", i, "looked up and grabbed of", n_taxids,"records total", sep = " "))
  }, error=function(e){cat("ERROR with taxid",i,":",conditionMessage(e), "\n")})
}
#tryCatch() prints an error and then allows loop to continue - e.g. when taxonomy ID doesn't return any search results  
endtime <- date()
#checked manually that the taxids that threw errors in above loop we ok to be missing - all were metagenomes, microbiomes, or only id-ed as "fungi"
#the error -> "Error in names(temp.lineage2) <- paste(names(temp.lineage2), "taxid",  : 
#'names' attribute [1] must be the same length as the vector [0]"
#is thrown when the returned taxonomy metadata contains only lineage levels with "no rank"
write.csv(df.building, file = paste0(workdir, "/9-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-TAXONOMYIDresultsfor_BioSamples-WITHseqdataofinterest.csv"))



#clean up taxonomy download --------------
rm(list=ls())
gc()
df <- read.csv(paste0(workdir, "/9-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-TAXONOMYIDresultsfor_BioSamples-WITHseqdataofinterest.csv"), header = T) %>% dplyr::select(-X)

df <- df %>% filter(is.na(taxid_passed_to_loop)==F) %>% dplyr::select(-iteration)

#check if taxid that we looked up in NCBI taxonomy database returned the same taxid
df %>% mutate(check = ifelse(taxid == taxid_passed_to_loop, "taxid match", "taxid don't match")) %>% filter(check == "taxid don't match")
  #all match

#add taxonomy suffix to all col names so we can keep track of where data is coming from
colnames(df) <- paste(colnames(df), "taxonomy", sep = "_")

unique(df) %>% dim()

write.csv(df, paste0(workdir, "/9.5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-TAXONOMYIDresultsfor_BioSamples-WITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD ASSOCIATED LINEAGE DATA USING TAXIZE PACKAGE AND NCBI TAXONOMY DATABASE  -------------------
rm(list=ls())
gc()

#this section looks up lineage info by species names and returns species names from NCBI taxonomy database
#we need to look this up bc we want to be able to compare apples to apples between 
#our NCBI species and tips in time tree that we looked up using this same method 
#(aka any resource outside NCBI doesn't have NCBI taxIDs as labels so we need species,family,etc names)


#get list of NCBI BioProjects/SRA that have some sort of (public) sequence data that we might want
df <- read.csv(paste0(workdir, "/8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOSAMPLEresultsforSRA-WITHseqdataofinterest.csv"), header = T, stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)

#extract out list of unique organism names - using organism name from biosamp because NCBI staff told us organism names in biosamp are the most currated/checked by NCBI staff
df <- unique(df$organism_biosamp) %>% as.data.frame() %>% rename(organism_biosamp = ".")

n_organismnames = nrow(df)
n_organismnames  #row of dataframe to end loop on
start_row = 1 #row of dataframe to start loop on
organismnamescolno <- which( colnames(df)=="organism_biosamp" ) #get col. number of PRJ BioProject number to use below in search

#build dummy dataframe to store output of loop in
df.building <- data.frame("V1"=NA, "V2"=NA, "V3"=NA, "iteration"=NA, "organismname_passed_to_loop"=NA)

#look up each species in our dataset in taxize
httr::set_config(httr::config(http_version = 0))

starttime <- date()
for(i in start_row:n_organismnames) {
  
  tryCatch({  
    #get organism name  
    query <- df[i, organismnamescolno]
    query <- as.character(query)
    print(query)
    
    #look up using taxize package, databases = NCBI
    result <- taxize::classification(query, db = 'ncbi', api_key = "3137dd8892dc5b2e9be454fe5e42eab93309", callopts=list(http_version = 0L))
    #Sys.sleep(time=5)
    temp <- t(matrix(unlist(result), byrow = T, nrow = 3)) %>% as.data.frame(.) %>% 
      mutate(iteration = i) %>% mutate(organismname_passed_to_loop = query)
    
    #bind onto running dataframe
    df.building <- rbind(df.building,temp)
    
    print(paste("Organism name", i, "looked up and grabbed of", n_organismnames,"records total", sep = " "))
    
  }, error=function(e){cat("ERROR with taxid",i,":", conditionMessage(e), "\n")}, finally={Sys.sleep(time=2)})
}  
endtime <- date()
system("say I stopped running loop")

write.csv(df.building, file = paste0(workdir, "/10-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-TAXONOMYNAMESfor_BioSamples-WITHseqdataofinterest.csv"))



# CHECK WHICH RECORDS IF ANY WE FAILED TO LOOK UP/ARE MISSING -------------
rm(list=ls())
gc()
#get full list of unique organism_biosamp names that we want taxonomy for
df.want <- read.csv(paste0(workdir, "/8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOSAMPLEresultsforSRA-WITHseqdataofinterest.csv"), header = T, stringsAsFactors = FALSE) %>% 
  dplyr::select(organism_biosamp) %>% distinct() %>% as.data.frame() %>% mutate(x = "full list we want")
head(df.want)

#get list of unique organism_biosamp that we already have downloaded taxonomy for
df.have <- read.csv(paste0(workdir, "/10-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-TAXONOMYNAMESfor_BioSamples-WITHseqdataofinterest.csv")) %>% 
  filter(is.na(organismname_passed_to_loop)==F) %>% dplyr::select(organismname_passed_to_loop) %>% distinct() %>% as.data.frame() %>% 
  mutate(y = "we have this record") %>% rename("organism_biosamp"="organismname_passed_to_loop")
head(df.have)

#make list of which biosamps we still need to get
df.get <- merge(df.want, df.have, by = "organism_biosamp", all.x = T, all.y = T)
df.get <- df.get %>% filter(x == "full list we want") %>% filter(is.na(y) == T)

df <- df.get %>% dplyr::select(organism_biosamp)

#go back to start of section to run df that we just created through the taxize lookup loop



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CLEAN UP taxize NCBI Taxonomy DOWNLOAD -----------------

#clean up taxize NCBI Taxonomy results some and convert long df to wide df
rm(list=ls())
gc()
df <- read.csv(paste0(workdir, "/10-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-TAXONOMYNAMESfor_BioSamples-WITHseqdataofinterest.csv")) %>% dplyr::select(-X)
head(df)
df <- df %>% filter(is.na(organismname_passed_to_loop)==F)

#filter out the no-rank groups (because we don't have a good way to handle them in a wide format)
#also filter out clade, because many organism_biosamps have many different/multiple clades listed
df <- df %>% filter(V2 != "no rank") %>% filter(V2 != "clade")

#check if any other organisms that we looked up have multiple entries for a taxonomic level
#n here is the number of "names" that organism has for that taxonomic level (e.g. an organism where 2 genera were returned)
df %>% group_by(organismname_passed_to_loop,V2) %>% summarise(n=n()) %>% filter(n>1)

#remove 3 bateria that have multiple classes, family, etc.
df <- df %>% filter(!organismname_passed_to_loop %in% 
                      c("Anatoma sp. MNHN-IM-2013-42003",
                        "Bathysciadiidae sp. MNHN-IM-2013-40843",
                        "Pseudococculinidae sp. MNHN-IM-2013-40847"))

df %>% group_by(organismname_passed_to_loop,V2) %>% summarise(n=n()) %>% filter(n>1)
  #this needs to be 0 for below long to wide to work

#get rid of a few ranks with barely any (nonhelpful) info
df %>% filter(V2 == "genotype")
df <- df %>% filter(V2 != "genotype")
df %>% filter(V2 == "isolate")
df <- df %>% filter(V2 != "isolate")
df %>% filter(V2 == "series")
df <- df %>% filter(V2 != "series")
df %>% filter(V2 == "species subgroup")
df <- df %>% filter(V2 != "species subgroup")
df %>% filter(V2 == "biotype")
df <- df %>% filter(V2 != "biotype")
df %>% filter(V2 == "subsection")
df <- df %>% filter(V2 != "subsection")
df %>% filter(V2 == "subvariety")
df <- df %>% filter(V2 != "subvariety")
df %>% filter(V2 == "forma")
df <- df %>% filter(V2 != "forma")
df %>% filter(V2 == "forma specialis")
df <- df %>% filter(V2 != "forma specialis")
df %>% filter(V2 == "section")
df <- df %>% filter(V2 != "section")
df %>% filter(V2 == "strain")
df <- df %>% filter(V2 != "strain")

#convert long to wide
#(it's ok that some species have more classification levels returned than others - code just fills missing ranks with NA)
df <- df %>% select(-V3,-iteration) %>% unique(.) %>% spread(.,V2,V1)

#add bioprj suffix to all col names so we can keep track of where data is coming from
colnames(df) <- paste(colnames(df), "ncbitaxonomy-via-taxize", sep = "_")

write.csv(df, paste0(workdir, "/10.5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-TAXONOMYNAMESfor_BioSamples-WITHseqdataofinterest.csv"))




# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD ASSOCIATED FULL RECORDS FROM NCBI BIOPROJECTS (so we have publication info ID) --------------
rm(list=ls())
gc()

#read in list of all bioprjs we are working with
df <- read.csv(paste0(workdir, "/6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BioProjectANDSRArecords-WITHseqdataofinterest.csv")) %>% 
  dplyr::select(-X)

#get list of unique bioprojIDs to fetch full entries for from NCBI bioprj database 
bioprjIDlist <- df$project_acc_bioprj %>% unique()

#look up each bioprj ID and fetch full record from NCBI
df.building <- data.frame("project_acc_bioprj" = NA, "full_fetched_record_bioprj" = NA)

for( bioprjID in bioprjIDlist ) {
  tryCatch({
    
      #look up ID for each bioprj
      search_results <- entrez_search(db="bioproject", term=bioprjID, use_history = T)
      #get full record from NCBI
      temp <- entrez_fetch(db = "bioproject", web_history=search_results$web_history, rettype = "xml", parsed = F)
  
      temp <- temp %>% as.data.frame() %>% rename("full_fetched_record_bioprj" = ".")
      temp <- temp %>% mutate(project_acc_bioprj = bioprjID) %>% dplyr::select(project_acc_bioprj,full_fetched_record_bioprj)
  
      df.building <- rbind(df.building,temp)
  
      which(bioprjIDlist[bioprjID] == 2585)
      print(paste("I am done getting record for bioprj ID",bioprjID, "iteration", which(bioprjIDlist[] == bioprjID), "of", length(bioprjIDlist), sep = " "))
      }, error=function(e){cat("ERROR with taxid", which(bioprjIDlist[] == bioprjID),":", conditionMessage(e), "\n")}, finally={Sys.sleep(time=2)})
}

write.csv(df.building, paste0(workdir, "/10.6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOPRJFULLRECORDS-WITHseqdataofinterest.csv"), row.names = FALSE)



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# CLEAN UP NCBI BioProject FULL RECORD DOWNLOAD -----------------
#some bioprj records don't have a publication ID (to look up in pubmed NCBI database) but they do have citation info in the full record,
#so we are extracting that citation info when available here
rm(list=ls())
gc()

#get output from previous step
df <- read.csv(paste0(workdir, "/10.6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-BIOPRJFULLRECORDS-WITHseqdataofinterest.csv")) %>% dplyr::select(-X)

#remove NA line from dummy df loop
df <- df %>% filter(is.na(project_acc_bioprj)==F)

#get NCBI publication ID (altho sometimes seems to be a DOI?)
df <- df %>% 
  separate(., full_fetched_record_bioprj, into = c("garbage","temp"), sep = "<Publication\ id=", remove = F) %>% 
  dplyr::select(-garbage) %>% separate(.,temp,into = c("garbage","publication_ID_bioprj"), sep = "\"", remove = T) %>% 
  dplyr::select(-garbage)

#get publication title
df <- df %>% separate(., full_fetched_record_bioprj, into = c("garbage","temp"), sep = "<StructuredCitation>", remove = F) %>%
  dplyr::select(-garbage) %>% separate(., temp, into = c("garbage","publication_title_bioprj"), sep = "Title>", remove = F) %>% 
  mutate(publication_title_bioprj = gsub("</", "", publication_title_bioprj)) %>% dplyr::select(-garbage)
#get publication journal info
df <- df %>% separate(., temp, into = c("garbage","temp1"), sep = "Journal>", remove = F) %>% dplyr::select(-garbage)
#get journal name
df <- df %>% separate(., temp1, into = c("garbage","publication_journalname_bioprj"), sep = "JournalTitle>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(publication_journalname_bioprj = gsub("</", "", publication_journalname_bioprj))
#get journal year
df <- df %>% separate(., temp1, into = c("garbage","publication_year_bioprj"), sep = "Year>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(publication_year_bioprj = gsub("</", "", publication_year_bioprj))
#get journal volume
df <- df %>% separate(., temp1, into = c("garbage","publication_volume_bioprj"), sep = "Volume>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(publication_volume_bioprj = gsub("</", "", publication_volume_bioprj))
#get journal issue
df <- df %>% separate(., temp1, into = c("garbage","publication_issue_bioprj"), sep = "Issue>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(publication_issue_bioprj = gsub("</", "", publication_issue_bioprj))
#get journal pages
df <- df %>% separate(., temp1, into = c("garbage","pagesstart"), sep = "PagesFrom>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(pagesstart = gsub("</", "", pagesstart)) %>% 
  separate(., temp1, into = c("garbage","pagesend"), sep = "PagesTo>", remove = F) %>% dplyr::select(-garbage) %>% 
  mutate(pagesend = gsub("</", "", pagesend)) %>% 
  mutate(publication_pages_bioprj = paste(pagesstart,pagesend, sep = "-")) %>% dplyr::select(-temp1)
#get first author name
df <- df %>% separate(., temp, into = c("garbage","tempname"), sep = "Author>", remove = F) %>% dplyr::select(-garbage) %>% 
  separate(., tempname, into = c("garbage","tempname1"), sep = "Name>", remove = T) %>% 
  separate(., tempname1, into = c("garbage","namefirst"), sep = "First>", remove = F) %>% 
  mutate(namefirst = gsub("</", "", namefirst)) %>% 
  separate(., tempname1, into = c("garbage","namelast"), sep = "Last>", remove = F) %>% 
  mutate(namelast = gsub("</", "", namelast)) %>% 
  mutate(publication_firstauthor_bioprj = paste(namelast,namefirst, sep = ", "))

#clean up two cols. so they have true NA not NA,NA as character
df <- df %>% mutate(publication_firstauthor_bioprj = gsub("NA, NA", NA, publication_firstauthor_bioprj)) %>% 
  mutate(publication_pages_bioprj = gsub("NA-NA", NA, publication_pages_bioprj))

#get final cols we want
df <- df %>% dplyr::select(project_acc_bioprj, full_fetched_record_bioprj,publication_ID_bioprj,
                                                                 publication_firstauthor_bioprj, publication_year_bioprj,
                                                                 publication_title_bioprj, publication_journalname_bioprj, 
                                                                 publication_volume_bioprj, publication_issue_bioprj, publication_pages_bioprj)

#save
write.csv(df, paste0(workdir, "/10.7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPRJFULLRECORDS-WITHseqdataofinterest.csv"), row.names = FALSE)



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# DOWNLOAD ASSOCIATED PUBMED INFO FROM NCBI (so we have publication info)--------------
#some bioprj full records downloaded in previous section are a pubmed NCBI ID and other are the DOI for the paper themselves 
#when pubmed NCBI ID is present in bioprj full record, citation info often isn't present
#so, for those bioprjs where found a pubmed NCBI ID from the full record, 
#we will pass those pubmed NCBI IDs to the pubmed NCBI database to get associated citation info 
#note - pubmed entrez_search will accept many types of IDs (including DOIs)
rm(list=ls())
gc()

#get bioprj full records info
df <- read.csv(paste0(workdir, "/10.7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPRJFULLRECORDS-WITHseqdataofinterest.csv"))

#get just rows that have a publication ID (found using previous step/section)
df <- df %>% filter(is.na(publication_ID_bioprj)==F)

#get list of unique publication IDs
pubmedIDlist <- df$publication_ID_bioprj %>% unique()

#strip off any url addresses from DOIs - they will work in a search without e.g. "https://doi.org/" but not with
pubmedIDlist <- pubmedIDlist %>% as.data.frame(.) %>% rename("publication_ID_bioprj" = ".") %>% 
  mutate(pubmedID_for_lookup = gsub("https://doi.org/","",publication_ID_bioprj)) %>% 
  mutate(pubmedID_for_lookup = gsub("http://dx.doi.org/","",pubmedID_for_lookup)) %>% 
  mutate(pubmedID_for_lookup = gsub("doi: ","",pubmedID_for_lookup)) %>% 
  mutate(pubmedID_for_lookup = gsub("dx.doi.org/","",pubmedID_for_lookup)) %>% 
  mutate(pubmedID_for_lookup = gsub("doi:","",pubmedID_for_lookup))

#look up publication in NCBI pubmed database to get more info and better formatted data

df.building <- data.frame("uid" = NA, "pubdate" = NA, "epubdate" = NA, "source" = NA, "lastauthor" = NA, "title" = NA, 
                          "sorttitle" = NA, "volume" = NA, "issue" = NA, "pages" = NA, "lang" = NA, "nlmuniqueid" = NA, 
                          "issn" = NA, "essn" = NA, "pubtype" = NA, "recordstatus" = NA, "pubstatus" = NA, "attributes" = NA, 
                          "pmcrefcount" = NA, "fulljournalname" = NA, "elocationid" = NA, "doctype" = NA, "booktitle" = NA, 
                          "medium" = NA, "edition" = NA, "publisherlocation" = NA, "publishername" = NA, "srcdate" = NA, 
                          "reportnumber" = NA, "availablefromurl" = NA, "locationlabel" = NA, "docdate" = NA, 
                          "bookname" = NA, "chapter" = NA, "sortpubdate" = NA, "sortfirstauthor" = NA, "vernaculartitle" = NA, 
                          "DOI" = NA, "publication_ID_bioprj" = NA)

for( pubmedID in pubmedIDlist$pubmedID_for_lookup) {
  tryCatch({ 
    pubmedID.og <- pubmedIDlist %>% filter(pubmedID_for_lookup == pubmedID)
    publication_ID_bioprj <- pubmedID.og$publication_ID_bioprj
    
    search_results <- entrez_search(db="pubmed", term=pubmedID, use_history = T)
    temp <- entrez_summary(db = "pubmed", web_history=search_results$web_history)
    
    #extract a dataframe of the specified elements/attributes present in each pubmed entry
    temp <- as.data.frame(t(extract_from_esummary(temp, elements = c("uid", "pubdate", "epubdate", "source", "authors", "lastauthor", 
                                                                     "title", "sorttitle", "volume", "issue", "pages", "lang", "nlmuniqueid", 
                                                                     "issn", "essn", "pubtype", "recordstatus", "pubstatus", "articleids", 
                                                                     "history", "references", "attributes", "pmcrefcount", "fulljournalname", 
                                                                     "elocationid", "doctype", "srccontriblist", "booktitle", "medium", 
                                                                     "edition", "publisherlocation", "publishername", "srcdate", "reportnumber", 
                                                                     "availablefromurl", "locationlabel", "doccontriblist", "docdate", 
                                                                     "bookname", "chapter", "sortpubdate", "sortfirstauthor", "vernaculartitle"))))
    
    #extract DOI out specifically
    temp.identifiers <- temp[[19]][[1]]
    temp.identifiers <- temp.identifiers %>% filter(idtype == "doi")
    doi <- temp.identifiers$value
    temp <- temp %>% mutate(DOI = doi)
    
    #remove cols that are lists so we can write out a 2D csv file
    temp <- temp %>% dplyr::select(-doccontriblist, -references, -srccontriblist, -authors, -history, -articleids)
    
    #convert list of lists to one dataframe
    temp <- data.frame(t(matrix(unlist(temp), nrow=length(temp), byrow=T)))
    
    #reassign col names
    names(temp) <- c("uid", "pubdate", "epubdate", "source", "lastauthor", "title", 
                     "sorttitle", "volume", "issue", "pages", "lang", "nlmuniqueid", 
                     "issn", "essn", "pubtype", "recordstatus", "pubstatus", "attributes", 
                     "pmcrefcount", "fulljournalname", "elocationid", "doctype", "booktitle", 
                     "medium", "edition", "publisherlocation", "publishername", "srcdate", 
                     "reportnumber", "availablefromurl", "locationlabel", "docdate", 
                     "bookname", "chapter", "sortpubdate", "sortfirstauthor", "vernaculartitle", 
                     "DOI")
    
    #save just the original full publication_ID_bioprj and the one we looked up (with http: etc. stripped off)
    temp <- temp %>% mutate(publication_ID_bioprj = publication_ID_bioprj)
    
    df.building <- rbind(df.building, temp)
    
    print(paste("I am done getting record for pubmed ID",pubmedID, "iteration", which(pubmedIDlist[] == pubmedID)[1], "of", nrow(pubmedIDlist), sep = " "))
    
  }, error=function(e){cat("Error for pubmedID",pubmedID,":",conditionMessage(e), "\n")}, finally={Sys.sleep(time=2)})
}
#note - tryCatch above allows loop to keep running if a pubmed ID throws an error

df.building.OG.BACKUP <- df.building

#reorder cols
df.building <- df.building %>% dplyr::select(publication_ID_bioprj, DOI, everything()) %>% filter(is.na(publication_ID_bioprj)==F)

#remove rows that have a DOI that doesn't start with "10."
df.building <- df.building %>% filter(grepl('^10.',DOI))

#inspect any duplicates
#df.building %>% group_by(publication_ID_bioprj) %>% filter(n()>1) %>% View()

#remove corrigendum and erratum, simplyfing our life and keeping one paper per bioprj for now
#df.building <- df.building %>% filter(sorttitle != "erratum") %>% filter(sorttitle != "corrigendum")

#keep unique rows based on publication_ID_bioprj and DOI
#df.building <- df.building %>% group_by(publication_ID_bioprj,DOI) %>% arrange(DOI,desc(pubdate)) %>% mutate(temp = 1:n()) %>% 
#  filter(temp == 1) %>% dplyr::select(-temp) %>% ungroup(.)

#check that we just have one row per publication_ID_bioprj ID now (should return 0)
#df.building %>% group_by(publication_ID_bioprj) %>% filter(n()>1) %>% nrow()

#add pubmed suffix to all col names so we can keep track of where data is coming from
colnames(df.building) <- paste(colnames(df.building), "pubmed", sep = "_")
df.building <- df.building %>% rename("publication_ID_bioprj" = "publication_ID_bioprj_pubmed")
names(df.building)

write.csv(df.building, paste0(workdir, "/10.8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-PUBMEDDOIS-WITHseqdataofinterest.csv"), row.names = FALSE)



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# JOIN BIOSAMP, TAXONOMY, and PUBLICATION DATA to BIOPROJECTS/SRA with WGS or RAD-Seq  -------------------
rm(list=ls())
gc()

#get list of bioprojects/sra with sequence data that we might want (WGS, Rad-Seq, or other)
df.bp.sra <- read.csv(paste0(workdir, "/6-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BioProjectANDSRArecords-WITHseqdataofinterest.csv")) %>% dplyr::select(-X, -our_unique_ID)

#get biosample data
df.biosamp <- read.csv(paste0(workdir, "/8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOSAMPLEresultsforSRA-WITHseqdataofinterest.csv")) %>% dplyr::select(-X)
        #note - some biosamp IDs returned multiple rows back, to see them run - 
        df.biosamp %>% group_by(biosamp_id_passed_to_loop_biosamp) %>% mutate(check = 1:n()) %>% filter(check>1) %>% dplyr::select(biosamp_id_passed_to_loop_biosamp,check,everything()) %>% View()

#merge
df <- merge(df.bp.sra, df.biosamp, by.x = "biosample_acc_sra", by.y = "biosamp_id_passed_to_loop_biosamp", all.x = T)

#get taxonomy id info
df.taxids <- read.csv(paste0(workdir, "/9.5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-TAXONOMYIDresultsfor_BioSamples-WITHseqdataofinterest.csv", header = T)) %>% dplyr::select(-X)

#merge
df <- merge(df, df.taxids, by.x = "taxonomy_biosamp", by.y = "taxid_passed_to_loop_taxonomy", all.x = T)

#get taxonomy names info
df.taxnames <- read.csv(paste0(workdir, "/10.5-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-TAXONOMYNAMESfor_BioSamples-WITHseqdataofinterest.csv", header = T)) %>% dplyr::select(-X)

#merge
df <- merge(df, df.taxnames, by.x = "organism_biosamp", by.y = "organismname_passed_to_loop_ncbitaxonomy.via.taxize", all.x = T)

#get publication DOIs from bioproject full records
#decided for now to keep just DOI - simpler - and can come back to get above additional publication info/details later if we want
df.pubinfo <- read.csv(paste0(workdir, "/10.7-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-BIOPRJFULLRECORDS-WITHseqdataofinterest.csv")) %>% dplyr::select(project_acc_bioprj,publication_ID_bioprj) %>% distinct()
#get additional publication DOIs from pubmed via publication ID in bioprj full records
df.pubmeddois <- read.csv(paste0(workdir, "/10.8-NCBI_bioprojects__not-human_bacteria_viral_metagenome__CLEAN-PUBMEDDOIS-WITHseqdataofinterest.csv")) %>% dplyr::select(publication_ID_bioprj,DOI_pubmed) %>% distinct()
#merge these together
df.pub <- merge(df.pubinfo, df.pubmeddois, by = "publication_ID_bioprj", all = T) %>% 
  dplyr::select(project_acc_bioprj,publication_ID_bioprj,DOI_pubmed)

#merge
df <- merge(df, df.pub, by = "project_acc_bioprj", all.x = T) 

#save
write.csv(df, paste0(workdir, "/10.9-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# FINAL PASS OF CREATING COLUMNS WE WANT FOR DOWNSTREAM STUFF -------------------
rm(list=ls())
gc()

#read in output from previous step
df <- read.csv(paste0(workdir, "/10.9-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv")) %>% select(-X)

#pull out single vs paired end info into clean col.
df <- df %>% separate(.,expxml_sra,into = c("garbage","temp1"), sep = "LIBRARY_LAYOUT", remove = FALSE) %>% dplyr::select(-garbage) %>% 
  separate(.,temp1,into = c("garbage","temp2"), sep = "<", remove = TRUE) %>% dplyr::select(-garbage) %>% 
  mutate(temp3 = gsub(" ","/",temp2)) %>% separate(.,temp3,into = c("read_type_sra","garbage"), sep = "/", remove = TRUE) %>% 
  dplyr::select(-temp2,-garbage)

#NOTE - decided not to pull out sequence file size (in bytes) from runs_sra column here bc it doesn't seem to match
#       file sizes we actually get when we download test files, total_number_reads_sra do match between test downloads and what this metadata says

#pull out collection date
df <- df %>% separate(., sampledata_biosamp, into = c("garbage","temp"), sep = "collection date", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(., temp, into = c("collection_date_biosamp","garbage"), sep = "<", remove = TRUE) %>% 
  mutate(collection_date_biosamp = gsub(">","",collection_date_biosamp)) %>% 
  mutate(collection_date_biosamp = gsub("\"","",collection_date_biosamp)) %>% 
  dplyr::select(-garbage)

#pull out sampleID type info about samples in a few different cols.
#note some helpful info might be in infraspecies_biosamp, but we don't need to extract this further
df <- df %>% separate(.,sampledata_biosamp, into = c("garbage","temp"), sep = "sample name", remove = FALSE) %>% 
  dplyr::select(-garbage) %>% separate(., temp, into = c("samplename_sampledata_biosamp","garbage"), sep = "<", remove = TRUE) %>% 
  mutate(samplename_sampledata_biosamp = gsub(">","",samplename_sampledata_biosamp)) %>% 
  mutate(samplename_sampledata_biosamp = gsub("\"","",samplename_sampledata_biosamp)) %>% 
  dplyr::select(-garbage)

df <- df %>% separate(.,identifiers_biosamp, into = c("garbage","temp1"), sep = "Sample name: ", remove = FALSE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp1, into = c("samplename_identifiers_biosamp","garbage"), sep = ";", remove = TRUE) %>% 
  dplyr::select(-garbage)

df <- df %>% separate(.,expxml_sra, into = c("garbage","temp1"), sep = "LIBRARY_NAME", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(.,temp1, into = c("library_name_expxml_sra","garbage"), sep = "<", remove = TRUE) %>% 
  dplyr::select(-garbage) %>% mutate(library_name_expxml_sra = gsub(">","",library_name_expxml_sra)) %>% 
  mutate(library_name_expxml_sra = gsub("/",NA,library_name_expxml_sra))

#pull out contact info
#pull out center name (e.g. University of Wisconsin-Madison) from SRA database col.
df <- df %>% separate(.,expxml_sra,into = c("garbage","temp1"), sep = "center_name=\"", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("center_name_sra","garbage"), sep = "\"", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out contact name from SRA database col.
df <- df %>% separate(.,expxml_sra,into = c("garbage","temp1"), sep = "contact_name=\"", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("contact_name_sra","garbage"), sep = "\"", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out lab name from SRA database col.
df <- df %>% separate(.,expxml_sra,into = c("garbage","temp1"), sep = "lab_name=\"", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("lab_name_sra","garbage"), sep = "\"", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out organization info from biosample database col.
df <- df %>% separate(.,sampledata_biosamp,into = c("garbage","temp1"), sep = "Owner>", remove = FALSE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("garbage","messy_organization_name_biosamp"), sep = "<Name url=\"", remove = TRUE) %>% 
  dplyr::select(-garbage) %>% separate(.,messy_organization_name_biosamp,into = c("organization_name_messy_biosamp","garbage"), sep = "Name>", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out contact email from biosample database col.
df <- df %>% separate(.,sampledata_biosamp,into = c("garbage","temp1"), sep = "Contact email=\"", remove = FALSE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("contact_email_biosamp","garbage"), sep = "\"", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out contact first name and last name from biosample database col.
df <- df %>% separate(.,sampledata_biosamp,into = c("garbage","temp1"), sep = "<Name>", remove = FALSE) %>% 
  separate(.,temp1,into = c("temp2","garbage"), sep = "Name>", remove = TRUE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp2,into = c("garbage","temp3"), sep = "<First>", remove = FALSE) %>%
  separate(.,temp3,into = c("contact_first_name_biosamp","garbage"), sep = "<", remove = TRUE) %>%
  dplyr::select(-garbage) %>% separate(.,temp2,into = c("garbage","temp4"), sep = "<Last>", remove = TRUE) %>%
  separate(.,temp4,into = c("contact_last_name_biosamp","garbage"), sep = "<", remove = TRUE) %>% 
  dplyr::select(-garbage)

#pull out collected by info from biosample database col.
df <- df %>% separate(.,sampledata_biosamp,into = c("garbage","temp1"), sep = "display_name=\"collected by\">", remove = FALSE) %>%
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("collectors_info_biosamp","garbage"), sep = "<", remove = TRUE) %>%
  dplyr::select(-garbage)

#NOTE!!! - put correct date range in here
#add info on about what date the metadata in this file was downloaded from NCBI
df <- df %>% mutate(date_this_metadata_downloaded = "10-19-2020 - 11-6-2020")


#work towards pulling out accession seq run info --------
#note - most returns from SRA have 1 run_acc_sra aka sequence ID per row
#but, some have more than 1 sequence per row aka multiple run_acc_sras, read counts, etc. per row
#these should all be for the same individual (e.g. an individual sequenced multiple times)
#since we care about unique organisms/individuals (not unique sequences) for purposes of sample size filtering
#we will extract the number of run acc IDs in runs_sra and the number of runs there should be according to NCBI and save these in columns for our records
#then we will extract the first run_acc_sra for every row (so some rows will have more buried run_acc_sras that we can go back and get later if we want)

#count ourselves how many times Run acc appears in runs_sra megacolumn
df <- df %>% mutate(run_acc_sra_count.us = str_count(runs_sra, pattern = "Run acc"))
#and also extract sequence count from NCBI records
df <- df %>% separate(.,expxml_sra,into = c("garbage","temp1") , sep = "total_runs=", remove = FALSE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp1,into = c("temp2","garbage") , sep = " ", remove = TRUE) %>% 
  dplyr::select(-garbage) %>% separate(.,temp2,into = c("run_acc_sra_count.ncbi","garbage") , sep = "/", remove = TRUE) %>%
  dplyr::select(-garbage) %>% mutate(run_acc_sra_count.ncbi = as.numeric(gsub(pattern = "\"", replacement = "", run_acc_sra_count.ncbi))) 

df %>% dplyr::select(expxml_sra,runs_sra,run_acc_sra_count.ncbi,run_acc_sra_count.us,project_acc_bioprj) %>% 
  mutate(check = ifelse(run_acc_sra_count.ncbi == run_acc_sra_count.us, "match","not a match")) %>% 
  group_by(check) %>% summarise(n=n())

#pull out FIRST run accession #, read counts, total number of base pairs, per row into clean col. 
df <- df %>% separate(.,runs_sra,into = c("empty","temp") , sep = "<", remove = FALSE) %>% dplyr::select(-empty) %>%
  separate(.,temp,into = c("x1","x2","x3","x4") , sep = " ", remove = TRUE) %>% dplyr::select(-x1) %>% 
  rename("run_acc_sra" = "x2") %>% 
  rename("total_number_reads_sra" = "x3") %>% 
  rename("total_number_base_pairs_sra" = "x4")
  
df <- df %>% mutate(run_acc_sra = gsub(pattern = "acc", replacement = "", df$run_acc_sra)) 
df <- df %>% mutate(run_acc_sra = gsub(pattern = "\\=", replacement = "", df$run_acc_sra))
df <- df %>% mutate(run_acc_sra = gsub(pattern = "\"", replacement = "", df$run_acc_sra))
df <- df %>% mutate(total_number_reads_sra = gsub(pattern = "total", replacement = "", df$total_number_reads_sra))
df <- df %>% mutate(total_number_reads_sra = gsub(pattern = "_spots", replacement = "", df$total_number_reads_sra))
df <- df %>% mutate(total_number_reads_sra = gsub(pattern = "\\=", replacement = "", df$total_number_reads_sra))
df <- df %>% mutate(total_number_reads_sra = gsub(pattern = "\"", replacement = "", df$total_number_reads_sra))
df <- df %>% mutate(total_number_base_pairs_sra = gsub(pattern = "total", replacement = "", df$total_number_base_pairs_sra))
df <- df %>% mutate(total_number_base_pairs_sra = gsub(pattern = "_bases", replacement = "", df$total_number_base_pairs_sra))
df <- df %>% mutate(total_number_base_pairs_sra = gsub(pattern = "\\=", replacement = "", df$total_number_base_pairs_sra))
df <- df %>% mutate(total_number_base_pairs_sra = gsub(pattern = "\"", replacement = "", df$total_number_base_pairs_sra))
df <- df %>% mutate(total_number_reads_sra = as.numeric(total_number_reads_sra))
df <- df %>% mutate(total_number_base_pairs_sra = as.numeric(total_number_base_pairs_sra))


#CREATE "LINK" COLUMN
#this is used extensively in downstream analyses to join diff. dataframes and to keep track of/label datasets
#this is the unit/level that we care about for our projects
#filling species name spaces with "-" so that we can use to e.g label directories and code won't be crabby about spaces in file/dir names
df <- df %>% mutate(link = paste(project_acc_bioprj, organism_biosamp, sep = "_")) %>% 
        mutate(link = gsub(pattern = " ", replacement = "-", link))

#visually check things look right in these new cols. we just created
df %>% select(link,run_acc_sra,total_number_reads_sra,total_number_base_pairs_sra,read_type_sra) %>% head()

#move some cols we care about/will use a lot up front
df <- df %>% select(link,project_acc_bioprj,organism_biosamp,run_acc_sra,everything())

#save temp file before we do the comparison of previous records to these records
write.csv(df, paste0(workdir, "/10.95-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# EXPLORE AND MARK WHICH RECORDS ARE NEW SINCE LAST TIME -------------------
rm(list=ls())
gc()
#get previous records
df.previous.og <- read.csv(paste0(workdir, "/firstgrab-11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest-pre11-9-2020.csv")) 
df.previous <- df.previous.og %>% dplyr::select(run_acc_sra,project_acc_bioprj,runs_sra) %>% distinct() %>% mutate(x = "had this seq id already") %>% filter(is.na(run_acc_sra)==F) %>% mutate(runs_sra = as.character(runs_sra))

#get records from this run
df.have.og <- read.csv(paste0(workdir, "/10.95-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"))
df.have <- df.have.og %>% dplyr::select(run_acc_sra,project_acc_bioprj,runs_sra) %>% distinct() %>% mutate(y = "got this seq id this time") %>% filter(is.na(run_acc_sra)==F) %>% mutate(runs_sra = as.character(runs_sra))

#compare run_acc_sras (finest level of ID in NCBI - one ID per sequence file) from last run and this run
#extract out all run_acc_sra IDs in long format (this extract all run IDs in runs_sra mega column, even if there are more than 1)
df.have <- df.have %>% mutate(temp = strsplit(runs_sra, "Run acc=")) %>% 
  unnest(temp) %>% separate(.,temp,into = c(NA,"run_acc_sra"), sep = "\"", remove = TRUE) %>% 
  filter(is.na(run_acc_sra)==F) %>% dplyr::select(-runs_sra) %>% distinct()
head(df.have)
df.previous <- df.previous %>% mutate(temp = strsplit(runs_sra, "Run acc=")) %>% 
  unnest(temp) %>% separate(.,temp,into = c(NA,"run_acc_sra"), sep = "\"", remove = TRUE) %>% 
  filter(is.na(run_acc_sra)==F) %>% dplyr::select(-runs_sra) %>% distinct()
head(df.previous)
#merge
status <- merge(df.previous, df.have, by = "run_acc_sra", all = T) %>% 
  mutate(project_acc_bioprj.x = as.character(project_acc_bioprj.x)) %>% 
  mutate(project_acc_bioprj.y = as.character(project_acc_bioprj.y))
head(status)
#assign some labels
status$y[is.na(status$y)==T] = "missing"
status$x[is.na(status$x)==T] = "missing"
status$project_acc_bioprj.x[is.na(status$project_acc_bioprj.x)==T] = "missing"
status$project_acc_bioprj.y[is.na(status$project_acc_bioprj.y)==T] = "missing"
status$check[status$x=="had this seq id already" & status$y=="missing"] = "lost"
status$check[status$x=="missing" & status$y=="got this seq id this time"] = "gained"
status$check[status$x=="had this seq id already" & status$y=="got this seq id this time"] = "no_change"
status %>% group_by(x,y,check) %>% summarise(n=n())
status.w <- status %>% group_by(project_acc_bioprj.x,project_acc_bioprj.y,check) %>% summarise(n=n()) %>% ungroup()
status.wx <- status.w %>% dplyr::select(project_acc_bioprj.x,check,n) %>% distinct() %>% rename("project_acc_bioprj" = "project_acc_bioprj.x") %>% filter(project_acc_bioprj != "missing")
status.wy <- status.w %>% dplyr::select(project_acc_bioprj.y,check,n) %>% distinct() %>% rename("project_acc_bioprj" = "project_acc_bioprj.y") %>% filter(project_acc_bioprj != "missing")
status.w.all <- rbind(status.wx,status.wy) %>% distinct() 
status.w.all <- status.w.all %>% spread(.,"check","n")
nrow(status.w.all)
rbind(df.have %>% distinct(project_acc_bioprj), df.previous %>% distinct(project_acc_bioprj)) %>% distinct() %>% dim()
  #above and status.w.all should have same number of rows (if there is a row for every unique bioprj found previously and/or now)  
status.w.all$lost[is.na(status.w.all$lost)==T] = 0 
status.w.all$gained[is.na(status.w.all$gained)==T] = 0 
status.w.all$no_change[is.na(status.w.all$no_change)==T] = 0 
status.w.all$status[status.w.all$lost == 0 & status.w.all$gained == 0] = "ok, had all seq records in this bioprj before"
status.w.all$status[status.w.all$lost == 0 & status.w.all$no_change == 0] = "ok, all seq records are in a new bioprj"
status.w.all$status[status.w.all$gained == 0 & status.w.all$no_change == 0] = "problem, lost all seq records in the bioprj"
status.w.all$status[status.w.all$lost == 0 & status.w.all$gained >= 1 & status.w.all$no_change >= 1] = "ok, added new seq records and already had some in this bioprj"
status.w.all$status[is.na(status.w.all$status)==T] = "problem, some combo of lost,gained,nochange seq records in the bioprj"


status.w.all %>% group_by(status) %>% summarise(n=n())

#get df of problem bioprjs - to investigate further manually
lost <- status.w.all %>% filter(grepl("problem",status))

#check if all run_acc_sras are unique in our check/investigating list (should return 0 if so)
status %>% nrow() - status %>% distinct(run_acc_sra) %>% nrow()

#check if all run_acc_sras in our df from download this time around are unique (should return 0 if so)
df.have.og %>% nrow() - df.have.og %>% distinct(run_acc_sra) %>% nrow()
#pull out run_acc_sras that are duplicated in our current df
dups <- df.have.og %>% group_by(run_acc_sra) %>% summarise(n=n()) %>% filter(n>1) %>% dplyr::select(run_acc_sra)
#this should return 0 if no duplciated run_acc_sra IDs are in multiple categories
merge(dups, status, by = "run_acc_sra", all.x = T) %>% dplyr::select(run_acc_sra,check) %>% 
  group_by(run_acc_sra,check) %>% summarise(n=n()) %>% group_by(run_acc_sra) %>% summarise(n=n()) %>% filter(n>1) %>% nrow()
#and double check that all run_acc_sra IDs are only in one category (e.g. no_change or gained but not both) (should also return 0)
merge(df.have.og, status, by = "run_acc_sra", all.x = T) %>% dplyr::select(run_acc_sra,check) %>% 
  group_by(run_acc_sra,check) %>% summarise(n=n()) %>% group_by(run_acc_sra) %>% summarise(n=n()) %>% filter(n>1) %>% nrow()
#finally, check that all run_acc_sras in our current df have a status (in our investigative work)
merge(df.have.og %>% filter(is.na(run_acc_sra)==F) %>% dplyr::select(run_acc_sra) %>% mutate(x = "in our current df"),
      status %>% dplyr::select(run_acc_sra,check),
      by = "run_acc_sra", all.x = T) %>% 
  group_by(x,check) %>% summarise(n=n())
  #records with NA for check are in our current df but not the status investigative list

#make key of when each run_acc_sra was added to our records
datekey <- status %>% dplyr::select(run_acc_sra,check) %>% filter(check != "lost") %>% distinct()
datekey$date_this_sequence_first_added_to_our_records[datekey$check == "no_change"] = "10-28-2019 - 11-7-2019"
datekey$date_this_sequence_first_added_to_our_records[datekey$check == "gained"] = "10-19-2020 - 11-6-2020"
datekey <- datekey %>% dplyr::select(-check)
datekey %>% group_by(date_this_sequence_first_added_to_our_records) %>% summarise(n=n())
head(datekey)

#add date that each run_acc_sra ID was first added to our records
dim(df.have.og)
merge(df.have.og, datekey, by = "run_acc_sra", all.x = T) %>% dim()
df <- merge(df.have.og, datekey, by = "run_acc_sra", all.x = T)

#save this "final" file!
write.csv(df, paste0(outdir, "/11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"))




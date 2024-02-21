#idea: run list of unique species in NCBI records we have gathered through WoRMS (marine taxonomy database) to find the marine species

#load libraries ------
library(rentrez)  #searching NCBI database
#library(worms)   #WoRMs database access, to scrape species names for marine species 
#NOTE: we need worms package but it masks a lot of dplyr functions with plyr functions so only calling it by :: and not loading it
library(dplyr)    #data handling
library(tidyr)    #data handling
library(stringr)  #counting number of characters in a string
library(ggplot2)  #graphing

# ************************************************************************************************************************************
# ************************************************************************************************************************************
#USE WoRMS TO SCRAPE NCBI BIOSAMPLE ORGANISM COL. FOR MARINE SPECIES -----------------

rm(list=ls())
gc()

#define final NCBI metadata ("all" species)
ncbi.allsp = "/Users/rachel/divdiv/data/NCBI_metadata/11-NCBI_bioprojects__not-human_bacteria_viral_metagenome__BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"

#define dir for intermediate files
workdir = "/Users/rachel/divdiv/data/NCBI_metadata/input_and_working"

#define dir for final output file
outdir = "/Users/rachel/divdiv/data/NCBI_metadata"


#get list of NCBI BioProjects and SRA that have some sort of (public) sequence data that we might want (all species right now minus humans, bacterial, viral, metagenome ...?)
df <- read.csv(ncbi.allsp, header = T, stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)



#working towards getting unique species list as character vector
dim(df)
df.u <- df %>% group_by(organism_name_bioprj,scientific_name_sra,scientific_name_taxonomy,organism_biosamp) %>% dplyr::summarise(n=n()) %>% ungroup()

df.u %>% mutate(check = ifelse(scientific_name_taxonomy == scientific_name_sra, "exact match", "not a match")) %>% group_by(check) %>% dplyr::summarise(n=n())
df.u %>% mutate(check = ifelse(scientific_name_taxonomy == organism_biosamp, "exact match", "not a match")) %>% group_by(check) %>% dplyr::summarise(n=n())
df.u %>% mutate(check = ifelse(scientific_name_sra == organism_biosamp, "exact match", "not a match")) %>% group_by(check) %>% dplyr::summarise(n=n())
  #no cols with species info appear to be exact matches/duplicates across all records

#NCBI puts most effort into biosample names
# from NCBI biosample staff - SRA and GenBank now get their organism information directly from BioSamples, 
# so those should always be the same, but SRA formerly was fairly lax in taxonomy vetting and there are legacy cases 
# where the organism may not be exactly the same, as you saw, though it won't be outright wrong. 
# SRA doesn't consider the taxonomy critical, so not much effort has been made to go back and clean those up.
# how many (SRA) records don't have a biosample species name?
df.u %>% filter(is.na(organism_biosamp) == T) %>% dim()
# view them (only 13 on 11-19-19; 15 on 11-13-2020)
df %>% filter(is.na(organism_biosamp) == T) %>% View()
  # one is inbred lines of fruit flies, some are Sorghum bicolor and Asclepias Sonoran Desert Clade
  # so far we are not missing anything important by focusing only on biosample taxonomy 

# now look at those biosample species names that don't match SRA species info
df.u %>% mutate(check = ifelse(scientific_name_sra == organism_biosamp, "exact match", "not a match")) %>% filter(check == "not a match") %>% 
  unique(.) %>% View()
#biosamples with names that don't match often seem to be from where biosample ID in SRA database that was passed to BioSample database is not a standard SAM#
#let's look at places where biosample ID is not an SAM number and biosamp species name doesn't match sra species name
df %>% filter(!grepl("SAM",biosample_acc_sra)) %>% 
  mutate(check = ifelse(scientific_name_sra == organism_biosamp, "exact match", "not a match")) %>% filter(check == "not a match") %>% 
  select(biosample_acc_sra,contains("organism"),contains("scientific"),everything()) %>% View()
  #as of 11-20-19 they are all sorghum bicolor samples with weird biosamp IDs from SRA - chuck these
#now let's check the rest of the biosamps that don't match but have normal SRA numbers 
df %>% filter(grepl("SAM",biosample_acc_sra)) %>% 
  mutate(check = ifelse(scientific_name_sra == organism_biosamp, "exact match", "not a match")) %>% filter(check == "not a match") %>% 
  select(biosample_acc_sra,contains("organism"),contains("scientific"),everything()) %>% View()
  #as of 11-20-19 most are plants and given names are only at genus level or higher
#let's go ahead with using only the biosamp species info (best currated according to NCBI contact)

#doing a quick test to check if old taxonomic names will be found with WORMS exact matching 
#new name (as of 2016)
matchf.new <- worms::wormsbynames(taxon_names="Chthamalus alani", match=FALSE, marine_only = "true")
matcht.new <- worms::wormsbynames(taxon_names="Chthamalus alani", match=TRUE, marine_only = "true")
#old name
matchf.old <- worms::wormsbynames(taxon_names="Chthamalus southwardorum", match=FALSE, marine_only = "true")
matcht.old <- worms::wormsbynames(taxon_names="Chthamalus southwardorum", match=TRUE, marine_only = "true")
  # yes - both old and new names are found by exact and fuzzy matching

#get unique species name vector 
species.list <- df %>% dplyr::select(organism_biosamp) %>% distinct() 
    #13,600 unique "species" 11-20-19; 18,556 on 11-13-2020
species.list <- as.character(species.list$organism_biosamp) 
str(species.list) #wormsbynames() requires character vector

#don't load worms package here bc it uses plyr which can mess other pkgs we are using up - just call using ::
#devtools::install_github("janhoo/worms") (the older version that R installs has bugs)
#library(worrms) 

#pass all of the "species" names in the organism col. from NCBI through WoRMs database to find the marine species
#note! - this takes a long time to run (minutes to hours)
#note! - match=T and match=F in wormsbynames() bring back diff results
#note! - per above test, taxonomic name changes/synonyms are caught/returned by exact matches, so we don't need fuzzy matches for taxonomy changes
#       - also, heard from NCBI Biosample staffer John Anderson that submission with names that don't match NCBI taxonomic database are flagged and 
#       -typos are corrected and/or new taxIDs assigned for new species
# for these reasons, decided not to use fuzzy matching anymore, as fuzzy matching seems to return mainly garbage - e.g. grass genera spelled similar to marine genus

#pass all of the species names from biosample species info (NCBI) through WoRMS database to find the marine species
#use exact matching only (not fuzzy name matching) - match=FALSE
starttime <- date()
marine.species.matchF <- worms::wormsbynames(taxon_names=species.list, match=FALSE, chunksize = 100, marine_only = "true", sleep_btw_chunks_in_sec = 0.1)
endtime <- date()

#write out results (insurance for time costly step)
write.csv(marine.species.matchF, paste0(workdir, "/NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-marineSPECIESLISTreturnedbyWoRM-matchF-biosamplespeicesnameonly.csv"))


#bring in and clean/filter WoRMS results ---------------
df.marine.sps.matchF <- read.csv(paste0(workdir, "/NCBI_bioprojects__not-human_bacteria_viral_metagenome__RAW-marineSPECIESLISTreturnedbyWoRM-matchF-biosamplespeicesnameonly.csv"), header = T) %>% 
  dplyr::select(-X) %>% mutate(iteration = 1:nrow(.))

#filter out NAs (biosample was not found in WoRMS)
dim(df.marine.sps.matchF)
exact <- df.marine.sps.matchF %>% filter(is.na(AphiaID)==F) %>% mutate(lookuptype = "exact name (no fuzzy matching)")
dim(exact) #1228 on 11-20-19; 1744 on 11-19-2020

#quick side track - get some summary stats ----------
#how many species are not at all marine?
exact %>% filter(isMarine == "0" | is.na(isMarine)==T) %>% dim()
#how many are only marine?
exact %>% filter(isMarine == "1") %>% filter(isBrackish == "0" | is.na(isMarine)==T) %>% dim()
#how many are marine/brackish?
exact %>% filter(isMarine == "1") %>% filter(isBrackish == "1") %>% dim()
#and how many are brackish but not marine?
exact %>% filter(isBrackish == "1") %>% filter(isMarine == "0" | is.na(isMarine)==T) %>% dim()

#here is everything that had an exact name match
marine.sp.matched <- exact

#make key to match scientific names returned by WoRMS to NCBI organism name that we submitted
species_list <- as.data.frame(species.list) %>% mutate(iteration = 1:nrow(.))
names(species_list) <- c("organism_name_passed_to_worms","iteration")
dim(species_list)

#add on organism names (so we can add no matches on next and cols. will all match)
marine.species <- merge(species_list,marine.sp.matched,by = "iteration",all.x = T, all.y = T)

#to check everything worked, grab species directly that weren't found at all in WoRMS
nomatch <- species_list %>% filter(!iteration %in% marine.sp.matched$iteration) %>% mutate(AphiaID=NA, url=NA, scientificname=NA)
#check that total number of species names we sent through WoRMS (nrow(species_list)) match sum of species with no matches and species with exact matches
nrow(species_list)
nrow(nomatch) + nrow(exact)
nrow(marine.species)
#these three above lines should give them same number of rows if everything worked

#write out a list of clean marine species results from WoRMS
write.csv(marine.species, paste0(workdir, "/NCBI_bioprojects__not-human_bacteria_viral_metgenome__CLEAN-marineSPECIESLISTreturnedbyWoRMEXACTmatchesonlyusingbiosampspeciesnames.csv"))



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# FILTER NCBI RECORDS BY "MARINE" SPECIES -------------------
rm(list=ls())
gc()

#bring full NCBI data in (filtered somewhat for sequence data but not yet lat/long)
df <- read.csv(ncbi.allsp, header = T, stringsAsFactors = F) %>% dplyr::select(-X)

#bring in marine species list from WoRMS - exact matches only
df.marine.sps <- read.csv(paste0(workdir, "/NCBI_bioprojects__not-human_bacteria_viral_metgenome__CLEAN-marineSPECIESLISTreturnedbyWoRMEXACTmatchesonlyusingbiosampspeciesnames.csv"), 
                          header = T, stringsAsFactors = F) %>% dplyr::select(-X)

#get list of organisms that were found in WoRMS
dim(df.marine.sps)
df.marine.sps.matched <- df.marine.sps %>% filter(is.na(AphiaID)==F)
dim(df.marine.sps.matched)

#filter our NCBI records to only include marine
df.marine <- df %>% filter(organism_biosamp %in% df.marine.sps.matched$organism_name_passed_to_worms)
df.marine <- unique(df.marine)

#join WoRMS info onto dataframe
#add worms suffix to all col names so we can keep track of where data is coming from
colnames(df.marine.sps.matched) <- paste(colnames(df.marine.sps.matched), "worms", sep = "_")

df.marine <- merge(df.marine,df.marine.sps.matched, by.x = "organism_biosamp", by.y = "organism_name_passed_to_worms_worms", all.y = T)


#CLEAN UP FULL MARINE DATAFRAME A LITTLE ------
#look at extinct species
df.marine %>% filter(isExtinct_worms == 1) %>% select(organism_biosamp,project_description_bioprj,contains("worms")) %>% View()
  #most are wrong genus (e.g. grapes), remove these
df.marine <- df.marine %>% filter(organism_biosamp != "Silene") %>% filter(organism_biosamp != "Oryza") %>% filter(organism_biosamp != "Musa") %>% filter(organism_biosamp != "Vitis") %>% filter(organism_biosamp != "Patu")
df.marine %>% filter(isExtinct_worms == 1) %>% select(organism_biosamp,project_description_bioprj,contains("worms")) %>% View()
#remove some cols we don't care about
df.marine <- df.marine %>% select(-contains("relevance")) %>% select(-isExtinct_worms) %>% select(-lsid_worms)
#check if two cols are redundant
df.marine %>% mutate(check = ifelse(valid_authority_worms == authority_worms, "exact match", "not a match")) %>% filter(check == "not a match") %>% dim()
  #distinct cols, keep both
#rename and reorder some cols with better names
df.marine <- df.marine %>% rename(download_batch_bioprj = X_bioprj) %>% rename(type_of_name_match_asked_for_worms = lookuptype_worms) %>% 
    rename(type_of_name_match_returned_worms = match_type_worms) %>% rename(date_modified_worms = modified_worms) %>% 
  mutate(organism_name_passed_to_worms = organism_biosamp) %>% select(contains("bioprj"),everything()) %>% select(-iteration_worms,-X.1)
#check if values are all the same in a col.
df.marine %>% group_by(project_type_bioprj) %>% summarise(n=n())
  #they are, remove this col.
df.marine <- df.marine %>% select(-project_type_bioprj)
#check if this col. has anything useful
df.marine %>% group_by(extlinks_sra) %>% summarise(n=n())
  #it doesn't, remove
df.marine <- df.marine %>% select(-extlinks_sra)
#more reordering
df.marine <- df.marine %>% select(link,contains("bioprj"),contains("biosamp"),contains("sra"),contains("taxonomy"),contains("worms"),everything())

#change NA's in lat/long to something else so they don't disappear when filtering later
df.marine <- df.marine %>% mutate(lat_long_biosamp = ifelse(is.na(lat_long_biosamp)==F,lat_long_biosamp,"not given in NCBI"))
df.marine %>% group_by(lat_long_biosamp) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head()


#write out list of bioproject info for bioprojects that have marine species and sequence data that we want
write.csv(df.marine, paste0(outdir, "/NCBI_bioprojects__not-human_bacteria_viral_metgenome__-MARINEexactmatchesbiosampspeciesnames-BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"))




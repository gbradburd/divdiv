#idea: filter to list of marine datasets to look at by hand that we actually want to include

#load libraries ------
library(dplyr)    #data handling
library(tidyr)    #data handling
library(ggplot2)  #graphing

rm(list=ls())
gc()

#define list of NCBI metadata that matched a marine species (in WoRMS)
ncbi.marine = "/Users/rachel/divdiv/data/NCBI_metadata/NCBI_bioprojects__not-human_bacteria_viral_metgenome__-MARINEexactmatchesbiosampspeciesnames-BioProjectsANDSRAANDBioSamplesANDTaxonomy-WITHseqdataofinterest.csv"

#define out dir
outdir = "/Users/rachel/divdiv/data/NCBI_metadata"

#get some overall stats quick before filtering
df <- read.csv(ncbi.marine, header = T, stringsAsFactors = F) %>% dplyr::select(-X)
df %>% group_by(lat_long_biosamp) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head()

#how many bioprojects?
df %>% group_by(project_acc_bioprj) %>% summarise(n=n()) %>% nrow()
  #930 11-21-19; 1336 12-2-2020

#how many species?
df %>% group_by(organism_biosamp) %>% summarise(n=n()) %>% nrow()
  #1223 11-21-19; 1738 12-2-2020

#how many sequence records?
df %>% group_by(uid_sra) %>% summarise(n=n()) %>% nrow()
  #72,973 unique SRA ids 11-21-19; 114,865 12-2-2020

#total records with 10+ indivs (regardless of lat/long)
df %>% group_by(biosample_acc_sra,project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  group_by(project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  filter(n > 9) %>% nrow()
  #474 11-20-19 bioprj/species combos (regardless of lat/long status); 749 12-2-2020

df %>% group_by(biosample_acc_sra,project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  group_by(project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  filter(n > 9) %>%
  group_by(organism_biosamp) %>% summarise(n=n()) %>%
  nrow()
  #251 species 11-20-19; 400 12-2-2-2020



# ************************************************************************************************************************************
# ************************************************************************************************************************************
# GET LIST OF POTENTIAL DATASETS - lat/long in NCBI or not, 10+ indivs, 2+ locations per bioprj/species, marine -----------

# filter to 10+ indivs per dataset -----------------------------------------------------------------------------------------------
#group at level of biorpj/species and get at least 10 indivs per combo
filterlist.10indivs <- df %>% filter(lat_long_biosamp != "unknown") %>%
  group_by(biosample_acc_sra,project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  group_by(project_acc_bioprj,organism_biosamp) %>% summarise(n=n()) %>% 
  filter(n > 9) %>% mutate(link = paste(project_acc_bioprj, organism_biosamp, sep = "_")) %>% 
  mutate(link = gsub(pattern = " ", replacement = "-", link))
nrow(filterlist.10indivs)
  #703 bioprj/species combos with at least 10 indivs

#get only above list from full records
df.filtered.10indivs <- df %>% filter(link %in% filterlist.10indivs$link) %>% 
  filter(lat_long_biosamp != "unknown")
nrow(df.filtered.10indivs)
  #99,065 12-2-2020 sequence records
df.filtered.10indivs %>% group_by(biosample_acc_sra) %>% summarise(n=n()) %>% nrow()
  #92,673 12-2-2020 indivs



# now filter based on number of unique lat/long coords per paper ------------------------------------------------------------------

#get number of unique lat/long coords per bioprj/species combo ( n col ), all of which combos have at least 10 indivs

#check that df contains only entries with lat/long = "not given in NCBI" or numerical values at this point (should be 0 if so)
df.filtered.10indivs %>% filter(lat_long_biosamp == "unknown") %>% nrow()

#if above is true, then assign "numerical" vs "MIA" to each row for lat/long status
#keep datasets where all samples are:
#missing lat/long in NCBI -OR-
#where there is 1 lat/long location given and the rest are missing -OR-
#where all lat/long are given and there are at least 2 distinct locations
filterlist.10indivsAND2latlong <- df.filtered.10indivs %>% 
  mutate(check = ifelse(lat_long_biosamp == "not given in NCBI", "MIA", "numerical")) %>% 
  group_by(biosample_acc_sra,project_acc_bioprj,organism_biosamp,lat_long_biosamp,check) %>% summarise(n=n()) %>% 
  group_by(project_acc_bioprj,organism_biosamp,lat_long_biosamp,check) %>% summarise(n=n()) %>% 
  group_by(project_acc_bioprj,organism_biosamp,check) %>% summarise(n.locales=n()) %>% 
  ungroup() %>% group_by(project_acc_bioprj,organism_biosamp) %>% 
  mutate(check2 = n())

#only records we want to toss are those where all samples in dataset have a lat/long and it is the same lat/long for all samples
#those cases are numerical (samples have lat/longs), 1 (there is only 1 location in the dataset), 1 (all samples have the same lat/long status aka within the dataset all samples are missing lat/long or all samples have lat/long)
filterlist.10indivsAND2latlong <- filterlist.10indivsAND2latlong %>% 
  mutate(status = ifelse(check == "numerical" & n.locales < 2 & check2 < 2, "toss","keep")) %>% 
  filter(status == "keep") %>% 
  mutate(link = paste(project_acc_bioprj, organism_biosamp, sep = "_")) %>% 
  mutate(link = gsub(pattern = " ", replacement = "-", link))
nrow(filterlist.10indivsAND2latlong)
  #690 bioprj/species combos with at least 10 indivs and 2+ lat/long coords. in NCBI or missing coords. in NCBI that might be in a published paper

filterlist.10indivsAND2latlong %>% ggplot() +
  geom_histogram(aes(x=n.locales), binwidth = 1) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x = "N locations", y = "Frequency",
       title = paste("Number of unique lat/long locations per bioprj/species combo with > 9 indivs."))

# get only above list from full records
df.filtered.10indivsAND2latlong <- df %>% filter(link %in% filterlist.10indivsAND2latlong$link) %>% 
  filter(lat_long_biosamp != "unknown")
nrow(df.filtered.10indivsAND2latlong)
  #96,760 sequences for 10+ indivs and 2+ lat/long (or missing lat/longs)

#how many individuals?
df.filtered.10indivsAND2latlong %>% group_by(biosample_acc_sra) %>% summarise(n=n()) %>% nrow()
  #90,381 indivs 12-2-2020



# EXPLORE SPECIES - and chuck some we def don't want -----------------------------------------------------------------------

#filter out some species we def don't want
df.filtered.10indivsAND2latlong.goodspecies <- df.filtered.10indivsAND2latlong %>% 
  filter(!grepl("Saccharomyces cerevisiae",organism_biosamp)) %>%
  filter(!grepl("Brassica",organism_biosamp)) %>% 
  filter(!grepl("Malvaceae",organism_biosamp)) %>% 
  filter(!grepl("Lycium",organism_biosamp)) %>% 
  filter(!grepl("Carex",organism_biosamp)) %>% 
  filter(!grepl("Polygonaceae",organism_biosamp)) %>% 
  filter(!grepl("Apocynaceae",organism_biosamp)) %>%
  filter(!grepl("Fungi",organism_biosamp)) %>% 
  filter(!grepl("Junco",organism_biosamp)) %>% 
  filter(!grepl("Limosa lapponica baueri",organism_biosamp))
  
nrow(df.filtered.10indivsAND2latlong)
nrow(df.filtered.10indivsAND2latlong.goodspecies)
  #80230 sequences 12-4-2020  

#which species do we have left?
#df.filtered.10indivsAND2latlong.goodspecies %>% group_by(organism_biosamp) %>% summarise(n=n()) %>% View()



# format for use on Google Drive as master dataset list --------------------------------------

#keep columns we want in same order as master divdiv dataset list
records.clean <- df.filtered.10indivsAND2latlong.goodspecies %>% dplyr::select(link,biosample_acc_sra,lat_long_biosamp,
                                                                               scientificname_worms,isMarine_worms,isBrackish_worms,isFreshwater_worms,
                                                                               isTerrestrial_worms,kingdom_worms,phylum_worms,class_worms,
                                                                               order_worms,family_worms,genus_worms,
                                                                               project_acc_bioprj,project_title_bioprj,project_description_bioprj,organism_biosamp,
                                                                               date_this_sequence_first_added_to_our_records) %>% 
  rename("date_dataset_first_added_to_this_list" = "date_this_sequence_first_added_to_our_records")

#count number of unique biosamps per bioprj/species combo
records.clean <- records.clean %>% group_by(link) %>% mutate(n_biosample_accs = length(unique(biosample_acc_sra))) %>% 
  ungroup() %>% dplyr::select(-biosample_acc_sra) %>% distinct()

#add lat/long status column
records.clean %>% filter(lat_long_biosamp == "unknown") %>% nrow() #should be 0
records.clean <- records.clean %>% mutate(has_latlong_in_NCBI.temp = ifelse(lat_long_biosamp == "not given in NCBI", "n", "y")) %>%
  group_by(link,has_latlong_in_NCBI.temp) %>% mutate(has_latlong_in_NCBI.temp2 = length(unique(has_latlong_in_NCBI.temp))) %>% 
  ungroup() %>% 
  mutate(has_latlong_in_NCBI = ifelse(has_latlong_in_NCBI.temp2 == 2, "mixed", has_latlong_in_NCBI.temp)) %>% 
  dplyr::select(-has_latlong_in_NCBI.temp, -has_latlong_in_NCBI.temp2, -lat_long_biosamp) %>% distinct()
records.clean %>% dplyr::distinct(link) %>% nrow()
  #538 datasets as of 12-4-2020

#add on some empty columns needed for master spreadsheet
records.clean <- records.clean %>% mutate(is_salmonoid = "", commenter = "", comments = "", commenter_spatialdata = "", 
                                          comments_about_spatial_data_location = "", keepinrunning_YN = "NA",
                                          link_to_published_paper = "", datathon_project_index = "",
                                          `Common Name` = "",  Category = "", )

#final ordering of cols
records.clean <- records.clean %>% dplyr::select(link,project_acc_bioprj,project_title_bioprj,project_description_bioprj,organism_biosamp,
                                                 n_biosample_accs,date_dataset_first_added_to_this_list,has_latlong_in_NCBI,
                                                 is_salmonoid,commenter,comments,commenter_spatialdata,comments_about_spatial_data_location,
                                                 keepinrunning_YN,link_to_published_paper,datathon_project_index,scientificname_worms,
                                                 `Common Name`,Category,isMarine_worms,isBrackish_worms,isFreshwater_worms,isTerrestrial_worms,
                                                 kingdom_worms,phylum_worms,class_worms,order_worms,family_worms,genus_worms)

#write out final list of new datasets to add to master divdiv dataset list/spreadsheet
write.csv(records.clean, paste0(outdir, "/working_list_marine_projects_with_10indivs-12-4-2020.csv"), row.names = FALSE)

#worked with final file on Google Drive for rest of project to keep track of datasets in and out, collect trait data, etc.
#note this script doesn't perfectly recreate "working_list_marine_projects_with_10indivs-12-4-2020.csv" bc we were already working with this file on google drive
#and then updated after pulling possible datasets from NCBI a second time





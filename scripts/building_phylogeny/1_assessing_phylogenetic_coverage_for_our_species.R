#idea - explore how many species in divdiv match TimeTree

#load libraries ------
library(dplyr)      #data handling
library(tidyr)      #data handling
library(ggplot2)    #graphing



#library(ape)        #phylogenetic trees
#library(ggtree)     #phylogenetic trees (devtools::install_github("YuLab-SMU/ggtree"))
#library(treeio)     #reading in phylogenetic trees (devtools::install_github("GuangchuangYu/treeio")
#library(tidytree)   #getting summary info about trees
#library(taxize)     #looking up linnaean info for tip labels

rm(list=ls())
gc()


# READ IN PHYLO TREE ---------------------------------------------------
#read in tree - from http://timetree.org/book (Timetree of Life V2015)
#with taxize lineage info tacked on, and no bacteria or virus tips
load(file="/Users/rachel/Desktop/DivDiv/divdiv_collecting_genetic_data/plotting_on_phylogeny/tree_files/TimeTree_140K_NCBI_Taxon_IDs-withtaxizeinfo.Robj")

#view summary about tree
tree

#get list of all lineage info in tree
intree <- tree@extraInfo %>% dplyr::select(node, species, everything())
intree <- intree %>% filter(is.na(species)==F) %>% mutate(in_timetree = "yes") %>% mutate(label = species)
str(intree)


# READ IN OUR DATA -----------------------------------------------------
#get master spreadsheet from Drive of species we want / that are "in" analyses
ourspp <- googledrive::shared_drive_find(pattern = "^divdiv$")
ourspp <- googledrive::drive_ls(path = ourspp, pattern = "working_datasheets", recursive = FALSE)
ourspp <- googledrive::drive_ls(path = ourspp, pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
ourspp$name
ourspp <- googlesheets4::range_read(ourspp, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#keep cols we want and only datasets that are "in"
ourspp <- ourspp %>% 
  filter(grepl("yes|Yes",keepinrunning_YN)) %>%
  dplyr::select(organism_biosamp, kingdom_worms, phylum_worms, class_worms, order_worms, family_worms, genus_worms) %>%
  mutate(organism_biosamp = unlist(organism_biosamp))
str(ourspp)

#grab just the linnaean rank columns for datasets that we want to use 
ourspp <- ourspp %>% select(organism_biosamp,kingdom_worms,phylum_worms,class_worms,order_worms,family_worms,genus_worms) %>% 
  unique(.) %>% mutate(in_our_NCBI_data = "yes") %>% mutate(label = organism_biosamp)
dim(ourspp)
str(ourspp)
#change one subspecies to species
ourspp$label[ourspp$organism_biosamp == "Halichoerus grypus atlantica"] = "Halichoerus grypus"


# JOIN OUR SPECIES LIST TO TREE ----------------------------------------------
#join our labels/data onto tip labels
df <- merge(intree, ourlabels, by = "label", all = T) %>% 
  mutate(in_timetree = ifelse(is.na(in_timetree)==T, "no", in_timetree)) %>% 
  mutate(in_our_NCBI_data = ifelse(is.na(in_our_NCBI_data)==T, "no", in_our_NCBI_data)) %>% 
  dplyr::select(label,in_timetree,in_our_NCBI_data,everything())


# EXPLORE MISSING TIPS -------------------------
#look at which species didn't find a match in tree
#(in_timetree == "no")
df %>% group_by(in_timetree,in_our_NCBI_data) %>% summarise(n=n())

df %>% filter(in_our_NCBI_data == "yes") %>% dplyr::select(label,in_timetree,in_our_NCBI_data,contains("worms")) %>% 
  arrange(in_our_NCBI_data,in_timetree) %>% View()


#looked at these species online in TimeTree tool directly and in lit. to understand who they were related to and 
#how we could sub them appropriately into Time Tree. 
#recorded these decisions in decisions in /data/phylo/missing_species_to_add_in.csv




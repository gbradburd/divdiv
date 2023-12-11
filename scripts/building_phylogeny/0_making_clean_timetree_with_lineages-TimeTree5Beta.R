#idea: read in Time Tree phylogeny, look up taxonomic lineage info (e.g. class, order, family, etc.) for each tip, 
#       remove all tips that aren't plants or animals, save, for use in collapse_phy_and_plot_species_of_interest.R

# timetree 2015 version originally downloaded from: http://timetree.org/book near bottom of page
# taxonomic lineage info came from looking tip labels up using taxize (which accesses NCBI taxonomic database)


#load libraries ------
library(dplyr)      #data handling
library(tidyr)      #data handling
library(ggplot2)    #graphing
library(ape)        #phylogenetic trees
library(ggtree)     #phylogenetic trees (devtools::install_github("YuLab-SMU/ggtree"))
library(treeio)     #reading in phylogenetic trees (devtools::install_github("GuangchuangYu/treeio")
library(tidytree)   #getting summary info about trees
library(taxize)     #look up lineage info from species names
#library(phytools)  #phylogenetic trees (devtools::install_github("liamrevell/phytools"))


# *****************************************
#### SECTION 1 - look up lineage info for tip labels in TimeTree ####
# *****************************************
rm(list=ls())
gc()

#read in time tree ---------------------------------------------------------------------------------------------------------------------
tree <- phytools::read.newick("data/phylo/TimeTree_140K_NCBI_Taxon_IDs.nwk")
tree

#convert to treedata object 
tree <- tidytree::as.treedata(tree)
tree

#get tip labels
tiplabels <- tree@phylo$tip.label %>% as.data.frame(.) %>% rename("tip.label" = ".") %>% mutate(tip.label = gsub("_", " ", tip.label))
head(tiplabels)
#get just vector of tip labels to send through taxize - so we can get lineage info for each tip
tiplabelsFORtaxize <- tiplabels$tip.label
str(tiplabelsFORtaxize)


#look up lineage info for each tip in taxize --------------------------------------------------------------------------------------------------
#need to look this up bc if we have only one species in a family and that species isn't in tree, 
#then we don't have a family from a species that IS in the tree to assign it to, unless we know families for/associated with all tips of tree

#function to take in NCBI taxID (as character) and return long list of taxonomic names from NCBI or a line with NAs if/when it breaks
lookuptiplabels <- function(taxID) {
  
  result <- tryCatch(
    {
      taxID = tiplabelsFORtaxize[loop.iter]
      result <- classification(taxID, db = 'ncbi', api_key = "3137dd8892dc5b2e9be454fe5e42eab93309", callopts=list(http_version = 0L))
      result <- t(matrix(unlist(result), byrow = T, nrow = 3)) %>% as.data.frame(.) %>% mutate(taxon = taxID)
    },
    error=function(cond) {
      print(paste0("Error for taxID: ", taxID, "; loop iter ", loop.iter))
      cat(message(cond), "\n")
      result <- data.frame(V1 = NA, V2 = NA, V3 = NA, taxon = taxID)
      return(result)
      },
    finally={
      Sys.sleep(time=0.2)
    }
  )    
  return(result)
  
}

#make blank df 
output <- data.frame(V1 = NA, V2 = NA, V3 = NA, taxon = NA)
httr::set_config(httr::config(http_version = 0))
#use above function to look up tip labels and make master df
for(loop.iter in 1:length(tiplabelsFORtaxize)) {
  taxID = tiplabelsFORtaxize[loop.iter]
  result <- lookuptiplabels(taxID)
  output <- rbind(output,result)
  print(paste0("done with taxID ", loop.iter, " of ", length(tiplabelsFORtaxize)))
}

#check if we have all the values
output %>% filter(is.na(taxon)==F) %>% filter(is.na(V1)==T & is.na(V2)==T & is.na(V3)==T) %>% View()
redo <- output %>% filter(is.na(taxon)==F) %>% filter(is.na(V1)==T & is.na(V2)==T & is.na(V3)==T) %>% dplyr::select(taxon)
tiplabelsFORtaxize <- redo$taxon
#rerun those we missed
output <- data.frame(V1 = NA, V2 = NA, V3 = NA, taxon = NA)
#use above function to look up tip labels and make master df
for(loop.iter in 1:length(tiplabelsFORtaxize)) {
  taxID = tiplabelsFORtaxize[loop.iter]
  result <- lookuptiplabels(taxID)
  output <- rbind(output,result)
  print(paste0("done with taxID ", loop.iter, " of ", length(tiplabelsFORtaxize)))
}

#get ready to tack onto tree tips --------------
#read in output from taxize
output <- output %>% mutate(taxon = as.character(taxon))

#filter out the no-rank groups and clade (bc clade is used multiple times with different values for some species)
output <- output %>% filter(V2 != "no rank") %>% filter(V2 != "clade")
output %>% group_by(taxon) %>% summarise(n.lineage.groups=n()) %>% group_by(n.lineage.groups) %>% summarise(n=n()) %>% as.data.frame()

#convert long to wide
#(it's ok that some species have more classification levels returned than others - this will just fill missing things with NA)
output.w <- output %>% dplyr::select(-V3) %>% unique(.) %>% pivot_wider(., names_from = V2, values_from = V1)

write.csv(output.w, "data/phylo/TimeTree_140K_NCBI_Taxon_IDs_linean_ranks_for_tips_from_taxize-wide.csv")



# *****************************************
#### SECTION 2 - join lineage info onto tree ####
# *****************************************
rm(list=ls())
gc()

#read in time tree
tree <- phytools::read.newick("data/phylo/TimeTree_140K_NCBI_Taxon_IDs.nwk")
tree
#convert to treedata object 
tree <- tidytree::as.treedata(tree)
tree

#get tip labels from timetree
tiplabels <- tree@phylo$tip.label %>% as.data.frame(.) %>% rename("tip.label" = ".") %>% mutate(tip.label = gsub("_", " ", tip.label))
head(tiplabels)

#get lineage info that we looked up for timetree
lineageinfo <- read.csv("data/phylo/TimeTree_140K_NCBI_Taxon_IDs_linean_ranks_for_tips_from_taxize-wide.csv") %>% select(-X)
head(lineageinfo)
#add lineage info on to tip labels
dim(tiplabels)
dim(lineageinfo)
lbs <- merge(tiplabels, lineageinfo, by.x = "tip.label", by.y = "taxon", all.x = T) %>% rename("label" = "tip.label")
dim(lbs)

#formatting of tip labels in tree
tree@phylo$tip.label %>% as.data.frame(.) %>% rename("tip.label" = ".") %>% head()
#formatting of our species from NCBI list
head(lbs$label)
og.tiplabels <- tree@phylo$tip.label %>% as.data.frame(.) %>% rename("label" = ".")
head(og.tiplabels)
#do a test merge of our info onto tip labels (not the actual tree obj. yet)
#if everything is correct, both of below should have 0 rows, meaning all tip labels identically match between tree and our df
merge(og.tiplabels %>% mutate(check = "in.og.tip.labels"), lbs, by = "label", all.y = T) %>% filter(is.na(check)==T)
merge(og.tiplabels, lbs %>% mutate(check = "in.new.df"), by = "label", all.x = T) %>% filter(is.na(check)==T)

dim(lbs)
dim(og.tiplabels)

#convert our df to tibble so it will join onto tree
lbs <- as_tibble(lbs)
head(lbs)

#join our labels/info onto tree
tree <- full_join(tree, lbs, by = "label")
tree

#save tree WITH taxon info tacked on - so we can use it later
save(tree, file = "data/phylo/TimeTree_140K_NCBI_Taxon_IDs-withtaxizeinfo.Robj")



# *****************************************
#### SECTION 3 - drop phy tips to get just metazoa and viridiplantae ####
# *****************************************

rm(list = setdiff(ls(), lsf.str()))
gc()

# read in Time Tree with taxonomic rank info saved in @extraInfo
# note, tax info looked up using taxize::classification(x = timetree_tip_labels, db = 'ncbi') - making_phylogeny.R
load("plotting_on_phylogeny/tree_files/TimeTree_140K_NCBI_Taxon_IDs-withtaxizeinfo.Robj")



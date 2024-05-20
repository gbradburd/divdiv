#idea: use TimeTree with df tacked on of taxonomic rank info/names (using taxize),
#list of our divdiv species, and list of missing species and how to sub them into tree,
#and build a subsetted TimeTree with all divdiv species in it
#for species with only one species to sub in that taxonomic level (e.g. an Acropora is missing from TimeTree but there are 
#other Acropora in TT and only one Acropora sp. in divdiv, pick a random species tip and sub it in)
#for species in divdiv and missing in TT that are at same tax level as another species in divdiv,
#use lit to find closest/sister species that is in TT and sub in there

#original code written/inspired by Bruce Martin

#load libraries
library(phytools)
library(tidytree)
library(dplyr)
library(ggtree)
library(ggplot2)


rm(list=ls())
gc()
set.seed(123)

setwd("/Users/rachel/divdiv")


# load in data -----------
load('data/phylo/inputs_and_working/TimeTree_140K_NCBI_Taxon_IDs-withtaxizeinfo.Robj') #TimeTree plus tacked on df of looked up rank info/names from taxize
missing.spp <- read.csv('data/phylo/inputs_and_working/missing_species_to_add_in.csv') #list of species not in TimeTree that we need to add in
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
  #120 species as of last time ran this script
#update name of Exaiptasia pallida to Exaiptasia diaphana (now current/accepted name)
ourspp$organism_biosamp[ourspp$organism_biosamp == "Exaiptasia pallida"] = "Exaiptasia diaphana" 



#tree data is a bit unusual in that it's actual labels are numeric codes, and species names are stored in a separate dataset
#the general approach here will be to find labels corresponding to tips we want to swap with,
##then editing the species names in the dataset accordingly
#afterwards, will bind an extra tip on for row 19 of missing.spp (where we want to make a polytomy), which is kinda a special case

# reformat tree a bit -----------
#grab phylo object
phy <- tree@phylo
#grab tip dataset
tip.seq <- match(seq_len(Ntip(phy)),tree@extraInfo[['node']]) #just in case things are out of order
dat <- tree@extraInfo[tip.seq,]
#add tip labels to tip dataset
dat$label <- phy$tip.label
#above is df of tip label (label), node id, and all the taxonomic rank info that we looked up and tacked onto tips earlier using taxize



# cheat and add some of our species names into tree so we find matches that we know we want to use later
#assign our species to the single tip with genus label (no other sp in divdiv on/in this tip)
dat$species[dat$node == 102902] = "Zenarchopterus dunckeri"
#assign our species to the single tip with genus label (no other sp in divdiv on/in this tip)
dat$species[dat$node == 69093] = "Anthopleura elegantissima"
#assign our species to the single tip with genus label (no other sp in divdiv on/in this tip)
dat$species[dat$node == 69091] = "Nematostella vectensis"
#assign our species to the tip with this family label (no other sp in divdiv on/in this tip)
dat$species[dat$node == 97363] = "Pyroteuthis margaritifera"
#change one species with subspecies to just species so it matches TT
dat$species[dat$species == "Halichoerus grypus"] = "Halichoerus grypus atlantica"



# first swap tips
#let's first drop the "_worms" suffix for easier index handling
colnames(missing.spp) <- gsub('_worms$','',colnames(missing.spp))
#add on a species col for the ones we already picked specific subs out for
missing.spp$species <- missing.spp$species_subbed
#and where we already picked a specific species to use as the sub, switch phy_level_to_sub_at col to species 
#(or for where we want to handle sp. ourselves aka bind in a tip somewhere manually - usually to make polytomy)
missing.spp$species_subbed[missing.spp$species_subbed == ""] = NA
missing.spp$Phy_level_to_sub_in_at[!is.na(missing.spp$species_subbed)]<-'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Acropora prolifera"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Bathyraja panthera"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Bathyraja aleutica"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Bartholomea annulata"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Exaiptasia brasiliensis"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Exaiptasia diaphana"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Systellaspis debilis"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Acanthephyra purpurea"] = 'species'
missing.spp$Phy_level_to_sub_in_at[missing.spp$label == "Cranchia scabra"] = 'species'

for(i in 1:nrow(missing.spp)){
  #get the Linnaean rank we're using for the swap
  lin.rank<-missing.spp$Phy_level_to_sub_in_at[i]
  #find the tips within that rank
  indices <- dat[[lin.rank]]==missing.spp[i,lin.rank]&
    !is.na(dat[[lin.rank]])
  tips <- dat$node[indices]
  #note that the Linnaean rank might have been para/polyphyletic
  #if we don't care about monophyly, we can stop here
  # missing.spp$species_to_sub[i]<-sample(dat$label[tips],1)
  #but if we DO care about monophyly, we first find the most recent common ancestor...
  if(length(tips)>1){
    mrca <- getMRCA(phy,tips)
    #get ALL descendants of that common ancestor...
    des <- getDescendants(phy,mrca)
    des <- des[des<Ntip(phy)] # get rid of internal nodes and just keep tip "nodes"
    #we can use this to store what "unexpected" tips end up in the group, just so we can see how extensive this potential issue is
    missing.spp$monophy_tips[i] <- list(des[!(des%in%tips)])
    #then just do what we would've done before with the FULL monophyletic group
    #pick a random tip within the group
    missing.spp$species_subbed[i] <- sample(dat$label[des],1)
  }else if(length(tips)==1){
    #note that we already have mrca if there's only 1 tip!
    missing.spp$species_subbed[i] <- dat$label[tips]
  }
}
missing.spp$species_subbed
cbind(missing.spp[,1:4],lengths(missing.spp$monophy_tips)) #so it looks like there was some monophyly issues with the Bathymodiolus genus
dat$species[missing.spp$monophy_tips[[5]]] #in this case, it seems like the taxonomic annotations for some tips were just incomplete
#swap out some names in dat with our divdiv names
matches <- match(missing.spp$species_subbed,dat$label)
nas <- is.na(matches) #exclude the special case for now (building polytomy of Acroporas)...
dat$species[matches[!nas]] <- missing.spp$label[!nas]



# special case: binding in new tips ------------

#phytool's bind.tip() unfortunately alters the node indexing in the tree
#this makes phy and dat no longer match up correctly
#so made a custom wrapper to alter dat properly
#just make sure the tip.label is a UNIQUE name!!!
#also assumes dat is correctly ordered from node 1 to # of tips (we did this above in lines 22-23)
custom.bind.tip <- function(tree,
                          dat,
                          tip.label, #the label of the new tip
                          sp.name=tip.label, #in case you want the tip label to be different from the actual species name in dat
                          edge.length=NULL, #length of new tip; defaults to making tip ultramteric if rest of tree is
                          where, #the node to connect the edge to
                          position=0){ #how far to shift the base of the edge "rootward" (0 if you want a polytomy)
  #first actually bind the tip
  tree <- bind.tip(tree,tip.label,edge.length,where,position)
  #find the node index of the new tip using it's label (this is why the label must be unique)
  new.tip <- which(tree$tip.label==tip.label)
  #duplicate the row corresponding to the new tip in dat
  reps <- rep(1,nrow(dat))
  reps[new.tip] <- 2
  inds <- rep(seq_len(nrow(dat)),reps)
  dat <- dat[inds,]
  #alter the new tip's row with appropriate name/species label
  dat[new.tip,-1] <- NA
  dat$label[new.tip] <- tip.label
  dat$species[new.tip] <- sp.name
  #increment the nodes below the new tip's row by 1
  nodes.to.inc <- seq_len(nrow(dat))>new.tip
  dat$node[nodes.to.inc] <- dat$node[nodes.to.inc]+1
  list(tree=tree,dat=dat)
}

#first find tips you want to make polytomy with
#then find their MRCA
#then use custom function to add a branch in 

#bind Acropora prolifera in
tips <- dat$node[dat$species%in%c('Acropora millepora','Acropora palmata')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Acropora prolifera",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Bathyraja panthera in
tips <- dat$node[dat$species%in%c('Bathyraja parmifera','Pavoraja nitida')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Bathyraja panthera",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Bathyraja aleutica in
tips <- dat$node[dat$species%in%c('Bathyraja parmifera','Pavoraja nitida')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Bathyraja aleutica",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Bartholomea annulata in
tips <- dat$node[dat$family%in%c('Aiptasiidae','Edwardsiidae')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Bartholomea annulata",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Exaiptasia brasiliensis in
tips <- dat$node[dat$family%in%c('Aiptasiidae','Edwardsiidae')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Exaiptasia brasiliensis",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Exaiptasia diaphana (pallida) in
# !!! NOTE !!! - using Exaiptasia diaphana instead of Exaiptasia pallida bc E. diaphana is accepted/current
tips <- dat$node[dat$family%in%c('Aiptasiidae','Edwardsiidae')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Exaiptasia diaphana",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Systellaspis debilis in
tips <- dat$node[dat$species%in%c('Pandalus montagui','Procaris ascensionis')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Systellaspis debilis",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Acanthephyra purpurea in
tips <- dat$node[dat$species%in%c('Pandalus montagui','Procaris ascensionis')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Acanthephyra purpurea",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat

#bind Cranchia scabra in
tips <- dat$node[dat$family%in%c('Pyroteuthidae','Bathyteuthidae')]
mrca <- getMRCA(phy,tips)
tmp <- custom.bind.tip(phy,
                     dat,
                     tip.label="Cranchia scabra",
                     where=mrca)
phy <- tmp$tree
dat <- tmp$dat


#if you didn't want to make a polytomy, just find where the mrca of the clade you want the new tip to be SISTER to
#then play around with the position=, which will slide the base of the new tip down the stem of that sister clade


# final steps -----------
#now we have a dataset with all the species we want, which gives the tip labels each species corresponds to in phy!
#this works if you have a single tip per each species you wish to keep
labels.to.keep<-match(ourspp$organism_biosamp,dat$species)
#otherwise, you wanna go with this
# labels.to.keep<-lapply(ourspp$organism_biosamp,function(ii) which(dat$species==ii))
# lengths(labels.to.keep) #thankfully, each species you want to keep correspond to either 0 or 1 labels
#which divdiv species are still not in the tree
ourspp$organism_biosamp[is.na(labels.to.keep)] #should return 0 aka they should all be in timetree as a tip now
#now we're ready to subset the tree using "keep.tip"
labels.to.keep <- dat$label[labels.to.keep]
labels.to.keep <- labels.to.keep[!is.na(labels.to.keep)]
phy <- keep.tip(phy,labels.to.keep)
#for better "readability"
phy$tip.label <- dat$species[match(phy$tip.label,dat$label)]

#view
plot(phy,cex=0.6)


# build df of label info for tree plotting aesthetics and order it same as phy tips
finaltreetips <- phy$tip.label %>% as.data.frame() %>% rename("tip.label"=".") %>% mutate(tiporder = 1:n())
lbs <- merge(ourspp, finaltreetips, by.x = "organism_biosamp", by.y = "tip.label", all.y = T) %>% 
  merge(., missing.spp, by.x = "organism_biosamp", by.y = "label", all.x = T) %>% 
  distinct() %>% mutate(sub_category = ifelse(is.na(sub_category), "sp in TT directly", sub_category)) %>% 
  arrange(tiporder)

p <- ggtree(phy, size=0.2) #make object and set line thickness
p %<+% lbs + #leftjoin labels for aesthetics onto phy
  geom_tiplab(aes(fill = sub_category),
              size = 1.6, #text size
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.03, "lines"), # amount of padding around the labels
              label.size = 0,
              fontface = "bold") + # size of label border
  scale_x_continuous(limits = c(0,1800), breaks = seq(0,max(phy$edge.length)+200,200), expand = c(0, 0)) + #change limits max to make space for tip labels
  scale_fill_manual(name = "Tip status", values = c("#E5446D","#E8871E","#49D49D","#49D49D","#49D49D","#6CBEED")) + 
  theme(legend.position = c(0.3,0.85),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5, ),
        axis.line.x = element_line(size = 0.2), #build scale bar
        axis.ticks.x = element_line(), #build scale bar
        axis.text.x = element_text(size = 2)) #plot ranges between 0 and 1

#have to save manually using export -> Save as PDF, point and click
#8.5 x 6 inches, "/figures/phy_115_divdiv_species.pdf"

#save phy as Robj ---------------
save(phy, file="data/phylo/divdiv_phy_from_timetreebeta5.Robj")





#sandbox to visualize what's going on ----------
tips<-dat$node[grepl("Bathyraja",dat$species)]
tips<-dat$node[dat$species%in%c('Bathyraja parmifera','Sphyrna tiburo')]
tips<-dat$node[grepl("Bathyraja|Okamejei",dat$genus)] #
tips<-dat$node[grepl("Arhynchobatidae",dat$family)]
tips<-dat$node[dat$species%in%c('Bathyraja parmifera','Pavoraja nitida')]
tips<-dat$node[grepl("Acropora",dat$species)]
tips<-dat$node[grepl("Actiniaria",dat$order)]
tips<-dat$node[grepl("Nematostella",dat$genus)]
tips<-dat$node[grepl("Aiptasiidae",dat$family)]
tips<-dat$node[grepl("Actiniaria",dat$order)]
tips<-dat$node[grepl("Bartholomea",dat$genus)]
tips<-dat$node[grepl("Anthozoa",dat$class)]
tips<-dat$node[grepl("Edwardsiidae|Actiniidae|Aiptasiidae|Parazoanthidae",dat$family)]
tips<-dat$node[grepl("Cephalopoda",dat$class)]
tips<-dat$node[grepl("Caridea|Polychelida",dat$infraorder)]
tips<-dat$node[grepl("Pleocyemata",dat$suborder)]
tips<-dat$node[grepl("Pandalus montagui|Polycheles aculeatus|Procaris ascensionis",dat$species)]
tips<-dat$node[grepl("Oegopsina",dat$suborder)]
tips<-dat$node[grepl("Teuthida",dat$order)] #plot at suborder
mrca <- getMRCA(phy,tips)
des <- getDescendants(phy,mrca)
toplot <- dat[dat$node %in% des,]
toplot <- keep.tip(phy,toplot$label)
toplot$tip.label <- dat$suborder[match(toplot$tip.label,dat$label)]
plot(toplot,cex=0.4)
#end sandbox


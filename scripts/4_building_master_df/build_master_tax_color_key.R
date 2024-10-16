#idea: build master tax group color key for divdiv species

#load libraries 
library(dplyr)
library(tidyr)

rm(list=ls())
gc()

#get phy so we can get species in each clade
load("data/phylo/divdiv_phy_from_timetreebeta5.Robj")

#get list of divdiv species in each clade
#bony fishes
des <- ape::getMRCA(phy,c("Anguilla rostrata","Sebastes diaconus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)] # get rid of internal nodes and just keep tip "nodes"
df1 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Bony fishes")
#sauropsida
des <- ape::getMRCA(phy,c("Pygoscelis papua","Caretta caretta")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df3 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Birds and reptiles")
#mammals
des <- ape::getMRCA(phy,c("Phocoena sinus","Halichoerus grypus atlantica")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df4 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Mammals")
#chondrichthyes
des <- ape::getMRCA(phy,c("Bathyraja parmifera","Sphyrna tiburo")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df2 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Cartilaginous fishes")
#echinoderms
des <- ape::getMRCA(phy,c("Apostichopus californicus","Pisaster ochraceus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df5 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Echinoderms")
#molluscs
des <- ape::getMRCA(phy,c("Bathymodiolus platifrons","Nautilus pompilius")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df6 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Molluscs")
#crustacea
des <- ape::getMRCA(phy,c("Callinectes sapidus","Isocladus armatus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df7 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Crustaceans")
#cnidarians
des <- ape::getMRCA(phy,c("Galaxea horrescens","Ectopleura larynx")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df8 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Cnidarians")
#vascular plants
des <- ape::getMRCA(phy,c("Rhizophora mangle","Laguncularia racemosa")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df9 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "Vascular plants")
#decided to make ochrophyta into other (only 1 divdiv sp in clade in end)
#ochrophyta
#des <- ape::getMRCA(phy,c("Fucus vesiculosus","Sargassum muticum")) %>% phytools::getDescendants(phy,.)
#des <- des[des<ape::Ntip(phy)]
#df10 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "ochrophyta")

#merge and tack on colors
taxcladelbs <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9)
taxcladelbs <- merge(taxcladelbs,
                     phy$tip.label %>% as.data.frame() %>% dplyr::rename("species"="."),
                     by = "species", all.x = T, all.y = T)
#fill in NA's with other
taxcladelbs$taxclade[is.na(taxcladelbs$taxclade) == T] = "Other"

taxcladelbs <- merge(taxcladelbs, 
                     data.frame(taxclade = c("Bony fishes", "Birds and reptiles", "Mammals",
                                             "Cartilaginous fishes", "Echinoderms", 
                                             "Molluscs", "Crustaceans", "Cnidarians", "Vascular plants", 
                                             "Other"),
                                cladecolor = c("#A6CEE3","#1F78B4","#B2DF8A",
                                               "#33A02C","#FB9A99",
                                               "#E31A1C","#FDBF6F","#FF7F00","#772E25",
                                               "#a4afb0")),
                     by = "taxclade", all.x = T)

taxcladelbs %>% group_by(taxclade, cladecolor) %>% summarise(n=n())

taxcladelbs <- taxcladelbs %>% dplyr::select(species, taxclade, everything())

write.csv(taxcladelbs, "data/master_tax_color_key.csv", row.names = FALSE, quote = FALSE)


#idea: build master tax group color key for divdiv species

#load libraries 
library(dplyr)
library(tidyr)

rm(list=ls())
gc()

#set taxonomic/clade colors (to match phy: figure-phylogeny_of_species.R)
cladeCols <- RColorBrewer::brewer.pal(10,"Paired")

#get phy so we can get species in each clade
load("data/phylo/divdiv_phy_from_timetreebeta5.Robj")

#get list of divdiv species in each clade
#bony fishes
des <- ape::getMRCA(phy,c("Anguilla rostrata","Sebastes diaconus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)] # get rid of internal nodes and just keep tip "nodes"
df1 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "bony fishes")
#chondrichthyes
des <- ape::getMRCA(phy,c("Bathyraja parmifera","Sphyrna tiburo")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df2 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "chondrichthyes")
#sauropsida
des <- ape::getMRCA(phy,c("Pygoscelis papua","Caretta caretta")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df3 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "sauropsida")
#mammals
des <- ape::getMRCA(phy,c("Phocoena sinus","Halichoerus grypus atlantica")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df4 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "mammals")
#echinoderms
des <- ape::getMRCA(phy,c("Apostichopus californicus","Pisaster ochraceus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df5 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "echinoderms")
#molluscs
des <- ape::getMRCA(phy,c("Bathymodiolus platifrons","Nautilus pompilius")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df6 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "molluscs")
#crustacea
des <- ape::getMRCA(phy,c("Callinectes sapidus","Isocladus armatus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df7 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "crustacea")
#cnidarians
des <- ape::getMRCA(phy,c("Galaxea horrescens","Ectopleura larynx")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df8 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "cnidarians")
#vascular plants
des <- ape::getMRCA(phy,c("Rhizophora mangle","Laguncularia racemosa")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df9 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "vascular plants")
#ochrophyta
des <- ape::getMRCA(phy,c("Fucus vesiculosus","Sargassum muticum")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df10 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "ochrophyta")

#merge and tack on colors
taxcladelbs <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)
taxcladelbs <- merge(taxcladelbs, 
                     data.frame(taxclade = c("bony fishes", "chondrichthyes", "sauropsida", "mammals", "echinoderms", 
                                             "molluscs", "crustacea", "cnidarians", "vascular plants", "ochrophyta"),
                                cladecolor = c("#A6CEE3","#33A02C","#1F78B4","#B2DF8A","#FB9A99",
                                               "#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")),
                     by = "taxclade", all.x = T)


taxcladelbs <- taxcladelbs %>% dplyr::select(species, taxclade, everything())

write.csv(taxcladelbs, "data/master_tax_color_key.csv", row.names = FALSE)


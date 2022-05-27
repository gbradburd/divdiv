#idea: use TimeTree to plot phylogeny of divdiv species

#load libraries
library(phytools)
library(tidytree)
library(dplyr)
library(ggtree)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
gc()

# load phylogeny (built from /DivDiv/divdiv_collecting_genetic_data/marine/scripts/building_marine_tree_from_TimeTree.R)
load(file="data/phylo/divdiv_phy_from_timetreebeta5.Robj")

# ! NOTE ! - don't remember the best way to subset tree to just species we have popgen data so far
# but here is list of which species those are
sps <- read.csv("data/master_df.csv") %>% filter(is.na(s.wish)==FALSE) %>% dplyr::select(species) %>% mutate(species = gsub("_"," ",species))


# ! NOTE ! - right now below builds a phy for all sps in divdiv

# build df of label info for tree plotting aesthetics and order it same as phy tips
finaltreetips <- phy$tip.label %>% as.data.frame() %>% rename("tip.label"=".") %>% mutate(tiporder = 1:n())

cladeCols <- brewer.pal(10,"Paired")

p <- ggtree(phy, size=0.2,layout="circular") #make object and set line thickness
p %<+% finaltreetips + #leftjoin labels for aesthetics onto phy
  geom_tiplab2(size = 2, #text size
              color = "black", # color for label font
              #geom = "label",  # labels not text
              #label.padding = unit(0.03, "lines"), # amount of padding around the labels
              #label.size = 0,
              fontface = "bold") + # size of label border
  scale_x_continuous(limits = c(0,1800), breaks = seq(0,max(phy$edge.length)+200,200), expand = c(0, 0)) + #change limits max to make space for tip labels
  theme(legend.position = c(0.3,0.85),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5, ),
        #axis.line.x = element_line(size = 0.2), #build scale bar
        axis.ticks.x = element_line(), #build scale bar
        axis.text.x = element_text(size = 2)) + #plot ranges between 0 and 1
	geom_hilight(node=MRCA(phy,"Anguilla rostrata","Sebastes diaconus"),fill=cladeCols[1],alpha=0.3) + #bony fishes
	geom_hilight(node=MRCA(phy,"Pygoscelis papua","Caretta caretta"),fill=cladeCols[2],alpha=0.3) + #sauropsida
	geom_hilight(node=MRCA(phy,"Phocoena sinus","Halichoerus grypus atlantica"),fill=cladeCols[3],alpha=0.3) + #mammals
	geom_hilight(node=MRCA(phy,"Bathyraja parmifera","Sphyrna tiburo"),fill=cladeCols[4],alpha=0.3) +  #chondrichthyes
	geom_hilight(node=MRCA(phy,"Apostichopus californicus","Pisaster ochraceus"),fill=cladeCols[5],alpha=0.3) + #echinoderms
	geom_hilight(node=MRCA(phy,"Bathymodiolus platifrons","Nautilus pompilius"),fill=cladeCols[6],alpha=0.3) + #molluscs
	geom_hilight(node=MRCA(phy,"Callinectes sapidus","Isocladus armatus"),fill=cladeCols[7],alpha=0.3) + #crustacea
	geom_hilight(node=MRCA(phy,"Galaxea horrescens","Ectopleura larynx"),fill=cladeCols[8],alpha=0.3) + #cnidarians
	geom_hilight(node=MRCA(phy,"Rhizophora mangle","Laguncularia racemosa"),fill=cladeCols[9],alpha=0.3) + #vascular plants
	geom_hilight(node=MRCA(phy,"Fucus vesiculosus","Sargassum muticum"),fill=cladeCols[10],alpha=0.3) #ochrophyta

ggsave(filename = "figures/phylo_distn.pdf", width = 15, height = 15, units = c("cm"))
	
#have to save manually using export point and click
#8.18 x 5.80 inches



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
load(file="../data/phylo/divdiv_phy_from_timetreebeta5.Robj")

havewishart <- read.csv("../data/master_df.csv") %>% mutate(havewishart = ifelse(is.na(s)==F, "yes", "no")) %>% 
  mutate(tip.label = gsub("_"," ",species)) %>% dplyr::select(tip.label, havewishart)
  
# build df of label info for tree plotting aesthetics and order it same as phy tips
phy <- drop.tip(phy,phy$tip.label[!phy$tip.label %in% havewishart$tip.label])
finaltreetips <- phy$tip.label %>% as.data.frame() %>% rename("tip.label"=".") %>% mutate(tiporder = 1:n())

# finaltreetips <- merge(finaltreetips, havewishart, by = "tip.label", all.x = T) %>% 
  # mutate(havewishart = ifelse(is.na(havewishart)==T, "no", havewishart))

cladeCols <- RColorBrewer::brewer.pal(10,"Paired")

p <- ggtree(phy, size=0.2,layout="circular") #make object and set line thickness
p %<+% finaltreetips + #leftjoin labels for aesthetics onto phy
  geom_tiplab(size = 2, #text size
              #color = "black", # color for label font
#              aes(color = havewishart),
              #geom = "label",  # labels not text
              #alpha = 0.5,
              #label.padding = unit(0.03, "lines"), # amount of padding around the labels
              #label.size = 0,
              fontface = "bold") + # size of label border
  scale_x_continuous(limits = c(0,1800), breaks = seq(0,max(phy$edge.length)+200,200), expand = c(0, 0)) + #change limits max to make space for tip labels
  scale_color_manual(values = c("black","red")) + #set text colors for if we have wishart or not
  theme(legend.position = "none",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        #axis.line.x = element_line(size = 0.2), #build scale bar
        axis.ticks.x = element_line(), #build scale bar
        axis.text.x = element_text(size = 2)) + #plot ranges between 0 and 1
	geom_hilight(node=MRCA(phy,"Engraulis encrasicolus","Sebastes diaconus"),fill=cladeCols[1],alpha=0.3) + #bony fishes
	geom_hilight(node=MRCA(phy,"Pygoscelis papua","Aythya marila"),fill=cladeCols[2],alpha=0.3) + #sauropsida
	geom_hilight(node=MRCA(phy,"Phocoena sinus","Halichoerus grypus atlantica"),fill=cladeCols[3],alpha=0.3) + #mammals
	geom_hilight(node=MRCA(phy,"Bathyraja panthera","Sphyrna tiburo"),fill=cladeCols[4],alpha=0.3) +  #chondrichthyes
	geom_hilight(node=MRCA(phy,"Paracentrotus lividus","Pisaster ochraceus"),fill=cladeCols[5],alpha=0.3) + #echinoderms
	geom_hilight(node=MRCA(phy,"Chlorostoma funebralis","Pyroteuthis margaritifera"),fill=cladeCols[6],alpha=0.3) + #molluscs
	geom_hilight(node=MRCA(phy,"Callinectes sapidus","Penaeus duorarum"),fill=cladeCols[7],alpha=0.3) + #crustacea
	geom_hilight(node=MRCA(phy,"Acropora palmata","Ectopleura larynx"),fill=cladeCols[8],alpha=0.3) + #cnidarians
	geom_hilight(node=MRCA(phy,"Rhizophora mangle","Laguncularia racemosa"),fill=cladeCols[9],alpha=0.3) + #vascular plants
#	geom_hilight(node=MRCA(phy,"Sargassum muticum","Sargassum muticum"),fill=cladeCols[10],alpha=0.3) #ochrophyta + 
	theme_void() #remove box
ggsave(filename = "../figures/phylo_distn.pdf", width = 15, height = 15, units = c("cm"))



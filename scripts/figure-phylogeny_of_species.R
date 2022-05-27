#idea: use TimeTree to plot phylogeny of divdiv species

#load libraries
library(phytools)
library(tidytree)
library(dplyr)
library(ggtree)
library(ggplot2)


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

p <- ggtree(phy, size=0.2) #make object and set line thickness
p %<+% finaltreetips + #leftjoin labels for aesthetics onto phy
  geom_tiplab(size = 1.6, #text size
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.03, "lines"), # amount of padding around the labels
              label.size = 0,
              fontface = "bold") + # size of label border
  scale_x_continuous(limits = c(0,1800), breaks = seq(0,max(phy$edge.length)+200,200), expand = c(0, 0)) + #change limits max to make space for tip labels
  theme(legend.position = c(0.3,0.85),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5, ),
        axis.line.x = element_line(size = 0.2), #build scale bar
        axis.ticks.x = element_line(), #build scale bar
        axis.text.x = element_text(size = 2)) #plot ranges between 0 and 1

#have to save manually using export point and click
#8.18 x 5.80 inches



#idea: create suppmat figures to show biotic trait coverage 

#load libraries -------
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list=ls())
gc()


setwd("/Users/rachel/divdiv")

#define dir to save figs in
figdir = "figures"

# get final, clean trait data we used
df <- read.csv("data/biotic/cleaned_numeric_biotic_traits.csv") %>% 
  dplyr::select(organism_biosamp, isBenthic, Body_Size, Generational_Structure, 
                Fecundity_EggSize, PLD_point2, ReturnToSpawningGround,
                Larval_feeding, Spawning_mode, isPlanktonic_atanypoint)



# SUPP FIGURE X - trait coverage by trait ------

#calc trait coverage by trait
total = nrow(df)
trtcov <- df %>% dplyr::select(-organism_biosamp) %>%
  summarise(across(everything(), ~ ((total - sum(is.na(.x)))/total)*100)) %>% 
  pivot_longer(., names_to = "trait", values_to = "percentage_coverage", cols = 1:ncol(.))

names <- data.frame(trait = trtcov$trait, trait.clean = c("Benthicity","Body size",
                                                          "Generational structure",
                                                          "Egg size","Pelagic larval duration",
                                                          "Philopatry","Larval provisioining",
                                                          "Spawning mode","Planktonicity"))

trtcov <- merge(trtcov, names, by = "trait", all.x = T, all.y = T)

#plot
trtcov %>%
  ggplot() + 
  geom_hline(yintercept = 50, colour = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 80, colour = "blue", linetype = "dashed", size = 0.5) +
  geom_point(aes(x = reorder(trait.clean, desc(trait.clean)), y = percentage_coverage), 
             shape = 21, colour = "black", size = 5, fill = "black") +
  labs(y = "Trait coverage across all species (%)",
       x = "Trait") +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(size=9.5),
        axis.title = element_text(size=11),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste0(figdir, "/biotic_trait_coverage-bytrait.pdf"), width = 8, height = 5, units = c("in"))



# SUPP FIGURE X - trait coverage by species ------------

#calc trait cov by species
spcov <- df %>%
  rowwise(organism_biosamp) %>% 
  summarise(spcov = (((ncol(.)-1) - (sum(is.na(cur_data()))))/(ncol(.)-1))*100) %>% 
  ungroup()

#plot
spcov %>%
  ggplot() +
  geom_hline(yintercept = 50, colour = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 80, colour = "blue", linetype = "dashed", size = 0.5) +
  geom_point(aes(x = reorder(organism_biosamp, desc(organism_biosamp)), y = spcov), 
             shape = 21, colour = "black", fill = "black", size = 2) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 3.5),
        axis.text.x = element_text(size=9.5),
        axis.title = element_text(size=11),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.3)) +
  labs(y = "Trait coverage per species (%)",
       x = "Species")

ggsave(paste0(figdir, "/biotic_trait_coverage-byspecies.pdf"), width = 8, height = 5, units = c("in"))




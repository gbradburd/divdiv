
# copied from another script - need to review/build out to be stand alone


# calculate trait and species coverage ------

#trait coverage
total = nrow(df)
trtcov <- df %>% dplyr::select(-organism_biosamp,-PLD_point,-PLD_Min, -PLD_Max) %>%
  summarise(across(everything(), ~ ((total - sum(is.na(.x)))/total)*100)) %>% 
  pivot_longer(., names_to = "trait", values_to = "percentage_coverage", cols = 1:ncol(.)) %>%
  mutate(trait_type = ifelse(trait %in% c("PLD_point2", "Fecundity_EggSize", 
                                          "Body_Size"), 
                             "green", "additional")) %>% 
  mutate(trait_type = ifelse(trait %in% c("isPlanktonic_atanypoint", "Spawning_mode", 
                                          "Larval_feeding", "ReturnToSpawningGround",
                                          "Generational_Structure", "isBenthic", "Large_Adult_Range"),
                             "yellow", trait_type)) %>% 
  mutate(trait_type = factor(trait_type, levels = c("green", "yellow", "additional")))

#species coverage
spcov.all <- df %>%
  rowwise(organism_biosamp) %>% summarise(spcov.all = (((ncol(.)-1) - (sum(is.na(cur_data()))))/(ncol(.)-1))*100) %>% ungroup()

spcov <- df %>% 
  dplyr::select(organism_biosamp, PLD_point2, isPlanktonic_atanypoint, 
                Spawning_mode, Larval_feeding, ReturnToSpawningGround, Fecundity_EggSize, Body_Size, Generational_Structure,
                isBenthic, Large_Adult_Range) %>%
  rowwise(organism_biosamp) %>% summarise(spcov.greenandyellow = (((ncol(.)-1) - (sum(is.na(cur_data()))))/(ncol(.)-1))*100) %>% ungroup()

#plot
trtcov %>% 
  arrange(desc(trait_type), desc(trait)) %>% mutate(order = 1:n()) %>% 
  mutate(trait = reorder(trait, order)) %>%  
  ggplot() + 
  geom_hline(yintercept = 50, colour = "red", alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 80, colour = "red", alpha = 0.5, linetype = "dashed") +
  geom_point(aes(x = trait, y = percentage_coverage, fill = trait_type), 
             shape = 21, colour = "black", size = 3) +
  scale_fill_manual(values = c("green","yellow","black")) +
  labs(fill = "Trait priority",
       y = "Trait coverage (%)",
       x = "Trait") +
  coord_flip() +
  theme_bw()

ggsave(paste0(figdir, "/biotic_trait_coverage-bytrait.pdf"), width = 8, height = 5, units = c("in"))



spcov %>%
  ggplot() +
  geom_hline(yintercept = 50, colour = "red", alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 80, colour = "red", alpha = 0.5, linetype = "dashed") +
  geom_point(aes(x = organism_biosamp, y = spcov.greenandyellow, fill = "turquoise2"), shape = 21, colour = "black", size = 2.5) +
  scale_fill_identity(guide = "legend", name = "Traits included",
                      breaks = c("black", "green", "yellow", "turquoise2"),
                      labels = c("all", "green", "yellow", "green + yellow")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 3.5)) +
  labs(y = "Trait coverage per species (%)",
       x = "Species")

ggsave(paste0(figdir, "/biotic_trait_coverage-byspecies.pdf"), width = 8, height = 5, units = c("in"))

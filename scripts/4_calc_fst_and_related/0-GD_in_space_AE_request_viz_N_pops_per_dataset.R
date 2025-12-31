#idea: explore number of sampling locations per dataset in response to AE comments

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()



# read in lat/long data for each dataset
alllatlongs <- read.csv("/Users/rachel/divdiv/data/final_genetic_latlongs.csv", header = TRUE) %>% 
  mutate(run_name = paste0("bioprj_",link))

alllatlongs %>% group_by(run_name) %>% summarise(n=n()) %>% nrow() #sanity check

n.locales <- alllatlongs %>% group_by(run_name, lat, lon) %>% summarise(n = n()) %>% 
  group_by(run_name) %>% summarise(n.locales = n())

nrow(alllatlongs)
alllatlongs %>% distinct(run_acc_sra) %>% nrow()

n.indivs <- alllatlongs %>% group_by(run_name, run_acc_sra) %>% summarise(n = n()) %>% 
  group_by(run_name) %>% summarise(n.indivs = n())

out <- merge(n.locales, n.indivs, by = "run_name", all.x = TRUE, all.y = TRUE)

out %>% ggplot() + geom_point(aes(x = n.locales, y = n.indivs)) + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() +
  theme_bw() +
  labs(title = "93 final datasets")
ggsave("/Users/rachel/divdiv/fst_stuff-N_indivs_x_locales.pdf", width = 8, height = 10.5, units = c("in"))

out %>% ggplot() + 
  geom_histogram(aes(x = n.locales), binwidth = 2, colour = "black", fill = "blue", alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,120,10)) +
  theme_bw() +
  labs(title = "93 final datasets")
ggsave("/Users/rachel/divdiv/fst_stuff-N_locales_histogram.pdf", width = 8, height = 10.5, units = c("in"))



# calc range of N indivs per location for each dataset ----------

df.ranges <- alllatlongs %>% group_by(run_name, lat, lon) %>% summarise(n = n()) %>% 
  ungroup() %>% group_by(run_name) %>% mutate(min.N.indivs.per.locale = min(n),
                                              max.N.indivs.per.locale = max(n),
                                              mean.N.indivs.per.locale = mean(n)) %>% 
  dplyr::select(run_name, min.N.indivs.per.locale, max.N.indivs.per.locale, mean.N.indivs.per.locale) %>% 
  distinct()

out <- merge(out, df.ranges, by = "run_name", all.x = TRUE, all.y = TRUE)

out %>% ggplot() +
  geom_segment(aes(x=min.N.indivs.per.locale, y=run_name, 
                   xend=max.N.indivs.per.locale, yend=run_name),
               linewidth = 1) +
  geom_point(aes(x=mean.N.indivs.per.locale, y=run_name), size = 2.5, colour = "black") +
  geom_point(aes(x=n.indivs, y=run_name), size = 2.5, colour = "orange", shape = 21, stroke = 1) +
  geom_vline(xintercept = 1, colour = "blue") +
  scale_x_continuous(trans='log10', breaks = c(1,5,10,20,50,100,150)) +
  labs(x = "Range and mean of N. indivs per locale; total indivs") +
  theme_bw() +
  theme(axis.text = element_text(size = 7.5))
ggsave("/Users/rachel/divdiv/fst_stuff-range_of_indivs_per_locale.pdf", width = 8, height = 10.5, units = c("in"))

df.ranges <- alllatlongs %>% group_by(run_name, lat, lon) %>% summarise(n = n()) %>% 
  ungroup() %>% group_by(run_name) %>% mutate(min.N.indivs.per.locale = min(n),
                                              max.N.indivs.per.locale = max(n),
                                              mean.N.indivs.per.locale = mean(n)) %>% 
  dplyr::select(run_name, min.N.indivs.per.locale, max.N.indivs.per.locale, mean.N.indivs.per.locale) %>% 
  distinct()

out <- merge(out, df.ranges, by = "run_name", all.x = TRUE, all.y = TRUE)

# species names
out %>% separate(., run_name, into = c("a","b","sptemp"), sep = "_", remove = F) %>% 
  dplyr::select(-a,-b) %>% mutate(species = gsub("-"," ",sptemp)) %>% 
  ggplot() +
  geom_segment(aes(x=min.N.indivs.per.locale, y=species, 
                   xend=max.N.indivs.per.locale, yend=species),
               linewidth = 1) +
  geom_point(aes(x=mean.N.indivs.per.locale, y=species), size = 2.5, colour = "black") +
  geom_point(aes(x=n.indivs, y=species), size = 2.5, colour = "orange", shape = 21, stroke = 1) +
  geom_vline(xintercept = 1, colour = "blue") +
  scale_x_continuous(trans='log10', breaks = c(1,5,10,20,50,100,150)) +
  scale_y_discrete(limits=rev) +
  labs(x = "Range and mean of N. indivs per locale; total indivs") +
  theme_bw() +
  theme(axis.text = element_text(size = 7.5))
ggsave("/Users/rachel/divdiv/fst_stuff-range_of_indivs_per_locale.pdf", width = 8, height = 10.5, units = c("in"))



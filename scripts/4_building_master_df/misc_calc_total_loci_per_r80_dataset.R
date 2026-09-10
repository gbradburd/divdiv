
#idea: get number of total loci per dataset 
#(calculated it for all stacks runs for all datasets but didn't 
#pull out into clean df for r80 V of each dataset)

#load libraries -------
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)


rm(list=ls())
gc()

setwd("/Users/rachel/divdiv")

master <- read.csv("data/master_df.csv")
df <- master %>% dplyr::select(run_name, stacksparams, bigm, littlem, n)

alllocistats.list <- list.files(path = "data/methodological/input_and_working/all_locisummary_csvs",
                                pattern = "allparams-locisummary_stats*",
                                full.names = TRUE)

alllocistats <- data.frame(run_name = NA, stacks_params = NA, n.polyloci = NA, n.nonpoly = NA, n.totalloci = NA)
for (i in 1:length(alllocistats.list)) {
  
  filename <- alllocistats.list[i]
  tmp <- read.csv(filename) %>% dplyr::select(run_name, stacks_params, n.polyloci, n.nonpoly, n.totalloci)
  alllocistats <- rbind(alllocistats, tmp)

  }
alllocistats <- alllocistats %>% filter(is.na(run_name)==FALSE)

#get just the values for the datasets we used aka r80 V

df <- df %>% mutate(link = paste0(run_name,".",stacksparams))

alllocistats <- alllocistats %>% 
  mutate(stacksparams = gsub("stacks_","",stacks_params)) %>% 
  mutate(link = paste0(run_name,".",stacksparams))

out <- merge(df, alllocistats, by = "link", all.x = TRUE)

mean(out$n.totalloci)

#86 minimum to 187,685 loci, mean = 37,463 (poly and nonpoly)



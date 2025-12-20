
library(dplyr)
library(tidyr)
library(stringr)

rm(list = ls())
gc()

setwd("/Users/rachel/divdiv/")

#get full df of FST and pi stats - to explore how GD is distributed across space
df <- read.csv("data/popgen/r80_FSTetc_stats-wide.csv", header = TRUE)

#note that there are the below character values in df:

# "not_calculated_HPCCtimedout" = 
#some datasets ran for more than 7 days 
#on the step to calculate the full pw FST matrix (between all pairs of pops)
#so we dropped this stat so job would finish in less than 7 day max time on HPCC
# applies to data columns: "mean.rowwisefst" and "sd.rowwisefst"

# "only_two_pops_total" = 
#some datasets only have 2 populations total so 
#there is only 1 pw FST value, and thus not a SD for this metric
# applies to data column: "sd.rowwisefst"

# "only_1indiv_per_every_pop" = 
#some datasets only have 1 indiv per every location (aka dummy pop)
#sampled, heirfstat returns 0s for global WC FST and pw WC FST, 
#regardless of what the genotypes between individuals are
#dropping these non-sensical values

#there should be no NAs in the dataset right now, check this, should return 0
sum(is.na(df))

#replace all of above notes with NA, so we can do math
df <- df %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "not_calculated_HPCCtimedout", "NA"))) %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "only_two_pops_total", "NA"))) %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "only_1indiv_per_every_pop", "NA")))
#convert to numbers
for(i in 3:ncol(df)) {
  df[,i] = as.numeric(df[,i])
}
#should be 18 NAs now in df
sum(is.na(df))

#calculate CV for all metrics
df <- df %>% mutate(cv.rowwisefst = sd.rowwisefst/mean.rowwisefst,
                    cv.1toallfst = sd.1toallfst/mean.1toallfst, 
                    cv.pwp = sd.pwp/mean.pwp)








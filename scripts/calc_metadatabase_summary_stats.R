#idea: calc some general summary stats about full meta database for paper

#some notes for reference:
# L plugged into WM is stored in below file type
# L = BPstats$nLoci
# BPstats$nLoci = df %>% dplyr::distinct(clocus) %>% nrow()
# where df is stacksFAfile
# so nLoci is total number of unique loci (of whatever length, about 100bp each) across dataset, after dropping low cov indivs and loci
#load("/Users/rachel/Downloads/WMfitwishart-bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_out.Robj")

#total number of unique base pairs in dataset
#BPstats$Nbp
#Nbp = df %>% dplyr::select(clocus,n.bp) %>% 
#  dplyr::distinct() %>% 
#  dplyr::summarise(n=sum(n.bp))
#Nbp = Nbp$n
#load("/Users/rachel/Downloads/bpstats.0.5.bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_BPstats.Robj")



#load libraries ---------
library(dplyr)
library(tidyr)

rm(list = ls())
gc()

setwd("/Users/rachel/divdiv")
bpstatspath = "/Users/rachel/divdiv/data/methodological/input_and_working/all_r80_bpstats"


# get master df
df <- read.csv("data/master_df.csv", header = TRUE)



# samples per dataset
#mean
mean(df$n_samples)
#range
min(df$n_samples)
max(df$n_samples)



# total raw reads
#mean
mean(df$mean_raw_read_cnt)
#range
min(df$mean_raw_read_cnt)
max(df$mean_raw_read_cnt)



# total aligned bps
bpstatslist <- list.files(path = bpstatspath, 
                          pattern="bpstats.0.5.", 
                          full.names = TRUE)

df.nbp <- data.frame(run_name = NA, Nbp = NA, nLoci = NA, stacksparams = NA)
for (i in 1:length(bpstatslist)) {
  
  bpFile <- bpstatslist[i]
  dataset = bpFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("bpstats.0.5.","",.) %>% gsub("_BPstats.Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  stacks_params = dataset %>% strsplit(., "_") %>% unlist() %>% .[5:10] %>%
    paste(., sep="", collapse="_")
  
  load(bpFile)
  
  df.i <- data.frame(run_name = run_name, Nbp = BPstats$Nbp, nLoci = BPstats$nLoci, stacksparams = stacks_params)
  
  df.nbp <- rbind(df.nbp, df.i)
  
}
df.nbp <- df.nbp %>% filter(is.na(run_name)==F)

df.nbp.kept <- df.nbp %>% filter(run_name %in% df$run_name)

#quick sanity check that stacks params match between BPstats Robjs we pulled and master df
merge(df.nbp.kept, 
      df %>% dplyr::select(run_name,stacksparams),
      by = "run_name", all.x = TRUE, all.y = TRUE) %>% 
  mutate(check = ifelse(stacksparams.x == stacksparams.y, "match", "not a match")) %>% 
  group_by(check) %>% summarise(n=n())

merge(df.nbp.kept, 
      df %>% dplyr::select(run_name,stacksparams,mean_raw_read_cnt,read_length),
      by = "run_name", all.x = TRUE, all.y = TRUE) %>% 
  mutate(check = ifelse(stacksparams.x == stacksparams.y, "match", "not a match")) %>% 
  mutate(est.nbp = nLoci*read_length) %>%
  mutate(offby = est.nbp - Nbp) %>% 
  dplyr::select(run_name, nLoci, read_length, Nbp, est.nbp, offby, everything()) %>%
  View()

#mean
mean(df.nbp.kept$Nbp)
#range
min(df.nbp.kept$Nbp)
max(df.nbp.kept$Nbp)

sum(df.nbp.kept$Nbp)
#344,759,017 total bps analyzed in divdiv


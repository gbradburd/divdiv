#idea: build df of methodological predictors

#load libraries
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)


rm(list=ls())
gc()


#get samp count and mean raw read count
sampkeys <- list.files(path = "data/all_samplenamekeys/",
                       pattern = "samplenamekey*",
                       full.names = TRUE)


df <- data.frame(run_name = NA, n_samples = NA, mean_raw_read_cnt = NA)
for (i in 1:length(sampkeys)) {
  
  filename <- sampkeys[i]
  file <- read.delim(filename)
  
  n_samples = nrow(file)
  mean_raw_read_cnt = mean(file$total_reads_written_to_final_fastq)
  run_name = strsplit(filename, split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub(".txt","",.) %>% gsub("samplenamekey-","",.)
  df.i <- data.frame(run_name = run_name, n_samples = n_samples, mean_raw_read_cnt = mean_raw_read_cnt)
  df <- rbind(df, df.i)
  
}
df <- df %>% filter(is.na(run_name)==F)


#get read length
readl <- read.csv("data/master_bookkeeping_sheet-preStacks.csv") %>% 
  dplyr::select(run_name, trimlength) %>% 
  rename("read_length" = "trimlength")

df <- merge(df, readl, by = "run_name", all.x = T)


#get mean locus depth, r80 stacks params (data/all_r80_gstacksoutlogs), and N total SNPs (data/all_locisummary_csvs)
#note N total SNPs is for r80 dataset, and is total number of SNPs in the dataset
#and each SNP is scored in at least 50% of indivs, 
#and this is indivs remaining after dropping indivs with cogenotyped bps <30% of mean for dataset
gstacklogs <- list.files(path = "data/all_r80_gstacksoutlogs",
                       pattern = "*.out",
                       full.names = TRUE)

locusd <- data.frame(run_name = NA, mean_locus_depth = NA, littlem = NA, bigm = NA, n = NA, n.totalsnps = NA)
for (i in 1:length(gstacklogs)) {
  
  filename <- gstacklogs[i]
  
  run_name = filename %>% strsplit("/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("run-gstacks_|_0.out","",.) %>% strsplit(., split = "_little") %>% 
    as.data.frame() %>% .[1,] %>% 
    stri_split_fixed(., pattern = "_", n = 2) %>% as.data.frame() %>% 
    .[nrow(.),]
  
  #get r80 stacks params
  stacks_params = filename %>% strsplit("/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("run-gstacks_|_0.out","",.) %>% strsplit(., split = "_little") %>% 
    as.data.frame() %>% .[nrow(.),] %>%
    as.data.frame() %>% rename("x"=".") %>% 
    separate(., x, into = c("a","littlem","c","bigm","n"), sep = "_") %>% 
    dplyr::select(littlem, bigm, n)
  
  #build stacks prefix of r80 params that matches formatting in loci summary .csv so we can pull out total n snps for r80 from loci summary later
  if(stacks_params$n == "nis1lessthanM"){stacks_params.n = as.numeric(stacks_params$bigm) - 1}
  if(stacks_params$n == "nequalM"){stacks_params.n = stacks_params$bigm}
  if(stacks_params$n == "nis1morethanM"){stacks_params.n = as.numeric(stacks_params$bigm) + 1}
  stacks_params_prefix = paste0("stacks_littlem_",stacks_params$littlem,"_bigm_",stacks_params$bigm,"_n",stacks_params.n,"_",stacks_params$n)
  print(stacks_params_prefix)
  
  file <- read.delim(filename, header = FALSE)
  
  #get mean locus depth
  locusd.temp <- file %>% filter(grepl("effective per-sample coverage",V1)) %>%
    as.character() %>% 
    strsplit(., split = " ") %>% as.data.frame() %>% 
    filter(grepl("mean",.[,1])) %>% 
    .[1,]
  locusd.temp <- locusd.temp %>% gsub("mean=|x|,","",.) %>% as.numeric()
  locusd.temp <- data.frame(run_name = run_name, mean_locus_depth = locusd.temp)
  locusd.temp <- cbind(locusd.temp, stacks_params)
  
  #get total number of SNPs in r80 V of dataset
  nsnps <- read.csv(paste0("data/all_locisummary_csvs/allparams-locisummary_stats-",run_name,".csv"))
  nsnp <- nsnps %>% filter(stacks_params == stacks_params_prefix)
  nsnps <- nsnp$n.totalsnps
  locusd.temp <- locusd.temp %>% mutate(n.totalsnps = nsnps)
  
  locusd <- rbind(locusd, locusd.temp)

}
locusd <- locusd %>% filter(is.na(run_name)==F)

#merge
df <- merge(df, locusd, by = "run_name", all.x = T, all.y = T)

#save
write.csv(df, "data/methodological/methodological_predictors-wide.csv", row.names = FALSE)


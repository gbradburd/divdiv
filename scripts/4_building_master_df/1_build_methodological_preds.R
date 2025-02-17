#idea: build df of methodological predictors

#load libraries -------
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)


rm(list=ls())
gc()

setwd("/Users/rachel/divdiv")
sampkeypath = "/Users/rachel/divdiv/data/all_samplenamekeys"
wmpath = "/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_gendiv_data/"

#get samp count and mean raw read count
sampkeys <- list.files(path = sampkeypath,
                       pattern = "samplenamekey*",
                       full.names = TRUE)
#get list of samples retained for popgen/gendiv calcs (so we only calc mean_raw_read_cnt for samples we used/kept in end)
wm <- list.files(path = wmpath, 
                 pattern="initPars", 
                 full.names = TRUE)

df <- data.frame(run_name = NA, n_samples = NA, mean_raw_read_cnt = NA, read_type = NA)
for (i in 1:length(wm)) {
  
  wmFile <- wm[i]
  load(wmFile)
  if (dim(initPars$parHom)[1] != dim(initPars$parHom)[2]) {"error - matrix dims are not the same"}
  n_samples = dim(initPars$parHom)[1]
  samps <- rownames(initPars$parHom)
  
  dataset = wmFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("WMfitwishart-","",.) %>% gsub("_initPars.Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  
  sampkeyFile <- paste0(sampkeypath,"/samplenamekey-",run_name,".txt")
  sampkeyFile <- read.delim(sampkeyFile)
  sampkeyFile <- sampkeyFile %>% filter(sampid_assigned_for_bioinf %in% samps)
  mean_raw_read_cnt = mean(sampkeyFile$total_reads_written_to_final_fastq)
  read_type_sra = unique(sampkeyFile$read_type_sra)
  df.i <- data.frame(run_name = run_name, n_samples = n_samples, 
                     mean_raw_read_cnt = mean_raw_read_cnt, 
                     read_type = read_type_sra)
  df <- rbind(df, df.i)
  
}
df <- df %>% filter(is.na(run_name)==F)

#check for datasets with mixed read type
df %>% dplyr::select(run_name, read_type) %>% distinct() %>% group_by(run_name) %>%  
  mutate(n=1:n()) %>% mutate(max = max(n)) %>% 
  mutate(check = ifelse(max > 1, "mixed","not mixed")) %>%
  filter(check == "mixed")
#assign these mixed read type
mixed <- df %>% dplyr::select(run_name, read_type) %>% distinct() %>% group_by(run_name) %>%  
  mutate(n=1:n()) %>% mutate(max = max(n)) %>% 
  mutate(check = ifelse(max > 1, "mixed","not mixed")) %>%
  filter(check == "mixed") %>% dplyr::select(run_name) %>% distinct()

df <- df %>% mutate(read_type = ifelse(run_name %in% mixed$run_name, "SINGLE AND PAIRED", read_type)) %>% 
  distinct()



#get read length
readl <- read.csv("data/methodological/input_and_working/master_bookkeeping_sheet-preStacks.csv") %>% 
  dplyr::select(run_name, trimlength, adapter_names_to_remove) %>% 
  rename("read_length" = "trimlength") %>% 
  rename("name_of_adapter_removed" = "adapter_names_to_remove")

df <- merge(df, readl, by = "run_name", all.x = T)


#get mean locus depth, r80 stacks params (data/all_r80_gstacksoutlogs), and N total SNPs (data/all_locisummary_csvs)
#note N total SNPs is for r80 dataset, and is total number of SNPs in the dataset
#and each SNP is scored in at least 50% of indivs, 
#and this is indivs remaining after dropping indivs with cogenotyped bps <30% of mean for dataset
gstacklogs <- list.files(path = "data/methodological/input_and_working/all_r80_gstacksoutlogs",
                       pattern = "*.out",
                       full.names = TRUE)

locusd <- data.frame(run_name = NA, mean_locus_depth = NA, littlem = NA, bigm = NA, n = NA, mean_depth_snp_filtered = NA, n.totalsnps = NA)
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
  print(paste0(stacks_params_prefix,", ",i," of ",length(gstacklogs)))
  
  #get mean locus depth (for all loci, no filtering, from Stacks gstacks module)
  gstackslog <- read.delim(filename, header = FALSE)
  locusd.temp <- gstackslog %>% filter(grepl("effective per-sample coverage",V1)) %>%
    as.character() %>% 
    strsplit(., split = " ") %>% as.data.frame() %>% 
    filter(grepl("mean",.[,1])) %>% 
    .[1,]
  locusd.temp <- locusd.temp %>% gsub("mean=|x|,","",.) %>% as.numeric()
  locusd.temp <- data.frame(run_name = run_name, mean_locus_depth = locusd.temp)
  locusd.temp <- cbind(locusd.temp, stacks_params)
  
  #get get mean locus depth for final filtered SNPs (with low cov samps and SNPs removed and for R80 params)
  depth <- read.csv(paste0("data/methodological/input_and_working/all_readdepth_csvs/allparams-readdepth_stats-",run_name,".csv"))
  depth <- depth %>% filter(stacks_params == stacks_params_prefix) %>% filter(level == "read_depth_per_snp")
  mean_depth_snp_filtered <- mean(depth$mean_readdepth)
  locusd.temp <- locusd.temp %>% mutate(mean_depth_snp_filtered = mean_depth_snp_filtered)
  
  #get total number of SNPs in r80 V of dataset
  nsnps <- read.csv(paste0("data/methodological/input_and_working/all_locisummary_csvs/allparams-locisummary_stats-",run_name,".csv"))
  nsnp <- nsnps %>% filter(stacks_params == stacks_params_prefix)
  nsnps <- nsnp$n.totalsnps
  locusd.temp <- locusd.temp %>% mutate(n.totalsnps = nsnps)
  
  locusd <- rbind(locusd, locusd.temp)

}
locusd <- locusd %>% filter(is.na(run_name)==F)

#merge
df <- merge(df, locusd, by = "run_name", all.x = T, all.y = T)
df <- df %>% na.omit()

#save
write.csv(df, "data/methodological/methodological_predictors-wide.csv", row.names = FALSE)




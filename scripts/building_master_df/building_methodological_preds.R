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
readl <- read.csv("/Volumes/mnemo2/desktop/DivDiv/divdiv_collecting_genetic_data/marine/sequence_level_data/master_bookkeeping_sheet-preStacks.csv") %>% 
  dplyr::select(run_name, trimlength) %>% 
  rename("read_length" = "trimlength")

df <- merge(df, readl, by = "run_name", all.x = T, all.y = T)


#get mean locus depth
gstacklogs <- list.files(path = "/Volumes/mnemo2/desktop/DivDiv/all_r80_gstacksoutlogs",
                       pattern = "*.out",
                       full.names = TRUE)

locusd <- data.frame(run_name = NA, mean_locus_depth = NA, littlem = NA, bigm = NA, n = NA)
for (i in 1:length(gstacklogs)) {
  
  filename <- gstacklogs[i]
  
  run_name = filename %>% strsplit("/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("run-gstacks_|_0.out","",.) %>% strsplit(., split = "_little") %>% 
    as.data.frame() %>% .[1,] %>% 
    stri_split_fixed(., pattern = "_", n = 2) %>% as.data.frame() %>% 
    .[nrow(.),]
  
  stacks_params = filename %>% strsplit("/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("run-gstacks_|_0.out","",.) %>% strsplit(., split = "_little") %>% 
    as.data.frame() %>% .[nrow(.),] %>%
    as.data.frame() %>% rename("x"=".") %>% 
    separate(., x, into = c("a","littlem","c","bigm","n"), sep = "_") %>% 
    dplyr::select(littlem, bigm, n)
  
  file <- read.delim(filename, header = FALSE)
  
  locusd.temp <- file %>% filter(grepl("effective per-sample coverage",V1)) %>%
    as.character() %>% 
    strsplit(., split = " ") %>% as.data.frame() %>% 
    filter(grepl("mean",.[,1])) %>% 
    .[1,]
  locusd.temp <- locusd.temp %>% gsub("mean=|x|,","",.) %>% as.numeric()
  locusd.temp <- data.frame(run_name = run_name, mean_locus_depth = locusd.temp)
  locusd.temp <- cbind(locusd.temp, stacks_params)
  locusd <- rbind(locusd, locusd.temp)

}
locusd <- locusd %>% filter(is.na(run_name)==F)

#merge
df <- merge(df, locusd, by = "run_name", all.x = T, all.y = T)

#save
write.csv(df, "data/methodological/methodological_predictors-wide.csv", row.names = FALSE)



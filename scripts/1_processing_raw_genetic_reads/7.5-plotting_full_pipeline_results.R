#idea: visualize the results of all bioinf steps pre-ustacks

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(scales)

rm(list = ls())
gc()

#set path to write plots out to
input_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results"
output_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results/figures"
#specify path to text list of datasets to process
listtokeep = "/Users/rachel/divdiv/scripts/master_keys/list-finalroundof77.txt"


#get a list of datasets that have made it all the way thru read cleaning/filtering pipeline to the very last step
list_of_completed_datasets <- list.files(input_file_path, pattern="summary_of_readcounts_posttrim", full.names = T) %>% as.data.frame() %>% 
  rename("x" = ".") %>% separate(., x, into = c("garbage","link.temp"), sep = "bioprj_", remove = F) %>% 
  separate(., link.temp, into = c("link"), sep = "\\.", remove = T) %>%
  mutate(run_name = paste("bioprj_",link, sep = "")) %>% 
  dplyr::select(link, run_name)
#optional - keep just a subset
keep <- read.delim(listtokeep, sep = " ", header = F) %>% mutate(keep = gsub(".txt","",V1))
list_of_completed_datasets <- list_of_completed_datasets %>% filter(run_name %in% keep$keep)


#**************************************************************
# ---- NUMBER OF READS DOWNLOADED AND POST LOWQUAL REMOVAL ----

#read data -------------
filenames <- list.files(input_file_path, pattern="summary_of_lowqual_read", full.names = T)
filenames

for( file in filenames ) {
  
  #read in file
  x <- read.csv(file)
  #get name
  name <- file %>% as.data.frame() %>% rename("x" = ".") %>% separate(., x, into = c("garbage","name"), sep = "bioprj_", remove = T)
  name <- paste("bioprj_", name$name, sep = "")
  run_name <- name %>% as.data.frame() %>% head() %>% rename("x" = ".") %>% separate(., x, into = c("name","garbage"), sep = "\\.", remove = T)
  run_name <- run_name$name
  #tack on run_name
  x <- x %>% mutate(run_name = run_name)
  #write out dataframe of results
  assign(paste0(name, sep=""), x) #give df a unique name for each vcf file
  
}

rm(list = c("name","file","x","filenames","run_name"))

#bind all dataframes together
df.postlowqual <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE))) %>% rename("sample_name" = "File_ID", 
                                                                                                "reads_downloaded" = "Total",
                                                                                                "reads_postlowqual" = "Retained_Reads") %>% 
  dplyr::select(run_name,sample_name,reads_downloaded,reads_postlowqual)
rownames(df.postlowqual) <- NULL

#wide to long
df.postlowqual <- df.postlowqual %>% gather(.,"stage","read_count",3:4)
head(df.postlowqual)

#filter to only include datasets that made it to end of pipeline
dim(df.postlowqual)
df.postlowqual <- df.postlowqual %>% filter(run_name %in% list_of_completed_datasets$run_name)
dim(df.postlowqual)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))



#**************************************************************
# ---- NUMBER OF READS AFTER REMOVING ADAPTERS ----

#read in data -------------
filenames <- list.files(input_file_path, pattern="summary_of_readlengths", full.names = T)
filenames

for( file in filenames ) {
  
  #read in file
  x <- read.csv(file)
  #get name
  name <- file %>% as.data.frame() %>% rename("x" = ".") %>% separate(., x, into = c("garbage","name"), sep = "bioprj_", remove = T)
  name <- paste("bioprj_", name$name, sep = "")
  run_name <- name %>% as.data.frame() %>% head() %>% rename("x" = ".") %>% separate(., x, into = c("name","garbage"), sep = "\\.", remove = T)
  run_name <- run_name$name
  #tack on run_name
  x <- x %>% mutate(run_name = run_name)
  #write out dataframe of results
  assign(paste0(name, sep=""), x) #give df a unique name for each vcf file
  
}

rm(list = c("name","file","x","filenames","run_name"))

#bind all dataframes together
df.postadapter <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE))) %>% rename("reads_postadapter" = "number_reads") %>% 
  dplyr::select(run_name,sample_name,reads_postadapter)
rownames(df.postadapter) <- NULL

#wide to long
df.postadapter <- df.postadapter %>% gather(.,"stage","read_count",3)
head(df.postadapter)

#filter to only include datasets that made it to end of pipeline
dim(df.postadapter)
df.postadapter <- df.postadapter %>% filter(run_name %in% list_of_completed_datasets$run_name)
dim(df.postadapter)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))



#**************************************************************
# ---- NUMBER OF READS AFTER TRIMMING TO CONSISTENT LENGTH ----

#read in data -------------
filenames <- list.files(input_file_path, pattern="summary_of_readcounts_posttrim", full.names = T)
filenames

for( file in filenames ) {
  
  print(paste0("file is ", file))
  #read in file
  x <- read.csv(file)
  
  #write out dataframe of results
  name <- x$run_name %>% unique()
  assign(paste0(name, sep=""), x) #give df a unique name for each vcf file
  
}

rm(list = c("name","file","x","filenames"))

#bind all dataframes together
df.posttrim <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE))) %>% rename("read_count" = "total_reads_written_to_final_fastq") %>% 
  mutate(stage = "reads_posttrim") %>% dplyr::select(run_name,sample_name,stage,read_count)
rownames(df.posttrim) <- NULL
head(df.posttrim)

#filter to only include datasets that made it to end of pipeline
dim(df.posttrim)
df.posttrim <- df.posttrim %>% filter(run_name %in% list_of_completed_datasets$run_name)
dim(df.posttrim)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))



#**************************************************************
# ---- BIND ALL DATA TOGETHER ----

df <- rbind(df.postadapter,df.postlowqual,df.posttrim)

#standardize sample names to just original NCBI ID/name
df <- df %>% mutate(sample_name = gsub("noadaptersremoved_","",sample_name)) %>% 
  mutate(sample_name = gsub("_1.fastq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub("_1.fq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub(".fq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub(".fastq.gz","",sample_name)) %>% 
  mutate(sample_name = gsub("readyforreadlengthcheck_","",sample_name))
#order values based on order in which they were performed
df <- df %>% mutate(stage = factor(stage, levels = c("reads_downloaded", "reads_postlowqual", "reads_postadapter", "reads_posttrim")))

#check if we are missing any data (should return no entries if all present)
df %>% group_by(run_name,stage) %>% summarise(n=n()) %>% group_by(run_name) %>% summarise(n=n()) %>% filter(n < 4)



#**************************************************************
# ---- PLOT ----

list_of_datasets <- unique(df$run_name)
list_of_datasets


for (run_name_i in list_of_datasets) {  
  
  print(paste("Building PDF of plots for dataset ", run_name_i, sep = ""))
  
  #open pdf for this dataset
  pdf(paste(output_file_path,"/full_preustacks_results_",run_name_i,".pdf",sep = ""), width = 11, height = 8.5)
  
  #get just one dataset
  df.g <- df %>% filter(run_name == run_name_i)
  
  lowbound <- 0
  highbound1 <- max(df.g$read_count, na.rm = T) + 100
  
  plot1 <- df.g %>% ggplot() + 
    geom_point(aes(x=read_count, y=sample_name, colour=stage, shape=stage), fill = NA) + 
    scale_colour_manual(values = c("orangered2","deepskyblue","magenta4","royalblue2")) +
    scale_shape_manual(values = c(22,21,24,23)) +
    scale_x_continuous(breaks = seq(from = 0, to = highbound1, by = 500000), limits = c(0,highbound1)) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5, angle = 45),
          axis.text.y = element_text(size = 4)
    ) +
    labs(title = "Read counts from NCBI download thru final cleaning - for each sample in dataset", 
         subtitle = paste("run_name = ", run_name_i, sep = ""),
         x = "Number of reads per indiv. sample",
         y = "sample_name")
  
  print(plot1)
  
  
  #make plots of number of final, cleaned, reads for each dataset
  #are there any samples with super high/low read numbers within a dataset that we want to remove before assembly?
  #histogram
  
  highbound2 <- df.g %>% filter(stage == "reads_posttrim") %>% filter(read_count == max(read_count, na.rm = T)) %>% dplyr::select(read_count) %>%
    mutate(read_count = read_count + 100) %>% .[[1]]
  
  plot2 <- df.g %>% filter(stage == "reads_posttrim") %>% ggplot(aes(x = read_count)) +
    geom_histogram(bins = 50, fill = "lightgray", colour = "black") + 
    scale_x_continuous(breaks = seq(from = 0, to = highbound2, by = 500000), limits = c(0,highbound2)) +
    geom_vline(xintercept = 100000, colour = "peachpuff") +
    geom_vline(xintercept = 50000, colour = "peachpuff") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 5, angle = 45),
          panel.grid = element_blank(), 
          legend.position = "none") +
    labs(title = "Histogram of pre-assembly read counts", 
         subtitle = paste("run_name = ", run_name_i, sep = ""),
         x = "Number of reads per indiv. sample (at end of cleaning/filtering pipeline)", 
         y = "Frequency (number of .fq sample files)")
  
  print(plot2)
  
  #dot plot
  #reorder origin by ascending count
  plot3 <- df.g %>% filter(stage == "reads_posttrim") %>% mutate(sample_name = reorder(sample_name,read_count)) %>% 
    ggplot(aes(y = sample_name, x = read_count)) +
    geom_point(shape = 21, fill = NA, colour = "black") +
    geom_vline(xintercept = 100000, colour = "peachpuff") +
    geom_vline(xintercept = 50000, colour = "peachpuff") +
    scale_x_continuous(breaks = seq(from = 0, to = highbound2, by = 500000), limits = c(0,highbound2)) +
    theme_bw() +
    #scale_x_discrete(labels = sample_name) +
    theme(axis.text.x = element_text(size = 5, angle = 45),
          axis.text.y = element_text(size = 4),
          legend.position = "none"
    ) +
    labs(title = "Dot plot of pre-assembly read counts", 
         subtitle = paste("run_name = ", run_name_i, sep = ""),
         y = "sample_name", 
         x = "Number of reads per indiv. sample (at end of cleaning/filtering pipeline)")
  
  print(plot3)

  #percent reads lost vs. absolute read count
  plot4 <- df.g %>% spread(.,stage,read_count) %>% mutate(perclost = round(((100-(reads_posttrim/reads_downloaded)*100)),2)) %>% 
    ggplot(aes(x = perclost, y = reads_posttrim)) +
    geom_hline(yintercept = 100000, colour = "peachpuff") +
    geom_hline(yintercept = 50000, colour = "peachpuff") +
    geom_point(shape = 21, fill = NA, colour = "black") +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5)
    ) +
    labs(title = "Read Loss", 
         subtitle = paste("run_name = ", run_name_i, sep = ""),
         y = "Number of reads per indiv. sample (at end of cleaning/filtering pipeline)",
         x = "Percent of reads lost (comparing initial to final read count)")
  
  print(plot4)
  
  #percent reads lost by sample ID
  #reorder origin by ascending count
  plot5 <- df.g %>% spread(.,stage,read_count) %>% mutate(perclost = round(((100-(reads_posttrim/reads_downloaded)*100)),2)) %>% 
    mutate(sample_name = reorder(sample_name,perclost)) %>% 
    ggplot(aes(y = sample_name, x = perclost)) +
    geom_point() +
    scale_x_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          legend.position = "none"
    ) +
    labs(title = "Percent read loss, arranged by sample", 
         subtitle = paste("run_name = ", run_name_i, sep = ""),
         y = "sample_name", 
         x = "Percent of reads lost (comparing initial to final read count)")
  
  print(plot5)
  
  
  #close and save pdf for this datasets
  dev.off()
  print("Closed PDF")
  
}



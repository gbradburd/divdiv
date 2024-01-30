#idea: visualize the read lengths of the datasets

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(stringr)

rm(list = ls())
gc()

#set path to write plots out to
input_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results"
output_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results/figures"
#specify path to text list of datasets to process
listtokeep = "/Users/rachel/divdiv/scripts/master_keys/list-finalroundof77.txt"

#read in read length summary stats -------------
filenames <- list.files(input_file_path, pattern="summary_of_readlengths", full.names = T)
filenames
keep <- read.delim(listtokeep, sep = " ", header = F) %>% mutate(keep = gsub(".txt","",V1))
filenames <- filenames[str_detect(filenames, str_c(keep$keep, collapse="|"))]

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
df.stats <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE)))
rownames(df.stats) <- NULL
head(df.stats)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))

#convert wide to long
df.stats.l <- df.stats %>% gather(., "type", "value", 2:6)
head(df.stats.l)  
str(df.stats.l)
rm(df.stats)

#read in read length histograms -------------
filenames <- list.files(input_file_path, pattern="histogram_of_readlengths", full.names = T)
filenames
#next line optional
filenames <- filenames[str_detect(filenames, str_c(keep$keep, collapse="|"))]

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
df.hist <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE)))
rownames(df.hist) <- NULL
head(df.hist)

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))



#plot read length results ---------------
list_of_datasets <- unique(df.stats.l$run_name)
list_of_datasets

for (run_name_i in list_of_datasets) {  
  
  print(paste("Building PDF of plots for dataset ", run_name_i, sep = ""))
  
  #open pdf for this dataset
  pdf(paste(output_file_path,"/read_length_results_",run_name_i,".pdf",sep = ""), width = 11, height = 8.5)
  
  #get just one dataset - read stats
  df.g1 <- df.stats.l %>% filter(run_name == run_name_i) %>% filter(type != "number_reads")
  
  lowbound <- min(df.g1$value, na.rm = T) - 1
  highbound <- max(df.g1$value, na.rm = T) + 1

  plot1 <- df.g1 %>% ggplot() + 
    geom_point(aes(x=value, y=sample_name, colour = type, shape = type)) + 
    scale_colour_manual(values = c("royalblue2","deepskyblue","orangered2","magenta4")) +
    scale_shape_manual(values = c(16,17,15,18)) +
    scale_x_continuous(breaks = seq(from = lowbound, to = highbound, by = 2), limits = c(lowbound,highbound)) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 4)
          #panel.grid.major.y = element_blank(),
          #panel.grid.minor.y = element_line(colour = "black"),
          #panel.grid.major.x = element_line(colour= "grey60"),
          #panel.grid.minor.x = element_line(colour= "grey60")
          ) +
    labs(title = "Read lengths - for each sample in dataset", 
        subtitle = paste("run_name = ", run_name_i, sep = ""),
        x = "read lengths",
        y = "sample_name")
  
  print(plot1)
  
  
  #get just one dataset - read hist
  df.g2 <- df.hist %>% filter(run_name == run_name_i)
  
  lowbound <- min(df.g2$read_length, na.rm = T) - 5
  highbound <- max(df.g2$read_length, na.rm = T) + 5
    
  plot2 <- df.g2 %>% ggplot() + 
      geom_col(aes(x=read_length, y=number_reads, colour = sample_name), fill = NA) + 
      #xlim(lowbound,highbound) +
      scale_x_continuous(breaks = seq(from = lowbound, to = highbound, by = 2)) +
      theme_bw() + 
      theme(axis.text = element_text(size = 5),
            legend.position = "none") +
      labs(title = "Read length histogram - for each sample in dataset", 
          subtitle = paste("run_name = ", run_name_i, sep = ""),
          x = "read lengths",
          y = "number of reads")
  
  print(plot2)
  
  
  #close and save pdf for this datasets
  dev.off()
  print("Closed PDF")
}

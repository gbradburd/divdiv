#idea: visualize the results of searching for common adapters in NCBI sequence data

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

rm(list = ls())
gc()

#set directory paths
input_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results"
output_file_path = "/Users/rachel/divdiv/data/bioinformatics/summary_file_of_results/figures"
dir.create(file.path(output_file_path))

#choose which pipeline step to look at
# preadapterremoval, adapterremoval, postadapterremoval
pipelinestep = "postadapterremoval"



# ----------------------------------------------

#read in adapter removal summary stats
filenames <- list.files(input_file_path, pattern=paste("adapter_hits-", pipelinestep, sep = ""), full.names = T)
filenames

#optional
keep <- read.delim("/Users/rachel/divdiv/scripts/master_keys/list-finalroundof77.txt", sep = " ", header = F) %>% mutate(V1 = gsub(".txt","",V1)) %>% 
  mutate(id = paste0(input_file_path,"/summary_of_adapter_hits-",pipelinestep,"-",V1,".csv",sep=""))
filenames <- filenames[filenames %in% keep$id]
#end optional

for( file in filenames ) {
  
  #read in file
  x <- read.csv(file)
  #get name
  name <- file %>% as.data.frame() %>% rename("x" = ".") %>% separate(., x, into = c("garbage","name"), sep = "bioprj_", remove = T)
  name <- paste("bioprj_", name$name, sep = "")
  #get run name from file name (in future iterations can get rid of this and use run_name in .csv itself, but code was broken at one point earlier and didn't write out a run_name in the file so working around it here)
  run_name = name %>% as.data.frame() %>% rename("run_name" = ".") %>% separate(.,run_name,into=c("run_name"),sep="\\.",remove = F)
  run_name_to_fill = run_name$run_name
  x <- x %>% mutate(run_name = run_name_to_fill)
  #write out dataframe of results
  assign(paste0(name, sep=""), x) #give df a unique name for each vcf file
  print(paste("done with ", name, sep = ""))
  
}

rm(list = c("name","file","x","run_name","run_name_to_fill"))

#bind all dataframes together
df <- do.call(rbind, mget(grep("bioprj_",names(.GlobalEnv),value=TRUE)))

#remove individual dataframes
rm(list = ls(pattern = "bioprj_"))


#do some data wrangling
df <- df %>% mutate(adapter_name = gsub(" ", "", adapter_name))
df <- df %>% mutate(perc_total_reads_trimmed = (reads_with_any_adapters/total_reads_processed)*100)
df <- df %>% mutate(perc_reads_trimmed_per_adapter = (reads_with_this_adapter/total_reads_processed)*100)

df.means <- df %>% group_by(run_name,adapter_name) %>% 
  summarise(mean_perc_reads_trimmed_per_adapter = mean(perc_reads_trimmed_per_adapter)) %>% ungroup(.) %>% 
  mutate(mean_perc_reads_trimmed_per_adapter = round(mean_perc_reads_trimmed_per_adapter, 2)) %>% 
  mutate(textcolor = ifelse(mean_perc_reads_trimmed_per_adapter > 1, "bhigh", "alow"))
dim(df)
df <- merge(df, df.means, by = c("run_name","adapter_name"), all.x = T)
dim(df)
str(df)
#reorder dataset names so they are in alphabetical order
ordered <- df %>% distinct(run_name) %>% arrange()
ordered <- ordered$run_name
df$run_name_ordered <- factor(df$run_name, levels = ordered)
# ! NOTE ! - perc_total_reads_trimmed can be greater than 100% (but no other percs should be) bc we remove 
#   Illumina common adapter first and then pass all of those reads through to rest of adapters, then
#   we extract the stat for how many reads were trimmed for common adapter and how many reads were trimmed for other 
#   adapters and if some reads were trimmed more than once (for a common adapter and another adapter in second pass),
#   this stat will report that more than 100% of the reads had an adapter (bc of the two rounds of adapter searching)

#overall which datasets seem to have adapters in them?
df %>% ggplot() + 
  #geom_hline(yintercept = 50, colour = "maroon", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 50, colour = "deeppink3", linetype = "dashed") +
  geom_boxplot(aes(y=perc_total_reads_trimmed, x=run_name_ordered)) + 
  scale_y_continuous(limits = c(0,max(df$perc_total_reads_trimmed))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Which datasets might have adapters?", 
       subtitle = "Percent of reads per sample that matched any adapter\nEach data point is an individual (e.g. sample0.fq.gz)",
       x = "run_name aka Dataset",
       y = "Reads containing an adapter seq. (%)")

ggsave(paste(output_file_path,"/adapter_removal_results_summary-",pipelinestep,"-",Sys.Date(),".pdf", sep = ""), 
       width = 11, height = 8.5, units = c("in"))



#which adapters seem to be present
list_of_datasets <- unique(df$run_name)
list_of_datasets

for (run_name_i in list_of_datasets) {  
  
  print(paste("Building PDF of plots for dataset ", run_name_i, sep = ""))
  
  #open pdf for this dataset
  pdf(paste(output_file_path,"/adapter_removal_results_",run_name_i,"_",pipelinestep,".pdf",sep = ""), width = 11, height = 8.5)
  
  #get just one dataset
  df.g <- df %>% filter(run_name == run_name_i)

  #filter out adapters not found at all and get stats about how many there were
  df.g.f <- df.g %>% filter(reads_with_this_adapter != 0) %>% filter(mean_perc_reads_trimmed_per_adapter != 0)

  n_adapters_not_found <- df.g %>% filter(reads_with_this_adapter == 0 | mean_perc_reads_trimmed_per_adapter == 0) %>% 
    group_by(adapter_name) %>% summarise(n=n()) %>% nrow(.)
  n_adapters_searched_for <- df.g %>% group_by(adapter_name) %>% summarise(n=n()) %>% nrow(.)

  # curves, fixed scales, only adapters we found
  plot1.f <- ggplot(data = df.g.f) + 
    geom_smooth(aes(x = length_trimmed, y = times_trimmed), formula = 'y ~ x', method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25, size = 0.6) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), formula = 'y ~ x', method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange", size = 0.6) + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter, x = Inf, y = Inf, colour = textcolor), hjust = 2, vjust = 2, size = 2.75) +
    scale_colour_manual(values = c("gray60","blue")) +
    theme_bw() + 
    theme(strip.text = element_text(size = 4),
          axis.text = element_text(size = 6),
          title = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.305)) +
    labs(title = paste("Which adapters are present?     Key: orange=expected; blue=observed; inset=mean % of reads with a hit",
                       "\nRun name: ",run_name_i,"    ",
                       n_adapters_not_found, " of ", n_adapters_searched_for, " adapters had \"0\" hits",
                       sep = ""),
         x = "Length of adapter seq.",
         y = "Number of hits") + 
    guides(color = FALSE) +
    facet_wrap(~ adapter_name, ncol = 7)
  print(plot1.f)
  print("Finished plot 1 of 4")
  
  # curves, free scales, only adapters we found
  plot2.f <- ggplot(data = df.g.f) + 
    geom_smooth(aes(x = length_trimmed, y = times_trimmed), formula = 'y ~ x', method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25, size = 0.6) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), formula = 'y ~ x', method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange", size = 0.6) + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter, x = Inf, y = Inf, colour = textcolor), hjust = 2, vjust = 2, size = 2.75) +
    scale_colour_manual(values = c("gray60","blue")) +
    theme_bw() + 
    theme(strip.text = element_text(size = 4),
          axis.text = element_text(size = 6),
          title = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.305)) +
    labs(title = paste("Which adapters are present?     Key: orange=expected; blue=observed; inset=mean % of reads with a hit",
                       "\nRun name: ",run_name_i,"    ",
                       n_adapters_not_found, " of ", n_adapters_searched_for, " adapters had \"0\" hits",
                       sep = ""),
         x = "Length of adapter seq.",
         y = "Number of hits") + 
    guides(color = FALSE) +
    facet_wrap(~ adapter_name, ncol = 7, scales = "free")
  print(plot2.f)
  print("Finished plot 2 of 4")
  
  # points, free scales, only adapters we found
  plot3.f <- ggplot(data = df.g.f) + 
    geom_point(aes(x = length_trimmed, y = times_trimmed), colour = "blue", fill = "blue", alpha = 0.25, size = 1) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), formula = 'y ~ x', method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange", size = 0.6) + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter, x = Inf, y = Inf, colour = textcolor), hjust = 2, vjust = 2, size = 2.75) +
    theme_bw() + 
    scale_colour_manual(values = c("gray60","blue")) +
    theme(strip.text = element_text(size = 4),
          axis.text = element_text(size = 6),
          title = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.305)) +
    labs(title = paste("Which adapters are present?     Key: orange=expected; blue=observed; inset=mean % of reads with a hit",
                       "\nRun name: ",run_name_i,"    ",
                       n_adapters_not_found, " of ", n_adapters_searched_for, " adapters had \"0\" hits",
                       sep = ""),
         x = "Length of adapter seq.",
         y = "Number of hits") + 
    guides(colour = FALSE) +
    facet_wrap(~ adapter_name, ncol = 7, scales = "free")
  print(plot3.f)
  print("Finished plot 3 of 4")
  
  # histogram behind inset values in previous graphs, free y scale, only adapters we found
  plot4.f <- ggplot(data = df.g.f %>% dplyr::select(adapter_name,sample_name,perc_reads_trimmed_per_adapter,mean_perc_reads_trimmed_per_adapter,textcolor) %>% unique() %>% as.data.frame()) + 
    geom_histogram(aes(x = perc_reads_trimmed_per_adapter), colour = NA, fill = "blue", alpha = 0.5, binwidth = 5) + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter, x = Inf, y = Inf, colour = textcolor), hjust = 2, vjust = 2, size = 2.75) +
    theme_bw() + 
    scale_colour_manual(values = c("gray60","blue")) +
    theme(strip.text = element_text(size = 4),
          axis.text = element_text(size = 6),
          title = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.305)) +
    xlim(0,100) +
    labs(title = paste("Distrib. behind inset value     Key: blue=observed; inset=mean % of reads with a hit across all samples",
                       "\nRun name: ",run_name_i,"    ",
                       n_adapters_not_found, " of ", n_adapters_searched_for, " adapters had \"0\" hits",
                       sep = ""),
         x = "Percent of reads that matched adapter, within ea. sample",
         y = "Frequency (should sum to total number of samples in dataset") + 
    guides(colour = FALSE) +
    facet_wrap(~ adapter_name, ncol = 7, scales = "free_y")
  print(plot4.f)
  print("Finished plot 4 of 4")
  
  dev.off()
  print("Closed PDF")
  
}
#  


  
  
#GRAVEYARD -----------------------------------------------------------  
  
  # points and curves, free scales, only adapters we found
  plot4.f <-  ggplot(data = df.g.f) + 
    geom_point(aes(x = length_trimmed, y = times_trimmed), colour = "black", fill = NA, alpha = 0.1, size = 1) + 
    geom_smooth(aes(x = length_trimmed, y = times_trimmed), method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25, size = 0.6) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange", size = 0.6) + 
    theme_bw() + 
    theme(strip.text = element_text(size = 4),
          axis.text = element_text(size = 7)) +
    labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,"    ",
                       n_adapters_not_found, " of ", n_adapters_searched_for, " adapters had \"0\" hits",
                       sep = ""),
         x = "Length of adapter seq.",
         y = "Number of hits") + 
    facet_wrap(~ adapter_name, ncol = 8, scales = "free")
  print(plot4.f)
  
  # Create labels
  grob <- grobTree(textGrob("Testing", x=0.45,  y=0.95, hjust=0,
                            gp=gpar(col="blue", fontsize=6)))
  # Plot
  plot3.f + annotation_custom(grob)
  
  
  # points, free scales, only adapters we found
plot4 <-  ggplot(data = df.g.f) + 
  geom_point(aes(x = length_trimmed, y = times_trimmed), colour = "blue", fill = "blue", alpha = 0.25) + 
  geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
  #geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter), x = xpos, y = ypos) +
  theme_bw() + 
  theme(strip.text = element_text(size = 4)) +
  labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,sep = ""),
       x = "Length of adapter trimmed",
       y = "Times adapter trimmed") + 
  facet_wrap(~ adapter_name, ncol = 8, scales = "free")
print(plot4)

# curves, fixed scale 
plot1 <-  ggplot(data = df.g) + 
  geom_smooth(aes(x = length_trimmed, y = times_trimmed), method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25) + 
  geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
  geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter), x = xpos, y = ypos, size = 1.75) +
  theme_bw() + 
  theme(strip.text = element_text(size = 4)) +
  labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,sep = ""),
       x = "Length of adapter trimmed",
       y = "Times adapter trimmed") + 
  facet_wrap(~ adapter_name, ncol = 10)
print(plot1)


# points, fixed scale 
plot2 <-  ggplot(data = df.g) + 
    geom_point(aes(x = length_trimmed, y = times_trimmed), colour = "blue", fill = "blue", alpha = 0.25) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter), x = xpos, y = ypos) +
    theme_bw() + 
    theme(strip.text = element_text(size = 4)) +
    labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,sep = ""),
         x = "Length of adapter trimmed",
         y = "Times adapter trimmed") + 
    facet_wrap(~ adapter_name, ncol = 10)
print(plot2)

# curves, free scales
plot3 <-  ggplot(data = df.g) + 
    geom_smooth(aes(x = length_trimmed, y = times_trimmed), method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
    theme_bw() + 
    theme(strip.text = element_text(size = 4)) +
    labs(title = paste(run_name_i)) + 
    facet_wrap(~ adapter_name, ncol = 10, scales = "free")
print(plot3)

# points, free scales
plot4 <-  ggplot(data = df.g) + 
    geom_point(aes(x = length_trimmed, y = times_trimmed), colour = "blue", fill = "blue", alpha = 0.25) + 
    geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
    geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter), x = xpos, y = ypos) +
    theme_bw() + 
    theme(strip.text = element_text(size = 4)) +
    labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,sep = ""),
         x = "Length of adapter trimmed",
         y = "Times adapter trimmed") + 
    facet_wrap(~ adapter_name, ncol = 3, scales = "free")
  

plot1 <-  ggplot(data = df.g) + 
  geom_smooth(aes(x = length_trimmed, y = times_trimmed), method = 'loess', colour = "blue", linetype = "solid", fill = "blue", alpha = 0.25) + 
  geom_smooth(aes(x = length_trimmed, y = times_expected), method = 'loess', colour = "darkorange2", linetype = "dashed", fill = "orange") + 
  geom_text(mapping = aes(label = mean_perc_reads_trimmed_per_adapter), x = xpos, y = ypos, size = 1.75) +
  theme_bw() + 
  theme(strip.text = element_text(size = 4)) +
  labs(title = paste("Which adapters are present?","\nRun name: ",run_name_i,sep = ""),
       x = "Length of adapter trimmed",
       y = "Times adapter trimmed") + 
  facet_wrap(~ adapter_name, ncol = 10)
print(plot1)

#idea: import raw/original .vcf file, read it into a matrix, 
#       get the genotype frequencies per SNP
#       get proportion of individuals that are heterozygous for each SNP
#       to be used for assessing Stacks assembly params (M/n)
#     assemble/compile all summary type stats into dataframes and save (.csv)
#     read .csv's back in and plot; choose optimal Stacks assembly parameters by assessing plots



#*****************************************************************************************
#*****************************************************************************************
#### Set up ####
#load libraries
suppressMessages(library(vcfR))             #read in vcf
library(tidyr)                              #dataframe manipulation
library(dplyr, warn.conflicts = FALSE)      #dataframe manipulation
library(ggplot2)                            #graphing
library(ggbeeswarm)                         #jittering points

options(dplyr.summarise.inform = FALSE)     #so that we don't get the summarise() grouped by message every time

#print session info and args for HPCC metadata/log files
print("STARTING TO RUN COMPARING_ASSEMBLIES_AND_IDING_r80_PARAMS.R SCRIPT")
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we just opened R and loaded all libraries")
print("session info:")
sessionInfo()
#get and print a vector of all arguments passed from "-export" in submit file
args <- commandArgs(trailingOnly = TRUE)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we are printing the contents of the vector args to show all variables passed from htcondor to R enviro")
print("arguments passed:")
cat(args, sep = "\n")

#define some variables (pulled in from bash script)
run_name = args[1] #dataset name
indir = args[3] #indir
outdir = args[4] #outdir
minPropIndivsScoredin = as.numeric(args[5]) #percent of indivs that locus must be scored in to save
workdir = args[6] #directory on execute node where work is being done

#source our functions
source(paste0(workdir,"/parsing.R"))
source(paste0(workdir,"/stats.R"))

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()




#get list of all gt files to process
gt_files = list.files(outdir, pattern=paste0("gt.",minPropIndivsScoredin,"."), full.names = TRUE)
print("gt_file_list contains:")
gt_files


#loop through all gt files and calculate summary statistics for each parameter set
print("CALCING SUMMARY STATS ACROSS PARAM COMBOS")
cat(paste("\nCALCING SUMMARY STATS ACROSS PARAM COMBOS\n", sep = " "), file = stderr())
for (gtFile in gt_files){
  
  load(gtFile)
  
  print(paste("processing gt file",gtFile, sep = " "))
  cat(paste("\nprocessing gt file", gtFile, "\n", sep = " "), file = stderr())
  #get file prefix to use in naming output
  tempfileprefix <- gtFile %>% strsplit(., split = "/") %>% as.data.frame()
  fileprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_gt.Robj","",.) %>% gsub(paste0("gt.",minPropIndivsScoredin,"."),"",.)
  stacks_params = fileprefix %>% gsub(paste0(run_name,"_"),"",.)
  rm(tempfileprefix)
  #change dataframe from wide matrix to long list
  Nindivs <- nrow(gt)
  gt.long <- gt %>% as.data.frame()
  gt.long$sampid = rownames(gt.long) 
  gt.long <- gt.long %>% dplyr::select(sampid,everything())
  gt.long <- gt.long %>% tidyfast::dt_pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.long), factor_key = T)
  
  #******** calc distibution of number of SNPs per locus ********
  #group data together at SNP level
  gt.long.SNP <- gt.long %>% group_by(SNPid) %>% na.omit %>% summarise(N=n())
  #calc how many SNPs are at each locus 
  gt.long.locus <- gt.long.SNP %>% tidyfast::dt_separate(., SNPid, into = c("locus_ID","SNP_position"), sep = ":") %>% 
    group_by(locus_ID) %>% summarise(N.snps=n())
  dim(gt.long.locus)
  gt.snpsperloc.distrib <- gt.long.locus %>% group_by(N.snps) %>% summarise(N=n()) %>% 
    mutate(stacks_params = stacks_params) %>% 
    mutate(N.total.snps = sum(N.snps*N)) %>% 
    mutate(N.polyloci = nrow(gt.long.locus))
  #write out dataframe of results
  assign(paste("locistats.df", run_name, stacks_params, sep="."), gt.snpsperloc.distrib) #give df a unique name for each vcf file
  
  #******** calculating heterozygosity per SNP ********
  #get a dataframe of genotype frequencies grouped by SNP (across all pops)
  freq.perSNP <- gt.long %>% group_by(SNPid,genotype) %>% summarise(N=n()) %>% na.omit() %>% ungroup()
  freq.perSNP <- freq.perSNP %>% tidyfast::dt_pivot_wider(., names_from = genotype, values_from = N)
  #replace NA's with 0's so math works later (NA here means SNP called but 0 of that genotype present)
  freq.perSNP[is.na(freq.perSNP)] <- 0
  #make "n.indiv" col that is addition of all genotypes category freqs (equals number of indivs with that SNP called in the pop)
  #and calc proportion of geno calls that are het per SNP
  freq.perSNP <- freq.perSNP %>% mutate(n.indiv = `0` + `1` + `2`) %>% mutate(prophet = `1`/n.indiv)
  #calculate hetero values for table
  n.polyloci <- freq.perSNP %>% tidyfast::dt_separate(.,SNPid, into = c("locus_ID","SNP_pos"), sep = ":") %>% distinct(locus_ID) %>% nrow()
  percent.allcalls.het <- round(((sum(freq.perSNP$`1`)/sum(freq.perSNP$n.indiv))*100), digits = 2) # total number of heterozygous calls / all SNP calls
  percenthet <- percent.allcalls.het #so graphing code works below
  percent.allcalls.over50perchet <- freq.perSNP %>% filter(prophet>0.5) %>% summarise(total.allcalls.over.50perchet = n()) %>% 
    mutate(percent.allcalls.over.50perchet = (total.allcalls.over.50perchet/nrow(freq.perSNP))*100) #percent of SNPs more than 50% hetero
  n.totalsnps <- nrow(freq.perSNP)
  hetero.df <- as.data.frame(cbind(stacks_params, n.polyloci, n.totalsnps, percent.allcalls.het, percent.allcalls.over50perchet))  
  #write out results
  assign(paste("hetero.df", run_name, stacks_params, sep="."), hetero.df) #give df a unique name for each vcf file
  
  #******** get depth per locus at indiv level ********
  vcf <- vcfR::read.vcfR(paste0(indir,"/",run_name,"_",stacks_params,"_populations.snps.vcf"), verbose = FALSE)
  dp <- vcfR::extract.gt(vcf, "DP", as.numeric = T)
  #convert to dataframe and transpose 
  #zero in dataframe means SNP not called/locus not present in that indiv
  dp <- as.data.frame(t(dp), stringsAsFactors = FALSE)
  #make a variable for sampleID (sampid) so IDs are not lost in data rearranging
  dp$sampid <- row.names(dp)
  #move sampid to be first col
  dp <- dp[,c(ncol(dp),1:(ncol(dp)-1))]
  #change dataframe from wide matrix to long list
  dp.long <- dp %>% tidyfast::dt_pivot_longer(., names_to = "SNPid", values_to = "read_depth", cols = 2:ncol(dp))
  readdepth.persamp <- dp.long %>% group_by(sampid) %>% summarise(mean_readdepth = mean(read_depth, na.rm = T)) %>% 
    mutate(stacks_params = stacks_params, level = "read_depth_per_samp", ID = sampid) %>% 
    dplyr::select(stacks_params, level, ID, mean_readdepth)
  readdepth.persnp <- dp.long %>% group_by(SNPid) %>% summarise(mean_readdepth = mean(read_depth, na.rm = T)) %>% 
    mutate(stacks_params = stacks_params, level = "read_depth_per_snp", ID = SNPid) %>% 
    dplyr::select(stacks_params, level, ID, mean_readdepth)
  readdepth.df <- rbind(readdepth.persamp,readdepth.persnp)
  #write out results
  assign(paste("readdepth.df", run_name, stacks_params, sep="."), readdepth.df) #give df a unique name for each vcf file
  
  #******** get number of non-polymorphic loci ********
  load(paste0(outdir,"/bpstats.",minPropIndivsScoredin,".",run_name,"_",stacks_params,"_BPstats.Robj"))
  n.nonpoly.value <- data.frame(n.nonpoly = BPstats$nLoci - hetero.df$n.polyloci) %>% mutate(stacks_params = stacks_params)
  assign(paste("nonpoly.df", run_name, stacks_params, sep="."), n.nonpoly.value) #give df a unique name for each vcf file
  
  #******** r80 stats ********
  #get file prefix to use in naming output
  popgenfile = paste0(outdir,"/popgenstats.",minPropIndivsScoredin,".",fileprefix,"_popgenstats.Robj")
  load(popgenfile)
  tempfileprefix <- popgenfile %>% strsplit(., split = "/") %>% as.data.frame()
  popgenprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_popgenstats.Robj","",.)
  filter = minPropIndivsScoredin
  params = popgenprefix %>% gsub("popgenstats.","",.) %>% gsub(paste(".",run_name,sep=""),"",.) %>% gsub(".*_stacks_","",.)
  #make df of popgen stats
  df <- cbind(rbind(popgenstats$thetaW,popgenstats$globalPi,mean(popgenstats$het)) %>% as.data.frame() %>% rename("value" = "V1"),
              data.frame(stat = c("thetaW","globalPi","mean.het")))
  #count base pairs scored in at least 80% of indivs, using bpstats file
  N_r80_bps = sum(BPstats$lociDistn[round(BPstats$Nindivs*0.8,0):length(BPstats$lociDistn)])
  #count polymorphic loci scored in at least 80% of indivs, using gt file
  N_r80_polyloci <- gt %>% as.data.frame() %>% mutate(sampid = row.names(.)) %>% dplyr::select(sampid,everything()) %>% 
    tidyfast::dt_pivot_longer(., names_to = "snpid",values_to = "genotype", cols = 2:(ncol(gt)+1)) %>% 
    na.omit() %>% group_by(snpid,genotype) %>% summarise(n_indivs=n()) %>% 
    ungroup() %>% 
    pivot_wider(., names_from = "genotype", values_from = "n_indivs") %>% 
    mutate_all(~replace(., is.na(.), 0)) %>% 
    mutate(total = `1`+`2`+`0`) %>% rowwise() %>% mutate(max = max(`2`,`0`)) %>% 
    ungroup() %>% 
    mutate(poly = total - max) %>% 
    filter(poly > 0) %>% tidyfast::dt_separate(., snpid, into = c("locus","garbage"), sep = ":", remove = F) %>%
    mutate(propindivs = total/Nindivs) %>% filter(propindivs >= 0.8) %>% 
    distinct(locus) %>% nrow()
  #add r80 stats onto df
  df <- cbind(rbind(df, 
                    data.frame(value = c(N_r80_bps,N_r80_polyloci), stat = c("N_r80_bps","N_r80_polyloci"))),
              filter,
              params,
              run_name)
  #give df a unique name
  assign(paste("r80stats.df", run_name, stacks_params, sep="."), df)

  #clean up
  rm(dp, dp.long, freq.perSNP, gt.long, gt.long.locus, gt.long.SNP, gt.snpsperloc.distrib, 
     hetero.df, n.nonpoly.value, percent.allcalls.over50perchet, 
     readdepth.df, readdepth.persamp, readdepth.persnp, vcf,
     n.polyloci, n.totalsnps, percent.allcalls.het, percenthet, N_r80_bps, N_r80_polyloci,
     popgenfile, popgenprefix, params, stacks_params, filter, fileprefix, tempfileprefix,
     popgenstats, df, gt, BPstats)
  
}




#************************* Collect all data into one clean dataframe to plot - loci summary stats *************************
print("MERGING AND SAVING ALL SUMMARY STATS")
cat(paste("\nMERGING AND SAVING ALL SUMMARY STATS\n", sep = " "), file = stderr())
#bind all dataframes together
full.data <- do.call(rbind, mget(grep("hetero.df",names(.GlobalEnv),value=TRUE)))
full.data.1 <- do.call(rbind, mget(grep("nonpoly.df",names(.GlobalEnv),value=TRUE)))
full.data <- merge(full.data, full.data.1, by = "stacks_params", all = T)

#tidy up dataframe
full.data.graph <- full.data %>% 
  separate(., stacks_params, into = paste("X", 1:7, sep = ""), sep = "_", fill = "right", remove = F) %>% 
  mutate(littlem=`X3`, M=`X5`, n=gsub("n","",`X6`), type=`X7`, 
         n.totalloci = n.polyloci + as.numeric(n.nonpoly),
         run_name = run_name) %>% 
  dplyr::select(-contains("X")) %>% dplyr::select(run_name, stacks_params, everything())

write.csv(full.data.graph, paste0(outdir, "/allparams-locisummary_stats-", run_name, ".csv"), row.names = FALSE)
rm(list=ls(pattern="hetero.df"))
rm(list=ls(pattern="nonpoly.df"))
rm(full.data, full.data.1, full.data.graph)


#************************* Collect all data into one clean dataframe to plot  - snp distrib *************************
#bind all dataframes together
full.data <- do.call(rbind, mget(grep("locistats.df",names(.GlobalEnv),value=TRUE)))

#tidy up dataframe
full.data.graph <- full.data %>% mutate(percofloci = (N/N.polyloci)*100) %>% 
  dplyr::select(stacks_params,N.snps,percofloci) %>% 
  pivot_wider(., names_from = N.snps, values_from = percofloci) %>% 
  mutate(`>5 snps` = rowSums(.[,7:ncol(.)], na.rm = TRUE)) %>% 
  pivot_longer(., names_to = "N.snps", values_to = "percofloci", cols = 2:ncol(.)) %>% 
  separate(., stacks_params, into = paste("X", 1:7, sep = ""), sep = "_", fill = "right", remove = F) %>% 
  mutate(littlem=`X3`, M=`X5`, n=gsub("n","",`X6`), type=`X7`, 
         run_name = run_name) %>% 
  dplyr::select(-contains("X")) %>% dplyr::select(run_name, stacks_params, everything())

write.csv(full.data.graph, paste0(outdir, "/allparams-snpdistrib_stats-", run_name, ".csv"), row.names = FALSE)
rm(list=ls(pattern="locistats.df"))
rm(full.data, full.data.graph)


#************************* Collect all data into one clean dataframe to plot  - read depth *************************
#bind all dataframes together
full.data <- do.call(rbind, mget(grep("readdepth.df",names(.GlobalEnv),value=TRUE)))

#tidy up dataframe
full.data.graph <- full.data %>% tidyfast::dt_separate(., stacks_params, into = paste("X", 1:7, sep = ""), sep = "_", fill = "right", remove = F) %>% 
  mutate(littlem=`X3`, M=`X5`, n=gsub("n","",`X6`), type=`X7`, 
         run_name = run_name) %>% 
  dplyr::select(-contains("X")) %>% dplyr::select(run_name, stacks_params, everything())

write.csv(full.data.graph, paste0(outdir, "/allparams-readdepth_stats-", run_name, ".csv"), row.names = FALSE)
rm(list=ls(pattern="readdepth.df"))
rm(full.data, full.data.graph)


#************************* Collect all data into one clean dataframe to plot - r80 summary stats *************************
#bind all dataframes together
full.data <- do.call(rbind, mget(grep("r80stats.df.",names(.GlobalEnv),value=TRUE)))

#tidy up dataframe
full.data.graph <- full.data %>% rename("stacks_params" = "params") %>% 
  tidyfast::dt_separate(., stacks_params, into = paste("X", 1:6, sep = ""), sep = "_", fill = "right", remove = F) %>% 
  mutate(littlem=`X2`, M=`X4`, n=gsub("n","",`X5`), type=`X6`, 
         run_name = run_name) %>% 
  dplyr::select(-contains("X")) %>% dplyr::select(run_name, stacks_params, everything())
full.data.graph <- full.data.graph %>% group_by(filter,run_name,stat) %>% 
  mutate(maxval = max(value)) %>% 
  mutate(maxr80bps = ifelse(maxval == value & stat == "N_r80_bps", "r80bp","no")) %>% 
  mutate(maxr80polyloci = ifelse(maxval == value & stat == "N_r80_polyloci", "r80polyloci","no")) %>% 
  ungroup() %>% dplyr::select(-maxval)

write.csv(full.data.graph, paste0(outdir, "/allparams-r80_stats-", run_name, ".csv"), row.names = FALSE)
rm(list=ls(pattern="r80stats.df"))
rm(full.data, full.data.graph)



# --------------- PLOTTING ---------------
print("MAKING FINAL PLOTS COMPARING PARAMS")
cat(paste("\nMAKING FINAL PLOTS COMPARING PARAMS\n", sep = " "), file = stderr())
#open PDF for param picking plots
pdf(file = paste(outdir,"/plots.ALL.",run_name,".plots.pdf",sep=""), width = 11, height = 8.5)


# Combo plots with all popgen summary stats for filtered and non-filtered and all Stacks params
popgen_files <- list.files(outdir, pattern="_popgenstats.Robj", full.names = TRUE)

for (popgenfile in popgen_files){
  
  load(popgenfile)
  #get file prefix to use in naming output
  tempfileprefix <- popgenfile %>% strsplit(., split = "/") %>% as.data.frame()
  popgenprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_popgenstats.Robj","",.)
  filter = popgenprefix %>% gsub("popgenstats.","",.) %>% gsub(paste(".",run_name,sep=""),"",.) %>% gsub("_.*","",.)
  params = popgenprefix %>% gsub("popgenstats.","",.) %>% gsub(paste(".",run_name,sep=""),"",.) %>% gsub(".*_stacks_","",.)
  #make into df
  df <- cbind(rbind(popgenstats$thetaW,popgenstats$globalPi,mean(popgenstats$het)) %>% as.data.frame() %>% rename("value" = "V1"),
              data.frame(stat = c("thetaW","globalPi","mean.het")),
              filter,
              params)
  #give df a unique name for each file
  assign(paste(popgenprefix, sep=""), df) 
  rm(popgenstats,tempfileprefix,df,filter,params,popgenprefix)
  
}
rm(popgen_files)
rm(popgenfile)

#bind all dataframes together
full.data <- do.call(rbind, mget(grep("popgenstats.",names(.GlobalEnv),value=TRUE)))
row.names(full.data) <- NULL
rm(list=ls(pattern="popgenstats"))



#************************* Plot number of polymorphic loci added for each increase in M parameter ************************* ---------
full.data.graph.raw <- read.csv(paste0(outdir, "/allparams-locisummary_stats-", run_name, ".csv"))
full.data.graph <- full.data.graph.raw %>% dplyr::select(littlem, type, M, n, n.polyloci) %>% 
  mutate(varset = paste("M",M,"n",n,sep="")) %>% arrange(M,n) %>% 
  dplyr::select(-type) 

out <- list()
for (i in 1:nrow(full.data.graph)) {
  if (i == 1) {
    varset = full.data.graph[i,]$varset
    n.polyloci.added = full.data.graph[i,]$n.polyloci
  } else {
    varset = paste(full.data.graph[i,]$varset, full.data.graph[i-1,]$varset, sep = "\\")
    n.polyloci.added = full.data.graph[i,]$n.polyloci - full.data.graph[i-1,]$n.polyloci
  }
  littlem = full.data.graph[i,]$littlem
  M = full.data.graph[i,]$M
  n = full.data.graph[i,]$n
  out[[i]] <- cbind(varset, n.polyloci.added, littlem, M, n) %>% as.data.frame()
}
full.data.graph <- do.call("rbind",out) %>% mutate(n.polyloci.added = as.numeric(n.polyloci.added), littlem = as.numeric(littlem))

r80 <- read.csv(paste0(outdir, "/allparams-r80_stats-", run_name, ".csv")) %>% 
  filter(stat == "N_r80_polyloci") %>% dplyr::select(value,M,n,maxr80polyloci)

full.data.graph <- merge(full.data.graph, r80, by = c("M","n"), all = T)

#mini function to be able to plot positive and negative values in ggplot
signed_log <- scales::trans_new("signed_log",
                                transform=function(x) sign(x)*log(abs(x)),
                                inverse=function(x) sign(x)*exp(abs(x)))

#plot
plot1 <- full.data.graph %>% ggplot(aes(x=varset, y=n.polyloci.added)) +
  geom_hline(yintercept = 1, colour = "gray60", linetype = "dashed") +
  geom_line(aes(group = littlem), colour = "black") +
  geom_point(aes(fill = value, colour = maxr80polyloci), shape = 21, size = 4, stroke = 1) +
  scale_colour_manual(values = c("transparent","red")) +
  scale_y_continuous(trans = signed_log, breaks = c(-100,-10,-1,0,1,10,100,100,1000)) +
  labs(x = "Values of Stacks parameters being compared", 
       y = "Number of additional polymorphic loci returned\n(compared to previous params)") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.title = element_text(face="bold", colour="black", size=9.5))

print(plot1)
#ggsave(paste0(outdir,"/fig_polylociadded.",run_name,".pdf") , width = 12.25, height = 8, dpi = 600, units = c("in"))
rm(full.data.graph, signed_log, out, r80)



#************************* Plot popgen summary stats  ************************* --------
#plot facetted
combo.1 <- full.data %>% ggplot(aes(group= filter)) +
  geom_point(aes(x = stat, y = value, fill = params, shape = filter), colour = "black", size=2.75, position = position_dodge(width = 0.3)) +
  #geom_beeswarm(aes(x = stat, y = value, fill = params, shape = filter), priority='density', cex=5, groupOnX=TRUE, colour = "black", size=2.5, position = position_dodge(width = 0.1)) +
  scale_shape_manual(values = c(21,23)) +
  scale_fill_manual(values = c("#d7b5d8","#dd1c77","#980043","#b2e2e2","#2ca25f","#006d2c","#f7f7f7","#969696","#252525")) +
  theme_bw() +
  labs(x = "Popgen stat - free scales",
       y = "Value",
       shape = "Min. prop. indiv.",
       title = "Popgen Stats",
       subtitle = paste(run_name)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~stat, scales = "free")
print(combo.1)

#plot with 1 scale
combo.2 <- full.data %>% ggplot() +
  geom_beeswarm(aes(x = stat, y = value, fill = params, shape = filter), priority='density', cex=1.75, groupOnX=TRUE, colour = "black", size=2.75) +
  scale_shape_manual(values = c(21,23)) +
  scale_fill_manual(values = c("#d7b5d8","#dd1c77","#980043","#b2e2e2","#2ca25f","#006d2c","#f7f7f7","#969696","#252525")) +
  theme_bw() +
  labs(x = "Popgen stat - fixed scales",
       y = "Value",
       shape = "Min. prop. indiv.",
       title = "Popgen Stats",
       subtitle = paste(run_name)) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
print(combo.2)

rm(full.data,combo.1,combo.2)


#************************* Plot loci summary stats  ************************* --------
#(total number of loci; total number of polymorphic loci, total number of SNPs, 
#percent of SNPs that are heterozygous, percent of SNPs heterozygous in more than 50% of indivs)
fill.colors <- colors()[c(374,343,142)]

full.data.graph.raw <- read.csv(paste0(outdir, "/allparams-locisummary_stats-", run_name, ".csv"))
full.data.graph <- full.data.graph.raw %>% mutate(M = as.numeric(M)) %>%
  dplyr::select(run_name, stacks_params, littlem, M, n, type, everything()) %>% 
  pivot_longer(., names_to = "variable", values_to = "value", cols = n.polyloci:n.totalloci) %>% 
  filter(!variable == "total.allcalls.over.50perchet") %>% filter(!variable == "n.nonpoly") %>% 
  mutate(variable_ordered = factor(variable, levels=c("n.totalloci","n.polyloci","n.totalsnps","percent.allcalls.het","percent.allcalls.over.50perchet")))

strip_text <- c(n.totalloci = "Total loci",
                n.polyloci = "Total\npolymorphic\nloci", 
                n.totalsnps = "Total SNPs",
                percent.allcalls.het = "Percent of all\nSNP calls that\nare heterozygous", 
                percent.allcalls.over.50perchet = "Percent of SNPs\nheterozygous\nin >50% of indivs")

plot2 <- full.data.graph %>% ggplot() + 
  geom_point(aes(x=M, y=value, fill = type), shape = 21, size = 4, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = fill.colors, labels=c("n = M", "n = M - 1", "n = M + 1")) +
  scale_x_continuous(name = "M (maximum nucleotide distance allowed between stacks)", breaks = seq(1,9,1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face="bold", color="black", size=8, angle=0),
        axis.text.x = element_text(face="bold", color="black", size=8, angle=45),
        axis.title.y = element_text(face="bold", colour="black", size=9),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size = 9),
        strip.text=element_text(size = 9, face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_flip() +
  facet_wrap(~ variable_ordered, ncol = 5, strip.position="bottom", scales = "free_x", labeller = as_labeller(strip_text))

print(plot2)
#ggsave(paste0(outdir,"/fig_loci_summary_stats.",run_name,".pdf") , width = 15, height = 8, dpi = 600, units = c("in"))
rm(full.data.graph, full.data.graph.raw)



#************************* Plot pi vs N r80 loci ************************* ---------
full.data.graph <- read.csv(paste0(outdir, "/allparams-r80_stats-", run_name, ".csv"))

p <- full.data.graph$stacks_params[full.data.graph$maxr80polyloci == "r80polyloci"]
full.data.graph$maxr80polyloci[full.data.graph$stacks_params == p] = "r80polyloci"

plot3 <- full.data.graph %>% filter(stat == "globalPi" | stat == "N_r80_polyloci") %>% 
  pivot_wider(., names_from = stat, values_from  = value) %>% 
  mutate(params = paste("M",M,"n",n,sep="")) %>% 
  mutate(params = factor(params, levels = c("M3n2","M3n3","M3n4","M6n5","M6n6","M6n7","M9n8","M9n9","M9n10"))) %>% 
  ggplot() +
  geom_point(aes(x = params, y = globalPi, fill = N_r80_polyloci, colour = maxr80polyloci), shape = 21, size = 4, stroke = 1) +
  scale_colour_manual(values = c("white","red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Stacks params",
       y = "Global pi") +
  guides(fill = guide_legend(override.aes=list(shape=21)))

print(plot3)
rm(full.data.graph, p)


#************************* Plot distribution of the number of SNPs per locus ************************* -------
full.data.graph <- read.csv(paste0(outdir, "/allparams-snpdistrib_stats-", run_name, ".csv")) %>% 
  dplyr::filter(N.snps %in% c("1","2","3","4","5",">5 snps")) %>% 
  mutate(N.snps_ordered = factor(N.snps, levels=c("1", "2", "3", "4", "5", ">5 snps"))) %>% 
  mutate(params = paste("M",M,"n",n,sep="")) %>% 
  mutate(params = factor(params, levels = c("M3n2","M3n3","M3n4","M6n5","M6n6","M6n7","M9n8","M9n9","M9n10")))

fill.colors <- c("#fff7ec", "#fee8c8", "#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")

#plot
strip_text <- c(nequalM = "n = M", 
                nis1lessthanM = "n = M - 1", 
                nis1morethanM = "n = M + 1")

plot4 <- full.data.graph %>% ggplot() + 
  geom_bar(aes(x=N.snps_ordered, y=percofloci, fill = params), stat = "identity", position = position_dodge(), colour = "black") +
  scale_fill_manual(values = fill.colors) +
  labs(fill="M", 
       x = "Number of SNPs per locus", 
       y = "Percent of total polymorphic loci with X SNPs per locus") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.title = element_text(face="bold", colour="black", size=9.5),
        legend.title = element_text(face="bold", colour="black", size=9.5),
        legend.text = element_text(colour="black", size = 9.5),
        strip.text=element_text(size = 9.5, face = "bold")) +
  facet_wrap(~ type, ncol = 3, strip.position="top", labeller = as_labeller(strip_text)) 

print(plot4)
#ggsave(paste0(outdir,"/fig_snpsperlocusdistrib1.",run_name,".pdf") , width = 12.25, height = 8, dpi = 600, units = c("in"))


plot5 <- full.data.graph %>% ggplot() + 
  geom_bar(aes(x=N.snps_ordered, y=percofloci, fill = params), stat = "identity", position = position_dodge(), colour = "black") +
  scale_fill_manual(values = fill.colors) +
  labs(fill="M", 
       x = "Number of SNPs per locus", 
       y = "Percent of total polymorphic loci with X SNPs per locus") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.title = element_text(face="bold", colour="black", size=9.5),
        legend.title = element_text(face="bold", colour="black", size=9.5))

print(plot5)
#ggsave(paste0(outdir,"/fig_snpsperlocusdistrib2.",run_name,".pdf") , width = 12.25, height = 8, dpi = 600, units = c("in"))
rm(full.data.graph)



#************************* Plot read depth ************************* -------
#plot
fill.colors <- colors()[c(374,343,142)]

full.data.graph <- read.csv(paste0(outdir, "/allparams-readdepth_stats-", run_name, ".csv")) %>% mutate(M = as.factor(M))

strip_text <- c(read_depth_per_snp = "Average read depth per SNP", 
                read_depth_per_samp = "Average read depth per individual")

plot6 <- full.data.graph %>% ggplot() + 
  geom_boxplot(aes(x=M, y=mean_readdepth, fill = type), lwd = 0.65) +
  xlab(label = "M (maximum nucleotide distance allowed between stacks)") +
  scale_fill_manual(values = fill.colors, labels=c("n = M", "n = M - 1", "n = M + 1")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.text.x = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.title.y = element_text(face="bold", colour="black", size=9.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size = 9.5),
        strip.text=element_text(size = 9.5, face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_flip() +
  facet_wrap(~ level, ncol = 2, strip.position="bottom", scales = "free_x", labeller = as_labeller(strip_text))

print(plot6)
#ggsave(paste0(outdir,"/fig_readdepth.",run_name,".pdf"), width = 7, height = 8, dpi = 600, units = c("in"))

#zoom in - filtering out high read depths
plot7 <- full.data.graph %>% ggplot() + 
  lims(y = c(0,50)) +
  geom_boxplot(aes(x=M, y=mean_readdepth, fill = type), lwd = 0.65) +
  xlab(label = "M (maximum nucleotide distance allowed between stacks)") +
  scale_fill_manual(values = fill.colors, labels=c("n = M", "n = M - 1", "n = M + 1")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.text.x = element_text(face="bold", color="black", size=9.5, angle=0),
        axis.title.y = element_text(face="bold", colour="black", size=9.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size = 9.5),
        strip.text=element_text(size = 9.5, face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_flip() +
  facet_wrap(~ level, ncol = 2, strip.position="bottom", scales = "free_x", labeller = as_labeller(strip_text))

print(plot7)
#ggsave(paste0(outdir,"/fig_readdepthZOOM.",run_name,".pdf"), width = 7, height = 8, dpi = 600, units = c("in"))

dev.off()


# finally, get R80 dataset so we can copy it to a new folder
r80 <- read.csv(paste0(outdir, "/allparams-r80_stats-", run_name, ".csv")) %>% 
  filter(stat == "N_r80_polyloci") %>% filter(maxr80polyloci == "r80polyloci")
r80 <- r80$stacks_params
write.table(r80, paste0(outdir,"/r80params.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)


print("ALL DONE MAKING COMBO PLOTS")
cat(paste("\nALL DONE MAKING COMBO PLOTS\n", sep = " "), file = stderr())

print("--------- THIS IS THE END OF THE R SCRIPTS ---------")

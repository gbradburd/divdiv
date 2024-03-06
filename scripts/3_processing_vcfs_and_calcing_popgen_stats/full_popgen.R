#R script to calc popgen stats and make exploratory pots on MSU HPCC
#this script is just to find samples with few bps and help us decide how to deal with them


#*****************************************************************************************
#*****************************************************************************************
#### Set up ####
#load libraries
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(tidyfast)
suppressMessages(library(vcfR))
library(ggbeeswarm)

options(dplyr.summarise.inform = FALSE) #so that we don't get the summarise() grouped by message every time

#print session info and args for HPCC metadata/log files
print("STARTING TO RUN FULL_POPGEN.R SCRIPT")
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
storagenode = args[2] #head directory where final files are copied to/stored
indir = args[3] #indir
outdir = args[4] #outdir
keysdir = args[7] #keydir - where samplename and lat/long live
nPCs = as.numeric(args[5]) #n principal components to save
minPropIndivsScoredin = as.numeric(args[6]) #percent of indivs that locus must be scored in to save
workdir = args[8] #directory on execute node where work is being done
outdir_final = args[9] #directory on storage node where final files are stored
manualsampstodrop = args[10] #file of additional samples to drop (e.g. ancient samples, captive samples, etc. that were missed earlier)

#source our functions 
source(paste0(workdir,"/parsing.R"))
source(paste0(workdir,"/stats.R"))

print("here1")

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#bring in samples to drop file
manualsampstodrop <- read.csv(manualsampstodrop, header = T)
#get just samples to drop for this dataset, if any
manualsampstodrop <- manualsampstodrop %>% filter(run_name == run_name)
manualsampstodrop <- manualsampstodrop$sampid_assigned_for_bioinf

print("here2")

#*****************************************************************************************
#*****************************************************************************************
#### Generate and save genotype, read depth, and n loci objects ####

print("STARTING PARSING OF POPGEN FILES")

Sys.time()

#read in samples.fa files and save N loci matrix
stacksFA_files <- list.files(indir, pattern=".samples.fa", full.names = TRUE)

for (stacksFAfile in stacksFA_files){
  
  #***********************************
  #get file prefix to use in naming output
  tempfileprefix <- stacksFAfile %>% strsplit(., split = "/") %>% as.data.frame()
  fileprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_populations.samples.fa","",.)
  rm(tempfileprefix)
  
  #***********************************
  #NO FILTERING (indivs or loci, !!! NOTE !!! = this means rows in gt might be longer than all stats/objects in BPstats if badApples dropped in BPstats)
  #get N loci matrix
  print(paste("NO FILTERING - starting to process samples.fa file",stacksFAfile, sep = " "))
  getBPstats(stacksFAfile = stacksFAfile, minPropIndivsScoredin = 0,
             checkforlowcovsamps = TRUE,
             outPath = paste(outdir,"/bpstats.","0",".",fileprefix,sep=""))
  #get genotype matrix
  vcfFile = list.files(indir, pattern=paste(fileprefix,"_populations.snps.vcf",sep = ""), full.names = TRUE)
  print(paste("NO FILTERING - starting to process vcf file",vcfFile, sep = " "))
  gt <- vcf2R(vcfFile = vcfFile, readDepth = FALSE, minPropIndivsScoredin = 0)
  save(gt,file=paste(outdir,"/gt.","0",".",fileprefix,"_gt.Robj",sep=""))
  rm(gt)
  #get mean read depths
  dp.long <- getallReadDepths(vcfFile = vcfFile)
  dp.mean <- dp.long %>% stats::na.omit() %>% dplyr::group_by(sampid) %>% dplyr::summarise(mean.dp = mean(read_depth))
  write.table(dp.mean,file=paste(outdir,"/meanreaddepth.","0",".",fileprefix,"_readDepth.txt",sep=""), row.names = FALSE, col.names = c("sampid","meanReadDepth"))
  rm(dp.long, dp.mean)
  
  #***********************************
  #FILTERING (indivs with < 30% of median cogenoed bps, then bps in < 50% of remaining indivs)
  print(paste("FILTERING - starting to process samples.fa file",stacksFAfile, sep = " "))
  #load list of samples in < 30% of median cogeno bps to remove
  lowcovsamps = list.files(outdir, pattern=paste(fileprefix,"_ALERT",sep = ""), full.names = TRUE)
  if (length(lowcovsamps)==1) {
    
    lowcovsamps <- read.delim(lowcovsamps, header = F) %>% mutate(V1 = gsub(" ","",V1))
    
    if (length(which(grepl("co-genotyped", lowcovsamps$V1)))==0) {
      print("I did not make a list of samples to filter out based on results of BPstats but I may still drop samples manually")
    }  else {
      lowcovsamps <- unique(lowcovsamps[which(grepl("^sample", lowcovsamps$V1)),])
    }
   }
  #load list of samples to manually drop (defined by text file input to pipeline)
  if (length(manualsampstodrop)>0) {
  lowcovsamps <- unique(c(manualsampstodrop, lowcovsamps))
  }
  #print out what we're going to do wrt dropping samples
  if (length(lowcovsamps)>0) {
  	cat(paste("Removing sample: ",lowcovsamps,"\n",sep = ""))
  }   else {
     print("I did not make a list of samples to filter out from BPstat results or manual list")
  }
  # ********
  #get N loci matrix (filtered)
  getBPstats(stacksFAfile = stacksFAfile, minPropIndivsScoredin = minPropIndivsScoredin, sampstodrop = lowcovsamps,
             checkforlowcovsamps = FALSE,
             outPath = paste(outdir,"/bpstats.",minPropIndivsScoredin,".",fileprefix,sep=""))
  # ********
  #get genotype matrix (minPropIndivsScoredin should be 0 here - bc using coGeno lists to filter so everything matches for sure)
  print(paste("FILTERING - starting to process vcf file",vcfFile, sep = " "))
  gt <- vcf2R(vcfFile = vcfFile, readDepth = FALSE, minPropIndivsScoredin = 0)
  #toss samples in gt that we already tossed from coGeno (bc they had few cogenotyped bps)
  bpstats <- list.files(outdir, pattern=paste("bpstats.",minPropIndivsScoredin,".",fileprefix,"_BPstats.Robj",sep = ""), full.names = TRUE)
  load(bpstats)
  gt <- gt[which(rownames(gt) %in% BPstats$sampIDskept$sample), ]
  #toss SNPs in loci in gt that we already tossed from coGeno (bc loci were scored in few indivs)
  #note - even if minPropIndivsScoredin = 0 may lose some SNPs here - aka SNPs only scored in low cov samps that were dropped in preceding step
  gt <- gt[ ,stringr::str_detect(colnames(gt), stringr::str_c(BPstats$lociIDskept$clocus, collapse = "|"))]
  save(gt,file=paste(outdir,"/gt.",minPropIndivsScoredin,".",fileprefix,"_gt.Robj",sep=""))
  # ********
  #get mean read depths
  dp.long <- getallReadDepths(vcfFile = vcfFile)
  #toss samples that we already tossed from coGeno (bc they had few cogenotyped bps)
  dp.long <- dp.long %>% filter(sampid %in% BPstats$sampIDskept$sample)
  #toss SNPs that we already tossed from coGeno (bc they were scored in few indivs)
  locusidstokeep <- BPstats$lociIDskept %>% dplyr::mutate(clocus = gsub("\\^|:","",clocus))
  dp.long <- dp.long %>% tidyfast::dt_separate(., SNPid, into = c("locusid","garbage"), sep = ":", remove = F) %>% 
    dplyr::filter(locusid %in% locusidstokeep$clocus)
  #calc mean read depth per indiv
  dp.mean <- dp.long %>% stats::na.omit() %>% dplyr::group_by(sampid) %>% dplyr::summarise(mean.dp = mean(read_depth))
  write.table(dp.mean,file=paste(outdir,"/meanreaddepth.",minPropIndivsScoredin,".",fileprefix,"_readDepth.txt",sep=""), row.names = FALSE, col.names = c("sampid","meanReadDepth"))
  
  
  rm(gt, bpstats, BPstats, vcfFile, locusidstokeep, dp.long, dp.mean, lowcovsamps)
  
  print(paste("FINISHED parsing for",fileprefix, sep = " "))
  
  rm(fileprefix)

workdirfiles = list.files(indir, pattern="meanreaddepth|gt|bpstats", full.names = TRUE)
print("workdirfiles are")
workdirfiles
print(paste("copying outputs to", outdir_final, sep = " "))
file.copy(from=workdirfiles, to=paste0(outdir_final,"/"),
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
rm(workdirfiles)
    

}
rm(stacksFA_files, stacksFAfile)

Sys.time()

#copy files to final location, in case job breaks or runs out of time, so we don't lose everything
#will still try to do a final copy in the bash script after R script completely finishes to be safe
workdirfiles = list.files(indir, pattern="meanreaddepth|gt|bpstats", full.names = TRUE)
print("workdirfiles are")
workdirfiles
print(paste("copying outputs to", outdir_final, sep = " "))
file.copy(from=workdirfiles, to=paste0(outdir_final,"/"), 
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
rm(workdirfiles)



#*****************************************************************************************
#*****************************************************************************************
#### Generate and save popgen stats ####

print("STARTING CALCING OF POPGEN STATS")

gt_files <- list.files(outdir, pattern="_gt.Robj", full.names = TRUE)

for (gtfile in gt_files){
  
  print(paste("starting to process gt file",gtfile, sep = " "))
  cat(paste("\nstarting to process gt file",gtfile,"\n", sep = " "), file = stderr())
  
  #get file prefix to use in naming output
  tempfileprefix <- gtfile %>% strsplit(., split = "/") %>% as.data.frame()
  fileprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_gt.Robj","",.) %>% gsub("gt.","",.)
  rm(tempfileprefix)
  
  #load inputs (from above steps)
  load(gtfile)
  bpstats <- list.files(outdir, pattern=paste("bpstats.",fileprefix,"_BPstats.Robj",sep = ""), full.names = TRUE)
  load(bpstats)
  
  #make sure we have same samples in all objects/inputs
  
  #calc popgen stats
  thetaW <- calcThetaW(gt=gt,lociDistn=BPstats$lociDistn)
    pwp.gt = gt[rownames(BPstats$coGeno),]  
  pwpList <- freqs2pairwisePi(freqs=pwp.gt/2,coGeno=BPstats$coGeno,quiet=TRUE)
  pwp <- pwpList$pwp
  se <- pwpList$se
  globalPi <- mean(pwp[upper.tri(pwp,diag=TRUE)])
  if (length(which(BPstats$coGeno == 0)) > 0) {
    print("skipping PCA bc at least one instance of no cogenotyped bps btwn indivs")
    pcs = NULL
  } else {
    cat("\n", file = stderr())
    pcs = NULL
    cat("I am trying to run PCA\n", file = stderr())
    try(pcs <- doPCA(gt=gt,nPCs=nPCs))
    cat("I am done trying to run PCA\n", file = stderr())
  }
  het <- calcHet(gt=pwp.gt,nLoci=diag(BPstats$coGeno))
  
  popgenstats <- list("thetaW" = thetaW,
                      "pwp" = pwp,
                      "se" = se,
                      "globalPi" = globalPi,
                      "pcs" = pcs,
                      "het" = het)
  
  save(popgenstats,file=paste0(outdir,"/popgenstats.",fileprefix,"_popgenstats.Robj"))
  
  rm(popgenstats, thetaW, pwp, pwpList, se, pcs, globalPi, gt, BPstats, het, bpstats)
  
  print(paste("finished processing",fileprefix, sep = " "))
  rm(fileprefix)
  
}
rm(gt_files, gtfile)

Sys.time()

print(paste0("ALL DONE WITH PARSING AND POPGEN STATS FOR ", run_name))

#copy files to final location, in case job breaks or runs out of time, so we don't lose everything
#will still try to do a final copy in the bash script after R script completely finishes to be safe
workdirfiles = list.files(indir, pattern="popgenstats", full.names = TRUE)
file.copy(from=workdirfiles, to=paste0(outdir_final,"/"), 
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
rm(workdirfiles)



#*****************************************************************************************
#*****************************************************************************************
#### Generate plots ####

print("STARTING TO MAKE PLOTS FOR EACH PARAM COMBO")

vcf_files <- list.files(indir, pattern="populations.snps.vcf", full.names = TRUE)

for (vcfFile in vcf_files){
  
  print(paste("starting to make plots with file",vcfFile, sep = " "))
  cat(paste("\nstarting to make plots with file",vcfFile,"\n", sep = " "), file = stderr())
  
  #get file prefix to use in naming output
  tempfileprefix <- vcfFile %>% strsplit(., split = "/") %>% as.data.frame()
  fileprefix = tempfileprefix[nrow(tempfileprefix),] %>% gsub("_populations.snps.vcf","",.)
  rm(tempfileprefix)
  
  pdf(file = paste(outdir,"/plots.",fileprefix,".plots.pdf",sep=""), width = 11, height = 8.5)
  
  #load sample name key
  sampkey <- read.delim(paste(keysdir,"/samplenamekey.txt",sep=""))
  
  #get N genoed bps (for filtered data)
  bpstats <- list.files(outdir, pattern=paste("bpstats.",minPropIndivsScoredin,".",fileprefix,"_BPstats.Robj",sep = ""), full.names = TRUE)
  load(bpstats)
  nbp <- diag(BPstats$coGeno) %>% as.data.frame() %>% rename("n.bp.genoed" = ".") %>% mutate(sampid_assigned_for_bioinf = row.names(.)) %>% 
    mutate(prop.missing.bps = (1 - (n.bp.genoed/BPstats$Nbp)))
  rm(bpstats)
  
  #load gt (from above steps) - no filtering
  gtfile = paste(outdir,"/gt.","0",".",fileprefix,"_gt.Robj", sep = "")
  load(gtfile)
  gt.long <- gt %>% as.data.frame()
  gt.long$sampid = rownames(gt.long) 
  gt.long <- gt.long %>% dplyr::select(sampid,everything())
  gt.long <- gt.long %>% tidyfast::dt_pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.long), factor_key = T)
  Nindivs.nof <- gt.long %>% distinct(sampid) %>% nrow()
  
  #load gt - filtered to loci in X prop. of indivs
  gtfile = paste(outdir,"/gt.",minPropIndivsScoredin,".",fileprefix,"_gt.Robj", sep = "")
  load(gtfile)
  gt.long.f <- gt %>% as.data.frame()
  gt.long.f$sampid = rownames(gt.long.f) 
  gt.long.f <- gt.long.f %>% dplyr::select(sampid,everything())
  gt.long.f <- gt.long.f %>% tidyfast::dt_pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.long.f), factor_key = T)
  Nindivs <- gt.long.f %>% distinct(sampid) %>% nrow()
  
  # EXPLORE SAMPLE COVERAGE -----------------------------------------------------------------
  #filter out the NAs - aka uncalled SNPs, and make historgram of how many indivs each SNP is scored in
  #N indivs each SNP is scored in
  #no filtering
  cat("\nstarting explore sample coverage plots\n", file = stderr())
  snpdistrib.plot1 <- gt.long %>% filter(is.na(genotype)==F) %>% group_by(SNPid) %>% summarise(n=n()) %>% 
    mutate(nsamps = Nindivs.nof) %>% mutate(percscoredin = ((n/nsamps)*100)) %>% 
    ggplot() +
    geom_vline(xintercept = minPropIndivsScoredin*100, colour = "red", linetype = "dashed") +
    geom_histogram(aes(x=percscoredin), binwidth = 1, colour = "black", fill = "gray40") +
    scale_x_continuous(breaks = seq(0,100,5), limits = c(0,100)) +
    theme_bw() +
    labs(x = "% of indivs SNP is scored in",
         y = "Number of SNPs",
         title = "Percent of indivs. SNPs scored in",
         subtitle = paste(fileprefix,", line = min indivs threshold", sep = ""))
  print(snpdistrib.plot1)
  
  #filtered to loci in X% of indivs
  snpdistrib.plot2 <- gt.long.f %>% filter(is.na(genotype)==F) %>% group_by(SNPid) %>% summarise(n=n()) %>% 
    mutate(nsamps = Nindivs) %>% mutate(percscoredin = ((n/nsamps)*100)) %>% 
    ggplot() +
    geom_vline(xintercept = minPropIndivsScoredin*100, colour = "red", linetype = "dashed") +
    geom_histogram(aes(x=percscoredin), binwidth = 1, colour = "black", fill = "gray40") +
    scale_x_continuous(breaks = seq(0,100,5), limits = c(0,100)) +
    theme_bw() +
    labs(x = "% of indivs SNP is scored in",
         y = "Number of SNPs",
         title = "Percent of indivs. SNPs scored in",
         subtitle = paste(fileprefix,", line = min indivs threshold, ", "filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep = ""))
  print(snpdistrib.plot2)
  
  
  # HETEROZYGOSITY -----------------------------------------------------------------
  cat("\nstarting heterozygosity plot section\n", file = stderr())
  load(paste(outdir,"/popgenstats.",minPropIndivsScoredin,".",fileprefix,"_popgenstats.Robj",sep=""))
  het.means.f <- mean(popgenstats$het)
  het.f <- popgenstats$het %>% as.data.frame() %>% rename("het" = ".")
  load(paste(outdir,"/popgenstats.","0",".",fileprefix,"_popgenstats.Robj",sep=""))
  het.means <- mean(popgenstats$het)
  het <- popgenstats$het %>% as.data.frame() %>% rename("het" = ".")
  #merge on raw read counts and N genoed bps
  het.f <- merge(het.f %>% mutate(sampid_assigned_for_bioinf = row.names(.)), sampkey, by = "sampid_assigned_for_bioinf", all.x = T)
  het.f <- merge(het.f, nbp, by = "sampid_assigned_for_bioinf", all.x = T)
  #histogram of indiv. heterozygosity
  het.plot1 <- ggplot() +
    geom_histogram(data = het, aes(x=het), alpha = 0.25, colour = "black", fill = "orange", bins = 30) + 
    geom_vline(aes(xintercept=het.means), linetype="dashed", colour = "orange") +
    geom_histogram(data = het.f, aes(x=het), alpha = 0.25, colour = "black", fill = "blue", bins = 30) + 
    geom_vline(aes(xintercept=het.means.f), linetype="dashed", colour = "blue") +
    theme_bw() +
    labs(x = "Avg. heterozygosity (per indiv)",
         y = "Frequency",
         title = "Heterozygosity",
         subtitle = paste(fileprefix, ",\nblue = filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs",
                          "\norange = no filtering",sep=""))
  print(het.plot1)
  #indiv heterozygosity vs. raw read count
  het.plot2 <- ggplot() +
    geom_smooth(data = het.f, aes(x=total_reads_written_to_final_fastq, y=het), formula = y~x, method = "lm", colour = "blue") +
    geom_point(data = het.f, aes(x=total_reads_written_to_final_fastq, y=het), colour = "blue", size = 2) +
    theme_bw() + 
    labs(y = "Avg. heterozygosity (per indiv)",
         title = "Heterozygosity vs. Missingness",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep = ""))
  print(het.plot2)
  #indiv heterozygosity vs. missingness
  het.plot3 <- ggplot() +
    geom_smooth(data = het.f, aes(x=prop.missing.bps, y=het), formula = y~x, method = "lm", colour = "blue") +
    geom_point(data = het.f, aes(x=prop.missing.bps, y=het), colour = "blue", size = 2) +
    theme_bw() + 
    labs(y = "Avg. heterozygosity (per indiv)",
         title = "Heterozygosity vs. Missingness",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep = ""))
  print(het.plot3)
  #indiv heterozygosity, no filter vs. minPropIndivsScoredin SNP filter
  het.both <- merge(het %>% mutate(sampid_assigned_for_bioinf = row.names(.)),
                    het.f, 
                    by = "sampid_assigned_for_bioinf", all = T)
  het.plot4 <- ggplot() + 
    geom_smooth(data = het.both, aes(x=het.x, y=het.y), formula = y~x, method = "lm", colour = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "gray") +
    geom_point(data = het.both, aes(x=het.x, y=het.y), colour = "black", size = 2) +
    theme_bw() + 
    labs(x = paste("Avg. heterozygosity (per indiv) - no locus filter"),
         y = paste("Avg. heterozygosity (per indiv) - filtered to loci in at least ",minPropIndivsScoredin*100,"% of indivs",sep=""),
         title = "Heterozygosity SNPs filtered vs. unfiltered",
         subtitle = paste(fileprefix,", gray line = 1:1 line",sep=""))
  print(het.plot4)
  
  rm(het, het.f, het.means, het.means.f, het.both, gt, gtfile)
  
  # SFS -----------------------------------------------------------------
  cat("\nstarting SFS plot section\n", file = stderr())
  #get a dataframe of genotype frequencies grouped by SNP (across all indvs/pops)
  freq.perSNP.f <- gt.long.f %>% group_by(SNPid,genotype) %>% summarise(N=n()) %>% ungroup() %>% na.omit()
  freq.perSNP.f <- freq.perSNP.f %>% tidyfast::dt_pivot_wider(., names_from = genotype, values_from = N) %>% as.data.frame()
  #replace NA's with 0's so math works later (NA here means SNP called but 0 of that genotype present)
  freq.perSNP.f[is.na(freq.perSNP.f)==T] <- 0
  #make column N - how many indivs is SNP scored/sequenced in?
  freq.perSNP.f <- freq.perSNP.f %>% mutate(N=`0`+`1`+`2`)
  #calculate allele counts from genotype frequencies
  freq.perSNP.f <- freq.perSNP.f %>% mutate(count_A = (`1` + `0`*2),
                                            count_a = (`1` + `2`*2),
                                            count_minor = pmin(count_A, count_a),
                                            propindivshet = `1`/N,
                                            freq_A = count_A/(2*N),
                                            freq_a = count_a/(2*N),
                                            He = 2*freq_A*freq_a,
                                            FIS = ((He-propindivshet)/He))
  
  #not filtered (to loci only in X prop. of indivs)
  # Get genotype and allele counts for each locus
  #get a dataframe of genotype frequencies grouped by SNP (across all pops)
  freq.perSNP <- gt.long %>% group_by(SNPid,genotype) %>% summarise(N=n()) %>% na.omit()
  freq.perSNP <- tidyfast::dt_pivot_wider(freq.perSNP, names_from = genotype, values_from = N) %>% as.data.frame()
  #replace NA's with 0's so math works later (NA here means SNP called but 0 of that genotype present)
  freq.perSNP[is.na(freq.perSNP)] <- 0
  #make column N - how many indivs is SNP scored/sequenced in?
  freq.perSNP <- freq.perSNP %>% mutate(N=`0`+`1`+`2`)
  #calculate genotype frequencies
  freq.perSNP <- freq.perSNP %>% mutate(freq_AA = `0`/N, freq_Aa = `1`/N, freq_aa = `2`/N)
  #calculate allele frequencies and counts and find the smaller minor allele freq for each SNP 
  freq.perSNP <- freq.perSNP %>% mutate(freq_A = ((`0`*2)+`1`)/(N*2), freq_a = ((`2`*2)+`1`)/(N*2),
                                        minor = pmin(freq_A, freq_a)) %>% 
    mutate(count_A = ((`0`*2)+`1`), count_a = ((`2`*2)+`1`), 
           count_minor = pmin(count_A, count_a), count_total = N*2,
           propindivshet = `1`/N)
  
  sfs.plot <- freq.perSNP.f %>% ggplot(aes(x = count_minor)) +
    geom_histogram(binwidth = 1, color = "black", fill = "gray40") +
    scale_x_continuous(breaks = seq(0,Nindivs,5), limits = c(0,Nindivs)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Count of minor allele (singletons,doubletons,etc.)",
         y = "Frequency",
         title = "SFS",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep="")) +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"))
  print(sfs.plot)
  
  sfs.zoom.plot <- freq.perSNP.f %>% filter(count_minor > 5) %>% ggplot(aes(x = count_minor)) +
    geom_histogram(binwidth = 1, color = "black", fill = "gray40") +
    scale_x_continuous(breaks = seq(6,Nindivs,5), limits = c(6,Nindivs)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Count of minor allele (singletons,doubletons,etc.)",
         y = "Frequency",
         title = "SFS - ZOOMED to counts > 5",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep="")) +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"))
  print(sfs.zoom.plot)
  
  # HWE -----------------------------------------------------------------
  cat("\nstarting HWE plot section\n", file = stderr())
  #filtered to loci scored in X prop. of indivs
  hwe.f.plot <- freq.perSNP.f %>% ggplot(aes(x = count_minor, y = propindivshet)) +
    geom_point(shape = 21) +
    stat_function(fun = function(x) (-0.5/(Nindivs^2))*((x-Nindivs)^2)+0.5, colour = "blue") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Count of minor allele (singletons,doubletons,etc.)",
         y = "Proportion of indivds that are hets",
         title = "Hardy Weinberg Expectation (per SNP position)",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(hwe.f.plot)
  
  #no filtering (to only SNPs scored in X prop. of indivs)
  hwe.plot <- freq.perSNP %>% ggplot(aes(x = count_minor, y = propindivshet)) +
    geom_point(shape = 21) +
    stat_function(fun = function(x) (-0.5/(Nindivs.nof^2))*((x-Nindivs.nof)^2)+0.5, colour = "orange") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Count of minor allele (singletons,doubletons,etc.)",
         y = "Proportion of indivds that are hets",
         title = "Hardy Weinberg Expectation (per SNP position)",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(hwe.plot)
  
  #histogram of FIS (per SNP) - filtered 
  fis.f.plot <- freq.perSNP.f %>% ggplot(aes(x = FIS)) +
    geom_vline(xintercept = mean(freq.perSNP.f$FIS, na.rm = T), colour = "red", linetype = "dashed") +
    geom_histogram(binwidth = 0.05) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "FIS",
         y = "Frequency",
         title = "Histogram of FIS (per SNP position)",
         subtitle = paste(fileprefix, ", filtered to loci in at least ", minPropIndivsScoredin*100, "% of indivs
                          \nred line = mean FIS; ", round(mean(freq.perSNP.f$FIS, na.rm = T),4), sep=""))
  print(fis.f.plot)
  
  rm(gt.long.f, Nindivs, Nindivs.nof, popgenstats)
  
  # IBD and PCA -------------------------------------------------------------------------------
  cat("\nstarting IBD and PCA plot section\n", file = stderr())
  #using data filtered to positions scored in at least X prop. of indivs
  # Load data
  #get lat/long table
  latlong <- read.delim(paste(keysdir,"/lat_long_table-",run_name,".txt",sep="")) %>% dplyr::select(-link) %>% distinct()
  #get read depth
  dp.mean <- read.delim(paste(outdir,"/meanreaddepth.",minPropIndivsScoredin,".",fileprefix,"_readDepth.txt",sep=""), sep = " ")
  #get coGeno
  cogeno <- BPstats$coGeno
  #get pi and pcs if present
  load(paste(outdir,"/popgenstats.",minPropIndivsScoredin,".",fileprefix,"_popgenstats.Robj",sep=""))
  #merge
  df <- merge(dp.mean, sampkey %>% rename("sampid" = "sampid_assigned_for_bioinf"), by = "sampid", all.x = T)
  df <- merge(df, latlong, by = "run_acc_sra", all.x = T)
  df <- merge(df, nbp %>% rename("sampid" = "sampid_assigned_for_bioinf"), by = "sampid", all.x = T)
  if (is.null(popgenstats$pcs) == TRUE) {
    print("popgenstats$pcs is empty")
  } else {
    print("popgenstats$pcs is not empty")
    #get pcs
    pcs <- popgenstats$pcs %>% as.data.frame() %>% mutate(sampid = row.names(.))
    df <- merge(df, pcs, by = "sampid", all.x = T)
  }
  df <- df %>% dplyr::arrange(sampid)
  #reformat
  coords <- cbind(df$long,df$lat)
  geoDist <- fields::rdist.earth(coords,miles=FALSE)
  rownames(geoDist) <- df$sampid
  colnames(geoDist) <- df$sampid
  diag(geoDist) <- NA
  geoDist <- as.vector(geoDist)
  pwp <- popgenstats$pwp %>% as.matrix()
  diag(pwp) <- NA
  pwp <- as.vector(pwp)
  diag(cogeno) <- NA
  cogeno <- as.vector(cogeno)
  
  # Plot
  #pwp vs geog dist
  IBD.plot <- ggplot() + 
    geom_point(aes(x = geoDist, y = pwp, fill = cogeno), shape = 21, colour = "black", size = 2) +
    geom_smooth(aes(x = geoDist, y = pwp), formula = y~x, method = "lm") + 
    theme_bw() +
    labs(x = "Geographic distance (km)",
         y = "Pairwise pi",
         title = "IBD - Pairwise pi vs. geog. distance",
         subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(IBD.plot)
  #plot map colored by n genotyped bps
  map.ngeno.plot <- df %>% ggplot(aes(x = long, y = lat)) +
    geom_point(aes(fill = n.bp.genoed), size = 2, shape = 21, colour = "black") +
    scale_fill_gradientn(colours = topo.colors(10)) +
    #geom_jitter(aes(colour = n.bp.genoed)) +
    theme_bw() +
    labs(x = "Longitude",
         y = "Latitude",
         title = "Map colored by N. genotyped bps",
         subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(map.ngeno.plot)
  #plot map colored by n genotyped bps - JITTERED
  map.ngeno.5.plot <- df %>% ggplot(aes(x = long, y = lat)) +
    geom_beeswarm(priority='density',cex=1.25,groupOnX=TRUE,aes(fill = n.bp.genoed),shape = 21, colour = "black",size=2) +
    scale_fill_gradientn(colours = topo.colors(10)) +
    theme_bw() +
    labs(x = "Longitude",
         y = "Latitude",
         title = "Map colored by N. genotyped bps - JITTERED",
         subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(map.ngeno.5.plot)
  #plot geogDist vs. coGeno (are samples further away from each other genotyped at fewer shared loci?)
  coGeno.forplotting <- BPstats$coGeno
  diag(coGeno.forplotting) <- NA
  coGeno.forplotting <- as.vector(coGeno.forplotting)
  geogVgeno.plot <- ggplot() + 
    geom_point(aes(x = geoDist, y = coGeno.forplotting), shape = 21, fill = NA, colour = "black", size = 2) +
    geom_smooth(aes(x = geoDist, y = coGeno.forplotting), formula = y~x, method = "lm") + 
    theme_bw() +
    labs(x = "Geographic distance (km)",
         y = "N cogenotyped bps",
         title = "Geog. distance vs. N cogenotyped bps",
         subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(geogVgeno.plot)
  #plot raw read count vs. coGeno (do samples with few coGenotyped bps just have low read counts?)
  coGenolong <- BPstats$coGeno %>% as.data.frame()
  coGenolong$sampid <- row.names(coGenolong)
  coGenolong <- coGenolong %>% dplyr::select(sampid,everything())
  coGenolong <- coGenolong %>% tidyfast::dt_pivot_longer(.,names_to="sampid2",values_to="ncogeno",cols=2:ncol(coGenolong), factor_key = T)
  coGenolong <- coGenolong %>% filter(sampid != sampid2)
  coGenolong <- coGenolong %>% group_by(sampid) %>% summarise(meancogeno = mean(ncogeno))
  coGenolong <- merge(coGenolong, sampkey %>% dplyr::select(sampid_assigned_for_bioinf,total_reads_written_to_final_fastq) %>% rename("sampid" = "sampid_assigned_for_bioinf"),
                      by = "sampid", all.x = T)
  bpVgeno.plot <- ggplot() + 
    geom_point(data = coGenolong, aes(x = total_reads_written_to_final_fastq, y = meancogeno), shape = 21, fill = NA, colour = "black", size = 2) +
    geom_smooth(data = coGenolong, aes(x = total_reads_written_to_final_fastq, y = meancogeno), formula = y~x, method = "lm") + 
    theme_bw() +
    labs(y = "mean N cogenotyped bps",
         title = "N raw reads vs. mean N cogenotyped bps",
         subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
  print(bpVgeno.plot)
  #PCA plots if present
  if (is.null(pcs) == TRUE) {
    print("popgenstats$pcs is NULL so skipping PCA plots")
  } else {
    cat("\nstarting to build PCA plots\n", file = stderr())
    #plot PCA colored by lat
    pc1 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = lat), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by Lat.",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc1)
    #plot PCA colored by long
    pc2 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = long), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by Long.",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc2)
    #plot PCA colored by raw read count
    pc3 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = total_reads_written_to_final_fastq), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by raw read count",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc3)
    #plot PCA colored by n genotyped bps
    pc4 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = n.bp.genoed), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by N. genotyped bps",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc4)
    #plot PCA colored by proportion of total genotyped bps that are missing in each indiv
    pc5 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = prop.missing.bps), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by proportion of missing genotyped bps",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc5)
    #plot PCA colored by mean read depth (of genotyped bps)
    pc6 <- df %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(aes(fill = meanReadDepth), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = "PC1",
           y = "PC2",
           title = "PC1 and 2 colored by mean read depth at SNP positions",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""))
    print(pc6)
    #plot map colored by PC1
    pc7 <- df %>% ggplot(aes(x = long, y = lat)) +
      geom_point(aes(fill = V1), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      #geom_jitter(aes(colour = V1)) +
      theme_bw() +
      labs(x = "Longitude",
           y = "Latitude",
           title = "Map colored by PC1",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""),
           fill = "PC1")
    print(pc7)
    #plot map colored by PC1 - JITTERED
    pc7.5 <- df %>% ggplot(aes(x = long, y = lat)) +
      geom_beeswarm(priority='density',cex=1.25,groupOnX=TRUE,aes(fill = V1),shape = 21, colour = "black",size=2) +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      labs(x = "Longitude",
           y = "Latitude",
           title = "Map colored by PC1 - JITTERED",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""),
           fill = "PC1")
    print(pc7.5)
    #plot map colored by PC2
    pc8 <- df %>% ggplot(aes(x = long, y = lat)) +
      geom_point(aes(fill = V2), size = 2, shape = 21, colour = "black") +
      scale_fill_gradientn(colours = topo.colors(10)) +
      #geom_jitter(aes(colour = V2)) +
      theme_bw() +
      labs(x = "Longitude",
           y = "Latitude",
           title = "Map colored by PC2",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""),
           fill = "PC2")
    print(pc8)
    #plot map colored by PC2 - JITTERED
    pc8.5 <- df %>% ggplot(aes(x = long, y = lat)) +
      geom_beeswarm(priority='density',cex=1.25,groupOnX=TRUE,aes(fill = V2),shape = 21, colour = "black",size=2) +
      scale_fill_gradientn(colours = topo.colors(10)) +
      theme_bw() +
      labs(x = "Longitude",
           y = "Latitude",
           title = "Map colored by PC2 - JITTERED",
           subtitle = paste(fileprefix, ", filtered to SNPs in at least ", minPropIndivsScoredin*100, "% of indivs", sep=""),
           colour = "PC2")
    print(pc8.5)
    rm(pcs,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc7.5,pc8,pc8.5)
  }
  rm(geoDist, coGenolong, coGeno.forplotting, ncoGeno, pwp, df,
     BPstats, dp.mean, popgenstats, coords)
  
  # READ DEPTH VS ALLELE FREQS/SINGLETONS -----------------------------------------------------------------
  cat("\nstarting to make read depth vs allele freq plots\n", file = stderr())
  # Get read depths
  vcf <- vcfR::read.vcfR(paste(indir,"/",fileprefix,"_populations.snps.vcf",sep=""), verbose = FALSE)
  #get depth per locus at indiv level 
  dp <-  vcfR::extract.gt(vcf, "DP", as.numeric = T)
  #convert to dataframe and transpose, NA in dataframe means SNP not called/locus not present in that indiv
  dp <- as.data.frame(t(dp), stringsAsFactors = FALSE)
  #make a variable for sampleID (sampid) so IDs are not lost in data rearranging and make it the first col
  dp$sampid <- row.names(dp)
  dp <- dp %>% dplyr::select(sampid, everything())
  #change dataframe from wide matrix to long list
  dp.long <- tidyfast::dt_pivot_longer(dp, names_to = "SNPid", values_to = "read_depth", cols = 2:ncol(dp))
  rm(vcf, dp)
  
  # Merge read_depth dataframe to allele_freq dataframe
  df <- merge(dp.long, freq.perSNP, by="SNPid", all.x = T) %>% mutate(link=(paste(sampid,SNPid, sep="_")))
  #also want individual genotype info so we can pick out which specific indiv has the singleton
  df <- merge(df, gt.long %>% mutate(link=(paste(sampid,SNPid, sep="_"))) %>% select(link,genotype), by="link", all.x = T)
  #filter out rows with read_depth = NA (non scored SNPs)
  df <- df %>% filter(read_depth>0)
  
  # Find singletons
  #find the individual with the singleton for each SNP that is the singleton
  df <- df %>% mutate(singleton = ifelse(count_minor == 1 & genotype == 1, "yes","no"))
  #identify each SNP position as a singleton or not
  df <- df %>% mutate(SNP_is_a_singleton = ifelse(count_a == 1, "yes","no"))
  #identify each locus as having a singleton or not
  df <- df %>% tidyfast::dt_separate(., SNPid, into = c("locusid","garbage"), sep = ":", remove = F) %>% 
    dplyr::select(-garbage) %>% 
    group_by(locusid) %>% mutate(temp = min(count_minor)) %>% ungroup() %>% 
    mutate(locus_has_a_singleton = ifelse(temp == 1, "yes","no")) %>% dplyr::select(-temp) %>% as.data.frame()
  
  # Plot
  #mean read depth per indiv (across all scored loci) vs. number of singletons in indiv
  dp.singletons1.plot <- df %>% group_by(sampid,singleton) %>% mutate(n.singletons=n()) %>% ungroup() %>% 
    group_by(sampid) %>% mutate(mean.dp = mean(read_depth)) %>% ungroup() %>% 
    dplyr::select(sampid,mean.dp,n.singletons,singleton) %>% distinct() %>% filter(singleton == "yes") %>%
    ggplot(aes(x = mean.dp, y = n.singletons)) +
    geom_point(colour = "black") +
    geom_smooth(formula = y~x, method  = "lm") +
    theme_bw() +
    labs(x = "Mean read depth for indiv.",
         y = "Number of singletons in indiv.",
         title = "N. singletons vs. read depth",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(dp.singletons1.plot)
  
  #read depth at *loci* that contain singletons vs those that dont, depth averaged across scored individuals
  means.locus <- df %>% group_by(locusid,locus_has_a_singleton) %>% summarise(mean.dp = mean(read_depth)) %>% group_by(locus_has_a_singleton) %>% 
    summarise(mean.dp = mean(mean.dp))
  #no x-axis limit
  dp.singletons2.plot <- df %>% group_by(locusid,locus_has_a_singleton) %>% summarise(mean.dp = mean(read_depth)) %>%
    ggplot() +
    geom_histogram(aes(x = mean.dp, group = locus_has_a_singleton, fill = locus_has_a_singleton), colour = "black", alpha = 0.4, binwidth = 1) +
    geom_vline(data = means.locus, aes(xintercept = mean.dp, colour = locus_has_a_singleton), linetype="dashed", size=0.8) +
    theme_bw() +
    labs(x = "Read depth (per locus)",
         y = "Frequency",
         title = "Read depth at loci with singletons vs. without",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(dp.singletons2.plot)
  #zoomed in to max depth 50 on x-axis
  dp.singletons3.plot <- df %>% group_by(locusid,locus_has_a_singleton) %>% summarise(mean.dp = mean(read_depth)) %>%
    ggplot() +
    geom_histogram(aes(x = mean.dp, group = locus_has_a_singleton, fill = locus_has_a_singleton), colour = "black", alpha = 0.4, binwidth = 1) +
    scale_x_continuous(breaks = seq(from = 0, to = 50, by = 5), limits = c(0,50)) +
    geom_vline(data = means.locus, aes(xintercept = mean.dp, colour = locus_has_a_singleton), linetype="dashed", size=0.8) +
    theme_bw() +
    labs(x = "Read depth (per locus)",
         y = "Frequency",
         title = "Read depth at loci with singletons vs. without - ZOOMED to read depths < 50",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(dp.singletons3.plot)
  
  #raw reads in indivdual vs number of singletons in indiv
  #merge onto current df
  df <- merge(df, sampkey, by.x = "sampid", by.y = "sampid_assigned_for_bioinf", all.x = T)
  dp.singletons4.plot <- df %>% group_by(sampid,singleton) %>% mutate(n.singletons=n()) %>% ungroup() %>% 
    dplyr::select(sampid,n.singletons,total_number_reads_sra,singleton) %>% distinct() %>% filter(singleton == "yes") %>%
    ggplot(aes(x = total_number_reads_sra, y = n.singletons)) +
    geom_point() +
    geom_smooth(formula = y~x, method  = "lm") +
    theme_bw() +
    labs(x = "Number raw reads into Stacks pipeline (per indiv.)",
         y = "Number of singletons (in indiv.)",
         title = "Number of singletons vs. raw reads",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(dp.singletons4.plot)
  
  #alternate alleles (aka SFS) vs mean read_depth at locus
  dp.singletons5.plot <- df %>% group_by(SNPid) %>% mutate(mean.dp = mean(read_depth)) %>% 
    dplyr::select(SNPid,mean.dp,minor) %>% distinct() %>%
    ggplot(aes(x = mean.dp, y = minor)) +
    geom_point() +
    geom_smooth(method  = "gam", formula = y ~ s(x, bs = "cs")) +
    theme_bw() +
    labs(x = "Mean read depth (per SNP position)",
         y = "Frequency of minor allele",
         title = "SFS vs. mean read depth",
         subtitle = paste(fileprefix, ", no filtering", sep=""))
  print(dp.singletons5.plot)
  
  dev.off()
  
  print(paste("finished processing",fileprefix, sep = " "))
  
  rm(dp.singletons1.plot, dp.singletons2.plot, dp.singletons3.plot, 
     dp.singletons4.plot, dp.singletons5.plot, map.ngeno.plot, map.ngeno.5.plot)
  rm(bpVgeno.plot, df, dp.long, freq.perSNP, freq.perSNP.f, geogVgeno.plot, gt.long, 
     het.plot1, het.plot2, het.plot3, hwe.f.plot, hwe.plot, IBD.plot, latlong, means.locus, 
     nbp, sampkey, sfs.plot, sfs.zoom.plot, snpdistrib.plot, 
     fileprefix, vcfFile)
  
}
rm(vcf_files)

print(paste0("ALL DONE WITH FULL_POPGEN.R FOR ", run_name))
cat(paste("\nALL DONE WITH FULL_POPGEN.R FOR", run_name,"\n", sep = " "), file = stderr())

#END


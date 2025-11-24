#idea: to respond to AE request, calc some metrics to get at how GS is distributed among populations / across space
# subset pwp matrix to just pairs of indivs in same location and calc pi per location
# calc thetaW per location
# calc pairwise FST between all pairs of locations; also calc pw FST of each population to all other indivs

# ! NOTE !
#we get list of samples used in final WM models to calc GD using final_genetic_latlongs.csv - 
#which was built from WM model outputs
#this lat/long final file MUST BE UP TO DATE for this code to work aka if the samples that go into WM code
#change, need to rebuild final_genetic_latlongs.csv to reflect that or rewrite how the code here works



# load libraries and source relevant functions
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(hierfstat, quietly = TRUE, warn.conflicts = FALSE)
#print session info and args for HPCC metadata/log files
print("STARTING TO RUN calc_metrics_about_how_GD_distributed_over_space.R SCRIPT")
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
workdir = args[3] #indir
outdir = args[4] #outdir
minPropIndivsScoredin = as.numeric(args[5]) #percent of indivs that locus must be scored in to save

#for local testing
#run_name="bioprj_PRJNA294760_Amphiprion-bicinctus"
#workdir="/Users/rachel/divdiv/for_troubleshooting"
#outdir="/Users/rachel/troubleshooting"
#minPropIndivsScoredin = 0.5
#loop.iter=1

#source our functions
#source(paste0(workdir,"/stats.R"))


#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#print all warnings/errors as they occur
options(warn=1)
Sys.time()



#get master list of all lat/long locations (for all datasets)
alllatlongs <- read.csv(paste0(workdir,"/final_genetic_latlongs.csv"), header = TRUE) %>% 
  mutate(run_name = paste0("bioprj_",link))
cat("here1\n", file = stderr())


#get input GS and related files for the dataset --------------

popgenfiles_list <- list.files(path = paste0(workdir), 
                               pattern = paste0("popgenstats.",minPropIndivsScoredin), 
                               full.names = TRUE)
print("popgenfiles_list contains:")
print(popgenfiles_list)

#POPGEN FILE LOOP
for (loop.iter in 1:length(popgenfiles_list)) {
  
  #get bookkeeping variables
  popgenfile = popgenfiles_list[loop.iter]   
  stacksparams = popgenfile %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("_popgenstats|Robj","",.) %>% gsub(run_name,"",.) %>% gsub(minPropIndivsScoredin,"",.) %>% 
    gsub("\\.","",.) %>% gsub("popgenstats_stacks_","",.)
  print(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams))
  #also print this to .err so we can tell when errors are thrown if they are
  cat(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams, " now\n"), file = stderr())
  
  #filter lat/longs down to just this dataset
  i.latlongs <- alllatlongs %>% filter(run_name %in% !!run_name)
  
  #get pw pi aka genetic matrix
  load(popgenfile)
  pwp <- popgenstats$pwp
  
  #get vector of number of bp genotyped for each number of individuals from 1:N (lociDistn)
  #(after dropping low cov indivs and loci)
  #load(paste0(workdir, "/bpstats.", minPropIndivsScoredin, ".", run_name, "_stacks_", stacksparams, "_BPstats.Robj"))
  #lociDistn = BPstats$lociDistn
  
  #get genotype x indiv matrix (gt)
  load(paste0(workdir, "/gt.", minPropIndivsScoredin, ".", run_name, "_stacks_", stacksparams, "_gt.Robj"))
  
  #get sample names and filter down to just samples used in end
  sampkey <- read.delim(paste0(workdir,"/samplenamekey.txt")) %>% 
    filter(run_acc_sra %in% i.latlongs$run_acc_sra)
  cat("here2\n", file = stderr())
  
  #assign dummy pop IDs by unique lat/long combos ---------------
  mysamps <- i.latlongs %>% mutate(locale = paste0(lat,"_",lon)) %>% 
    group_by(run_name, lat, lon) %>% 
    mutate(dummy_pop_fst = paste0("dummypop_", cur_group_id())) %>% 
    ungroup() %>% 
    dplyr::select(run_acc_sra, link, dummy_pop_fst, lat, lon)
  #attach on bioinf sample0 sample names
  mysamps <- merge(mysamps,
                   sampkey %>% dplyr::select(run_acc_sra, sampid_assigned_for_bioinf),
                   by = "run_acc_sra", all.x = TRUE)
  cat("here3\n", file = stderr())
  
  #build gt matrix formatted for FST calcs --------------
  # recode alleles from 0,1,2  to 11,12,22 format
  gt.fst <- gt
  gt.fst[gt.fst == 0] <- 11
  gt.fst[gt.fst == 1] <- 12
  gt.fst[gt.fst == 2] <- 22
  # make sure cols are all integers
  for(i in 1:ncol(gt.fst)) {
    gt.fst[,i] = as.integer(gt.fst[,i])
  }
  #set sample IDs as row names, first column as pop IDs, no other non genotype cols present
  #(and rest of columns as numerical genotypes)
  gt.fst <- gt.fst %>% as.data.frame() %>% mutate(sampid_assigned_for_bioinf = row.names(.))
  gt.fst <- merge(mysamps %>% dplyr::select(sampid_assigned_for_bioinf,dummy_pop_fst), 
                  gt.fst, 
                  by = "sampid_assigned_for_bioinf", all.x = TRUE)
  row.names(gt.fst) <- gt.fst$sampid_assigned_for_bioinf
  gt.fst <- gt.fst %>% dplyr::select(-sampid_assigned_for_bioinf)
  cat("here4\n", file = stderr())
  
  
  
  #LOCATION LOOP
  #calc GD summary stats per location --------------
  unique.pop.list <- unique(mysamps$dummy_pop_fst)
  
  out <- data.frame("dummy_pop_fst"=NA, "n.indiv.i"=NA, "mean.pwp.pop.i"=NA, "fst.wc.1toall"=NA)
  for (loop.iter.pop in 1:length(unique.pop.list)) {

    pop.i = unique.pop.list[loop.iter.pop]
    
    #get list of indivs in this pop
    samps.in.pop.i <- mysamps %>% filter(dummy_pop_fst == pop.i)
    #get how many indivs are at location
    n.indiv.i <- nrow(samps.in.pop.i)
    
    #calc mean pi ---------
    #subset pwp matrix to just samps in this dummy population aka location
    pwp.i <- pwp[rownames(pwp) %in% samps.in.pop.i$sampid_assigned_for_bioinf, 
                 rownames(pwp) %in% samps.in.pop.i$sampid_assigned_for_bioinf]
    #calc mean pi for this location
    mean.pwp.pop.i = mean(pwp[upper.tri(pwp.i,diag=TRUE)])
    cat("here5\n", file = stderr())
    
    #calc FST of each population to all rest of samples ------------
    #recode gt and calc FST so focal population i is one pop and all other pops another pop
    gt.fst.i <- gt.fst %>% mutate(dummy_pop_1toall = 
                                    ifelse(rownames(.) %in% samps.in.pop.i$sampid_assigned_for_bioinf,
                                           pop.i, "dummy_all")) %>%
      dplyr::select(-dummy_pop_fst) %>% 
      dplyr::select(dummy_pop_1toall, everything())
    #calc pairwise Fst - SLOW STEP !
    fst.wc.1toall  <- as.data.frame(hierfstat::pairwise.WCfst(gt.fst.i, diploid=TRUE))
    fst.wc.1toall <- fst.wc.1toall$dummy_all[2]
    cat("here6\n", file = stderr())
    
    #save all the per location results -------------
    out.i <- data.frame(dummy_pop_fst = pop.i, n.indiv.i = n.indiv.i, 
                        mean.pwp.pop.i = mean.pwp.pop.i,
                        fst.wc.1toall = fst.wc.1toall)
    out <- rbind(out, out.i)
    rm(pop.i, samps.in.pop.i, n.indiv.i, pwp.i,
       mean.pwp.pop.i, gt.fst.i, fst.wc.1toall, out.i)
    cat("here7\n", file = stderr())
  }
  out <- out %>% filter(is.na(dummy_pop_fst)==FALSE)
  cat("here8\n", file = stderr())
  
  
  #calc full pairwise Fst pop x pop matrix -----------
  #SLOW STEP !
  fst.pw.wc  <- as.data.frame(hierfstat::pairwise.WCfst(gt.fst, diploid=TRUE))
  cat("here9\n", file = stderr())
  
  #calculate global WC FST -------------
  fst.wc <- hierfstat::wc(gt.fst, diploid = TRUE)
  cat("here10\n", file = stderr())
  
  out.mysamps <- mysamps %>% dplyr::select(-lat, -lon) %>% 
    dplyr::select(link, run_acc_sra, sampid_assigned_for_bioinf, 
                  dummy_pop_fst)
  
  #save everthing ---------
  fststats <- list("fst_samp_info" = out.mysamps,
                   "fst_summaries" = out,
                   "fst.pw.wc" = fst.pw.wc,
                   "fst.wc" = fst.wc)
  
  save(fststats, file=paste0(outdir,"/fststats.",run_name,"_",stacksparams,"_fststats.Robj"))
  
  print(paste("finished processing",run_name, stacksparams, sep = " "))
  rm(popgenfile, stacksparams, i.latlongs, popgenstats, pwp, gt, 
     sampkey, mysamps, gt.fst, unique.pop.list, out, 
     fst.pw.wc, fst.wc, out.mysamps, fststats)
  
}

Sys.time()

print(paste0("ALL DONE WITH CALCING FST STATS FOR ", run_name))
cat(paste("\nALL DONE WITH CALCING FST STATS FOR", run_name,"\n", sep = " "), file = stderr())

#END


  
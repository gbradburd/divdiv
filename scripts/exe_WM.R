################################################################
################################################################
#	running inference of wright-malecot model on SLiMulated data
################################################################
################################################################


# load libraries and source relevant functions
library(rstan, warn.conflicts = FALSE, quietly = TRUE)
library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
#print session info and args for HPCC metadata/log files
print("STARTING TO RUN exe_WM.R SCRIPT")
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

#for local testing
#indir="/path_to_input_files_and_scripts"
#outdir="/path_for_outputs"
#run_name="bioprj_PRJNA294760_Amphiprion-bicinctus"
#minPropIndivsScoredin = 0.5

#source our functions/load models
source(paste0(indir,"/wm_lib.R"))

stanFile <- "wm_hom_cmpPar_mod_block_scaled.R"
source(paste0(indir,"/",stanFile))
ibsMod <- stan_model(model_code=stanBlock)

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#print all warnings/errors as they occur
options(warn=1)


#get pwp and distance matrices --------------

popgenfiles_list <- list.files(path = paste0(indir), pattern = "popgenstats", full.names = TRUE)

for (loop.iter in 1:length(popgenfiles_list)) {
  
  #get sample names
  sampkey <- read.delim(paste0(indir,"/samplenamekey.txt"))
  
  #get pw distance matrix
  load(paste0(indir, "/max_and_pw_dists", run_name, ".Robj", sep = ""))
  #geoDist <- max_and_pw_dists$pw.gcd.genetic
  geoDist <- max_and_pw_dists$pw.seadist.genetic
  
  popgenfile = popgenfiles_list[loop.iter]   
  stacksparams = popgenfile %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("_popgenstats|Robj","",.) %>% gsub(run_name,"",.) %>% gsub(minPropIndivsScoredin,"",.) %>% 
    gsub("\\.","",.) %>% gsub("popgenstats_stacks_","",.)
  
  print(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams))
  #also print this to .err so we can tell when errors are thrown if they are
  cat(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams, " now\n"), file = stderr())
  
  load(popgenfile)
  pwp <- popgenstats$pwp
  hom <- 1-pwp
  # add inbreeding
  diag(hom) <- 1
  
  #get number of polymorphic loci
  load(paste0(indir, "/bpstats.", minPropIndivsScoredin, ".", run_name, "_stacks_", stacksparams, "_BPstats.Robj"))
  Npolyloci = BPstats$nLoci
  
  # rename geoDist, keep just samps in pwp, and order to match pwp
  #keep just distances for samples in pwp matrix
  geoDistnames <- colnames(geoDist) %>% as.data.frame() %>% rename("run_acc_sra" = ".") %>% 
    mutate(order = 1:n()) %>% 
    merge(., sampkey %>% dplyr::select(run_acc_sra,sampid_assigned_for_bioinf), by = "run_acc_sra", all.x = T) %>% 
    arrange(order)
  geoDistnames <- geoDistnames$sampid_assigned_for_bioinf
  colnames(geoDist) <- geoDistnames
  rownames(geoDist) <- geoDistnames
  geoDist <- geoDist[rownames(geoDist) %in% rownames(pwp), colnames(geoDist) %in% rownames(pwp)]
  #make sure all matrices are ordered the same
  order <- rownames(pwp)
  #geoDist <- graph4lg::reorder_mat(mat = geoDist, order = order)
  geoDist.ordered <- geoDist[order, order]
  geoDist <- geoDist.ordered
  if (identical(row.names(pwp),row.names(geoDist)) == TRUE) {
    print("names and order of rownames for pwp and geoDist match")
  } else {
    print("ERROR ! - names and/or order of rownames for pwp and geoDist do not match!")
  }
  
  # make dataBlock for stan
  # N = number of samples
  # L = number of loci
  # hom = pairwise homozygosity
  # k = geographic distance w/in which W-M breaks down, (should be ~2*sigma)
  # geoDist = pairwise geographic distance
  dataBlock <- list("N"=nrow(pwp),
                    "L" = Npolyloci,
                    "hom"=hom,
                    "k" = 0.25,
                    "geoDist"=geoDist)
  
  # run inference
  try(runWM(stanMod = ibsMod,
        dataBlock = dataBlock,
        nChains = 3,
        nIter = 4e3,
        prefix = paste0("WMfit-",run_name)))

}



#idea: compile all lat/long etc. metadata gathered during DIPnet datathon into one file that we can use
#students and Eric C looked up lat/long in published papers and filled into spreadsheets we will read in here by hand
#Eric C QC'd all of these spreadsheets

#11/7/2024
#note there were errors in lat/long for samples SAMN09212008 (SRR7189624) and SAMN09212013 (SRR7189621)
#contacted joseph.pfaller@noaa.gov to get corrected lat/long coords and updated in 
#divdiv/data/abiotic/input_and_working/QCd_latlong_andothermetadata_from_datathon/

#load libraries ------
library(dplyr)    #data handling
library(tidyr)    #data handling
library(readxl)   #read in excel sheets
library(gtools)   #smartbind instead of rbind

rm(list=ls())
gc()

indir = "/Users/rachel/divdiv/data/abiotic/input_and_working/QCd_latlong_andothermetadata_from_datathon/"
outdir = "/Users/rachel/divdiv/data/abiotic/input_and_working"


# ************************************************************************************************************************************
# ************************************************************************************************************************************
# READ IN ALL DATA -----------------

#get list of files to read in
file_list <- list.files(path=indir, pattern=".xlsx", full.names = F)

#read in sheet with data in it
for (file in file_list) {
  
  df <- read_excel(paste(indir, file, sep = ""), sheet = "Samples") %>% as.data.frame()
  assign(paste("df", file, sep = "_"), df)
  
}

#bind all df's together
rm(df)
df <- do.call(smartbind, mget(grep("df_",names(.GlobalEnv),value=TRUE)))

rm(list = ls(pattern = "df_"))

#write out
write.csv(df, paste0(outdir,"/datathon_metadata_all.csv"))




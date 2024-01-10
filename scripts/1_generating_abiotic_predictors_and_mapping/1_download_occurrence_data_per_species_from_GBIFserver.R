#from http://ropensci.org/tutorials/rgbif_tutorial.html
#if on linux and this fails because rgeos fails, see https://stackoverflow.com/questions/38924767/error-installing-r-package-for-linux
#also, make sure ggplot2 is installed, because gbifmap depends on it
#this script originally written by Jamie Pringle

library('rgbif')
library('tidyverse')
library(googledrive)

#load species name and occurrence list for genetic data
#get list of all datasets that we should have
#get master spreadsheet from Drive
using <- googledrive::shared_drive_find(pattern = "^divdiv$")
using <- googledrive::drive_ls(path = using, pattern = "working_datasheets", recursive = FALSE)
using <- googledrive::drive_ls(path = using, pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
using$name
using <- googlesheets4::range_read(using, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
#keep cols we want and only datasets that are "in"
using <- using %>% dplyr::select(organism_biosamp,run_name,link,keepinrunning_YN) %>%
  filter(grepl("yes|Yes",keepinrunning_YN))
speciesList <- using$organism_biosamp
#speciesList<-speciesList[1:2] #for debugging quickly

#this is the directory to save data in
datDir<-'/Users/rachel/divdiv/data/abiotic/input_and_working/speciesLocDat_Sep08_2021_withFilters/'


#this code plots occurrence data from GBIF and OBIS for some marine species
#what species do we care about
for (someSpecies in speciesList)
{
  if (lengths(strsplit(someSpecies,' '))>1) #only do if we have at least species and genus
  {
    #this is the file to save plot in
    datFile<-paste(datDir,'GBIFLocations_',someSpecies,'.csv',sep='')
    
    #now only try to download data and make file if plot file does not exist
    if (file.exists(datFile)) 
    {
      cat(datFile,'exists, skipping\n')
    }
    else
    {
      if (!(someSpecies=="Exaiptasia brasiliensis")) {
        #get the data from GBIF
        print(' ')
        cat('Getting',someSpecies,'from GBIF')
        
        #this is the old way, but fails when large numbers of occurrences
        #key <- name_backbone(name=someSpecies)$speciesKey
        #datGBIF <- occ_data(taxonKey=key, hasCoordinate=TRUE,limit=100000)
        #print(paste('   found',nrow(datGBIF$data)))
        #jnk<-datGBIF$data[,c("decimalLongitude","decimalLatitude")]
        #write.csv(jnk,datFile)
        
        #get taxon key
        taxonKey <- name_backbone(name=someSpecies)$speciesKey
        
        #how bad of an idea can I have?
        #hard code authentication info
        #set up authentication info
        user<-'jpringle' 
        password<- 'Prlx;2533'
        email<-'jpringle@unh.edu'
        
        #spin up search
        #old with fewer filters
        #searchRequest<-occ_download(pred('taxonKey',taxonKey),pred("hasCoordinate", TRUE),format='SIMPLE_CSV',
        #                            user=user,pwd=password,email=email)
        searchRequest<-occ_download(type='and',pred('taxonKey',taxonKey),
                                    pred("hasGeospatialIssue", FALSE),
                                    pred("hasCoordinate", TRUE),
                                    #pred_gte("year", 1900),
                                    pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
                                    pred_or(
                                      pred_not(pred("establishmentMeans","MANAGED")),
                                      pred_not(pred_notnull("establishmentMeans"))
                                    ),
                                    pred_or(
                                      pred("occurrenceStatus","PRESENT"),
                                      pred_not(pred_notnull("occurrenceStatus"))
                                    ),
                                    format='SIMPLE_CSV',
                                    user=user,pwd=password,email=email)
        
                
        #run to find status
        metaStatus=occ_download_meta(searchRequest)
        print(metaStatus)
        
        #now wait for it to be done
        print('waiting')
        occ_download_wait(searchRequest)
        print('done waiting')
        
        #now make file name to save to
        #whereSave<-paste('dataFromGbif',someSpecies,sep='')
        print(paste('getting data for',someSpecies))
        whatGot=occ_download_get(searchRequest,path='dataFromGbif',overwrite=TRUE)
        print('done getting data')
        
        print('importingData')
        datGBIF=occ_download_import(x=whatGot)
        
        jnk<-datGBIF[,c("decimalLongitude","decimalLatitude",'establishmentMeans','occurrenceStatus')]
        #this is the file to save plot in
        write.csv(jnk,datFile)
        paste('done with',someSpecies)
      }
    }
  }
}

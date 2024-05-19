#idea: calc ratio of genetic "sampled size" to GBIF "range size"

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)


rm(list = ls())
gc()

indir = "/Users/rachel/divdiv/data/abiotic/input_and_working/marmapdists-output"
outdir = "/Users/rachel/divdiv/data/abiotic"


#get estimates for params for each dataset -------------
marmap_fileslist <- list.files(path = indir, pattern = "max_and_pw", full.names = TRUE)
#filter to make sure we only use inputs we know are good/done (rather than all that are on google drive) by getting master status sheet from Google Drive
done <- googledrive::shared_drive_find(pattern = "^divdiv_datafiles")
done <- googledrive::drive_ls(path = done, pattern = "master_notes_file", recursive = F)
done$name
done <- googlesheets4::range_read(done, sheet = 2) %>% as.data.frame() %>% dplyr::select(run_name, done_with_marmap) %>% 
  filter(grepl("^bioprj",run_name)) %>% as.data.frame() %>% filter(done_with_marmap == "DONE")
marmap_fileslist <- marmap_fileslist %>% as.data.frame() %>% rename("x"=".") %>% 
  mutate(run_name = gsub("marmapdists-output/max_and_pw_dists|.Robj","",x)) #%>% filter(run_name %in% done$run_name)
marmap_fileslist <- marmap_fileslist$x

df <- data.frame(loop.iter=NA, run_name=NA, 
                 max.95.sea.gbif=NA, max.100.sea.gbif=NA, maxgenetic.sea=NA, 
                 max.95.gcd.gbif=NA, max.100.gcd.gbif=NA, maxgenetic.gcd=NA)

for (loop.iter in 1:length(marmap_fileslist)) {
  
  print(paste0("iteration ", loop.iter, " of ", length(marmap_fileslist)))
  marmapFile = marmap_fileslist[loop.iter]
  load(marmapFile)
  
  #get run_name and stacks params
  run_name = marmapFile %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("max_and_pw_dists","",.) %>% gsub(".Robj","",.)
  
  #extract distances
  max.95.sea.gbif <- max_and_pw_dists$max.95.sea.gbif
  max.100.sea.genetic <- max_and_pw_dists$max.100.sea.genetic
  max.95.gcd.gbif <- max_and_pw_dists$max.95.gcd.gbif
  max.100.gcd.genetic <- max_and_pw_dists$max.100.gcd.genetic
  max.100.sea.gbif <- max_and_pw_dists$max.100.sea.gbif
  max.100.gcd.gbif <- max_and_pw_dists$max.100.gcd.gbif
  
  out <- data.frame(loop.iter = loop.iter, run_name = run_name, 
                    max.95.sea.gbif = max.95.sea.gbif, max.100.sea.gbif = max.100.sea.gbif, maxgenetic.sea = max.100.sea.genetic,
                    max.95.gcd.gbif = max.95.gcd.gbif, max.100.gcd.gbif = max.100.gcd.gbif, maxgenetic.gcd = max.100.gcd.genetic)
  
  df <- rbind(df, out)
  rm(run_name,max.95.sea.gbif,max.100.sea.genetic,max.95.gcd.gbif,max.100.gcd.genetic,max.100.sea.gbif,max.100.gcd.gbif,out)
  
}
df <- df %>% filter(is.na(loop.iter)==F)
rownames(df) <- NULL

#filter to just datasets we are using (so we aren't debugging problems in datasets we've kicked out)
#get master spreadsheet from Drive
# using <- googledrive::shared_drive_find(pattern = "^divdiv$")
# using <- googledrive::drive_ls(path = using, pattern = "working_datasheets", recursive = F)
# using <- googledrive::drive_ls(path = using, pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = F)
# using$name
# using <- googlesheets4::range_read(using, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
# #keep cols we want and only datasets that are "in"
# using <- using %>% dplyr::select(organism_biosamp,run_name,link,keepinrunning_YN) %>%
#   filter(grepl("yes|Yes",keepinrunning_YN))

#check for errors (should return no rows if all good)
df %>% filter(max.100.sea.gbif > 30000 | maxgenetic.sea > 30000)

#drop issues for now (will correct upstream and rerun this script)
#df <- df %>% filter(max.100.sea.gbif < 30000 & maxgenetic.sea < 30000)

#calc ratio of sampled distance to range size distance
df <- df %>% 
  mutate(ratio.sea.95 = maxgenetic.sea/max.95.sea.gbif) %>% 
  mutate(ratio.gcd.95 = maxgenetic.gcd/max.95.gcd.gbif) %>% 
  mutate(ratio.sea.100 = maxgenetic.sea/max.100.sea.gbif) %>%
  mutate(ratio.gcd.100 = maxgenetic.gcd/max.100.gcd.gbif)

#explore ratios
df %>% dplyr::select(run_name,contains("ratio"),everything()) %>% 
  filter(ratio.sea.95 > 1 | ratio.gcd.95 > 1 | ratio.sea.100 > 1 | ratio.gcd.100 > 1) %>% 
  arrange(run_name) %>% 
  View()
  #checked visually, all make sense why ratios ended up > 1 from looking visually at output plots from marmap step

#force ratios > 1 to be 1
df <- df %>% mutate(ratio.sea.95 = ifelse(ratio.sea.95 > 1, 1, ratio.sea.95),
                    ratio.gcd.95 = ifelse(ratio.gcd.95 > 1, 1, ratio.gcd.95),
                    ratio.sea.100 = ifelse(ratio.sea.100 > 1, 1, ratio.sea.100),
                    ratio.gcd.100 = ifelse(ratio.gcd.100 > 1, 1, ratio.gcd.100))

#save for modeling
write.csv(df %>% dplyr::select(-loop.iter), paste0(outdir, "/range_lengths_and_ratios.csv"), row.names = FALSE)




# MISC. ------------------------------

#visualize
df %>% ggplot(aes(x = max.95.sea.gbif, y = max.100.sea.gbif)) + geom_point() + geom_abline(intercept = 0, slope = 1) #+ geom_text(aes(label = run_name))
df %>% ggplot(aes(x = max.95.gcd.gbif, y = max.100.gcd.gbif)) + geom_point() + geom_abline(intercept = 0, slope = 1)


df %>% ggplot() + geom_histogram(aes(x = ratio.sea.95), binwidth = 5000)
df %>% ggplot() + geom_histogram(aes(x = ratio.sea.100), binwidth = 0.1)
df %>% ggplot(aes(x = ratio.sea.95, y = ratio.sea.100)) + geom_point() + geom_abline(intercept = 0, slope = 1) + geom_text(aes(label = run_name))
df %>% filter(ratio.sea.95 < 100) %>% ggplot(aes(x = ratio.sea.95, y = ratio.sea.100)) + geom_point() + geom_abline(intercept = 0, slope = 1)

df %>% filter(ratio.sea.95 <= 1) %>% ggplot() + geom_histogram(aes(x = ratio.sea.95), binwidth = 0.05)

df %>% ggplot(aes(x = ratio.sea.95, y = ratio.gcd.95)) + geom_point() + geom_abline(slope = 1, intercept = 0)
df %>% ggplot(aes(x = ratio.sea.100, y = ratio.gcd.100)) + geom_point() + geom_abline(slope = 1, intercept = 0)


#view the datasets where genetic distance is greater than GBIF
df %>% filter(ratio.sea.95 > 1 | ratio.sea.100 > 1)



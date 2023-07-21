#idea pull in all the .Robjs to get genetic diversity extrapolated values, 
#to build master df for downstream modeling
#run script copy_genetic_outputs_for_building_master_df-local.sh prior to this script (to get all of the inputs for this script onto local computer)


#load libraries
library(dplyr)
library(tidyr)


rm(list = ls())
gc()



# collate response variables/data from WM models ----------
#load("/Users/rachel/WMfit-bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_pars.Robj")

#indir = "/Users/rachel/Desktop/DivDiv/divdiv_data_analysis/ALL_r80_gendiv_data"
#indir = "/Volumes/mnemo2/base_rachel/ALL_r80_gendiv_data"
indir = "/Users/rachel/ALL_r80_gendiv_data"


#get WM og aka Wishart
file_list <- list.files(indir, pattern="WMfit", full.names = TRUE)
file_list <- file_list[grepl("pars.Robj", file_list)]
  
df <- data.frame("species"=NA, "run_name"=NA, "s"=NA, "nugget"=NA, "m"=NA, "nbhd"=NA, "inDeme"=NA, "stacksparams"=NA)

for (loop.iter in 1:length(file_list)) {
  
  gendivFile <- file_list[loop.iter]
  dataset = gendivFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    strsplit(., "\\.") %>% as.data.frame() %>% .[1,]
  run_name = dataset %>% gsub("WMfitwishart-","",.) %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  print(paste0("run_name ", run_name,"; loop ", loop.iter, " of ", length(file_list )))
  stacksparams = dataset %>% gsub("WMfit-","",.) %>% strsplit(., "_") %>% unlist() %>% .[5:10] %>% 
    paste(., sep="", collapse="_")
  species = dataset %>% gsub("WMfitwishart-","",.) %>% strsplit(., "_") %>% unlist() %>% .[3] %>% gsub("-"," ",.)
  gendivFile <- load(gendivFile)
  gendivFile <- data.frame(species = species, run_name = run_name,
                           s = outPars$pt$s, nugget = outPars$pt$nugget, m = outPars$pt$m, nbhd = outPars$pt$nbhd, inDeme = outPars$pt$inDeme, 
                           stacksparams = stacksparams)
  df <- rbind(df, gendivFile)
  
}
df <- df %>% filter(is.na(run_name)==F)

#save
write.csv(df, "data/popgen/r80_WM_stats-wide.csv", row.names = FALSE)




# collate response variables/data from popgen stats ----------
#load("/Users/rachel/popgenstats.0.5.bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_popgenstats.Robj")

indir="/Users/rachel/ALL_r80_popgen_data"

#get WM og aka Wishart
file_list <- list.files(indir, pattern="popgenstats.0.5", full.names = TRUE)

df <- data.frame("species"=NA, "run_name"=NA, "thetaW"=NA, "globalPi"=NA, "meanHet"=NA, "stacksparams"=NA)

for (loop.iter in 1:length(file_list)) {
  
  popgenFile <- file_list[loop.iter]
  dataset = popgenFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    strsplit(., "\\.") %>% as.data.frame() %>% .[4,]
  run_name = dataset %>% gsub("_popgenstats","",.) %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  print(paste0("run_name ", run_name,"; loop ", loop.iter, " of ", length(file_list )))
  stacksparams = dataset %>% gsub("_popgenstats","",.) %>% strsplit(., "_") %>% unlist() %>% .[5:10] %>% 
    paste(., sep="", collapse="_")
  species = dataset %>% gsub("_popgenstats","",.) %>% strsplit(., "_") %>% unlist() %>% .[3] %>% gsub("-"," ",.)
  popgenFile <- load(popgenFile)
  popgenFile <- data.frame(species = species, run_name = run_name,
                           thetaW = popgenstats$thetaW, globalPi = popgenstats$globalPi, meanHet = mean(popgenstats$het), 
                           stacksparams = stacksparams)
  df <- rbind(df, popgenFile)
  
}
df <- df %>% filter(is.na(run_name)==F)

#save
write.csv(df, "data/popgen/r80_popgen_stats-wide.csv", row.names = FALSE)







#graveyard -----

# #from when we were doing wishart and cmplkl - and comparing outputs
# 
# df.1 <- df %>% filter(is.na(species)==F) %>% mutate(model = "ogWM")
# 
# 
# #get WM cpmLnl
# file_list <- list.files(indir, pattern="WMfitcmpLnl-bioprj", full.names = TRUE)
# file_list <- file_list[grepl("out.Robj", file_list)]
# 
# df <- data.frame("species"=NA, "run_name"=NA, "s"=NA, "m"=NA, "nbhd"=NA, "inDeme"=NA, "stacksparams"=NA)
# 
# for (loop.iter in 1:length(file_list)) {
#   
#   print(paste0("loop ", loop.iter, " of ", length(file_list )))
#   gendivFile <- file_list[loop.iter]
#   dataset = gendivFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
#     strsplit(., "\\.") %>% as.data.frame() %>% .[1,]
#   run_name = dataset %>% gsub("WMfitcmpLnl-","",.) %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
#     paste(., sep="", collapse="_")
#   stacksparams = dataset %>% gsub("WMfitcmpLnl-","",.) %>% strsplit(., "_") %>% unlist() %>% .[5:10] %>% 
#     paste(., sep="", collapse="_")
#   species = dataset %>% gsub("WMfitcmpLnl-","",.) %>% strsplit(., "_") %>% unlist() %>% .[3] %>% gsub("-"," ",.)
#   gendivFile <- load(gendivFile)
#   s <- rstan::extract(out$fit,"s",inc_warmup=FALSE)
#   m <- rstan::extract(out$fit,"m",inc_warmup=FALSE)
#   nbhd <- rstan::extract(out$fit,"nbhd",inc_warmup=FALSE)
#   inDeme <- rstan::extract(out$fit,"inDeme",inc_warmup=FALSE)
#   pt <- list("s" = mean(s[[1]]),
#              "m" = mean(m[[1]]),
#              "nbhd" = mean(nbhd[[1]]),
#              "inDeme" = mean(inDeme[[1]]))
#   gendivFile <- data.frame(species = species, run_name = run_name,
#                            s = pt$s, m = pt$m, nbhd = pt$nbhd, inDeme = pt$inDeme, 
#                            stacksparams = stacksparams)
#   df <- rbind(df, gendivFile)
#   
# }
# df.2 <- df %>% filter(is.na(species)==F) %>% mutate(model = "WMcmpLnl")
# 
# 
# df <- merge(df.1 %>% rename("s.wish"="s","m.wish"="m","nbhd.wish"="nbhd","inDeme.wish"="inDeme") %>% dplyr::select(-model), 
#             df.2 %>% rename("s.cmpl"="s","m.cmpl"="m","nbhd.cmpl"="nbhd","inDeme.cmpl"="inDeme") %>% dplyr::select(-model),
#             by = c("species","run_name","stacksparams"),
#             all.x = T, all.y = T)
# 
# df %>% ggplot(aes(x = s.wish, y = s.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1)
# df %>% filter(s.wish > 0.85) %>% ggplot(aes(x = s.wish, y = s.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1) + geom_text(aes(label = run_name))
# df %>% filter(s.wish > 0.85) %>% ggplot(aes(x = s.wish, y = s.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1)
# 
# df %>% ggplot(aes(x = m.wish, y = m.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1)
# 
# df %>% ggplot(aes(x = nbhd.wish, y = nbhd.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1)
# 
# df %>% ggplot(aes(x = inDeme.wish, y = inDeme.cmpl)) + geom_point() + geom_abline(intercept = 0, slope = 1)
# 
# 
# 
# 

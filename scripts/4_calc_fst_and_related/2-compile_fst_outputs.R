#idea pull in all the .Robjs to get genetic diversity extrapolated values, 
#to build master df for downstream modeling
#run script copy_genetic_outputs_for_building_master_df-local.sh prior to this script (to get all of the inputs for this script onto local computer)


#load libraries
library(dplyr)
library(tidyr)


rm(list = ls())
gc()

setwd("/Users/rachel/divdiv")

# collate FST Robjs ----------

indir = "data/popgen/input_and_working/ALL_fst"
file_list <- list.files(indir, pattern="fststats", full.names = TRUE)

#iterate thru each dataset to get summary stats related to spatial distribution of genetic div
out <- data.frame("species"=NA, "run_name"=NA, "fst.global.wc"=NA, 
                  "mean.rowwisefst"=NA, "sd.rowwisefst"=NA, 
                  "mean.1toallfst"=NA, "sd.1toallfst"=NA, 
                  "mean.pwp"=NA, "sd.pwp"=NA)

for (loop.iter in 1:length(file_list)) {
  
  fstFile <- file_list[loop.iter]
  dataset = fstFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    as.data.frame() %>% .[1,] %>% gsub("fststats.","",.) %>% 
    gsub("_Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_") %>% gsub("","",.)
  print(paste0("run_name ", run_name,"; loop ", loop.iter, " of ", length(file_list )))
  stacksparams = dataset %>% gsub(paste0(run_name,"_"),"",.)
  species = dataset %>% strsplit(., "_") %>% unlist() %>% .[3]
  
  load(fstFile)
  
  #global WC FST (value per dataset)
  fst.global.wc <- fststats$fst.wc$FST
  
  #row-wise FST (value per pop)
  if(is.character(fststats$fst.pw.wc)){
    mean.rowwisefst <- "not_calculated_HPCCtimedout"
    sd.rowwisefst <- "not_calculated_HPCCtimedout"
  } else {
    fst.pw.rowwise.mean <- fststats$fst.pw.wc %>% 
      mutate(dummy.pop.id = row.names(.)) %>% rowwise() %>% 
      mutate(rowwise_mean_fst = mean(c_across(2:nrow(.)), na.rm = TRUE)) %>%
      ungroup() %>% dplyr::select(dummy.pop.id,rowwise_mean_fst)
    
    #if only two pops we'll get an NA here for the diagonal for the second, bc only one value total
    #deal with dropping this erroneous NA
    if(nrow(fst.pw.rowwise.mean) == 2){
      mean.rowwisefst <- mean(fst.pw.rowwise.mean[1,2]$rowwise_mean_fst)
      sd.rowwisefst <- "only_two_pops_total"
    } else {
      #mean of rowwise FST (value per dataset)
      mean.rowwisefst <- mean(fst.pw.rowwise.mean$rowwise_mean_fst)
      sd.rowwisefst <- sd(fst.pw.rowwise.mean$rowwise_mean_fst) 
    }
  }
  
  #one to all FST (value per pop)
  mean.1toallfst <- mean(fststats$fst_summaries$fst.wc.1toall)
  sd.1toallfst <- sd(fststats$fst_summaries$fst.wc.1toall)
  
  #pw pi (value per pop)
  mean.pwp <- mean(fststats$fst_summaries$mean.pwp.pop.i)
  sd.pwp <- sd(fststats$fst_summaries$mean.pwp.pop.i)
  
  out.i <- data.frame(species = species, run_name = run_name, 
                      fst.global.wc = fst.global.wc, 
                      mean.rowwisefst = mean.rowwisefst,
                      sd.rowwisefst = sd.rowwisefst,
                      mean.1toallfst = mean.1toallfst,
                      sd.1toallfst = sd.1toallfst,
                      mean.pwp = mean.pwp, sd.pwp = sd.pwp)
  
  out <- rbind(out, out.i)
 
}
out <- out %>% filter(is.na(run_name)==F)

#3 datasets have only 1 indiv per every unique location aka dummy pop
#hierfstat returns 0 for global and pw WC FST when comparing 1 indiv to 1 indiv (regardless of their genotypes)
#recoding these 0s to notes so we know what's up
out$fst.global.wc[out$run_name=="bioprj_PRJNA559677_Caretta-caretta"] = "only_1indiv_per_every_pop"
out$mean.rowwisefst[out$run_name=="bioprj_PRJNA559677_Caretta-caretta"] = "only_1indiv_per_every_pop"
out$sd.rowwisefst[out$run_name=="bioprj_PRJNA559677_Caretta-caretta"] = "only_1indiv_per_every_pop"
out$fst.global.wc[out$run_name=="bioprj_PRJNA480308_Rhizophora-mangle"] = "only_1indiv_per_every_pop"
out$mean.rowwisefst[out$run_name=="bioprj_PRJNA480308_Rhizophora-mangle"] = "only_1indiv_per_every_pop"
out$sd.rowwisefst[out$run_name=="bioprj_PRJNA480308_Rhizophora-mangle"] = "only_1indiv_per_every_pop"
out$fst.global.wc[out$run_name=="bioprj_PRJNA379028_Porites-astreoides"] = "only_1indiv_per_every_pop"
out$mean.rowwisefst[out$run_name=="bioprj_PRJNA379028_Porites-astreoides"] = "only_1indiv_per_every_pop"
out$sd.rowwisefst[out$run_name=="bioprj_PRJNA379028_Porites-astreoides"] = "only_1indiv_per_every_pop"

#save
write.csv(out, "data/popgen/r80_FSTetc_stats-wide.csv", row.names = FALSE)



#graveyard ----------

#calc CVs
out <- out %>% mutate(cv.rowwisefst = sd.rowwisefst/mean.rowwisefst,
                      cv.1toallfst = sd.1toallfst/mean.1toallfst, 
                      cv.pwp = sd.pwp/mean.pwp)

# viz exploration
out %>% ggplot(aes(x=mean.1toallfst, y=mean.rowwisefst)) + 
  geom_point() +
  geom_smooth(formula = "y~x", method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()
ggsave("figures/fst_mean_comparions1.pdf", 
       width = 6, height = 6, units = c("in"), dpi = 600)

out %>% ggplot(aes(x=fst.global.wc, y=mean.1toallfst)) + 
  geom_point() +
  geom_smooth(formula = "y~x", method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()
ggsave("figures/fst_mean_comparions2.pdf", 
       width = 6, height = 6, units = c("in"), dpi = 600)

out %>% ggplot(aes(x=fst.global.wc, y=mean.rowwisefst)) + 
  geom_point() +
  geom_smooth(formula = "y~x", method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()
ggsave("figures/fst_mean_comparions3.pdf", 
       width = 6, height = 6, units = c("in"), dpi = 600)





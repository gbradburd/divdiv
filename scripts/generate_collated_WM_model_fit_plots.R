#idea: generate model fit plots for all datasets, collated onto pages, for WM gendiv extrapolating model


#load libraries
library(dplyr)
library(tidyr)

rm(list = ls())
gc()

#for testing
#load(paste0(indir,"WMfit-bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_initPars.Robj"))
#load(paste0(indir,"WMfit-bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_out.Robj"))
#load(paste0(indir,"WMfit-bioprj_PRJNA294760_Amphiprion-bicinctus_stacks_littlem_3_bigm_3_n2_nis1lessthanM_pars.Robj"))

indir = "/Users/rachel/Desktop/DivDiv/divdiv_data_analysis/ALL_r80_gendiv_data/"
file_list <- list.files(path = indir, pattern = "_out.Robj", full.names = TRUE)

#load a few functions
getPhom <- function(model.fit,chain.no,N){
  par.cov <- array(NA,dim=c(model.fit@sim$n_save[chain.no],N,N))
  for(i in 1:N){
    for(j in 1:N){
      my.par <- sprintf("pHom[%s,%s]",i,j)
      par.cov[,i,j] <- rstan::extract(model.fit,pars=my.par,inc_warmup=FALSE,permuted=FALSE)[,chain.no,]
    }
  }
  return(par.cov)
}

plotFit <- function(out,pHom,chainCol,run_name,stacksparams){
  ut <- upper.tri(out$dataBlock$geoDist,diag=TRUE)
  plot(out$dataBlock$geoDist[ut],out$dataBlock$hom[ut],
       ylim=range(c(pHom,out$dataBlock$hom)),
       xlab="pairwise distance",ylab="pairwise homozygosity",
       pch=19,col=adjustcolor(1,0.05),
       main = paste0(run_name,"; ",stacksparams))
  invisible(
    lapply(seq(1,250,length.out=25),function(i){
      points(out$dataBlock$geoDist[ut],pHom[i,,][ut],pch=20,col=adjustcolor(chainCol,0.05))
    }))
}


#open PDF to save all plots to
pdf(file = "figures/WM_collated_model_fits-chain1.pdf", width = 14, height = 10)

#generate model fit plots
for ( wmOutfile in file_list ) {
  
  load(wmOutfile)
  
  #get labels
  dataset = wmOutfile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("WMfit-","",.) %>% gsub(".Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  stacksparams = dataset %>% strsplit(., "_") %>% unlist() %>% .[5:10] %>% 
    paste(., sep="", collapse="_")
  print(paste0("processing ", run_name," with params ", stacksparams))
  
  #get inputs for plot
  try(post <- rstan::get_logposterior(out$fit,inc_warmup=FALSE))
  try(pHom <- invisible(lapply(1:length(post),
                           function(i){
                             getPhom(out$fit,i,out$dataBlock$N)})))
  nChains <- length(post)
  chainCols <- c("blue","goldenrod1","red","forestgreen","purple","black")[1:nChains]
  
  #make plot
  #change below to 1:length(post) to make a plot for each chain
  for(i in 1:1){
    par(mfrow=c(1,1))
    try(plotFit(out,pHom[[i]],chainCols[i],run_name,stacksparams))
  }
  
}

#save pdf
dev.off()




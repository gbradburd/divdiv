#idea: download trait data from Google drive,
#calculate and explore trait coverage - by trait and by species,
#explore trait-trait correlations

#load libraries
library(dplyr)
library(tidyr)
library(googledrive)
library(googlesheets4)
library(stringr)
library(ggplot2)


rm(list = ls())
gc()

#define path to save outputs to
workingdir = "/Users/rachel/divdiv/data/biotic/inputs_and_working"
outdir = "/Users/rachel/divdiv/data/biotic"



# get trait data -------

#get master spreadsheet from Drive
df.raw  <- googledrive::shared_drive_find(pattern = "^divdiv$")
df.raw  <- googledrive::drive_ls(path = df.raw , pattern = "working_datasheets", recursive = FALSE)
df.raw  <- googledrive::drive_ls(path = df.raw , pattern = "working_list_marine_projects_with_10indivs-12-4-2020", recursive = FALSE)
df.raw $name
df.raw <- googlesheets4::range_read(df.raw , sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))

# this is older way we were building the trait df
# #keep traits we care about (one row per unique entry, species level)
# df <- df.raw %>% filter(grepl("yes|Yes",keepinrunning_YN))
# df <- df %>% dplyr::select(organism_biosamp, somePopsFWorBrackish, isInvasiveRange, isMaybeHybrid, 
#                                DiffbyDepthSuspected, isAdultBenthicMarine, isNektonSometimes, 
#                                PLD_point, PLD_Min, PLD_Max, PLD_point2,
#                                isPlanktonic_atanypoint, Spawning_mode, Larval_feeding,
#                                Large_Adult_Range, ReturnToSpawningGround,
#                                Fecundity_EggsFemaleSpawn, Fecundity_SpawnFrequency, Fecundity_EggSize, 
#                                Body_Size, Reproductive_Age, Generational_Structure,
#                                Asexual_Stage, Census_Size,
#                                Genome_Size_MB, Genome_Size_picoGrams, Coded_longevity, LogEgg) %>% 
#   distinct()

# new/current way to build trait df, including only traits and species that were curated in final trait curation effort
df <- df.raw %>% filter(final_keepers_with_popgen_data == "TRUE" | final_keepers_with_popgen_data == "TOBEADDED")
nrow(df)
df <- df %>% dplyr::select(organism_biosamp, 
                           PLD_point, PLD_Min, PLD_Max, PLD_point2,
                           isPlanktonic_atanypoint, 
                           Spawning_mode, Larval_feeding, ReturnToSpawningGround,
                           Fecundity_EggSize, Body_Size, Generational_Structure,
                           isBenthic, Large_Adult_Range,
                           isInvert, isFish, isAmniote, isPlant, isCoral)

#check if there are any species with non-identical info
#should return 0 if no duplicates
df %>% group_by(organism_biosamp) %>% mutate(nonid = 1:n()) %>% filter(nonid > 1) %>% nrow()
# get a list of the non-identical species
dups <- df %>% group_by(organism_biosamp) %>% mutate(nonid = 1:n()) %>% mutate(max = max(nonid)) %>% filter(max > 1)
nrow(dups)
dups

#check coding of variables
for ( trait in names(df) ) {
  print(df %>% group_by(!!as.symbol(trait)) %>% summarise(n=n()))
}

#replace values that aren't actually data with true NAs - NULL, ?, TBD, NA
df <- df %>% 
  mutate(across(everything(), gsub, pattern = "unknown", replacement = NA)) %>% 
  mutate(across(everything(), gsub, pattern = "NULL", replacement = NA)) %>% 
  mutate(across(everything(), gsub, pattern = "NA", replacement = NA)) %>% 
  mutate(across(everything(), gsub, pattern = "\\?", replacement = NA)) %>% 
  mutate(across(everything(), gsub, pattern = "TBD", replacement = NA))

#check coding of variables again
for ( trait in names(df) ) {
  print(df %>% group_by(!!as.symbol(trait)) %>% summarise(n=n()))
}

# #change ranges to midpoints
# range_midpoint <- function(x){
#   split <- str_split(x, "–")
#   midpoint <- mean(as.numeric(split[[1]]))
#   return(midpoint)
# }
# 
# df$Body_Size <- purrr::map_dbl(df$Body_Size, range_midpoint)
# df$Reproductive_Age <- purrr::map_dbl(df$Body_Size, range_midpoint)
# df$Fecundity_EggsFemaleSpawn <- purrr::map_dbl(df$Fecundity_EggsFemaleSpawn, range_midpoint)
# df$Fecundity_SpawnFrequency <- purrr::map_dbl(df$Fecundity_SpawnFrequency , range_midpoint)

write.csv(df, file = paste0(workingdir,"/marinerds_traits_05-19-2024.csv"), row.names = FALSE)



# convert trait data to numerically coded for analysis/modeling ---------------
rm(list = ls() %>% stringr::str_subset(., c("workingdir|outdir"), negate = T))
gc()

df <- read.csv(paste0(workingdir,"/marinerds_traits_05-19-2024.csv"))

#recode some variables
#for dispersal related traits, larger values mean higher dispersal
df <- df %>% 
  mutate(isPlanktonic_atanypoint = ifelse(isPlanktonic_atanypoint == "TRUE", 1, 0)) %>% 
  mutate(ReturnToSpawningGround = ifelse(ReturnToSpawningGround == "TRUE", 1, 0)) %>% 
  mutate(Large_Adult_Range = ifelse(Large_Adult_Range == "TRUE", 1, 0)) %>% 
  mutate(isInvert = ifelse(isInvert == "TRUE", 1, 0)) %>%
  mutate(isFish = ifelse(isFish == "TRUE", 1, 0)) %>%
  mutate(isAmniote = ifelse(isAmniote == "TRUE", 1, 0)) %>%
  mutate(isPlant = ifelse(isPlant == "TRUE", 1, 0))
df$Spawning_mode[df$Spawning_mode == "I"] = 0
df$Spawning_mode[df$Spawning_mode == "N"] = 1
df$Spawning_mode[df$Spawning_mode == "F"] = 2
df$Larval_feeding[df$Larval_feeding == "N"] = 0
df$Larval_feeding[df$Larval_feeding == "L"] = 1
df$Larval_feeding[df$Larval_feeding == "P"] = 2
df$Generational_Structure[df$Generational_Structure == "I"] = 0
df$Generational_Structure[df$Generational_Structure == "S"] = 1
df$isBenthic[df$isBenthic == "A"] = 0
df$isBenthic[df$isBenthic == "S"] = 1
df$isBenthic[df$isBenthic == "N"] = 2
#make sure numbers are numbers
df <- df %>% 
  mutate(Spawning_mode = as.numeric(Spawning_mode), 
         Larval_feeding = as.numeric(Larval_feeding),
         Generational_Structure = as.numeric(Generational_Structure),
         isBenthic = as.numeric(isBenthic))

#check coding of variables again
for ( trait in names(df) ) {
  print(df %>% group_by(!!as.symbol(trait)) %>% summarise(n=n()))
}
str(df)

#save cleaned, numerically coded traits
write.csv(df, paste0(outdir,"/cleaned_numeric_biotic_traits.csv"), row.names = FALSE)




# trait-trait correlation matrix ------

#function to make correlation plot with pretty colors and ellipses
my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"),
                         lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"),
                         digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "",
                         asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
  # this is a modified version of the plotcorr function from the ellipse package
  # modified by Esteban Buz
  if (!require('ellipse', quietly = TRUE, character = TRUE)) {
    stop("Need the ellipse library")
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}

#can only put numerically coded traits into correlation table
df.numeric <- df %>% dplyr::select(PLD_point2, isPlanktonic_atanypoint, Spawning_mode, 
                                   Larval_feeding, ReturnToSpawningGround, Fecundity_EggSize, 
                                   Body_Size, Generational_Structure, isBenthic, Large_Adult_Range)
for ( trait in names(df.numeric) ) {
  print(df %>% group_by(!!as.symbol(trait)) %>% summarise(n=n()))
}
str(df.numeric)

#cor table using pairwise obs
#aka for each cell of table, use all indivs that have those data for those two trts, repeat for each cell
trtcor.pairwise <- df.numeric %>% cor(., use = "pairwise.complete.obs", method = "pearson")
#note that log() handles NA values and just carries over/prints NA whenever a value passed to log() is missing
trtcor.pairwise

# set colors
# colors are selected from a list with colors recycled
# Thus to map correlations to colors we need to make a list of suitable colors
# To start, pick the end (and mid) points of a scale, here a red to white to blue for neg to none to pos correlation
colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
# Build a ramp function to interpolate along the scale, I've opted for the Lab interpolation rather than the default rgb,
#check the documentation about the differences
colramp = colorRampPalette(colsc, space='Lab')
# make a scale with 100 points
colors = colramp(100)
# then pick colors along this 100 point scale given the correlation value * 100 rounded down to the nearest integer
# to do that we need to move the correlation range from [-1, 1] to [0, 100]
# col=colors[((corr.mtcars + 1)/2) * 100]

#plot
my.plotcorr(trtcor.pairwise,
            upper.panel="number", lower.panel = "ellipse", diag = "none",
            mar=c(0,0,0,0),
            outline = F,
            col=colors[((trtcor.pairwise + 1)/2) * 100],
            axes = FALSE, xlab = "", ylab = "",
            cex.lab = 0.6, cex = 0.6)
# !!! have to save manually above figure (ggsave() doesn't work) !!!
#pred_correlation_matrix-biotic_traits.pdf , 14x14 in

# end


#idea: combine all predictors and response dataframes into one master df for input to models

#load libraries
library(dplyr)
library(tidyr)

rm(list = ls())
gc()


#read in data
lat <- read.csv("data/abiotic/mean_min_max_latitude_stats-wide.csv") %>% dplyr::select(-X) %>% mutate(link = paste0("bioprj_",link))
ecor <- read.csv("data/abiotic/number_of_ecoregions_stats-wide.csv") %>% mutate(species = gsub(" ","_",species))
rangesize <- read.csv("data/abiotic/range_lengths_and_ratios.csv") %>% dplyr::rename("link"="run_name")

biotic <- read.csv("data/biotic/cleaned_numeric_biotic_traits.csv") %>% dplyr::rename("species"="organism_biosamp") %>% 
  mutate(species = gsub(" ","_",species))

popg <- read.csv("data/popgen/r80_popgen_WM_stats-wide.csv") %>% dplyr::rename("link"="run_name") %>% dplyr::select(-species)


#keep just traits in hypoth bingo
biotic <- biotic %>% dplyr::select(species, Body_Size, Fecundity_EggSize, Generational_Structure, 
                                   ReturnToSpawningGround, Spawning_mode, Larval_feeding,
                                   PLD_point2, isPlanktonic_atanypoint)

#make all latitude values positive
lat %>% ggplot() + geom_histogram(aes(x = meanlat.gbif)) + scale_x_continuous(limits = c(-90,90))
lat <- lat %>% mutate(meanlat.gbif = abs(meanlat.gbif))
lat %>% ggplot() + geom_histogram(aes(x = meanlat.gbif)) + scale_x_continuous(limits = c(-90,90))

#sync some species names
biotic$species[biotic$species == "Exaiptasia_pallida"] = "Exaiptasia_diaphana" 
biotic$species[biotic$species == "Seriola_dorsalis"] = "Seriola_lalandi_dorsalis" 

#name issues to correct/sync
#Exaiptasia diaphana in gbif ecor
#Exaiptasia diaphana in gbif mean/min/max latitudes
#Exaiptasia pallida in full list

#Seriola dorsalis in full list
#Seriola lalandi dorsalis in gbif ecor
#Seriola lalandi dorsalis in gbif mean/min/max latitudes


#do fake merge first to check that name mismatches are corrected
df <- merge(lat, ecor, by = "species", all = T)
df <- merge(df, rangesize, by = "link", all = T)
df <- merge(df, biotic, by = "species", all = T)
df <- merge(df, popg, by = "link", all = T)
df %>% filter(grepl("Exaiptasia",link)) %>% dplyr::select(link,species)
df %>% filter(grepl("Exaiptasia",species)) %>% dplyr::select(link,species)
df %>% filter(grepl("Seriola",link)) %>% dplyr::select(link,species)
df %>% filter(grepl("Seriola",species)) %>% dplyr::select(link,species)


#do real merge
df <- merge(lat, ecor, by = "species", all = F)
df <- merge(df, rangesize, by = "link", all = F)
df <- merge(df, biotic, by = "species", all = F)
df.preds <- df
df <- merge(df, popg, by = "link", all = F)

names(df)


#save
#write.csv(df, paste0("data/master_df-",Sys.Date(),".csv"), row.names = FALSE)
write.csv(df, "data/master_df.csv", row.names = FALSE)




#**********************************************
#look at predictor correlations -------

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

#correlations for datapoints in model aka with popgen data right now --------------
#grab our predictors (all need to be numeric for code to work)
df.using <- df %>% dplyr::select(meanlat.gbif, n_ECOREGIONS.all, maxgbif.sea_km, ratio.sea,
                                 Body_Size, Fecundity_EggSize, Generational_Structure, ReturnToSpawningGround,
                                   Spawning_mode, Larval_feeding, PLD_point2, isPlanktonic_atanypoint)
str(df.using)

# !!! JUST FOR NOW - drop one outlier sea ratio point !!! ------
df.using$ratio.sea[df.using$ratio.sea > 100000] = NA

#cor table using pairwise obs
#aka for each cell of table, use all indivs that have those data for those two trts, repeat for each cell
trtcor.pairwise <- df.using %>% cor(., use = "pairwise.complete.obs", method = "pearson")
trtcor.pairwise

# set colors
colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

#plot
my.plotcorr(trtcor.pairwise,
            upper.panel="number", lower.panel = "ellipse", diag = "none",
            mar=c(0,0,0,0),
            outline = F,
            col=colors[((trtcor.pairwise + 1)/2) * 100],
            axes = FALSE, xlab = "", ylab = "",
            cex.lab = 0.6, cex = 0.6)

#have to save manually above figure (ggsave() doesn't work)



#correlations for all potential datapoints --------------
#grab our predictors (all need to be numeric for code to work)
df.using <- df.preds %>% dplyr::select(meanlat.gbif, n_ECOREGIONS.all, maxgbif.sea_km, ratio.sea,
                                 Body_Size, Fecundity_EggSize, Generational_Structure, ReturnToSpawningGround,
                                 Spawning_mode, Larval_feeding, PLD_point2, isPlanktonic_atanypoint)
str(df.using)

# !!! JUST FOR NOW - drop one outlier sea ratio point !!! ------
df.using$ratio.sea[df.using$ratio.sea > 100000] = NA

#cor table using pairwise obs
#aka for each cell of table, use all indivs that have those data for those two trts, repeat for each cell
trtcor.pairwise <- df.using %>% cor(., use = "pairwise.complete.obs", method = "pearson")
trtcor.pairwise

# set colors
colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

#plot
my.plotcorr(trtcor.pairwise,
            upper.panel="number", lower.panel = "ellipse", diag = "none",
            mar=c(0,0,0,0),
            outline = F,
            col=colors[((trtcor.pairwise + 1)/2) * 100],
            axes = FALSE, xlab = "", ylab = "",
            cex.lab = 0.6, cex = 0.6)

#have to save manually above figure (ggsave() doesn't work)

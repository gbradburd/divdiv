#idea: read in all ecoregion files per species and process
#to get total number of ecoregions each species falls in (GBIF points)

#load libraries -----
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()


indir="/Users/rachel/divdiv/data/abiotic/input_and_working/ecoregions-output"
outdir="/Users/rachel/divdiv/data/abiotic"

#read in and collate data --------
#path to where individual ecoregion files live
# ! NOTE ! - /ecoregions-output directory too big to live in git repo itself (~600MB), lives on Google Drive (divdiv_datafiles:ecoregions-output)

ecor_files_list <- list.files(indir, pattern="number_of_ecoregions", full.names = TRUE)

df <- data.frame("point_i"=NA, "FID"=NA, "ECO_CODE"=NA, "ECOREGION"=NA, "PROV_CODE"=NA, "PROVINCE"=NA, 
                   "RLM_CODE"=NA, "REALM"=NA, "ALT_CODE"=NA, "ECO_CODE_X"=NA, "Lat_Zone"=NA, "SHAPE_Leng"=NA, 
                   "SHAPE_Area"=NA, "Realm_Costello"=NA, "species"=NA, "dataset"=NA)
for (loop.iter in 1:length(ecor_files_list)) {
  
  print(paste0("loop ", loop.iter, " of ", length(ecor_files_list)))
  ecorFile <- ecor_files_list[loop.iter]
  dataset = ecorFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    strsplit(., "\\.") %>% as.data.frame() %>% .[2,]
  ecorFile <- read.csv(ecorFile) %>% mutate(dataset = dataset)
  if (length(names(ecorFile) %>% grep("link", .)) != 0) {
    ecorFile <- ecorFile %>% rename("species" = "link")
  }
  df <- rbind(df, ecorFile)
  
}
df <- df %>% filter(is.na(point_i)==F)

#match sp names between GBIF and genetic outputs
df <- df %>% dplyr::rename("sp" = "species") %>% tidyfast::dt_separate(sp, into = c("x","y","sp.temp"), sep = "_", remove = F) %>% 
  mutate(sp.temp = gsub("-", " ", sp.temp)) %>% mutate(species = ifelse(is.na(sp.temp)==F, sp.temp, sp)) %>% 
  dplyr::select(-sp.temp,-sp, -x, -y)

#calc number of regions each species falls in (number of unique regions across/including genetic and GBIF points)
eco.all <- df %>% group_by(species,ECO_CODE) %>% summarise(n=n()) %>% group_by(species) %>% summarise(n_ECOREGIONS.all = n())
prov.all <- df %>% group_by(species,PROV_CODE) %>% summarise(n=n()) %>% group_by(species) %>% summarise(n_PROVINCES.all = n())
realm.all <- df %>% group_by(species,REALM) %>% summarise(n=n()) %>% group_by(species) %>% summarise(n_REALMS.all = n())
realm.c.all <- df %>% group_by(species,Realm_Costello) %>% summarise(n=n()) %>% group_by(species) %>% summarise(n_Realm_Costello.all = n())

#calc number of regions each species falls in - just GBIF points and just genetic points
eco <- df %>% group_by(species,ECO_CODE,dataset) %>% summarise(n=n()) %>% group_by(species,dataset) %>% summarise(n_ECOREGIONS = n()) %>%
  pivot_wider(., names_from = dataset, values_from = n_ECOREGIONS, names_prefix = "n_ECOREGIONS.") 
prov <- df %>% group_by(species,PROV_CODE,dataset) %>% summarise(n=n()) %>% group_by(species,dataset) %>% summarise(n_PROVINCES = n()) %>% 
  pivot_wider(., names_from = dataset, values_from = n_PROVINCES, names_prefix = "n_PROVINCES.")
realm <- df %>% group_by(species,REALM,dataset) %>% summarise(n=n()) %>% group_by(species,dataset) %>% summarise(n_REALMS = n()) %>% 
  pivot_wider(., names_from = dataset, values_from = n_REALMS, names_prefix = "n_REALMS.")
realm.c <- df %>% group_by(species,Realm_Costello,dataset) %>% summarise(n=n()) %>% group_by(species,dataset) %>% summarise(n_Realm_Costello = n()) %>% 
  pivot_wider(., names_from = dataset, values_from = n_Realm_Costello, names_prefix = "n_Realm_Costello.")

#merge all stats together
df.summary <- full_join(eco.all, prov.all, by = "species") %>% 
  full_join(., realm.all, by = "species") %>% full_join(., realm.c.all, by = "species") %>% 
  full_join(., eco, by = "species") %>% full_join(., prov, by = "species") %>% 
  full_join(., realm, by = "species") %>% full_join(., realm.c, by = "species") %>% 
  pivot_longer(., names_to = "level", values_to = "n_regions", cols = 2:ncol(.))

#save
#write.csv(df.summary, paste0(outdir,"/number_of_ecoregions_stats-long.csv", row.names = FALSE))

df.summary.w <- df.summary %>% pivot_wider(., names_from = level, values_from = n_regions)
write.csv(df.summary.w, paste0(outdir,"/number_of_ecoregions_stats-wide.csv"), row.names = FALSE)





# MISC. -------------------


#see which ones are missing
df.summary %>% filter(is.na(n_regions))

# plot -------
df.summary %>% ggplot(aes(x = n_regions)) + 
  geom_histogram(binwidth = 1, colour = "black", fill = "gray70") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(~level, scales = "free") + 
  theme(strip.text = element_text(size = 4.5),
        axis.text.x = element_text(angle = 45))

df.summary %>% filter(level %in% c("n_Realm_Costello.all","n_REALMS.all")) %>%
  pivot_wider(., names_from = "level", values_from = "n_regions") %>%
  ggplot(aes(x = n_REALMS.all, y = n_Realm_Costello.all)) +
  geom_point() + geom_smooth(formula = y~x, method = "lm") + 
  ggpmisc::stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)


#plot correlation matrix

#function to make correlation plot with pretty colors and ellipses
my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"),
                         lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"),
                         digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "",
                         asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
  # this is a modified version of the plotcorr function from the ellipse package
  # this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
  # diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
  # digits specifies the number of digits after the . to round to
  # unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
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

#usage of my.plotcorr
#much like the my.plotcorr function, this is modified from the plotcorr documentation
#this function requires the ellipse library, though, once installed you don't need to load it - it is loaded in the function
#install.packages(c('ellipse'))
#library(ellipse)

#cor table using pairwise obs
#aka for each cell of table, use all indivs that have those data for those two trts, repeat for each cell
trtcor.pairwise <- df.summary %>% 
  pivot_wider(., names_from = level, values_from = n_regions) %>% dplyr::select(-species) %>% 
  cor(., use = "pairwise.complete.obs", method = "pearson")
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
            cex.lab = 0.7, cex = 0.7)

#have to save manually above figure (ggsave() doesn't work)
#saved as: /divdiv/figures/CORMATRIX-number_of_ecoregions_stats.png



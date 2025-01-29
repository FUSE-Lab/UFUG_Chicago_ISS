
library('lidR')
library('leafR')
library('stringr')
library('dplyr')
library("viridis")
library("sf")
library("sp")
library("terra")
library("doParallel")
library("data.table")
library("doSNOW")
library("future")
library("lwgeom")
library("raster")

#####################################################################################
#### Structural Metric calculations
# These are setup in two functions and then several from the leafR package in a for loop
##################################################################
######################################################
##################### Functions - From Jianmin
##################### Functions 
##################### Functions 
{
  ### metrics based on chms
  ### z0 is considered as the cut-off height of shrub/grass
  chm_metrics <- function(grid_30m, z0){
    grid_30m@data$Z[grid_30m@data$Z < z0] <- 0
    #if then statment to filter out any bad plots (no points or no points > 0)
    if(sum(grid_30m@data$Z) > 0) {
      chm <- lidR::grid_canopy(grid_30m, res = 1, p2r()) # this uses all the points in the chm
      #metrics from a chm (grid_canopy)
      rumple <- rumple_index(chm)
      mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE)
      max.canopy.ht <- max(chm@data@values, na.rm=TRUE)
      top.rugosity <- sd(chm@data@values, na.rm = TRUE)
      #deep gap fraction (fraction of total rasterized area (1 m2 cells) w/points that are
      #zero height (> 3 m))
      cells <- length(chm@data@values)
      chm.0 <- chm
      chm.0[is.na(chm.0)] <- 0
      deepgaps <- sum(chm.0@data@values == 0)
      deepgap.fraction <- deepgaps/cells
      out.chm <- c(rumple, top.rugosity, mean.max.canopy.ht, max.canopy.ht, deepgap.fraction)
    } else {
      #set rumple to 0, gap fraction to 1 if there were not points > .5 m height or no points
      # out.chm <- c(0, 0, 0, 0, 1)
      out.chm <- rep(NA, 5)
    }
    names(out.chm) <- c("rumple", "top.rugosity", "mean.max.canopy.ht",   "max.canopy.ht", "deepgap.fraction")
    return(out.chm)
  }
  structural_diversity_metrics <- function(grid_30m, z0) {
    ##### added by JW
    grid_30m@data$Z[grid_30m@data$Z < z0] <- NA
    Zs <- grid_30m@data$Z
    Zs <- Zs[!is.na(Zs)]
    if(length(Zs)  > 1 ) {
      #metrics from cloud_metrics
      vert.sd <- sd(Zs, na.rm = TRUE)  # <- corrected
      meanH <- mean(Zs, na.rm = TRUE)   # <- corrected
      vertCV <- vert.sd / meanH
      #metrics from grid_metrics
      sd.9m2 <- lidR::grid_metrics(grid_30m, ~sd(Z, na.rm = T), 3) #9 m2 grid for sd.sd, na.rm =T and add ~ JW
      sd.sd <- sd(sd.9m2@data@values, na.rm = TRUE)
      mean.sd <- mean(sd.9m2@data@values, na.rm = TRUE)
      # vertq <- quantile(las_file$Z, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE)  # <- corrected
      ### Gini coefficient
      # Change GC function to load laz file from the envrionment, not from directory
      n <- length(Zs)
      x <- sort(Zs)
      G <- 2 * sum(x * 1L:n)/sum(x) - (n + 1L)
      GC <- G/(n - 1L)
      out.plot <- c(meanH, vert.sd, vertCV, mean.sd, sd.sd, GC) ##, vertq
    } else{
      out.plot <- rep(NA, 6)
    }
    names(out.plot) <- c("meanH", "vert.sd",  "vertCV", "mean.sd", "sd.sd", "GC")
    ## ,    "q0", "q1", "q5", "q25", "q50", "q75", "q95", 'q99', "q100"
    return(out.plot)
  }
  LAI_lidR <- function(grid_30m, k, z0){
    # grid_30m <- las_file
    #metrics from a Z vector
    Zs <- grid_30m@data$Z
    Zs <- Zs[!is.na(Zs)]
    if(length(Zs) > 1) {
      gap_frac <- gap_fraction_profile(Zs, dz = 1, z0 = z0) #ignore points < 3 m
      GFP <- mean(gap_frac$gf) ###the proportion of accumulated number of points in neighboring layers (z0 to x + dz) / (z0 to x)
      # LADen <- lidR::LAD(Zs, dz = 1, k=k, z0= z0) #ignore points < 3 m
      # ##### LADen <- LAD2(Zs, dz = 1, k=k, z0= z0) #ignore points < 3 m
      # VAI <- sum(LADen$lad, na.rm=TRUE)
      Zs <- Zs[Zs>=z0]  #### add by Jianmin
      VCI <- lidR::VCI(Zs, by = 1, zmax=100)
      out.plot <- c(GFP, VCI)  ###VAI,
    } else{
      out.plot <- c(NA, NA) ##3NA,
    }
    names(out.plot) <- c('GFP', 'VCI') ###'LAI.lidR',
    return(out.plot)
  }
  # Change lad.voxels function to load laz file from the envrionment, not from directory
  lad.voxels2 <- function (normlas.file, grain.size = 3, k = 1) {
    # normlas.file <- las_file_38
    # normlas.file <- las_file; grain.size = 30; k = 0.5
    LAD_VOXELS <- list()
    Z <- NA
    .las <- normlas.file
    .las@data$Z[.las@data$Z < 0] <- 0
    maxZ <- floor(max(.las@data$Z))
    func <- formula(paste0("~leafR::pointsByZSlice(Z, ", maxZ, ")"))
    t.binneds <- lidR::grid_metrics(.las, func, res = grain.size,
                                    start = c(min(.las@data$X), max(.las@data$Y)))
    # t.binneds <- rasterToPoints(t.binneds, na.rm =F )
    ### If it only has 1 band, then raster::values will not record names.
    ### It it has multiple bands, then it will record the names in the matrix or array.
    ### So add the names to the 1 band
    nlayer <- raster::nlayers(t.binneds)
    if(nlayer <= 1 ) return(NULL)
    # if(nlayer == 1)  colname <- names(t.binneds)
    t.binneds <- data.frame(sp::coordinates(t.binneds), raster::values(t.binneds))
    names(t.binneds)[1:2] <- c("X", "Y")
    # if(nlayer == 1) names(t.binneds)[3]  <- colname
    ### Difference than rasterToPoints is that it includes NA
    pulses.profile.dz1 <- t.binneds[, length(t.binneds):3]
    total.pulses.matrix.dz1 = matrix(apply(pulses.profile.dz1,
                                           1, sum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1))
    cumsum.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, cumsum),
                               ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1),
                               byrow = TRUE)
    rm(pulses.profile.dz1)
    pulse.out.dz1 = total.pulses.matrix.dz1 - cumsum.matrix.dz1
    if (nrow(t.binneds) != 1) {
      pulse.in.dz1 <- cbind(total.pulses.matrix.dz1[, 1],
                            pulse.out.dz1[, -c(ncol(pulse.out.dz1))])
    } else {
      pulse.in.dz1 <- c(total.pulses.matrix.dz1[, 1], pulse.out.dz1[,
                                                                    -c(ncol(pulse.out.dz1))])
    }
    rm(total.pulses.matrix.dz1, cumsum.matrix.dz1)
    dz = 1
    LAD.dz1 = log(pulse.in.dz1/pulse.out.dz1) * 1/k * 1/dz
    rm(pulse.in.dz1, pulse.out.dz1)
    LAD.dz1[is.infinite(LAD.dz1)] <- NA
    LAD.dz1[is.nan(LAD.dz1)] <- NA
    LAD.dz1 <- LAD.dz1[, -c(ncol(LAD.dz1))] ### This step will cause ncol-1 columns and will make
    LAD_VOXELS[["LAD"]] <-  matrix(LAD.dz1, ncol = nlayer - 1)
    LAD_VOXELS[["coordenates"]] = t.binneds[, c("X", "Y")]
    rm(LAD.dz1, t.binneds)
    return(LAD_VOXELS)
  }
  LAI_leafR <- function(grid_30m, grain.size, z0, k){
    VOXELS_LAD <- lad.voxels2(grid_30m, grain.size = grain.size, k = k)  ### try i <- 38, 60
    if(!is.null(VOXELS_LAD)){
      lad_profile <- leafR::lad.profile(VOXELS_LAD, relative = T) #uses relative total LAI
      lad_profile <- lad_profile %>% dplyr::filter(height >= z0)
      fhd <- leafR::FHD(lad_profile, evenness = FALSE, LAD.threshold = -1) #default
      lad_profile_lai <- leafR::lad.profile(VOXELS_LAD, relative = F) #doesn't use relative LAI
      LAI <- leafR::lai(lad_profile_lai, min = z0, max = 100)
      # LAI_subcanopy <- leafR::lai(lad_profile_lai, min = 1, max = 5)
      out.plot <- c(fhd, LAI)
    } else {
      out.plot <- c(NA, NA)
    }
    names(out.plot) <- c('FHD', 'LAI')
    return(out.plot)
  }
  LAIsubcanopy_leafR <- function(grid_30m, grain.size, z0, max, k){
    VOXELS_LAD <- lad.voxels2(grid_30m, grain.size = grain.size, k = k)  ### try i <- 38, 60
    if(!is.null(VOXELS_LAD)){
      lad_profile <- leafR::lad.profile(VOXELS_LAD, relative = T) #uses relative total LAI
      lad_profile <- lad_profile %>% dplyr::filter(height >= z0)
      lad_profile_lai <- leafR::lad.profile(VOXELS_LAD, relative = F) #doesn't use relative LAI
      LAI_subcanopy <- leafR::lai(lad_profile_lai, min = z0, max = max)
      out.plot <- c(LAI_subcanopy)
    } else {
      out.plot <- c(NA)
    }
    names(out.plot) <- c('LAI_subcanopy')
    return(out.plot)
  }
  
  ##### Additional Metrics Section #####
  modeBinary <- function (vxl, bk = 1) {
    Zs <- vxl$Z
    vxlHist <- graphics::hist(Zs, breaks = seq(0, max(Zs), by = bk), plot = T)
    countsOfConcern <- vxlHist$counts[2:length(vxlHist$counts)]
    midsOfConcern <- vxlHist$mids[2:length(vxlHist$mids)]
    modeIndex <- max(countsOfConcern) == countsOfConcern
    maxs <- midsOfConcern[modeIndex]
    modeSmall <- 0
    indexM <- 1
    if (!(length(maxs) == 0)) {
      while (indexM <= length(maxs) & modeSmall == 0) {
        if (maxs[[indexM]] <= 10) {
          modeSmall <- 1
        }
        indexM <- indexM + 1
      }
    }
    return(modeSmall)
  }
  
  
  # Peaks Location
  peaks <- function (vxl, bk = 1) {
    Zs <- vxl$Z
    # if(length(Zs) <= 0) return(NA)
    vxlHist <- graphics::hist(Zs, breaks = seq(0, max(Zs), by = bk), plot = F)
    peakx <- vxlHist$breaks[which(diff(sign(diff(vxlHist$density)))==-2)]
    peaksNumb <- length(peakx)
    if (peaksNumb == 0){
      peaksVar <- 0
    } else if (peaksNumb == 1) {
      peaksVar <- 0
    } else if (peaksNumb == 2) {
      peaksVar <- 0
    } else {
      peaksVar <- var(peakx)
    }
    return(peaksVar)
  }
  # Ratio of bottom points over top points
  canopy_class = function(Z) {
    # if(length(Z) <= 0) return(NA)
    hRatio =  sum(Z < 10)/length(Z) * 100
    return(hRatio)
  }
  # All together now
  additional_metrics <- function(data.30ft) {
    if(length(data.30ft$Z) >0 ){
      data.vxl <- voxelize_points(data.30ft, res = 1)
      # Mode of height after voxelizing
      modeLoc <- modeBinary(data.vxl, bk = 1)
      peaksVar <- peaks(data.vxl, bk = 1)
      hRatio <- canopy_class(data.vxl$Z)
      # Get data in easy data table with metric as column and value based on plotID as row
      additional_metrics_i <- c(modeLoc, peaksVar, hRatio)
    } else{
      additional_metrics_i <- c(NA, NA, NA)
    }
    
    names(additional_metrics_i) <- c("ModeLocation", "PeaksVarience", "HeightRatio")
    return(additional_metrics_i)
  }
  # Mode of height
  
}

################################################################################
###########################    FSD calculation     ############################# # 1~63, 64~126
################################################################################
las <- list.files("C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/PlotLiDAR/01.AllPlots_Original_Density/Final_LiDAR", pattern = ".las$", full.names = T)
las.id <-list.files("C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/PlotLiDAR/01.AllPlots_Original_Density/Final_LiDAR", pattern = ".las$", full.names = F)
las.id <- str_replace(las.id, ".las", "")

ctg <- readLAScatalog(las)
z0 <- 0.5
beer.lambert.constant <- 1
grain.size <- 3

OUT <- data.frame()
for (k in 1:length(las)) {
  las_file <- readLAS(las[k])
  las_file <- filter_poi(classify_noise(las_file, ivf(res=3,n=0)),
                         Classification %in% c(1, 2, 3, 4, 5)) %>%
    filter_duplicates() %>%
    # filter_poi(Z >= 0.5 & Z < mean(las_file@data$Z)+ 6*sd(las_file@data$Z)) %>%
    filter_poi(Z >= 0) %>%
    filter_poi(Z < mean(las_file@data$Z)+ 6*sd(las_file@data$Z))
  
  if(any(las_file$Z>z0)) {
    out.i <- las.id[k]
    
    las_file2 <- las_file %>% filter_poi(Z >= z0)
    vertq <- quantile(las_file2$Z, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE)  # <- corrected
    names(vertq) <- c('q0', 'q1', 'q5', 'q25', 'q50', 'q75', 'q95', 'q99', 'q100')
    plot.area <- area(las_file2)
    den <- length(las_file2@data$Z)/plot.area
    
    chm.metrics <- chm_metrics(las_file, z0)
    #if there are no points b/c all veg is < 0.5 m keep loop running
    FSD.i <- structural_diversity_metrics(las_file, z0)
    LAI.lidR <- LAI_lidR(las_file, beer.lambert.constant, z0)
    
    out.z <- c(plot.area = plot.area, density = den, vertq, chm.metrics, FSD.i, LAI.lidR)
    
    LAI.leafR <- LAI_leafR(las_file, grain.size, z0, beer.lambert.constant)
    
    LAIsub.leafR <- LAIsubcanopy_leafR(las_file, grain.size = grain.size, z0 = z0, max = 5, k = 1)
    add_met <- additional_metrics(las_file)
    out.z <- c(out.z, LAI.leafR, LAIsub.leafR, add_met)
    out.i <- c(out.i, out.z)
    out.i <- as.data.frame(t(out.i))
  }
  colnames(out.i)[1] <- "plotID"
  OUT <- rbind(OUT, out.i)
} #215 error


write.csv(OUT, "X:/01.PROJECT/02.UrbanStudy/02.UrbanInvasive/0000000.Final_all/LiDAR_data/Final_LiDAR/FSD_plots.csv", row.names = F)

# head(OUT)
# colnames(OUT)[1:2] <- c("FileName", "gbifID")
# getwd()
# print(k)




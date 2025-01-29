################################################################################
#####################          Modelling             ###########################
################################################################################
###### Table Organization


###### Table Organization
################################################################################
#####################          Modelling             ###########################
################################################################################
###### Table Organization
library(sf)


plots <- st_read("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/PlotPolygon/plot_points.shp")
gbif_plots <- st_read("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/PlotPolygon/gbif_2017_2023_ISS_Forest_March_2024.shp")


plots <- plots[!grepl("D", plots$PlotID), ]

colnames(gbif_plots)[1] <- "PlotID"
gbif_plots$invaded <- 1
gbif_plots <- gbif_plots %>% dplyr::select(PlotID, invaded, geometry)


plots <- rbind(plots, gbif_plots)
plots

plots <- st_transform(plots, "EPSG:6455")
# mapview::mapview(gbif_selected.prj.NADft)
bufferedData <- st_buffer(plots, 60)
bufferedData_30ft <- st_buffer(plots, 30)


# st_write(bufferedData_30ft, "C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/Final_Publish/Last_plots.shp")
bufferedData_30ft_union <- st_union(bufferedData_30ft) %>% st_cast("POLYGON")


# bufferedData_30ft_union$ID <- 1:length(bufferedData_30ft_union)[1]

overalapped_polygons <- st_intersects(bufferedData_30ft_union, bufferedData_30ft)

which(unlist(lapply(overalapped_polygons, length)) != 1)
overalapped_polygons[which(unlist(lapply(overalapped_polygons, length)) != 1)]

a <- unlist(overalapped_polygons[which(unlist(lapply(overalapped_polygons, length)) == 1)])
b <- sapply(overalapped_polygons[which(unlist(lapply(overalapped_polygons, length)) != 1)],"[[",2)

PlotID <- bufferedData_30ft[c(a,b), ]$PlotID


###### Table Organization
library(dplyr)
library(sf)
# setwd("S:/Labs/HardimanLab/Isaac/Playing with Buckthorn Data/Playing with Buckthorn Data")
# setwd("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/Output")
setwd("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/Final_Publish")

# plotTable <- read.csv("metricsJoined2.csv")
# species_BA_list <- read.csv("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/Output/Lindsay_Table/plotBAForBuckthorn1.csv")
species_BA_list <- read.csv("C:/Users/star1/Documents/GitHub/Chicago-Invasives-Analysis/Output/Lindsay_Table/plotBAForBuckthorn1.csv")

library(ggplot2)
FSD <- read.csv("./FSD_last_March_02.csv")
NDVI <- read.csv("./plots_Chicago_NAIP2017_NDVI_leaf_on_30ft_last.csv")

colnames(FSD)[1] <- "PlotID"
FSD_NDVI <- merge(NDVI, FSD, by = "PlotID")
colnames(FSD_NDVI)[3] <- "value"

FSD_NDVI <- FSD_NDVI[which(FSD_NDVI$PlotID %in% PlotID), ]
sum(grepl("A", FSD_NDVI$PlotID))
sum(grepl("B", FSD_NDVI$PlotID))
sum(grepl("C", FSD_NDVI$PlotID))
sum(grepl("D", FSD_NDVI$PlotID))
sum(grepl("E", FSD_NDVI$PlotID))




FSD_NDVI <- FSD_NDVI[FSD_NDVI$plot.area > 260 * 0.7, ]
dim(FSD_NDVI)

# Final Plots
Final_plots <- FSD_NDVI$PlotID

bufferedData_30ft <- bufferedData_30ft[which(bufferedData_30ft$PlotID %in% Final_plots), ]
sum(grepl("A", bufferedData_30ft$PlotID))
sum(grepl("B", bufferedData_30ft$PlotID))
sum(grepl("C", bufferedData_30ft$PlotID))
sum(grepl("E", bufferedData_30ft$PlotID))
# 443 - 136-36-55
# st_write(bufferedData_30ft, "C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/Final_Publish/Plots.shp")
# FSD_NDVI <- FSD_NDVI[which(FSD_NDVI$PlotID %in% PlotID), ]
# dim(FSD_NDVI)
# FSD_NDVI <- FSD_NDVI[,c(1, 3:38)]
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(ranger)
library(dplyr)
library(ggplot2)
library(Boruta)
library(caret)
library(randomForest)
# FinalTable_filtered <- read.csv("S:/Labs/HardimanLAb/Isaac/05.StatTables/FinalTable.csv")
# FinalTable_filtered <- read.csv("Y:/Isaac/05.StatTables/FinalTable.csv")

# FinalTable_filtered$ModeIn10meter <- as.factor(FinalTable_filtered$ModeIn10meter)
# FinalTable_filtered$Perc_Shrub <- as.numeric(FinalTable_filtered$Perc_Shrub)
# FinalTable_filtered_scaled <- FinalTable_filtered |>
#   dplyr::select(!c("PlotID", "Rhamnus_ca", "Perc_Rhamn", "Perc_Shrub", "BAShrub", "shrubby")) |>
#   dplyr::mutate(across(where(is.numeric), datawizard::standardize))

################################################################################
# xgboost

#only lidar
FinalTable_filtered <- FSD_NDVI

aucc <- c()
model_list <- vector("list")
train_data <- c()
test_data <- c()
predicted <- c()
sss <- c()
# library(ROCR)
# library(plotROC)
library(pROC)

set.seed(4342)
rn <- sample(1:10000,9999)

cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$VCI)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$LAI)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$value)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$LAI_subcanopy)
cor(FinalTable_filtered$LAI, FinalTable_filtered$VCI)
cor(FinalTable_filtered$LAI, FinalTable_filtered$LAI_subcanopy)
cor(FinalTable_filtered$LAI, FinalTable_filtered$deepgap.fraction)

for (iter in 1:9999) {
  set.seed(rn[iter])
  ix <- sample(nrow(FinalTable_filtered), 0.7 * nrow(FinalTable_filtered))
  train_data <- FinalTable_filtered[ix,]
  test_data <- FinalTable_filtered[-ix,]
  gls.binominal <- glm(invaded ~ 
                         scale(density) + 
                         scale(mean.max.canopy.ht) +
                         scale(LAI_subcanopy) + 
                         scale(LAI) +
                         scale(value)+ 
                         scale(VCI)
                       ,
                       family = "binomial",
                       data = train_data)
  model_list <- gls.binominal
  # acc.train <- c(acc.train, mean(as.numeric(gls.binominal$fitted.values > 0.85) == train_data$invaded))
  predicted <- predict(gls.binominal, test_data, type="response")
  aucc <- c(aucc, auc(test_data$invaded, predicted))
}

mean(aucc)
max(aucc)

i <- which(aucc == max(aucc))
# i <- which.min(abs(aucc - mean(aucc)))

set.seed(rn[i])
ix <- sample(nrow(FinalTable_filtered), 0.7 * nrow(FinalTable_filtered))
train_data <- FinalTable_filtered[ix,]
test_data <- FinalTable_filtered[-ix,]

gls.binominal <- glm(invaded ~ 
                       scale(density) + 
                       scale(mean.max.canopy.ht) +
                       scale(LAI_subcanopy) + 
                       scale(LAI) +
                       scale(value)+ 
                       scale(VCI),
                       # scale(VCI):scale(LAI_subcanopy),
                     family = "binomial",
                     data = train_data)
predicted <- predict(gls.binominal, test_data, type="response")
auc(test_data$invaded, predicted)
performance::check_model(gls.binominal)


df_test <- data.frame(predicted, test_data$invaded)
colnames(df_test) <- c("ypred", "true")

library(ROCR)
library(plotROC)
basicplot <- ggplot(df_test, aes(d = true, m = ypred)) + geom_roc() 
basicplot + style_roc(theme = theme_grey) +
  theme(axis.text = element_text(colour = "blue")) +
  ggtitle("Invasive Shrub Species - Occurrence") + 
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))

hist(aucc)

threshold <- seq(0.2, 0.99, 0.001) # which is 1 or 0. 
sensitivity_th <- c()
specificity_th <- c()
sss_s <- c()

for (th in threshold) {
  predicted_values<-ifelse(predict(gls.binominal, test_data, type="response")>th, 1, 0)
  actual_values<-test_data$invaded
  conf_matrix<-table(predicted_values,actual_values)
  
  sensitivity_th <- c(sensitivity_th, sensitivity(conf_matrix))
  specificity_th <- c(specificity_th, specificity(conf_matrix))
  sss_s <- sensitivity_th + specificity_th
}

df_sss <- data.frame(threshold, sensitivity_th, specificity_th, sss_s)
df_sss_melted <- reshape2::melt(df_sss, id = "threshold")
ggplot(df_sss_melted, aes(threshold, value, col = variable)) + geom_line()

threshold_last <- df_sss[df_sss$sss_s == max(df_sss), ]$threshold
threshold_last

acc.train <- c()
for (th in threshold_last) {
  acc.train <- c(acc.train, mean(ifelse(predict(gls.binominal, test_data, type="response")>th, 1, 0)== test_data$invaded))
}

pred_inva <- as.vector(ifelse(predict(gls.binominal, test_data, type="response")>th, 1, 0))
test_data$invaded

library("MLmetrics")
F1_Score(pred_inva,test_data$invaded)
library("caret")

pred_inva <- as.factor(pred_inva)
test_data$invaded <- as.factor(test_data$invaded)
confusionMatrix(pred_inva, test_data$invaded,
                mode = "everything",
                positive="1")


library(MoMAColors)
col2 <- moma.colors(
  "Andri",
  2,
  type = c("discrete"),
  # direction = c(1),
)


Effect_plots <- sjPlot::plot_model(gls.binominal,
                   show.values = T,
                   type = "std2",
                   transform = NULL,
                   # axis.label = "Buffer Radius (m)",
                   title = "Presence of Invasive Shrub Species",
                   sort.est = T,
                   value.offset = 0.3,
                   value.size = 7) + theme_bw() + geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        text = element_text(size = 20)) + 
  scale_color_manual(values = c(col2[1], col2[2]))


Effect_plots$data <- Effect_plots$data[-1, ]
Effect_plots

################################################################################
# model compariosns

# wo NDVI
gls.binominal.wo_NDVI <- glm(invaded ~ 
                               scale(density) + 
                               scale(mean.max.canopy.ht) +
                               scale(LAI_subcanopy) +
                               scale(LAI) +
                               scale(VCI),
                               # scale(value),
                     family = "binomial",
                     data = train_data)
predicted <- predict(gls.binominal.wo_NDVI, test_data, type="response")
auc(test_data$invaded, predicted)

# wo Height
gls.binominal.wo_height <- glm(invaded ~ 
                                 scale(density) + 
                                 # scale(mean.max.canopy.ht) +
                                 scale(LAI_subcanopy) +
                                 scale(LAI) +
                                 scale(VCI) +
                                 scale(value),
                     # scale(VCI):scale(LAI_subcanopy),
                     family = "binomial",
                     data = train_data)
predicted <- predict(gls.binominal.wo_NDVI, test_data, type="response")
auc(test_data$invaded, predicted)

# wo LAI_sub
gls.binominal.wo_LAIsub <- glm(invaded ~ 
                                 scale(density) + 
                                 scale(mean.max.canopy.ht) +
                                 # scale(LAI_subcanopy) + 
                                 scale(LAI) +
                                 scale(VCI) +
                                 scale(value),
                             # scale(VCI):scale(LAI_subcanopy),
                             family = "binomial",
                             data = train_data)
predicted <- predict(gls.binominal.wo_LAIsub, test_data, type="response")
auc(test_data$invaded, predicted)

# wo LiDAR
gls.binominal.wo_LiDAR <- glm(invaded ~ 
                                # scale(density) + 
                                # scale(mean.max.canopy.ht) +
                                # scale(LAI_subcanopy) + 
                                # scale(LAI) +
                                # scale(VCI) + 
                                scale(value),
                               # scale(VCI):scale(LAI_subcanopy),
                               family = "binomial",
                               data = train_data)
predicted <- predict(gls.binominal.wo_LiDAR, test_data, type="response")
auc(test_data$invaded, predicted)

# w height, ndvi
gls.binominal.w_h_nd <- glm(invaded ~ 
                              scale(density) + 
                                scale(mean.max.canopy.ht) +
                                # scale(LAI_subcanopy) + 
                                # scale(LAI) +
                                # scale(VCI) + 
                                scale(value),
                              # scale(VCI):scale(LAI_subcanopy),
                              family = "binomial",
                              data = train_data)
predicted <- predict(gls.binominal.wo_LiDAR, test_data, type="response")
auc(test_data$invaded, predicted)


# wo density
gls.binominal.wo.density <- glm(invaded ~ 
                              # scale(density) + 
                              scale(mean.max.canopy.ht) +
                              scale(LAI_subcanopy) +
                              scale(LAI) +
                              scale(VCI) +
                              scale(value),
                            # scale(VCI):scale(LAI_subcanopy),
                            family = "binomial",
                            data = train_data)
predicted <- predict(gls.binominal.wo.density, test_data, type="response")
auc(test_data$invaded, predicted)



gls.binominal.w.density_ndvi <- glm(invaded ~ 
                                  scale(density) +
                                  # scale(mean.max.canopy.ht) +
                                  # scale(LAI_subcanopy) +
                                  # scale(LAI) +
                                  # scale(VCI) +
                                  scale(value),
                                # scale(VCI):scale(LAI_subcanopy),
                                family = "binomial",
                                data = train_data)
predicted <- predict(gls.binominal.wo.density, test_data, type="response")
auc(test_data$invaded, predicted)

sjPlot::tab_model(gls.binominal, gls.binominal.wo.density, gls.binominal.wo_NDVI, 
                  gls.binominal.wo_height, gls.binominal.wo_LAIsub, 
                  gls.binominal.wo_LiDAR, gls.binominal.w_h_nd,
                  gls.binominal.w.density_ndvi,
                  show.aic = T,
                  show.aicc = T)


acc.train


######## Correlogram 
library(GGally)

library(MoMAColors)
col2 <- moma.colors(
  "Andri",
  2,
  type = c("discrete"), return_hex = T
  # direction = c(1),
)
# FinalTable_filtered[,c(21, 31, 33, 34, 4)]

pm <- ggpairs(FinalTable_filtered, columns = c(5, 17, 27, 29, 30, 3), 
        ggplot2::aes(colour=factor(invaded), alpha = 0.2),
        columnLabels = c("Density", "MoM_CH", "VCI", "VAI", "VAI_subcanpoy", "NDVI"),
        upper = list(continuous = wrap("cor", size = 5)),
        diag = list(continuous = wrap("box_no_facet", alpha=0.5)), 
        lower = list(continuous = "smooth"))  + 
  scale_colour_manual(values = c("#15134b", "#f56455")) + 
  scale_fill_manual(values = c("#15134b", "#f56455")) + theme_bw() + 
  theme(strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 13))
pm
# p_ <- GGally::print_if_interactive
# p_(pm)

cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$VCI)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$LAI)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$value)
cor(FinalTable_filtered$mean.max.canopy.ht, FinalTable_filtered$LAI_subcanopy)
cor(FinalTable_filtered$LAI, FinalTable_filtered$VCI)
cor(FinalTable_filtered$LAI, FinalTable_filtered$LAI_subcanopy)
cor(FinalTable_filtered$LAI, FinalTable_filtered$deepgap.fraction)


################################################################################
# overall <- st_read("C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/Output/Total_last0222_final_old_new_Forest_orginal_density.shp")
# colnames(overall)[c(15, 28)] <- c("mean.max.canopy.ht", "LAI_subcanopy")
# overall$pred <- predict(gls.binominal, overall, type="response")
# 
# bbox <- st_bbox(overall)
# r <- terra::rast(res = 10, xmin = bbox[1], xmax = bbox[3], ymin = bbox[2], ymax= bbox[4], crs = "EPSG:32616")
# r <- raster::raster(r)
# r <- fasterize::fasterize(overall, r, field = "pred", fun="last")
# 
# plot(r)
# getwd()
# terra::writeRaster(r, "C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/Final_Publish/PredMap_fin_0222_24_MoMheight_original_density_fin_0302_2024_07cover.tif")
################################################################################
# th 0.587

r_1 <- r > 0
plot(r_1)
terra::writeRaster(r_1, "C:/Users/dhchoi/Documents/GitHub/Chicago-Invasives-Analysis/Output/forest_binary.tif")


################################################################################
# increase ÷ original number × 100

library(ggstatsplot)
ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  LAI_subcanopy,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,
  pairwise.comparisons = FALSE, 
  pairwise.display = "all",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg",
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1])) + 
  theme(text = element_text(size=20))

(0.65-0.16)/0.16 * 100 #306.25


ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  LAI,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg",
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1]))+
  theme(text = element_text(size=20))

(1.66-1.36)/1.36 * 100 #22

ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  VCI,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg", 
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1]))+
  theme(text = element_text(size=20))


ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  value,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg",
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1]))+
  theme(text = element_text(size=20))
(0.28-0.39)/0.39 * 100 #-28



ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  mean.max.canopy.ht,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg",
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1]))+
  theme(text = element_text(size=20))


ggbetweenstats(
  data = FinalTable_filtered,
  x = invaded,
  y =  density,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg",
  centrality.label.args = list(size = 8, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)) + ggplot2::scale_color_manual(values = c(col2[2], col2[1]))+
  theme(text = element_text(size=20))
##########

ggbetweenstats(
  data = dplyr::filter(grid, Year == 2017, invaded == 0),
  x = pointgroup,
  y = value,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)


ggbetweenstats(
  data = dplyr::filter(grid, Year == 2017, invaded == 0),
  x = pointgroup,
  y = value,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)

ggbetweenstats(
  data = dplyr::filter(grid, Year == 2017, pointgroup == "B"),
  x = invaded,
  y = value,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)


ggbetweenstats(
  data = dplyr::filter(grid, Year == 2017),
  x = invaded,
  y = value,
  violin.args = list(width = 0),
  type = "p",
  conf.level = 0.95,pairwise.comparisons = TRUE, pairwise.display = "significant",
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)


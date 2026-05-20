library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(reshape2)
library(RColorBrewer)
#scenerio0
#Trial1

sdwls_markers_stlsq <- read_csv("scenerio0-default/scenerio0_stlsq_trial1_2_positivew.dat")
spotlight_scenerio0 <- read_csv("scenerio0-default/scenerio0_spotlight.dat")
norm_weights_scenerio0 <- read.csv("scenerio0-default/scenerio0_rctd_results.dat")
norm_weights_scenerio0 <- dplyr::select(norm_weights_scenerio0, -1)
giotto_scenerio0 <- readRDS('scenerio0-default/scenerio0_spot_giotto_Results.rds')
giotto_scenerio0_2<- giotto_scenerio0@spatial_enrichment$DWLS
giotto_scenerio0_2 <- dplyr::select(giotto_scenerio0_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio0_2))))
giotto_scenerio0_2 <- as.matrix(giotto_scenerio0_2)
giotto_scenerio0_2_rear <- giotto_scenerio0_2[, s]

scenerio0_dwls <- readRDS('scenerio0-default/scenerio0_DWLS_Results.rds')
scenerio0_dwls <- t(scenerio0_dwls)
ss = order(unique(as.numeric(colnames(scenerio0_dwls))))

cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14) #rearrange columns based on the column order given in ground truth data
scenerio0_dwls_rear <- as.matrix(scenerio0_dwls[, cols])
head(scenerio0_dwls_rear)
head(scenerio0_dwls)


sdwls_markers_stlsq <- as.matrix(sdwls_markers_stlsq)
giotto_scenerio0_2_rear2 <- as.matrix(giotto_scenerio0_2_rear)

s = order(unique(as.numeric(colnames(spotlight_scenerio0))))

spotlight_scenerio0 <- spotlight_scenerio0[,s]

spotlight_scenerio0 <- as.matrix(spotlight_scenerio0)

sdwls_markers_stlsq_prop <- sdwls_markers_stlsq / rowSums(sdwls_markers_stlsq)
sdwls_markers_stlsq_prop[is.na(sdwls_markers_stlsq_prop)] = 0

stereoscope01 <- read.csv("scenerio0-default/stereoscope_results_trial1/counts.st-heart_with_fraction_0.1/W.2023-08-02180430.699341.tsv", sep = "\t")
stereoscope01 <- dplyr::select(stereoscope01, -1)

dim(sdwls_markers_stlsq); dim(norm_weights_scenerio0); dim(giotto_scenerio0_2_rear2); dim(scenerio0_dwls_rear); dim(spotlight_scenerio0); dim(stereoscope01)


ground_tr <- read.csv("scenerio0-default/proportions.st-heart_with_fraction_0.1.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(sdwls_markers_stlsq_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls_rear))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio0))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio0_2_rear2))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio0))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope01))


sdwls_markers_stlsq_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

#Trial2
sdwls_markers_stlsq2 <- read_csv("scenerio0-default/scenerio0_stlsq_trial2_2_positivew.dat")
spotlight_scenerio2 <- read_csv("scenerio0-default/scenerio0_spotlight_trial2.dat")
norm_weights_scenerio2 <- read.csv("scenerio0-default/scenerio0_rctd_results_trial2.dat")
norm_weights_scenerio2 <- dplyr::select(norm_weights_scenerio2, -1)
giotto_scenerio2 <- readRDS('scenerio0-default/scenerio0_spot_giotto_Results_trial2.rds')
giotto_scenerio2_2<- giotto_scenerio2@spatial_enrichment$DWLS
giotto_scenerio2_2 <- dplyr::select(giotto_scenerio2_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio2_2))))
giotto_scenerio2_2 <- as.matrix(giotto_scenerio2_2)
giotto_scenerio2_2_rear <- giotto_scenerio2_2[, s]
scenerio2_dwls <- readRDS('scenerio0-default/scenerio0_DWLS_Results_trial2.rds')
scenerio2_dwls <- t(scenerio2_dwls)
cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14)
scenerio2_dwls_rear <- as.matrix(scenerio2_dwls[, cols])
head(scenerio2_dwls_rear)
sdwls_markers_stlsq2 <- as.matrix(sdwls_markers_stlsq2)
giotto_scenerio2_2_rear2 <- as.matrix(giotto_scenerio2_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio2))))

spotlight_scenerio2 <- spotlight_scenerio2[,s]
spotlight_scenerio2 <- as.matrix(spotlight_scenerio2)

sdwls_markers_stlsq2_prop <- sdwls_markers_stlsq2 / rowSums(sdwls_markers_stlsq2)
sdwls_markers_stlsq2_prop[is.na(sdwls_markers_stlsq2_prop)] = 0

stereoscope02 <- read.csv("scenerio0-default/stereoscope_results_trial2/counts_trial2.st-heart_with_fraction_0.1_trial2_scenerio0_trial2/W.2023-08-03104112.827834.tsv", sep = "\t")
stereoscope02 <- dplyr::select(stereoscope02, -1)

dim(sdwls_markers_stlsq2); dim(norm_weights_scenerio2); dim(giotto_scenerio2_2_rear2); dim(scenerio2_dwls_rear); dim(spotlight_scenerio2); dim(stereoscope02)


ground_tr2 <- read.csv("scenerio0-default/proportions_trial2.st-heart_with_fraction_0.1_trial2_scenerio0_trial2.tsv", sep = "\t")
ground_tr2 <- dplyr::select(ground_tr2, -1)
ground_tr2 <- as.matrix(ground_tr2)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(sdwls_markers_stlsq2_prop))
dwls_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(scenerio2_dwls_rear))
rctd_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(norm_weights_scenerio2))
giotto_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(giotto_scenerio2_2_rear2))
spotlight_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(spotlight_scenerio2))
stereoscope_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(stereoscope02))


sdwls_markers_stlsq2_deneme$Name <- c("WISpR")
dwls_rmse2_deneme$Name <- c("DWLS")
rctd_rmse2_deneme$Name <- c("RCTD")
giotto_rmse2_deneme$Name <- c("S-DWLS")
spotlight_rmse2_deneme$Name <- c("SPOTlight")
stereoscope_rmse2_deneme$Name <- c("Stereoscope")

#Trial3
sdwls_markers_stlsq3 <- read_csv("scenerio0-default/scenerio0_stlsq_trial3_2_positivew.dat")
spotlight_scenerio3 <- read_csv("scenerio0-default/scenerio0_spotlight_trial3.dat")
norm_weights_scenerio3 <- read.csv("scenerio0-default/scenerio0_rctd_results_trial3.dat")
norm_weights_scenerio3 <- dplyr::select(norm_weights_scenerio3, -1)
giotto_scenerio3 <- readRDS('scenerio0-default/scenerio0_spot_giotto_Results_trial3.rds')
giotto_scenerio3_2<- giotto_scenerio3@spatial_enrichment$DWLS
giotto_scenerio3_2 <- dplyr::select(giotto_scenerio3_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio3_2))))
giotto_scenerio3_2 <- as.matrix(giotto_scenerio3_2)
giotto_scenerio3_2_rear <- giotto_scenerio3_2[, s]
scenerio3_dwls <- readRDS('scenerio0-default/scenerio0_DWLS_Results_trial3.rds')
scenerio3_dwls <- t(scenerio3_dwls)
#cols <- c(12,4,3,6,11,13,14,1,5,7,10,2,8,15,9)
cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14)
scenerio3_dwls_rear <- as.matrix(scenerio3_dwls[, cols])
sdwls_markers_stlsq3 <- as.matrix(sdwls_markers_stlsq3)
giotto_scenerio3_2_rear2 <- as.matrix(giotto_scenerio3_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio3))))

spotlight_scenerio3 <- spotlight_scenerio3[,s]
spotlight_scenerio3 <- as.matrix(spotlight_scenerio3)

sdwls_markers_stlsq3_prop <- sdwls_markers_stlsq3 / rowSums(sdwls_markers_stlsq3)
sdwls_markers_stlsq3_prop[is.na(sdwls_markers_stlsq3_prop)] = 0

stereoscope03 <- read.csv("scenerio0-default/stereoscope_results_trial3/counts_trial3.st-heart_with_fraction_0.1_trial3_scenerio0_trial3/W.2023-08-03232054.842727.tsv", sep = "\t")
stereoscope03 <- dplyr::select(stereoscope03, -1)
dim(sdwls_markers_stlsq3); dim(norm_weights_scenerio3); dim(giotto_scenerio3_2_rear2); dim(scenerio3_dwls_rear); dim(spotlight_scenerio3); dim(stereoscope03)


ground_tr3 <- read.csv("scenerio0-default/proportions_trial3.st-heart_with_fraction_0.1_trial3_scenerio0_trial3.tsv", sep = "\t")
ground_tr3 <- dplyr::select(ground_tr3, -1)
ground_tr3 <- as.matrix(ground_tr3)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(sdwls_markers_stlsq3_prop))
dwls_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(scenerio3_dwls_rear))
rctd_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(norm_weights_scenerio3))
giotto_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(giotto_scenerio3_2_rear2))
spotlight_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(spotlight_scenerio3))
stereoscope_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(stereoscope03))

sdwls_markers_stlsq3_deneme$Name <- c("WISpR")
dwls_rmse3_deneme$Name <- c("DWLS")
rctd_rmse3_deneme$Name <- c("RCTD")
giotto_rmse3_deneme$Name <- c("S-DWLS")
spotlight_rmse3_deneme$Name <- c("SPOTlight")
stereoscope_rmse3_deneme$Name <- c("Stereoscope")


all_merged3_scenerio0 <- rbind(sdwls_markers_stlsq_deneme, sdwls_markers_stlsq2_deneme, 
                               sdwls_markers_stlsq3_deneme, dwls_rmse_deneme, dwls_rmse2_deneme, 
                               dwls_rmse3_deneme, rctd_rmse_deneme, rctd_rmse2_deneme, 
                               rctd_rmse3_deneme, giotto_rmse_deneme, giotto_rmse2_deneme, 
                               giotto_rmse3_deneme, spotlight_rmse_deneme, spotlight_rmse2_deneme, 
                               spotlight_rmse3_deneme, stereoscope_rmse_deneme, stereoscope_rmse2_deneme, 
                               stereoscope_rmse3_deneme)

all_merged3_scenerio0$Version <- c("SCENERIO_0")
library(plyr)
df <- ddply(all_merged3_scenerio0, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3_scenerio0
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses


colnames(all_merged3_scenerio0) <- c("RMSE", "Name")
dim(all_merged3_scenerio0)

#Scenerio1
#Trial1

sdwls_markers_stlsq <- read_csv("scenerio1_visium/scenerio1_stlsq_trial1_filtered_positivew.dat")
spotlight_scenerio0 <- read_csv("scenerio1_visium/scenerio1_spotlight.dat")
norm_weights_scenerio0 <- read.csv("scenerio1_visium/scenerio1_rctd_results.dat")
norm_weights_scenerio0 <- dplyr::select(norm_weights_scenerio0, -1)
giotto_scenerio0 <- readRDS('scenerio1_visium/scenerio1_spot_giotto_Results.rds')
giotto_scenerio0_2<- giotto_scenerio0@spatial_enrichment$DWLS
giotto_scenerio0_2 <- dplyr::select(giotto_scenerio0_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio0_2))))
giotto_scenerio0_2 <- as.matrix(giotto_scenerio0_2)
giotto_scenerio0_2_rear <- giotto_scenerio0_2[, s]
scenerio0_dwls <- readRDS('scenerio1_visium/scenerio1_DWLS_Results.rds')
scenerio0_dwls <- t(scenerio0_dwls)

cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14)
scenerio0_dwls_rear <- as.matrix(scenerio0_dwls[, cols])

sdwls_markers_stlsq <- as.matrix(sdwls_markers_stlsq)
giotto_scenerio0_2_rear2 <- as.matrix(giotto_scenerio0_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio0))))
spotlight_scenerio0 <- spotlight_scenerio0[,s]

spotlight_scenerio0 <- as.matrix(spotlight_scenerio0)
dim(sdwls_markers_stlsq); dim(norm_weights_scenerio0); dim(giotto_scenerio0_2_rear2); dim(scenerio0_dwls_rear); dim(spotlight_scenerio0)

sdwls_markers_stlsq_prop <- sdwls_markers_stlsq / rowSums(sdwls_markers_stlsq)
sdwls_markers_stlsq_prop[is.na(sdwls_markers_stlsq_prop)] = 0
stereoscope11 <- read.csv("scenerio1_visium/stereoscope_results1/countsvisium.st-heart_with_fraction_0.1visium/W.2023-08-04130737.016781.tsv", sep = "\t")
stereoscope11 <- dplyr::select(stereoscope11, -1)

ground_tr <- read.csv("scenerio1_visium/proportionsvisium.st-heart_with_fraction_0.1visium.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

sdwls_markers_stlsq_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(sdwls_markers_stlsq_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls_rear))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio0))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio0_2_rear2))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio0))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope11))

sdwls_markers_stlsq_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

#Trial2
sdwls_markers_stlsq2 <- read_csv("scenerio1_visium/scenerio1_stlsq_trial2_filtered_positivew.dat")
spotlight_scenerio2 <- read_csv("scenerio1_visium/scenerio1_spotlight_trial2.dat")
norm_weights_scenerio2 <- read.csv("scenerio1_visium/scenerio1_rctd_results_trial2.dat")
norm_weights_scenerio2 <- dplyr::select(norm_weights_scenerio2, -1)
giotto_scenerio2 <- readRDS('scenerio1_visium/scenerio1_spot_giotto_Results_trial2.rds')
giotto_scenerio2_2<- giotto_scenerio2@spatial_enrichment$DWLS
giotto_scenerio2_2 <- dplyr::select(giotto_scenerio2_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio2_2))))
giotto_scenerio2_2 <- as.matrix(giotto_scenerio2_2)
giotto_scenerio2_2_rear <- giotto_scenerio2_2[, s]
scenerio2_dwls <- readRDS('scenerio1_visium/scenerio1_DWLS_Results_trial2.rds')
scenerio2_dwls <- t(scenerio2_dwls)
cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14)
scenerio2_dwls_rear <- as.matrix(scenerio2_dwls[, cols])

sdwls_markers_stlsq2 <- as.matrix(sdwls_markers_stlsq2)
giotto_scenerio2_2_rear2 <- as.matrix(giotto_scenerio2_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio2))))

spotlight_scenerio2 <- spotlight_scenerio2[,s]

spotlight_scenerio2 <- as.matrix(spotlight_scenerio2)
dim(sdwls_markers_stlsq2); dim(norm_weights_scenerio2); dim(giotto_scenerio2_2_rear2); dim(scenerio2_dwls_rear); dim(spotlight_scenerio2)

sdwls_markers_stlsq2_prop <- sdwls_markers_stlsq2 / rowSums(sdwls_markers_stlsq2)
sdwls_markers_stlsq2_prop[is.na(sdwls_markers_stlsq2_prop)] = 0

stereoscope12 <- read.csv("scenerio1_visium/stereoscope_results2/countstrial2_visium.st-heart_with_fraction_0.1_trial2_scenerio1trial2_visium/W.2023-08-05073948.076123.tsv", sep = "\t")
stereoscope12 <- dplyr::select(stereoscope12, -1)
head(stereoscope12)

ground_tr2 <- read.csv("scenerio1_visium/proportionstrial2_visium.st-heart_with_fraction_0.1_trial2_scenerio1trial2_visium.tsv", sep = "\t")
ground_tr2 <- dplyr::select(ground_tr2, -1)
ground_tr2 <- as.matrix(ground_tr2)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm2 <- ground_tr2/rowSums(ground_tr2)

sdwls_markers_stlsq2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(sdwls_markers_stlsq2_prop))
dwls_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(scenerio2_dwls_rear))
rctd_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(norm_weights_scenerio2))
giotto_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(giotto_scenerio2_2_rear2))
spotlight_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(spotlight_scenerio2))
stereoscope_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope12))

sdwls_markers_stlsq2_deneme$Name <- c("WISpR")
dwls_rmse2_deneme$Name <- c("DWLS")
rctd_rmse2_deneme$Name <- c("RCTD")
giotto_rmse2_deneme$Name <- c("S-DWLS")
spotlight_rmse2_deneme$Name <- c("SPOTlight")
stereoscope_rmse2_deneme$Name <- c("Stereoscope")

#Trial3
sdwls_markers_stlsq3 <- read_csv("scenerio1_visium/scenerio1_stlsq_trial3_filtered_positivew.dat")
spotlight_scenerio3 <- read_csv("scenerio1_visium/scenerio1_spotlight_trial3.dat")
norm_weights_scenerio3 <- read.csv("scenerio1_visium/scenerio1_rctd_results_trial3.dat")
norm_weights_scenerio3 <- dplyr::select(norm_weights_scenerio3, -1)
giotto_scenerio3 <- readRDS('scenerio1_visium/scenerio1_spot_giotto_Results_trial3.rds')
giotto_scenerio3_2<- giotto_scenerio3@spatial_enrichment$DWLS
giotto_scenerio3_2 <- dplyr::select(giotto_scenerio3_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio3_2))))
giotto_scenerio3_2 <- as.matrix(giotto_scenerio3_2)
giotto_scenerio3_2_rear <- giotto_scenerio3_2[, s]
scenerio3_dwls <- readRDS('scenerio1_visium/scenerio1_DWLS_Results_trial3.rds')
scenerio3_dwls <- t(scenerio3_dwls)
cols <- c(8,12,3,2,9,4,10,13,15,11,5,1,6,7,14)
scenerio3_dwls_rear <- as.matrix(scenerio3_dwls[, cols])
sdwls_markers_stlsq3 <- as.matrix(sdwls_markers_stlsq3)
giotto_scenerio3_2_rear2 <- as.matrix(giotto_scenerio3_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio3))))
spotlight_scenerio3 <- spotlight_scenerio3[,s]
spotlight_scenerio3 <- as.matrix(spotlight_scenerio3)
dim(sdwls_markers_stlsq3); dim(norm_weights_scenerio3); dim(giotto_scenerio3_2_rear2); dim(scenerio3_dwls_rear); dim(spotlight_scenerio3)

sdwls_markers_stlsq3_prop <- sdwls_markers_stlsq3 / rowSums(sdwls_markers_stlsq3)
sdwls_markers_stlsq3_prop[is.na(sdwls_markers_stlsq3_prop)] = 0

stereoscope13 <- read.csv("scenerio1_visium/stereoscope_results3/countstrial3_visium.st-heart_with_fraction_0.1_trial2_scenerio1trial3_visium/W.2023-08-06074025.620135.tsv", sep = "\t")
stereoscope13 <- dplyr::select(stereoscope13, -1)

ground_tr3 <- read.csv("scenerio1_visium/proportionstrial3_visium.st-heart_with_fraction_0.1_trial2_scenerio1trial3_visium.tsv", sep = "\t")
ground_tr3 <- dplyr::select(ground_tr3, -1)
ground_tr3 <- as.matrix(ground_tr3)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm3 <- ground_tr3/rowSums(ground_tr3)

sdwls_markers_stlsq3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(sdwls_markers_stlsq3_prop))
dwls_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(scenerio3_dwls_rear))
rctd_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(norm_weights_scenerio3))
giotto_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(giotto_scenerio3_2_rear2))
spotlight_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(spotlight_scenerio3))
stereoscope_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope13))

sdwls_markers_stlsq3_deneme$Name <- c("WISpR")
dwls_rmse3_deneme$Name <- c("DWLS")
rctd_rmse3_deneme$Name <- c("RCTD")
giotto_rmse3_deneme$Name <- c("S-DWLS")
spotlight_rmse3_deneme$Name <- c("SPOTlight")
stereoscope_rmse3_deneme$Name <- c("Stereoscope")

all_merged3_scenerio1 <- rbind(sdwls_markers_stlsq_deneme, sdwls_markers_stlsq2_deneme, 
                               sdwls_markers_stlsq3_deneme, dwls_rmse_deneme, dwls_rmse2_deneme, 
                               dwls_rmse3_deneme, rctd_rmse_deneme, rctd_rmse2_deneme, 
                               rctd_rmse3_deneme, giotto_rmse_deneme, giotto_rmse2_deneme, 
                               giotto_rmse3_deneme, spotlight_rmse_deneme, spotlight_rmse2_deneme, 
                               spotlight_rmse3_deneme, stereoscope_rmse_deneme, stereoscope_rmse2_deneme, 
                               stereoscope_rmse3_deneme)

dim(all_merged3_scenerio1)
all_merged3_scenerio1$Version <- c("SCENERIO_1")
library(plyr)
df <- ddply(all_merged3_scenerio1, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3_scenerio1

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
colnames(all_merged3_scenerio1) <- c("RMSE", "Name")

#scenerio2
#Trial1

sdwls_markers_stlsq <- read_csv("scenerio2-missingsc_celltype/scenerio2_stlsq_trial1_2_positivew.dat")
spotlight_scenerio0 <- read_csv("scenerio2-missingsc_celltype/scenerio2_spotlight.dat")
norm_weights_scenerio0 <- read.csv("scenerio2-missingsc_celltype/scenerio2_rctd_results.dat")
norm_weights_scenerio0 <- dplyr::select(norm_weights_scenerio0, -1)
giotto_scenerio0 <- readRDS('scenerio2-missingsc_celltype/scenerio2_spot_giotto_Results.rds')
giotto_scenerio0_2<- giotto_scenerio0@spatial_enrichment$DWLS
giotto_scenerio0_2 <- dplyr::select(giotto_scenerio0_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio0_2))))
giotto_scenerio0_2 <- as.matrix(giotto_scenerio0_2)
giotto_scenerio0_2_rear <- giotto_scenerio0_2[, s]
scenerio0_dwls <- readRDS('scenerio2-missingsc_celltype/scenerio2_DWLS_Results.rds')
cols <- c(8, 12, 3, 2, 9, 4, 10, 14,  11, 5, 1, 6, 7, 13)
scenerio0_dwls_rear <- as.matrix(scenerio0_dwls[, cols])
scenerio0_dwls <- t(scenerio0_dwls)
dim(scenerio0_dwls)

sdwls_markers_stlsq <- as.matrix(sdwls_markers_stlsq)
giotto_scenerio0_2_rear2 <- as.matrix(giotto_scenerio0_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio0))))
spotlight_scenerio0 <- spotlight_scenerio0[,s]

spotlight_scenerio0 <- as.matrix(spotlight_scenerio0)
dim(sdwls_markers_stlsq); dim(norm_weights_scenerio0); dim(giotto_scenerio0_2_rear2); dim(scenerio0_dwls_rear); dim(spotlight_scenerio0)

sdwls_markers_stlsq_prop <- sdwls_markers_stlsq / rowSums(sdwls_markers_stlsq)
sdwls_markers_stlsq_prop[is.na(sdwls_markers_stlsq_prop)] = 0

stereoscope21 <- read.csv("scenerio2-missingsc_celltype/stereoscope_results1/countsvisium.st-heart_with_fraction_0.1_scenerio2visium/W.2023-08-10192804.743760.tsv", sep = "\t")
stereoscope21 <- dplyr::select(stereoscope21, -1)
head(stereoscope21)

ground_tr <- read.csv("scenerio2-missingsc_celltype/proportionsvisium.st-heart_with_fraction_0.1_scenerio2visium.tsv", sep = "\t")

ground_tr <- dplyr::select(ground_tr, -1)
ground_tr = subset(ground_tr, select = -c(X7))
ground_tr <- as.matrix(ground_tr)
dim(ground_tr)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)
gr_norm0[is.na(gr_norm0)] = 0

sdwls_markers_stlsq_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(sdwls_markers_stlsq_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls_rear))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio0))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio0_2_rear2))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio0))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope21))

sdwls_markers_stlsq_deneme[is.na(sdwls_markers_stlsq_deneme)] = 0
dwls_rmse_deneme[is.na(dwls_rmse_deneme)] = 0
rctd_rmse_deneme[is.na(rctd_rmse_deneme)] = 0
giotto_rmse_deneme[is.na(giotto_rmse_deneme)] = 0
spotlight_rmse_deneme[is.na(spotlight_rmse_deneme)] = 0
stereoscope_rmse_deneme[is.na(stereoscope_rmse_deneme)] = 0

sdwls_markers_stlsq_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

#Trial2
sdwls_markers_stlsq2 <- read_csv("scenerio2-missingsc_celltype/scenerio2_stlsq_trial2_2_positivew.dat")
spotlight_scenerio2 <- read_csv("scenerio2-missingsc_celltype/scenerio2_spotlight_trial2.dat")
norm_weights_scenerio2 <- read.csv("scenerio2-missingsc_celltype/scenerio2_rctd_results_trial2.dat")
norm_weights_scenerio2 <- dplyr::select(norm_weights_scenerio2, -1)
giotto_scenerio2 <- readRDS('scenerio2-missingsc_celltype/scenerio2_spot_giotto_Results_trial2.rds')
giotto_scenerio2_2<- giotto_scenerio2@spatial_enrichment$DWLS
giotto_scenerio2_2 <- dplyr::select(giotto_scenerio2_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio2_2))))
giotto_scenerio2_2 <- as.matrix(giotto_scenerio2_2)
giotto_scenerio2_2_rear <- giotto_scenerio2_2[, s]
scenerio2_dwls <- readRDS('scenerio2-missingsc_celltype/scenerio2_DWLS_Results_trial2.rds')
scenerio2_dwls <- t(scenerio2_dwls)
cols <- c(8, 12, 3, 2, 9, 4, 10, 14,  11, 5, 1, 6, 7, 13)
scenerio2_dwls_rear <- as.matrix(scenerio2_dwls[, cols])

sdwls_markers_stlsq2 <- as.matrix(sdwls_markers_stlsq2)
giotto_scenerio2_2_rear2 <- as.matrix(giotto_scenerio2_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio2))))
spotlight_scenerio2 <- spotlight_scenerio2[,s]

spotlight_scenerio2 <- as.matrix(spotlight_scenerio2)
dim(sdwls_markers_stlsq2); dim(norm_weights_scenerio2); dim(giotto_scenerio2_2_rear2); dim(scenerio2_dwls_rear); dim(spotlight_scenerio2)

sdwls_markers_stlsq2_prop <- sdwls_markers_stlsq2 / rowSums(sdwls_markers_stlsq2)
sdwls_markers_stlsq2_prop[is.na(sdwls_markers_stlsq2_prop)] = 0

stereoscope22 <- read.csv("scenerio2-missingsc_celltype/stereoscope_results2/countstrial2_visium.st-heart_with_fraction_0.1_scenerio2_trial2trial2_visium/W.2023-08-11085953.695786.tsv", sep = "\t")
stereoscope22 <- dplyr::select(stereoscope22, -1)

ground_tr2 <- read.csv("scenerio2-missingsc_celltype/proportionstrial2_visium.st-heart_with_fraction_0.1_scenerio2_trial2trial2_visium.tsv", sep = "\t")
ground_tr2 <- dplyr::select(ground_tr2, -1)
ground_tr2 = subset(ground_tr2, select = -c(X7))
ground_tr2 <- as.matrix(ground_tr2)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm2 <- ground_tr2/rowSums(ground_tr2)
gr_norm2[is.na(gr_norm2)] = 0

sdwls_markers_stlsq2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(sdwls_markers_stlsq2_prop))
dwls_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(scenerio2_dwls_rear))
rctd_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(norm_weights_scenerio2))
giotto_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(giotto_scenerio2_2_rear2))
spotlight_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(spotlight_scenerio2))
stereoscope_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope22))

sdwls_markers_stlsq2_deneme[is.na(sdwls_markers_stlsq2_deneme)] = 0
dwls_rmse2_deneme[is.na(dwls_rmse2_deneme)] = 0
rctd_rmse2_deneme[is.na(rctd_rmse2_deneme)] = 0
giotto_rmse2_deneme[is.na(giotto_rmse2_deneme)] = 0
spotlight_rmse2_deneme[is.na(spotlight_rmse2_deneme)] = 0
stereoscope_rmse2_deneme[is.na(stereoscope_rmse2_deneme)] = 0


sdwls_markers_stlsq2_deneme$Name <- c("WISpR")
dwls_rmse2_deneme$Name <- c("DWLS")
rctd_rmse2_deneme$Name <- c("RCTD")
giotto_rmse2_deneme$Name <- c("S-DWLS")
spotlight_rmse2_deneme$Name <- c("SPOTlight")
stereoscope_rmse2_deneme$Name <- c("Stereoscope")


#Trial3
sdwls_markers_stlsq3 <- read_csv("scenerio2-missingsc_celltype/scenerio2_stlsq_trial3_2_positivew.dat")
spotlight_scenerio3 <- read_csv("scenerio2-missingsc_celltype/scenerio2_spotlight_trial3.dat")
norm_weights_scenerio3 <- read.csv("scenerio2-missingsc_celltype/scenerio2_rctd_results_trial3.dat")
norm_weights_scenerio3 <- dplyr::select(norm_weights_scenerio3, -1)
giotto_scenerio3 <- readRDS('scenerio2-missingsc_celltype/scenerio2_spot_giotto_Results_trial3.rds')
giotto_scenerio3_2<- giotto_scenerio3@spatial_enrichment$DWLS
giotto_scenerio3_2 <- dplyr::select(giotto_scenerio3_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio3_2))))
giotto_scenerio3_2 <- as.matrix(giotto_scenerio3_2)
giotto_scenerio3_2_rear <- giotto_scenerio3_2[, s]
scenerio3_dwls <- readRDS('scenerio2-missingsc_celltype/scenerio2_DWLS_Results_trial3.rds')
dim(scenerio3_dwls)
scenerio3_dwls <- t(scenerio3_dwls)
cols <- c(8, 12, 3, 2, 9, 4, 10, 14,  11, 5, 1, 6, 7, 13)
scenerio3_dwls_rear <- as.matrix(scenerio3_dwls[, cols])
sdwls_markers_stlsq3 <- as.matrix(sdwls_markers_stlsq3)
giotto_scenerio3_2_rear2 <- as.matrix(giotto_scenerio3_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio3))))
spotlight_scenerio3 <- spotlight_scenerio3[,s]
spotlight_scenerio3 <- as.matrix(spotlight_scenerio3)
dim(sdwls_markers_stlsq3); dim(norm_weights_scenerio3); dim(giotto_scenerio3_2_rear2); dim(scenerio3_dwls_rear); dim(spotlight_scenerio3)

sdwls_markers_stlsq3_prop <- sdwls_markers_stlsq3 / rowSums(sdwls_markers_stlsq3)
sdwls_markers_stlsq3_prop[is.na(sdwls_markers_stlsq3_prop)] = 0

stereoscope23 <- read.csv("scenerio2-missingsc_celltype/stereoscope_results3/countstrial3_visium.st-heart_with_fraction_0.1_scenerio2_trial2trial3_visium/W.2023-08-11184106.957224.tsv", sep = "\t")
stereoscope23 <- dplyr::select(stereoscope23, -1)

ground_tr3 <- read.csv("scenerio2-missingsc_celltype/proportionstrial3_visium.st-heart_with_fraction_0.1_scenerio2_trial2trial3_visium.tsv", sep = "\t")
ground_tr3 <- dplyr::select(ground_tr3, -1)
ground_tr3 = subset(ground_tr3, select = -c(X7))
ground_tr3 <- as.matrix(ground_tr3)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(sdwls_markers_stlsq3_prop))
dwls_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(scenerio3_dwls_rear))
rctd_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(norm_weights_scenerio3))
giotto_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(giotto_scenerio3_2_rear2))
spotlight_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(spotlight_scenerio3))
stereoscope_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope23))

sdwls_markers_stlsq3_deneme[is.na(sdwls_markers_stlsq3_deneme)] = 0
dwls_rmse3_deneme[is.na(dwls_rmse3_deneme)] = 0
rctd_rmse3_deneme[is.na(rctd_rmse3_deneme)] = 0
giotto_rmse3_deneme[is.na(giotto_rmse3_deneme)] = 0
spotlight_rmse3_deneme[is.na(spotlight_rmse3_deneme)] = 0
stereoscope_rmse3_deneme[is.na(stereoscope_rmse3_deneme)] = 0

sdwls_markers_stlsq3_deneme$Name <- c("WISpR")
dwls_rmse3_deneme$Name <- c("DWLS")
rctd_rmse3_deneme$Name <- c("RCTD")
giotto_rmse3_deneme$Name <- c("S-DWLS")
spotlight_rmse3_deneme$Name <- c("SPOTlight")
stereoscope_rmse3_deneme$Name <- c("Stereoscope")

all_merged3_scenerio2 <- rbind(sdwls_markers_stlsq_deneme, sdwls_markers_stlsq2_deneme, 
                               sdwls_markers_stlsq3_deneme, dwls_rmse_deneme, dwls_rmse2_deneme, 
                               dwls_rmse3_deneme, rctd_rmse_deneme, rctd_rmse2_deneme, 
                               rctd_rmse3_deneme, giotto_rmse_deneme, giotto_rmse2_deneme, 
                               giotto_rmse3_deneme, spotlight_rmse_deneme, spotlight_rmse2_deneme, 
                               spotlight_rmse3_deneme, stereoscope_rmse_deneme, stereoscope_rmse2_deneme, 
                               stereoscope_rmse3_deneme)

all_merged3_scenerio2$Version <- c("SCENERIO_2")
library(plyr)
df <- ddply(all_merged3_scenerio2, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3_scenerio2
all_merged3_2
cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
colnames(all_merged3_scenerio2) <- c("RMSE", "Name")

#scenerio3
#Trial1

sdwls_markers_stlsq <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_stlsq_trial1_2_positivew.dat")
spotlight_scenerio0 <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_spotlight.dat")
norm_weights_scenerio0 <- read.csv("scenerio3-mislabelledst_celltype/scenerio3_rctd_results.dat")
norm_weights_scenerio0 <- dplyr::select(norm_weights_scenerio0, -1)
giotto_scenerio0 <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_spot_giotto_Results.rds')
giotto_scenerio0_2<- giotto_scenerio0@spatial_enrichment$DWLS
giotto_scenerio0_2 <- dplyr::select(giotto_scenerio0_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio0_2))))
giotto_scenerio0_2 <- as.matrix(giotto_scenerio0_2)
giotto_scenerio0_2_rear <- giotto_scenerio0_2[, s]
scenerio0_dwls <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_DWLS_Results.rds')
scenerio0_dwls <- t(scenerio0_dwls)
dim(scenerio0_dwls)
cols <- c(7, 11, 2, 8, 3, 9, 12,14,10,  4,1, 5, 6, 13)
scenerio0_dwls_rear <- as.matrix(scenerio0_dwls[, cols])

sdwls_markers_stlsq <- as.matrix(sdwls_markers_stlsq)
giotto_scenerio0_2_rear2 <- as.matrix(giotto_scenerio0_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio0))))
spotlight_scenerio0 <- spotlight_scenerio0[,s]

spotlight_scenerio0 <- as.matrix(spotlight_scenerio0)
dim(sdwls_markers_stlsq); dim(norm_weights_scenerio0); dim(giotto_scenerio0_2_rear2); dim(scenerio0_dwls_rear); dim(spotlight_scenerio0)

sdwls_markers_stlsq_prop <- sdwls_markers_stlsq / rowSums(sdwls_markers_stlsq)
sdwls_markers_stlsq_prop[is.na(sdwls_markers_stlsq_prop)] = 0

stereoscope31 <- read.csv("scenerio3-mislabelledst_celltype/stereoscope_results1/countsvisium.st-heart_with_fraction_0.1_mislabelled_ct2to3visium/W.2023-08-12090904.632209.tsv", sep = "\t")
stereoscope31 <- dplyr::select(stereoscope31, -1)
head(stereoscope31)

ground_tr <- read.csv("scenerio3-mislabelledst_celltype/proportionsvisium.st-heart_with_fraction_0.1_mislabelled_ct2to3visium.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)
dim(ground_tr)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(sdwls_markers_stlsq_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls_rear))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio0))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio0_2_rear2))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio0))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope31))

sdwls_markers_stlsq_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

#Trial2
sdwls_markers_stlsq2 <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_stlsq_trial2_positivew.dat")
spotlight_scenerio2 <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_spotlight_trial2.dat")
norm_weights_scenerio2 <- read.csv("scenerio3-mislabelledst_celltype/scenerio3_rctd_results_trial2.dat")
norm_weights_scenerio2 <- dplyr::select(norm_weights_scenerio2, -1)
giotto_scenerio2 <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_spot_giotto_Results_trial2.rds')
giotto_scenerio2_2<- giotto_scenerio2@spatial_enrichment$DWLS
giotto_scenerio2_2 <- dplyr::select(giotto_scenerio2_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio2_2))))
giotto_scenerio2_2 <- as.matrix(giotto_scenerio2_2)
giotto_scenerio2_2_rear <- giotto_scenerio2_2[, s]
scenerio2_dwls <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_DWLS_Results_trial2.rds')
scenerio2_dwls <- t(scenerio2_dwls)
cols <- c(7, 11, 2, 8, 3, 9, 12,14,10,  4,1, 5, 6, 13)
scenerio2_dwls_rear <- as.matrix(scenerio2_dwls[, cols])

sdwls_markers_stlsq2 <- as.matrix(sdwls_markers_stlsq2)
giotto_scenerio2_2_rear2 <- as.matrix(giotto_scenerio2_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio2))))
spotlight_scenerio2 <- spotlight_scenerio2[,s]
spotlight_scenerio2 <- as.matrix(spotlight_scenerio2)
dim(sdwls_markers_stlsq2); dim(norm_weights_scenerio2); dim(giotto_scenerio2_2_rear2); dim(scenerio2_dwls_rear); dim(spotlight_scenerio2)

sdwls_markers_stlsq2_prop <- sdwls_markers_stlsq2 / rowSums(sdwls_markers_stlsq2)
sdwls_markers_stlsq2_prop[is.na(sdwls_markers_stlsq2_prop)] = 0

stereoscope32 <- read.csv("scenerio3-mislabelledst_celltype/stereoscope_results2/countstrial3_visium.st-heart_with_fraction_0.1_mislabelled_ct2to3_trial2trial3_visium/W.2023-08-12195045.474513.tsv", sep = "\t")
stereoscope32 <- dplyr::select(stereoscope32, -1)

ground_tr2 <- read.csv("scenerio3-mislabelledst_celltype/proportionstrial3_visium.st-heart_with_fraction_0.1_mislabelled_ct2to3_trial2trial3_visium.tsv", sep = "\t")
ground_tr2 <- dplyr::select(ground_tr2, -1)
ground_tr2 <- as.matrix(ground_tr2)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(sdwls_markers_stlsq2_prop))
dwls_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(scenerio2_dwls_rear))
rctd_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(norm_weights_scenerio2))
giotto_rmse2_deneme<- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(giotto_scenerio2_2_rear2))
spotlight_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr2/rowSums(ground_tr2)), as.matrix(spotlight_scenerio2))
stereoscope_rmse2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope32))

sdwls_markers_stlsq2_deneme$Name <- c("WISpR")
dwls_rmse2_deneme$Name <- c("DWLS")
rctd_rmse2_deneme$Name <- c("RCTD")
giotto_rmse2_deneme$Name <- c("S-DWLS")
spotlight_rmse2_deneme$Name <- c("SPOTlight")
stereoscope_rmse2_deneme$Name <- c("Stereoscope")

#Trial3
sdwls_markers_stlsq3 <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_stlsq_trial3_2_positivew.dat")
spotlight_scenerio3 <- read_csv("scenerio3-mislabelledst_celltype/scenerio3_spotlight_trial3.dat")
norm_weights_scenerio3 <- read.csv("scenerio3-mislabelledst_celltype/scenerio3_rctd_results_trial3.dat")
norm_weights_scenerio3 <- dplyr::select(norm_weights_scenerio3, -1)
giotto_scenerio3 <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_spot_giotto_Results_trial3.rds')
giotto_scenerio3_2<- giotto_scenerio3@spatial_enrichment$DWLS
giotto_scenerio3_2 <- dplyr::select(giotto_scenerio3_2, -1)
s = order(unique(as.numeric(colnames(giotto_scenerio3_2))))
giotto_scenerio3_2 <- as.matrix(giotto_scenerio3_2)
giotto_scenerio3_2_rear <- giotto_scenerio3_2[, s]
scenerio3_dwls <- readRDS('scenerio3-mislabelledst_celltype/scenerio3_DWLS_Results_trial3.rds')
scenerio3_dwls <- t(scenerio3_dwls)
cols <- c(7, 11, 2, 8, 3, 9, 12,14,10,  4,1, 5, 6, 13)
scenerio3_dwls_rear <- as.matrix(scenerio3_dwls[, cols])
sdwls_markers_stlsq3 <- as.matrix(sdwls_markers_stlsq3)
giotto_scenerio3_2_rear2 <- as.matrix(giotto_scenerio3_2_rear)
s = order(unique(as.numeric(colnames(spotlight_scenerio3))))
spotlight_scenerio3 <- spotlight_scenerio3[,s]
spotlight_scenerio3 <- as.matrix(spotlight_scenerio3)
dim(sdwls_markers_stlsq3); dim(norm_weights_scenerio3); dim(giotto_scenerio3_2_rear2); dim(scenerio3_dwls_rear); dim(spotlight_scenerio3)

sdwls_markers_stlsq3_prop <- sdwls_markers_stlsq3 / rowSums(sdwls_markers_stlsq3)
sdwls_markers_stlsq3_prop[is.na(sdwls_markers_stlsq3_prop)] = 0

stereoscope33 <- read.csv("scenerio3-mislabelledst_celltype/stereoscope_results3/countstrial3_visium.st-heart_with_fraction_0.1_mislabelled_ct2to3_trial3trial3_visium/W.2023-08-13080229.706961.tsv", sep = "\t")
stereoscope33 <- dplyr::select(stereoscope33, -1)


ground_tr3 <- read.csv("scenerio3-mislabelledst_celltype/proportionstrial3_visium.st-heart_with_fraction_0.1_mislabelled_ct2to3_trial3trial3_visium.tsv", sep = "\t")
ground_tr3 <- dplyr::select(ground_tr3, -1)
ground_tr3 <- as.matrix(ground_tr3)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

sdwls_markers_stlsq3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(sdwls_markers_stlsq3_prop))
dwls_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(scenerio3_dwls_rear))
rctd_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(norm_weights_scenerio3))
giotto_rmse3_deneme<- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(giotto_scenerio3_2_rear2))
spotlight_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr3/rowSums(ground_tr3)), as.matrix(spotlight_scenerio3))
stereoscope_rmse3_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope33))

sdwls_markers_stlsq3_deneme$Name <- c("WISpR")
dwls_rmse3_deneme$Name <- c("DWLS")
rctd_rmse3_deneme$Name <- c("RCTD")
giotto_rmse3_deneme$Name <- c("S-DWLS")
spotlight_rmse3_deneme$Name <- c("SPOTlight")
stereoscope_rmse3_deneme$Name <- c("Stereoscope")


all_merged3_scenerio3 <- rbind(sdwls_markers_stlsq_deneme, sdwls_markers_stlsq2_deneme, 
                               sdwls_markers_stlsq3_deneme, dwls_rmse_deneme, dwls_rmse2_deneme, 
                               dwls_rmse3_deneme, rctd_rmse_deneme, rctd_rmse2_deneme, 
                               rctd_rmse3_deneme, giotto_rmse_deneme, giotto_rmse2_deneme, 
                               giotto_rmse3_deneme, spotlight_rmse_deneme, spotlight_rmse2_deneme, 
                               spotlight_rmse3_deneme, stereoscope_rmse_deneme, stereoscope_rmse2_deneme, 
                               stereoscope_rmse3_deneme)
#all_merged_cbind <- cbind(nuray01_rmse_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme)
dim(all_merged3_scenerio3)
all_merged3_scenerio3$Version <- c("SCENERIO_3")
library(plyr)
library(RColorBrewer)
library(ggsci)
df <- ddply(all_merged3_scenerio3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3_scenerio3
all_merged3_2
cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
colnames(all_merged3_scenerio3) <- c("RMSE", "Name")

#50_0_50
#Trial1
readsc <- read.csv("mixture_50_0_50/Mix_50_0_50_1000spot_sc.dat", sep = ",")
head(readsc)

Wispr <- read.csv("mixture_50_0_50/WISpR_deconvoluted_deneme.Mix_50_0_50_1000spot_sc.dat", sep = "\t")
dim(readsc); dim(Wispr)
colnames(Wispr) <- colnames(readsc)
head(Wispr)
Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_0_50/mixture50_0_50_Spotlight_trial_1000spot_results.csv")
head(spotlight_scenerio0)
spotlight_scenerio0 <- dplyr::select(spotlight_scenerio0, -1)
norm_weights_scenerio02 <- read.csv("mixture_50_0_50/mixture_50_0_50_RCTD_trial_results_1000spots.csv")
head(norm_weights_scenerio02)
giotto_scenerio0 <- read.csv('mixture_50_0_50/mixture_50_0_50_giotto_Results_visium_1000.csv')
head(giotto_scenerio0)
scenerio0_dwls <- read.csv('mixture_50_0_50/DWLS_Results_visium_1000spots.dat')
head(scenerio0_dwls)
stereoscope50_0_50 <- read.csv("mixture_50_0_50/stereoscope_results/counts.synthetic_50_0_50_1000spots_visium/W.2023-07-26160940.733692.tsv", sep = "\t")
stereoscope50_0_50 <- dplyr::select(stereoscope50_0_50, -1)
colnames(stereoscope50_0_50)

dim(scenerio0_dwls); dim(giotto_scenerio0); dim(norm_weights_scenerio02); dim(spotlight_scenerio0); dim(Wispr) ;dim(stereoscope50_0_50)

scenerio0_dwls2 <- scenerio0_dwls %>% dplyr::select(c(1,5,6,4,3,8,12,2,9,13,10,11,7,15,14,"Astrocytes_14", 
                                                      "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                                      "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                                      "Neurons_26","Neurons_27", "Neurons_63",
                                                      "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                                      "Vascular_67", "Vascular_69"), everything())


giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(c(1,5,6,4,3,8,12,2,9,13,10,11,7,15,14,"Astrocytes_14", 
                                                          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                                          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                                          "Neurons_26","Neurons_27", "Neurons_63",
                                                          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                                          "Vascular_67", "Vascular_69"), everything())
norm_weights_scenerio022 <- norm_weights_scenerio02 %>% dplyr::select(c(1,2,8,9,10,11,12,13,14,15,3,4,5,6,7,"Astrocytes_14", 
                                                                        "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                                                        "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                                                        "Neurons_26","Neurons_27", "Neurons_63",
                                                                        "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                                                        "Vascular_67", "Vascular_69"), everything())
spotlight_scenerio02 <- spotlight_scenerio0 %>% dplyr::select(c(1,2,8,9,10,11,12,13,14,15,3,4,5,6,7,"Astrocytes_14", 
                                                                "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                                                "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                                                "Neurons_26","Neurons_27", "Neurons_63",
                                                                "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                                                "Vascular_67", "Vascular_69"), everything())

Wispr2 <- Wispr %>% dplyr::select(c(1,2,8,9,10,11,12,13,14,15,3,4,5,6,7,"Astrocytes_14", 
                                    "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                    "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                    "Neurons_26","Neurons_27", "Neurons_63",
                                    "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                    "Vascular_67", "Vascular_69"), everything())
stereoscope50_0_50_2 <- stereoscope50_0_50 %>% dplyr::select(c(1,2,8,9,10,11,12,13,14,15,3,4,5,6,7,"Astrocytes_14", 
                                                               "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
                                                               "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
                                                               "Neurons_26","Neurons_27", "Neurons_63",
                                                               "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
                                                               "Vascular_67", "Vascular_69"), everything())
head(Wispr2)
colnames(Wispr2) == colnames(scenerio0_dwls2)

cols <- c("Astrocytes_14", 
          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
          "Neurons_26","Neurons_27", "Neurons_63",
          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
          "Vascular_67", "Vascular_69")

not_cols <- c("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14")

Wispr3_r <- Wispr2[,cols]
norm_weights_scenerio03_r <- norm_weights_scenerio02[,cols]
giotto_scenerio03_r <- giotto_scenerio02[,cols]
scenerio0_dwls3_r <- scenerio0_dwls2[,cols]
spotlight_scenerio03_r <- spotlight_scenerio02[,cols]
stereoscope50_0_50_3_r <- stereoscope50_0_50_2[,cols]

Wispr3_prop <- Wispr3_r / rowSums(Wispr3_r)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03_r/rowSums(giotto_scenerio03_r)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03_r/rowSums(norm_weights_scenerio03_r)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3_r/rowSums(scenerio0_dwls3_r)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03_r/rowSums(spotlight_scenerio03_r)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_0_50_4 <- stereoscope50_0_50_3_r / rowSums(stereoscope50_0_50_3_r)
stereoscope50_0_50_4[is.na(stereoscope50_0_50_4)] = 0

ground_tr <- read.csv("mixture_50_0_50/proportions.synthetic_50_0_50_1000spots_visium.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope50_0_50_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_0_50_4))

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope50_0_50_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope50_0_50_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses

colnames(all_merged3) <- c("RMSE", "Name")

###
Wispr3 <- Wispr2[,not_cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,not_cols]
giotto_scenerio03 <- giotto_scenerio02[,not_cols]
giotto_scenerio03<- na.omit(giotto_scenerio03)
scenerio0_dwls3 <- scenerio0_dwls2[,not_cols]
scenerio0_dwls3<- na.omit(scenerio0_dwls3)
spotlight_scenerio03 <- spotlight_scenerio02[,not_cols]
stereoscope50_0_50_3 <- stereoscope50_0_50_2[,not_cols]

Wispr3_sum <- sqrt(rowSums(Wispr3)^2) %>% as.data.frame()
norm_weights_scenerio03_sum <- sqrt(rowSums(norm_weights_scenerio03)^2) %>% as.data.frame()
giotto_scenerio03_sum <- sqrt(rowSums(giotto_scenerio03)^2) %>% as.data.frame()
scenerio0_dwls3_sum <- sqrt(rowSums(scenerio0_dwls3)^2) %>% as.data.frame()
spotlight_scenerio03_sum <- sqrt(rowSums(spotlight_scenerio03)^2) %>% as.data.frame()
stereoscope50_0_50_3_sum <- sqrt(rowSums(stereoscope50_0_50_3)^2) %>% as.data.frame()

Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_0_50_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_0_50_3_sum)
df <- ddply(non_real_merged, c("Name"), summarize, Mean = mean(.), SD = sd(.))
colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")
names(all_merged3)
mix_50_0_50 <- rbind(all_merged3, non_real_merged)
dim(mix_50_0_50)
mix_50_0_50$Version <- c("MIX_50_0_50")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_0_50, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))

#50_10_40
#Trial1
readsc <- read.csv("mixture_50_10_40/Mix_50_10_40_1000spot_sc.dat", sep = ",")
head(readsc)

Wispr <- read.csv("mixture_50_10_40/WISpR_deconvoluted_deneme.Mix_50_10_40_1000spot_sc.dat", sep = "\t")
dim(readsc); dim(Wispr)
colnames(Wispr) <- colnames(readsc)
head(Wispr)
Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_10_40/mixture50_10_40_Spotlight_trial_1000spot_results.csv")
head(spotlight_scenerio0)
spotlight_scenerio0 <- dplyr::select(spotlight_scenerio0, -1)
norm_weights_scenerio02 <- read.csv("mixture_50_10_40/mixture_50_10_40_RCTD_trial_results_1000spots.csv")
head(norm_weights_scenerio02)
giotto_scenerio0 <- read.csv('mixture_50_10_40/mixture_50_10_40_giotto_Results_visium_1000.csv')
head(giotto_scenerio0)
scenerio0_dwls <- read.csv('mixture_50_10_40/DWLS_Results_visium_1000spots.dat')
head(scenerio0_dwls)
col <- colnames(Wispr)
colnames(Wispr) == colnames(giotto_scenerio02)

stereoscope50_10_40 <- read.csv("mixture_50_10_40/stereoscope_results/counts.synthetic_50_10_40_visium_1000spots/W.2023-07-28090423.208142.tsv", sep = "\t")
stereoscope50_10_40 <- dplyr::select(stereoscope50_10_40, -1)
colnames(stereoscope50_10_40)

scenerio0_dwls2 <- scenerio0_dwls %>% dplyr::select(col, everything())
giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(col, everything())
stereoscope50_10_40_2 <- stereoscope50_10_40 %>% dplyr::select(col, everything())


cols <- c("Astrocytes_14", 
          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
          "Neurons_26","Neurons_27", "Neurons_63",
          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
          "Vascular_67", "Vascular_69")

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0)
colnames(Wispr); colnames(ground_tr);colnames(norm_weights_scenerio02)
Wispr3 <- Wispr[,cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,cols]
giotto_scenerio03 <- giotto_scenerio02[,cols]
scenerio0_dwls3 <- scenerio0_dwls2[,cols]
spotlight_scenerio03 <- spotlight_scenerio02[,cols]
stereoscope50_10_40_3 <- stereoscope50_10_40_2[,cols]

Wispr3_prop <- Wispr3 / rowSums(Wispr3)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03/rowSums(giotto_scenerio03)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03/rowSums(norm_weights_scenerio03)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3/rowSums(scenerio0_dwls3)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03/rowSums(spotlight_scenerio03)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_10_40_4 <- stereoscope50_10_40_3/rowSums(stereoscope50_10_40_3)
stereoscope50_10_40_4[is.na(stereoscope50_10_40_4)] = 0

ground_tr <- read.csv("mixture_50_10_40/proportions.synthetic_50_10_40_visium_1000spots.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_10_40_4))

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses

colnames(all_merged3) <- c("RMSE", "Name")

####
dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0); dim(stereoscope_rmse_deneme)
colnames(scenerio0_dwls2)
non_cols <- c("X0","X1","X10","X11","X12","X13","X14","X2","X3","X4","X5","X6","X7","X8","X9",
              "Astrocytes_42", "Neurons_11", "Neurons_12","Neurons_14","Neurons_15",
              "Neurons_18","Neurons_21","Neurons_23","Neurons_24", "Neurons_48", 
              "Neurons_51","Neurons_52", "Oligos_5", "Oligos_53", "Vascular_68")

Wispr3_non <- Wispr[, non_cols]
norm_weights_scenerio03_non <- norm_weights_scenerio02[,non_cols]
giotto_scenerio03_non <- giotto_scenerio02[,non_cols]
scenerio0_dwls3_non <- scenerio0_dwls[,non_cols]
spotlight_scenerio03_non <- spotlight_scenerio0[,non_cols]
stereoscope50_10_40_3_non <- stereoscope50_10_40_2[,non_cols]

Wispr3_sum <- sqrt((rowSums(Wispr3_non))^2) %>% as.data.frame()
norm_weights_scenerio03_sum <- sqrt((rowSums(norm_weights_scenerio03_non))^2) %>% as.data.frame()
giotto_scenerio03_sum <- sqrt((rowSums(giotto_scenerio03_non))^2) %>% as.data.frame()
scenerio0_dwls3_sum <- sqrt((rowSums(scenerio0_dwls3_non))^2) %>% as.data.frame()
spotlight_scenerio03_sum <- sqrt((rowSums(spotlight_scenerio03_non))^2) %>% as.data.frame()
stereoscope50_10_40_3_sum <- sqrt((rowSums(stereoscope50_10_40_3_non))^2) %>% as.data.frame()

scenerio0_dwls3_non2 <- round(scenerio0_dwls3_non, digits = 10)

scenerio0_dwls3_sum <- round(scenerio0_dwls3_sum, digits = 10)
giotto_scenerio03_sum <- round(giotto_scenerio03_sum, digits = 10)
giotto_scenerio03_sum[is.na(giotto_scenerio03_sum)] = 0

Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_10_40_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_10_40_3_sum)

level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")
names(all_merged3)
mix_50_10_40 <- rbind(all_merged3, non_real_merged)
dim(mix_50_10_40)
mix_50_10_40$Version <- c("MIX_50_10_40")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_10_40, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))

#50_20_30
#Trial1
readsc <- read.csv("mixture_50_20_30/Mix_50_20_30_1000spot_sc.dat", sep = ",")

Wispr <- read.csv("mixture_50_20_30/WISpR_deconvoluted_deneme.Mix_50_20_30_1000spot_sc.dat", sep = "\t")
colnames(Wispr) <- colnames(readsc)
Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_20_30/mixture50_20_30_Spotlight_trial_1000spot_results.csv")
head(spotlight_scenerio0)
spotlight_scenerio0 <- dplyr::select(spotlight_scenerio0, -1)
norm_weights_scenerio02 <- read.csv("mixture_50_20_30/mixture_50_20_30_RCTD_trial_results_1000spots.csv")
head(norm_weights_scenerio02)
giotto_scenerio0 <- read.csv('mixture_50_20_30/mixture_50_20_30_giotto_Results_visium_1000.csv')
head(giotto_scenerio0)
scenerio0_dwls <- read.csv('mixture_50_20_30/DWLS_Results_visium_1000spots.dat')
head(scenerio0_dwls)
col <- colnames(Wispr)
colnames(Wispr) == colnames(spotlight_scenerio0)

stereoscope50_20_30 <- read.csv("mixture_50_20_30/stereoscope_results/counts.synthetic_50_20_30_visium_1000spots/W.2023-07-29204744.086255.tsv", sep = "\t")
stereoscope50_20_30 <- dplyr::select(stereoscope50_20_30, -1)
colnames(stereoscope50_20_30)

scenerio0_dwls2 <- scenerio0_dwls %>% dplyr::select(col, everything())
giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(col, everything())
stereoscope50_20_30_2 <- stereoscope50_20_30 %>% dplyr::select(col, everything())
colnames(giotto_scenerio02)

cols <- c("Astrocytes_14", 
          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
          "Neurons_26","Neurons_27", "Neurons_63",
          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
          "Vascular_67", "Vascular_69")

non_cols <- c("X0","X1","X10","X11","X12","X13","X14","X2","X3","X4","X5","X6","X7","X8","X9",
              "Astrocytes_42", "Neurons_11", "Neurons_12","Neurons_14","Neurons_15",
              "Neurons_18","Neurons_21","Neurons_23","Neurons_24", "Neurons_48", 
              "Neurons_51","Neurons_52", "Oligos_5", "Oligos_53", "Vascular_68")

Wispr3 <- Wispr[,cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,cols]
giotto_scenerio03 <- giotto_scenerio02[,cols]
scenerio0_dwls3 <- scenerio0_dwls2[,cols]
spotlight_scenerio03 <- spotlight_scenerio02[,cols]
stereoscope50_20_30_3 <- stereoscope50_20_30_2[,cols]

Wispr3_prop <- Wispr3 / rowSums(Wispr3)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03/rowSums(giotto_scenerio03)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03/rowSums(norm_weights_scenerio03)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3/rowSums(scenerio0_dwls3)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03/rowSums(spotlight_scenerio03)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_20_30_4 <- stereoscope50_20_30_3/rowSums(stereoscope50_20_30_3)
stereoscope50_20_30_4[is.na(stereoscope50_20_30_4)] = 0

ground_tr <- read.csv("mixture_50_20_30/proportions.synthetic_50_20_30_visium_1000spots.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_20_30_4))

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses


colnames(all_merged3) <- c("RMSE", "Name")

####

Wispr3_non <- Wispr[, non_cols]
norm_weights_scenerio03_non <- norm_weights_scenerio02[,non_cols]
giotto_scenerio03_non <- giotto_scenerio02[,non_cols]
scenerio0_dwls3_non <- scenerio0_dwls2[,non_cols]
spotlight_scenerio03_non <- spotlight_scenerio0[,non_cols]
stereoscope50_20_30_3_non <- stereoscope50_20_30_2[,non_cols]

scenerio0_dwls3_non[is.na(scenerio0_dwls3_non)] = 0
giotto_scenerio03_non[is.na(giotto_scenerio03_non)] = 0

giotto_scenerio03_non <- round(giotto_scenerio03_non, digits = 15)
scenerio0_dwls3_non <- round(scenerio0_dwls3_non, digits = 13)

Wispr3_sum <- (sqrt(rowSums(Wispr3_non))) %>% as.data.frame()
norm_weights_scenerio03_sum <- (sqrt(rowSums(norm_weights_scenerio03_non))) %>% as.data.frame()
giotto_scenerio03_sum <- (sqrt(rowSums(giotto_scenerio03_non))) %>% as.data.frame()
scenerio0_dwls3_sum <- (sqrt(rowSums(scenerio0_dwls3_non))) %>% as.data.frame()
spotlight_scenerio03_sum <- (sqrt(rowSums(spotlight_scenerio03_non))) %>% as.data.frame()
stereoscope50_20_30_3_sum <- (sqrt(rowSums(stereoscope50_20_30_3_non))) %>% as.data.frame()

Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_20_30_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_20_30_3_sum)

colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")

mix_50_20_30 <- rbind(all_merged3, non_real_merged)
mix_50_20_30$Version <- c("MIX_50_20_30")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_20_30, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))

#50_30_20
#Trial1
readsc <- read.csv("mixture_50_30_20/Mix_50_30_20_1000spot_sc.dat", sep = ",")
head(readsc)

Wispr <- read.csv("mixture_50_30_20/WISpR_deconvoluted_deneme.Mix_50_30_20_1000spot_sc.dat", sep = "\t")

colnames(Wispr) <- colnames(readsc)

Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_30_20/mixture50_30_20_Spotlight_trial_1000spot_results.csv")

spotlight_scenerio0 <- dplyr::select(spotlight_scenerio0, -1)
norm_weights_scenerio02 <- read.csv("mixture_50_30_20/mixture_50_30_20_RCTD_trial_results_1000spots.csv")

giotto_scenerio0 <- read.csv('mixture_50_30_20/mixture_50_30_20_giotto_Results_visium_1000.csv')

scenerio0_dwls <- read.csv('mixture_50_30_20/DWLS_Results_visium_1000spots.dat')

col <- colnames(Wispr)

stereoscope50_30_20 <- read.csv("mixture_50_30_20/stereoscope_results/counts.synthetic_50_30_20_visium_1000spots/W.2023-07-30150828.129298.tsv", sep = "\t")
stereoscope50_30_20 <- dplyr::select(stereoscope50_30_20, -1)


scenerio0_dwls2 <- scenerio0_dwls %>% dplyr::select(col, everything())
giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(col, everything())
stereoscope50_30_20_2 <- stereoscope50_30_20 %>% dplyr::select(col, everything())
colnames(Wispr) == colnames(scenerio0_dwls2)

cols <- c("Astrocytes_14", 
          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
          "Neurons_26","Neurons_27", "Neurons_63",
          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
          "Vascular_67", "Vascular_69")

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0); dim(stereoscope50_30_20_2)

Wispr3 <- Wispr[,cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,cols]
giotto_scenerio03 <- giotto_scenerio02[,cols]
scenerio0_dwls3 <- scenerio0_dwls2[,cols]
spotlight_scenerio03 <- spotlight_scenerio02[,cols]
stereoscope50_30_20_3 <- stereoscope50_30_20_2[,cols]

Wispr3_prop <- Wispr3 / rowSums(Wispr3)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03/rowSums(giotto_scenerio03)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03/rowSums(norm_weights_scenerio03)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3/rowSums(scenerio0_dwls3)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03/rowSums(spotlight_scenerio03)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_30_20_4 <- stereoscope50_30_20_3/ rowSums(stereoscope50_30_20_3)
stereoscope50_30_20_4[is.na(stereoscope50_30_20_4)] = 0

ground_tr <- read.csv("mixture_50_30_20/proportions.synthetic_50_30_20_visium_1000spots.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_30_20_4))

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))

all_merged3_2 <- all_merged3
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses


colnames(all_merged3) <- c("RMSE", "Name")

#######
non_cols <- c("X0","X1","X10","X11","X12","X13","X14","X2","X3","X4","X5","X6","X7","X8","X9",
              "Astrocytes_42", "Neurons_11", "Neurons_12","Neurons_14","Neurons_15",
              "Neurons_18","Neurons_21","Neurons_23","Neurons_24", "Neurons_48", 
              "Neurons_51","Neurons_52", "Oligos_5", "Oligos_53", "Vascular_68")

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0)


Wispr3_non <- Wispr[, non_cols]
norm_weights_scenerio03_non <- norm_weights_scenerio02[,non_cols]
giotto_scenerio03_non <- giotto_scenerio02[,non_cols]
scenerio0_dwls3_non <- scenerio0_dwls2[,non_cols]
spotlight_scenerio03_non <- spotlight_scenerio0[,non_cols]
stereoscope50_30_20_3_non <- stereoscope50_30_20_2[,non_cols]

#giotto_scenerio03_non[is.na(giotto_scenerio03_non)] = 0
giotto_scenerio03_non <- round(giotto_scenerio03_non, digits = 13)
scenerio0_dwls3_non <- round(scenerio0_dwls3_non, digits = 13)
giotto_scenerio03_non[is.na(giotto_scenerio03_non)] = 0


Wispr3_sum <- (sqrt(rowSums(Wispr3_non))) %>% as.data.frame()
norm_weights_scenerio03_sum <- (sqrt(rowSums(norm_weights_scenerio03_non))) %>% as.data.frame()
giotto_scenerio03_sum <- (sqrt(rowSums(giotto_scenerio03_non))) %>% as.data.frame()
scenerio0_dwls3_sum <- (sqrt(rowSums(scenerio0_dwls3_non))) %>% as.data.frame()
spotlight_scenerio03_sum <- (sqrt(rowSums(spotlight_scenerio03_non))) %>% as.data.frame()
stereoscope50_30_20_3_sum <- (sqrt(rowSums(stereoscope50_30_20_3_non))) %>% as.data.frame()

Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_30_20_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_30_20_3_sum)

level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")

mix_50_30_20 <- rbind(all_merged3, non_real_merged)

mix_50_30_20$Version <- c("MIX_50_30_20")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_30_20, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))

#50_40_10
#Trial1
readsc <- read.csv("mixture_50_40_10/Mix_50_40_10_1000spot_sc.dat", sep = ",")


Wispr <- read.csv("mixture_50_40_10/WISpR_deconvoluted_deneme.Mix_50_40_10_1000spot_sc.dat", sep = "\t")

colnames(Wispr) <- colnames(readsc)

Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_40_10/mixture50_40_10_Spotlight_trial_1000spot_results.csv")

spotlight_scenerio0 <- dplyr::select(spotlight_scenerio0, -1)
norm_weights_scenerio02 <- read.csv("mixture_50_40_10/mixture_50_40_10_RCTD_trial_results_1000spots.csv")

giotto_scenerio0 <- read.csv('mixture_50_40_10/mixture_50_40_10_giotto_Results_visium_1000.csv')

scenerio0_dwls <- read.csv('mixture_50_40_10/DWLS_Results_visium_1000spots.dat')

col <- colnames(Wispr)
colnames(Wispr) == colnames(spotlight_scenerio0)

stereoscope50_40_10 <- read.csv("mixture_50_40_10/stereoscope_results/counts.synthetic_50_40_10_visium_1000spots/W.2023-07-31102251.461536.tsv", sep = "\t")
stereoscope50_40_10 <- dplyr::select(stereoscope50_40_10, -1)

scenerio0_dwls2 <- scenerio0_dwls %>% dplyr::select(col, everything())
giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(col, everything())
stereoscope50_40_10_2 <- stereoscope50_40_10 %>% dplyr::select(col, everything())


cols <- c("Astrocytes_14", 
          "Astrocytes_40", "Astrocytes_41", "Blood_73", "Ependymal_47", 
          "Excluded_38", "Immune_32", "Immune_34", "Immune_35",  "Neurons_25",
          "Neurons_26","Neurons_27", "Neurons_63",
          "Oligos_0",  "Oligos_1", "Oligos_14", "Vascular_14",
          "Vascular_67", "Vascular_69")

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0); dim(stereoscope50_40_10_2)

Wispr3 <- Wispr[,cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,cols]
giotto_scenerio03 <- giotto_scenerio02[,cols]
scenerio0_dwls3 <- scenerio0_dwls2[,cols]
spotlight_scenerio03 <- spotlight_scenerio02[,cols]
stereoscope50_40_10_3 <- stereoscope50_40_10_2[,cols]

Wispr3_prop <- Wispr3 / rowSums(Wispr3)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03/rowSums(giotto_scenerio03)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03/rowSums(norm_weights_scenerio03)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3/rowSums(scenerio0_dwls3)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03/rowSums(spotlight_scenerio03)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_40_10_4 <- stereoscope50_40_10_3 / rowSums(stereoscope50_40_10_3)
stereoscope50_40_10_4[is.na(stereoscope50_40_10_4)] = 0

ground_tr <- read.csv("mixture_50_40_10/proportions.synthetic_50_40_10_visium_1000spots.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)

newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_40_10_4)) 

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope_rmse_deneme)
#all_merged_cbind <- cbind(nuray01_rmse_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
all_merged3_2 <- all_merged3
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses


colnames(all_merged3) <- c("RMSE", "Name")

######

non_cols <- c("X0","X1","X2","X3","X4","X5","X7",
              "Astrocytes_42", "Neurons_11", "Neurons_12","Neurons_14","Neurons_15",
              "Neurons_18","Neurons_21","Neurons_23","Neurons_24", "Neurons_48", 
              "Neurons_51","Neurons_52", "Oligos_5", "Oligos_53", "Vascular_68")

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls2); dim(spotlight_scenerio0)


Wispr3_non <- Wispr[, non_cols]
norm_weights_scenerio03_non <- norm_weights_scenerio02[,non_cols]
giotto_scenerio03_non <- giotto_scenerio02[,non_cols]
scenerio0_dwls3_non <- scenerio0_dwls2[,non_cols]
spotlight_scenerio03_non <- spotlight_scenerio0[,non_cols]
stereoscope50_40_10_3_non <- stereoscope50_40_10_2[,non_cols]

giotto_scenerio03_non[is.na(giotto_scenerio03_non)] = 0
giotto_scenerio03_non <- round(giotto_scenerio03_non, digits = 14)
scenerio0_dwls3_non <- round(scenerio0_dwls3_non, digits = 14)

Wispr3_sum <- (sqrt(rowSums(Wispr3_non))) %>% as.data.frame()
norm_weights_scenerio03_sum <- (sqrt(rowSums(norm_weights_scenerio03_non))) %>% as.data.frame()
giotto_scenerio03_sum <- (sqrt(rowSums(giotto_scenerio03_non))) %>% as.data.frame()
scenerio0_dwls3_sum <- (sqrt(rowSums(scenerio0_dwls3_non))) %>% as.data.frame()
spotlight_scenerio03_sum <- (sqrt(rowSums(spotlight_scenerio03_non))) %>% as.data.frame()
stereoscope50_40_10_3_sum <- (sqrt(rowSums(stereoscope50_40_10_3_non))) %>% as.data.frame()

Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_40_10_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_40_10_3_sum)

level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")

mix_50_40_10 <- rbind(all_merged3, non_real_merged)

mix_50_40_10$Version <- c("MIX_50_40_10")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_40_10, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))

#50_50_0
#Trial1
readsc <- read.csv("mixture_50_50_0/Mix_50_50_0_1000spot_sc.dat", sep = ",")


Wispr <- read.csv("mixture_50_50_0/WISpR_deconvoluted_deneme.Mix_50_50_0_1000spot_sc.dat", sep = "\t")

colnames(Wispr) <- colnames(readsc)

Wispr <- dplyr::select(Wispr, -1)

spotlight_scenerio0 <- read.csv("mixture_50_50_0/mixture50_50_0_Spotlight_trial_1000spot_results.csv")

norm_weights_scenerio02 <- read.csv("mixture_50_50_0/mixture_50_50_0_RCTD_trial_results_1000spots.csv")

giotto_scenerio0 <- read.csv('mixture_50_50_0/mixture_50_50_0_giotto_Results_visium_1000.csv')

scenerio0_dwls <- read.csv('mixture_50_50_0/DWLS_Results_visium_1000spots.dat')

colnames(Wispr) == colnames(scenerio0_dwls)
cols = colnames(Wispr)

stereoscope50_50_0 <- read.csv("mixture_50_50_0/stereoscope_results/counts.synthetic_50_50_0_visium_1000spot/W.2023-08-01100026.663212.tsv", sep = "\t")
stereoscope50_50_0 <- dplyr::select(stereoscope50_50_0, -1)


giotto_scenerio02 <- giotto_scenerio0 %>% dplyr::select(cols, everything())
colnames(Wispr) == colnames(norm_weights_scenerio02)

dim(Wispr); dim(norm_weights_scenerio02); dim(giotto_scenerio02); dim(scenerio0_dwls); dim(spotlight_scenerio0); dim(stereoscope50_50_0)

ground_tr <- read.csv("mixture_50_50_0/proportions.synthetic_50_50_0_visium_1000spot.tsv", sep = "\t")
ground_tr <- dplyr::select(ground_tr, -1)
ground_tr <- as.matrix(ground_tr)

cols <- colnames(ground_tr)

colnames(Wispr); colnames(ground_tr)

Wispr3 <- Wispr[,cols]
norm_weights_scenerio03 <- norm_weights_scenerio02[,cols]
giotto_scenerio03 <- giotto_scenerio02[,cols]
scenerio0_dwls3 <- scenerio0_dwls[,cols]
spotlight_scenerio03 <- spotlight_scenerio0[,cols]
stereoscope50_50_0_3 <- stereoscope50_50_0[,cols]


Wispr3_prop <- Wispr3 / rowSums(Wispr3)
Wispr3_prop[is.na(Wispr3_prop)] = 0

giotto_scenerio04 <- giotto_scenerio03/rowSums(giotto_scenerio03)
giotto_scenerio04[is.na(giotto_scenerio04)] = 0

norm_weights_scenerio04 <- norm_weights_scenerio03/rowSums(norm_weights_scenerio03)
norm_weights_scenerio04[is.na(norm_weights_scenerio04)] = 0

scenerio0_dwls4 <- scenerio0_dwls3/rowSums(scenerio0_dwls3)
scenerio0_dwls4[is.na(scenerio0_dwls4)] = 0

spotlight_scenerio04 <- spotlight_scenerio03/rowSums(spotlight_scenerio03)
spotlight_scenerio04[is.na(spotlight_scenerio04)] = 0

stereoscope50_50_0_4 <- stereoscope50_50_0_3 / rowSums(stereoscope50_50_0_3)
stereoscope50_50_0_4[is.na(stereoscope50_50_0_4)] = 0


dim(Wispr3_prop); dim(norm_weights_scenerio04); dim(giotto_scenerio04); dim(scenerio0_dwls4); dim(spotlight_scenerio04);dim(stereoscope50_50_0_4)


newRMSESpot <- function(sig, spot) {
  rmse_spot <- as.data.frame(matrix(0, nrow = ncol(spot), ncol = 1))
  cc <- seq(1, nrow(spot), by=1)
  for(i in cc)
    rmse_spot[i,]<-hydroGOF::rmse(sig[i,],spot[i,])
  print(i)
  return(rmse_spot)
}

gr_norm0 <- ground_tr/rowSums(ground_tr)

Wispr2_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(Wispr3_prop))
dwls_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(scenerio0_dwls4))
rctd_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(norm_weights_scenerio04))
giotto_rmse_deneme<- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(giotto_scenerio04))
spotlight_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(spotlight_scenerio04))
stereoscope_rmse_deneme <- newRMSESpot(as.matrix(ground_tr/rowSums(ground_tr)), as.matrix(stereoscope50_50_0_4))

Wispr2_deneme$Name <- c("WISpR")
dwls_rmse_deneme$Name <- c("DWLS")
rctd_rmse_deneme$Name <- c("RCTD")
giotto_rmse_deneme$Name <- c("S-DWLS")
spotlight_rmse_deneme$Name <- c("SPOTlight")
stereoscope_rmse_deneme$Name <- c("Stereoscope")

all_merged3 <- rbind(Wispr2_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme, stereoscope_rmse_deneme)
#all_merged_cbind <- cbind(nuray01_rmse_deneme, dwls_rmse_deneme, rctd_rmse_deneme, giotto_rmse_deneme, spotlight_rmse_deneme)

library(plyr)
df <- ddply(all_merged3, c("Name"), summarize, Mean = mean(V1), SD = sd(V1))
df
all_merged3_2 <- all_merged3
#all_merged3_2$Name <- relevel(all_merged3_2$Name, ref="SparSc")
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442")
library("ggsci")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses


colnames(all_merged3) <- c("RMSE", "Name")

###

non_cols <- c("Astrocytes_42", "Neurons_11", "Neurons_12","Neurons_14","Neurons_15",
              "Neurons_18","Neurons_21","Neurons_23","Neurons_24", "Neurons_48", 
              "Neurons_51","Neurons_52", "Oligos_5", "Oligos_53", "Vascular_68")

Wispr3_non <- Wispr[, non_cols]
norm_weights_scenerio03_non <- norm_weights_scenerio02[,non_cols]
giotto_scenerio03_non <- giotto_scenerio02[,non_cols]
scenerio0_dwls3_non <- scenerio0_dwls[,non_cols]
spotlight_scenerio03_non <- spotlight_scenerio0[,non_cols]
stereoscope50_50_0_3_non <- stereoscope50_50_0[,non_cols]


Wispr3_sum <- sqrt((rowSums(Wispr3_non))^2) %>% as.data.frame()
norm_weights_scenerio03_sum <- sqrt((rowSums(norm_weights_scenerio03_non))^2) %>% as.data.frame()
giotto_scenerio03_sum <- sqrt((rowSums(giotto_scenerio03_non))^2) %>% as.data.frame()
scenerio0_dwls3_sum <- sqrt((rowSums(scenerio0_dwls3_non))^2) %>% as.data.frame()
spotlight_scenerio03_sum <- sqrt((rowSums(spotlight_scenerio03_non))^2) %>% as.data.frame()
stereoscope50_50_0_3_sum <- sqrt((rowSums(stereoscope50_50_0_3_non))^2) %>% as.data.frame() 

scenerio0_dwls3_sum <- round(scenerio0_dwls3_sum, digits = 15)
giotto_scenerio03_sum <- round(giotto_scenerio03_sum, digits = 15)

giotto_scenerio03_sum[is.na(giotto_scenerio03_sum)] = 0


Wispr3_sum$Name <- c("WISpR")
norm_weights_scenerio03_sum$Name <- c("RCTD")
giotto_scenerio03_sum$Name <- c("S-DWLS")
scenerio0_dwls3_sum$Name <- c("DWLS")
spotlight_scenerio03_sum$Name <- c("SPOTlight")
stereoscope50_50_0_3_sum$Name <- c("Stereoscope")

non_real_merged <- rbind(Wispr3_sum, norm_weights_scenerio03_sum, giotto_scenerio03_sum, scenerio0_dwls3_sum, spotlight_scenerio03_sum, stereoscope50_50_0_3_sum)

level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
colnames(non_real_merged) <- c("RMSE", "Name")
colnames(all_merged3) <- c("RMSE", "Name")

mix_50_50_0 <- rbind(all_merged3, non_real_merged)

mix_50_50_0$Version <- c("MIX_50_50_0")
level_order <- c('WISpR', 'DWLS', 'RCTD', 'S-DWLS', 'SPOTlight', "Stereoscope") #this vector might be useful for other plots/analyses
df <- ddply(mix_50_50_0, c("Name"), summarize, Mean = mean(RMSE), SD = sd(RMSE))


#merge all results
final_mix <- rbind(mix_50_0_50, mix_50_10_40, mix_50_20_30, mix_50_30_20,
                   mix_50_40_10, mix_50_50_0)

colnames(all_merged3_scenerio0) <- c("RMSE","Name","Version")
colnames(all_merged3_scenerio1) <- c("RMSE","Name","Version")
colnames(all_merged3_scenerio2) <- c("RMSE","Name","Version")
colnames(all_merged3_scenerio3) <- c("RMSE","Name","Version")

final_scen <- rbind(all_merged3_scenerio0, all_merged3_scenerio1, all_merged3_scenerio2, 
                    all_merged3_scenerio3)

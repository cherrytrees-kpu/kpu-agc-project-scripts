#### This document is revised from RRS pilot (https://github.com/notothen/radpilot)

#### Source script containing restriction enzymes,
#### size windows and functions 

#### load packages required for these functions
#### install them first if you don't have them
library(here) # to shorten file paths
library(SimRAD) # for in silico digestion
library(bioanalyzeR) # to read bioanalyzer files
library(tidyverse) # to arrange data and plot
library(scales) # to calculate percentages
library(RColorBrewer) # for plotting colors

#### load restriction enzymes
#####
## SbfI : 5'--CCTGCA  GG--3
sbf1_5 <- "CCTGCA"
sbf1_3 <- "GG"
## EcoRI : 5'--G  AATTC--3
ecor1_5 <- "G"
ecor1_3 <- "AATTC"
## SphI : 5'--GCATG  C--3
sph1_5 <- "GCATG"
sph1_3 <- "C"
## PstI : 5'--CTGCA  G--3
pst1_5 <- "CTGCA"
pst1_3 <- "G"
## MspI : 5'--C CGG--3
msp1_5 <- "C"
msp1_3 <- "CGG"
## MseI : 5'--T  TAA--3
mse1_5 <- "T"
mse1_3 <- "TAA"

## prepare data frame
run <- c(1:21)
names <- c("sbf1", "ecor1", "sph1", "pst1", "msp1", "mse1",
           "sbf1-ecor1", "sbf1-sph1", "sbf1-pst1", "sbf1-msp1", "sbf1-mse1",
           "ecor1-sph1", "ecor1-pst1", "ecor1-msp1", "ecor1-mse1",
           "sph1-pst1", "sph1-msp1", "sph1-mse1",
           "pst1-msp1", "pst1-mse1", 
           "msp1-mse1")
forward <- c(sbf1_5, ecor1_5, sph1_5, pst1_5, msp1_5, mse1_5,
             sbf1_5, sbf1_5, sbf1_5, sbf1_5, sbf1_5, 
             ecor1_5, ecor1_5, ecor1_5, ecor1_5, 
             sph1_5, sph1_5, sph1_5, 
             pst1_5, pst1_5, 
             msp1_5)
reverse <- c(sbf1_3, ecor1_3, sph1_3, pst1_3, msp1_3, mse1_3,
             sbf1_3, sbf1_3, sbf1_3, sbf1_3, sbf1_3, 
             ecor1_3, ecor1_3, ecor1_3, ecor1_3, 
             sph1_3, sph1_3, sph1_3, 
             pst1_3, pst1_3, 
             msp1_3)
forward2 <- c((rep('', 6)),
              ecor1_5, sph1_5, pst1_5, msp1_5, mse1_5, sph1_5, pst1_5, msp1_5, mse1_5, pst1_5, msp1_5, mse1_5, msp1_5, mse1_5, mse1_5)
reverse2 <- c((rep('', 6)),
              ecor1_3, sph1_3, pst1_3, msp1_3, mse1_3, sph1_3, pst1_3, msp1_3, mse1_3, pst1_3, msp1_3, mse1_3, msp1_3, mse1_3, mse1_3)

## combine
recto_REs <- data.frame(run, names, forward, reverse, forward2, reverse2)



## add size windows
lower_size <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
upper_size <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550)


#### create digest function
#####
recto_digest <- function(genome, enzyme, minsize, maxsize, ratio){
    dig1 <- insilico.digest(genome, enzyme$forward[1], enzyme$reverse[1], verbose = F)
    dig1_1 <- (length(dig1)-1)*ratio
    dig1_2 <- (length(size.select(dig1, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig1_3 <- (length(size.select(dig1, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig1_4 <- (length(size.select(dig1, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig1_5 <- (length(size.select(dig1, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig1_6 <- (length(size.select(dig1, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig1_7 <- (length(size.select(dig1, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig1_8 <- (length(size.select(dig1, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig1_9 <- (length(size.select(dig1, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig1_10 <- (length(size.select(dig1, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig1_11 <- (length(size.select(dig1, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig1_12 <- (length(size.select(dig1, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig2 <- insilico.digest(genome, enzyme$forward[2], enzyme$reverse[2], verbose = F)
    dig2_1 <- (length(dig2)-1)*ratio
    dig2_2 <- (length(size.select(dig2, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig2_3 <- (length(size.select(dig2, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig2_4 <- (length(size.select(dig2, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig2_5 <- (length(size.select(dig2, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig2_6 <- (length(size.select(dig2, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig2_7 <- (length(size.select(dig2, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig2_8 <- (length(size.select(dig2, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig2_9 <- (length(size.select(dig2, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig2_10 <- (length(size.select(dig2, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig2_11 <- (length(size.select(dig2, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig2_12 <- (length(size.select(dig2, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig3 <- insilico.digest(genome, enzyme$forward[3], enzyme$reverse[3], verbose = F)
    dig3_1 <- (length(dig3)-1)*ratio
    dig3_2 <- (length(size.select(dig3, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig3_3 <- (length(size.select(dig3, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig3_4 <- (length(size.select(dig3, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig3_5 <- (length(size.select(dig3, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig3_6 <- (length(size.select(dig3, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig3_7 <- (length(size.select(dig3, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig3_8 <- (length(size.select(dig3, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig3_9 <- (length(size.select(dig3, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig3_10 <- (length(size.select(dig3, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig3_11 <- (length(size.select(dig3, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig3_12 <- (length(size.select(dig3, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig4 <- insilico.digest(genome, enzyme$forward[4], enzyme$reverse[4], verbose = F)
    dig4_1 <- (length(dig4)-1)*ratio
    dig4_2 <- (length(size.select(dig4, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig4_3 <- (length(size.select(dig4, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig4_4 <- (length(size.select(dig4, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig4_5 <- (length(size.select(dig4, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig4_6 <- (length(size.select(dig4, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig4_7 <- (length(size.select(dig4, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig4_8 <- (length(size.select(dig4, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig4_9 <- (length(size.select(dig4, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig4_10 <- (length(size.select(dig4, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig4_11 <- (length(size.select(dig4, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig4_12 <- (length(size.select(dig4, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig5 <- insilico.digest(genome, enzyme$forward[5], enzyme$reverse[5], verbose = F)
    dig5_1 <- (length(dig5)-1)*ratio
    dig5_2 <- (length(size.select(dig5, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig5_3 <- (length(size.select(dig5, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig5_4 <- (length(size.select(dig5, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig5_5 <- (length(size.select(dig5, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig5_6 <- (length(size.select(dig5, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig5_7 <- (length(size.select(dig5, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig5_8 <- (length(size.select(dig5, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig5_9 <- (length(size.select(dig5, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig5_10 <- (length(size.select(dig5, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig5_11 <- (length(size.select(dig5, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig5_12 <- (length(size.select(dig5, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig6 <- insilico.digest(genome, enzyme$forward[6], enzyme$reverse[6], verbose = F)
    dig6_1 <- (length(dig6)-1)*ratio
    dig6_2 <- (length(size.select(dig6, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig6_3 <- (length(size.select(dig6, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig6_4 <- (length(size.select(dig6, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig6_5 <- (length(size.select(dig6, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig6_6 <- (length(size.select(dig6, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig6_7 <- (length(size.select(dig6, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig6_8 <- (length(size.select(dig6, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig6_9 <- (length(size.select(dig6, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig6_10 <- (length(size.select(dig6, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig6_11 <- (length(size.select(dig6, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig6_12 <- (length(size.select(dig6, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
  
    dig7 <- insilico.digest(genome, enzyme$forward[7], enzyme$reverse[7], enzyme$forward2[7], enzyme$reverse2[7], verbose = F)
    dig7 <- adapt.select(dig7, type = "AB+BA", enzyme$forward[7], as.character(enzyme$reverse[7]), enzyme$forward2[7], as.character(enzyme$reverse2[7]))
    dig7_1 <- (length(dig7)-1)*ratio
    dig7_2 <- (length(size.select(dig7, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig7_3 <- (length(size.select(dig7, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig7_4 <- (length(size.select(dig7, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig7_5 <- (length(size.select(dig7, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig7_6 <- (length(size.select(dig7, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig7_7 <- (length(size.select(dig7, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig7_8 <- (length(size.select(dig7, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig7_9 <- (length(size.select(dig7, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig7_10 <- (length(size.select(dig7, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig7_11 <- (length(size.select(dig7, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig7_12 <- (length(size.select(dig7, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig8 <- insilico.digest(genome, enzyme$forward[8], enzyme$reverse[8], enzyme$forward2[8], enzyme$reverse2[8], verbose = F)
    dig8 <- adapt.select(dig8, type = "AB+BA", enzyme$forward[8], as.character(enzyme$reverse[8]), enzyme$forward2[8], as.character(enzyme$reverse2[8]))
    dig8_1 <- (length(dig8)-1)*ratio
    dig8_2 <- (length(size.select(dig8, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig8_3 <- (length(size.select(dig8, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig8_4 <- (length(size.select(dig8, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig8_5 <- (length(size.select(dig8, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig8_6 <- (length(size.select(dig8, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig8_7 <- (length(size.select(dig8, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig8_8 <- (length(size.select(dig8, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig8_9 <- (length(size.select(dig8, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig8_10 <- (length(size.select(dig8, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig8_11 <- (length(size.select(dig8, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig8_12 <- (length(size.select(dig8, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig9 <- insilico.digest(genome, enzyme$forward[9], enzyme$reverse[9], enzyme$forward2[9], enzyme$reverse2[9], verbose = F)
    dig9 <- adapt.select(dig9, type = "AB+BA", enzyme$forward[9], as.character(enzyme$reverse[9]), enzyme$forward2[9], as.character(enzyme$reverse2[9]))
    dig9_1 <- (length(dig9)-1)*ratio
    dig9_2 <- (length(size.select(dig9, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig9_3 <- (length(size.select(dig9, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig9_4 <- (length(size.select(dig9, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig9_5 <- (length(size.select(dig9, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig9_6 <- (length(size.select(dig9, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig9_7 <- (length(size.select(dig9, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig9_8 <- (length(size.select(dig9, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig9_9 <- (length(size.select(dig9, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig9_10 <- (length(size.select(dig9, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig9_11 <- (length(size.select(dig9, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig9_12 <- (length(size.select(dig9, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig10 <- insilico.digest(genome, enzyme$forward[10], enzyme$reverse[10], enzyme$forward2[10], enzyme$reverse2[10], verbose = F)
    dig10 <- adapt.select(dig10, type = "AB+BA", enzyme$forward[10], as.character(enzyme$reverse[10]), enzyme$forward2[10], as.character(enzyme$reverse2[10]))
    dig10_1 <- (length(dig10)-1)*ratio
    dig10_2 <- (length(size.select(dig10, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig10_3 <- (length(size.select(dig10, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig10_4 <- (length(size.select(dig10, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig10_5 <- (length(size.select(dig10, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig10_6 <- (length(size.select(dig10, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig10_7 <- (length(size.select(dig10, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig10_8 <- (length(size.select(dig10, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig10_9 <- (length(size.select(dig10, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig10_10 <- (length(size.select(dig10, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig10_11 <- (length(size.select(dig10, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig10_12 <- (length(size.select(dig10, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig11 <- insilico.digest(genome, enzyme$forward[11], enzyme$reverse[11], enzyme$forward2[11], enzyme$reverse2[11], verbose = F)
    dig11 <- adapt.select(dig11, type = "AB+BA", enzyme$forward[11], as.character(enzyme$reverse[11]), enzyme$forward2[11], as.character(enzyme$reverse2[11]))
    dig11_1 <- (length(dig11)-1)*ratio
    dig11_2 <- (length(size.select(dig11, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig11_3 <- (length(size.select(dig11, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig11_4 <- (length(size.select(dig11, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig11_5 <- (length(size.select(dig11, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig11_6 <- (length(size.select(dig11, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig11_7 <- (length(size.select(dig11, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig11_8 <- (length(size.select(dig11, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig11_9 <- (length(size.select(dig11, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig11_10 <- (length(size.select(dig11, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig11_11 <- (length(size.select(dig11, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig11_12 <- (length(size.select(dig11, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig12 <- insilico.digest(genome, enzyme$forward[12], enzyme$reverse[12], enzyme$forward2[12], enzyme$reverse2[12], verbose = F)
    dig12 <- adapt.select(dig12, type = "AB+BA", enzyme$forward[12], as.character(enzyme$reverse[12]), enzyme$forward2[12], as.character(enzyme$reverse2[12]))
    dig12_1 <- (length(dig12)-1)*ratio
    dig12_2 <- (length(size.select(dig12, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig12_3 <- (length(size.select(dig12, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig12_4 <- (length(size.select(dig12, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig12_5 <- (length(size.select(dig12, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig12_6 <- (length(size.select(dig12, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig12_7 <- (length(size.select(dig12, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig12_8 <- (length(size.select(dig12, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig12_9 <- (length(size.select(dig12, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig12_10 <- (length(size.select(dig12, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig12_11 <- (length(size.select(dig12, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig12_12 <- (length(size.select(dig12, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig13 <- insilico.digest(genome, enzyme$forward[13], enzyme$reverse[13], enzyme$forward2[13], enzyme$reverse2[13], verbose = F)
    dig13 <- adapt.select(dig13, type = "AB+BA", enzyme$forward[13], as.character(enzyme$reverse[13]), enzyme$forward2[13], as.character(enzyme$reverse2[13]))
    dig13_1 <- (length(dig13)-1)*ratio
    dig13_2 <- (length(size.select(dig13, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig13_3 <- (length(size.select(dig13, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig13_4 <- (length(size.select(dig13, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig13_5 <- (length(size.select(dig13, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig13_6 <- (length(size.select(dig13, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig13_7 <- (length(size.select(dig13, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig13_8 <- (length(size.select(dig13, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig13_9 <- (length(size.select(dig13, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig13_10 <- (length(size.select(dig13, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig13_11 <- (length(size.select(dig13, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig13_12 <- (length(size.select(dig13, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig14 <- insilico.digest(genome, enzyme$forward[14], enzyme$reverse[14], enzyme$forward2[14], enzyme$reverse2[14], verbose = F)
    dig14 <- adapt.select(dig14, type = "AB+BA", enzyme$forward[14], as.character(enzyme$reverse[14]), enzyme$forward2[14], as.character(enzyme$reverse2[14]))
    dig14_1 <- (length(dig14)-1)*ratio
    dig14_2 <- (length(size.select(dig14, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig14_3 <- (length(size.select(dig14, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig14_4 <- (length(size.select(dig14, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig14_5 <- (length(size.select(dig14, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig14_6 <- (length(size.select(dig14, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig14_7 <- (length(size.select(dig14, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig14_8 <- (length(size.select(dig14, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig14_9 <- (length(size.select(dig14, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig14_10 <- (length(size.select(dig14, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig14_11 <- (length(size.select(dig14, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig14_12 <- (length(size.select(dig14, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig15 <- insilico.digest(genome, enzyme$forward[15], enzyme$reverse[15], enzyme$forward2[15], enzyme$reverse2[15], verbose = F)
    dig15 <- adapt.select(dig15, type = "AB+BA", enzyme$forward[15], as.character(enzyme$reverse[15]), enzyme$forward2[15], as.character(enzyme$reverse2[15]))
    dig15_1 <- (length(dig15)-1)*ratio
    dig15_2 <- (length(size.select(dig15, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig15_3 <- (length(size.select(dig15, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig15_4 <- (length(size.select(dig15, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig15_5 <- (length(size.select(dig15, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig15_6 <- (length(size.select(dig15, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig15_7 <- (length(size.select(dig15, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig15_8 <- (length(size.select(dig15, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig15_9 <- (length(size.select(dig15, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig15_10 <- (length(size.select(dig15, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig15_11 <- (length(size.select(dig15, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig15_12 <- (length(size.select(dig15, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig16 <- insilico.digest(genome, enzyme$forward[16], enzyme$reverse[16], enzyme$forward2[16], enzyme$reverse2[16], verbose = F)
    dig16 <- adapt.select(dig16, type = "AB+BA", enzyme$forward[16], as.character(enzyme$reverse[16]), enzyme$forward2[16], as.character(enzyme$reverse2[16]))
    dig16_1 <- (length(dig16)-1)*ratio
    dig16_2 <- (length(size.select(dig16, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig16_3 <- (length(size.select(dig16, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig16_4 <- (length(size.select(dig16, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig16_5 <- (length(size.select(dig16, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig16_6 <- (length(size.select(dig16, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig16_7 <- (length(size.select(dig16, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig16_8 <- (length(size.select(dig16, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig16_9 <- (length(size.select(dig16, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig16_10 <- (length(size.select(dig16, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig16_11 <- (length(size.select(dig16, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig16_12 <- (length(size.select(dig16, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig17 <- insilico.digest(genome, enzyme$forward[17], enzyme$reverse[17], enzyme$forward2[17], enzyme$reverse2[17], verbose = F)
    dig17 <- adapt.select(dig17, type = "AB+BA", enzyme$forward[17], as.character(enzyme$reverse[17]), enzyme$forward2[17], as.character(enzyme$reverse2[17]))
    dig17_1 <- (length(dig17)-1)*ratio
    dig17_2 <- (length(size.select(dig17, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig17_3 <- (length(size.select(dig17, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig17_4 <- (length(size.select(dig17, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig17_5 <- (length(size.select(dig17, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig17_6 <- (length(size.select(dig17, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig17_7 <- (length(size.select(dig17, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig17_8 <- (length(size.select(dig17, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig17_9 <- (length(size.select(dig17, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig17_10 <- (length(size.select(dig17, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig17_11 <- (length(size.select(dig17, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig17_12 <- (length(size.select(dig17, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig18 <- insilico.digest(genome, enzyme$forward[18], enzyme$reverse[18], enzyme$forward2[18], enzyme$reverse2[18], verbose = F)
    dig18 <- adapt.select(dig18, type = "AB+BA", enzyme$forward[18], as.character(enzyme$reverse[18]), enzyme$forward2[18], as.character(enzyme$reverse2[18]))
    dig18_1 <- (length(dig18)-1)*ratio
    dig18_2 <- (length(size.select(dig18, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig18_3 <- (length(size.select(dig18, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig18_4 <- (length(size.select(dig18, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig18_5 <- (length(size.select(dig18, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig18_6 <- (length(size.select(dig18, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig18_7 <- (length(size.select(dig18, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig18_8 <- (length(size.select(dig18, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig18_9 <- (length(size.select(dig18, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig18_10 <- (length(size.select(dig18, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig18_11 <- (length(size.select(dig18, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig18_12 <- (length(size.select(dig18, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig19 <- insilico.digest(genome, enzyme$forward[19], enzyme$reverse[19], enzyme$forward2[19], enzyme$reverse2[19], verbose = F)
    dig19 <- adapt.select(dig19, type = "AB+BA", enzyme$forward[19], as.character(enzyme$reverse[19]), enzyme$forward2[19], as.character(enzyme$reverse2[19]))
    dig19_1 <- (length(dig19)-1)*ratio
    dig19_2 <- (length(size.select(dig19, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig19_3 <- (length(size.select(dig19, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig19_4 <- (length(size.select(dig19, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig19_5 <- (length(size.select(dig19, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig19_6 <- (length(size.select(dig19, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig19_7 <- (length(size.select(dig19, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig19_8 <- (length(size.select(dig19, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig19_9 <- (length(size.select(dig19, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig19_10 <- (length(size.select(dig19, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig19_11 <- (length(size.select(dig19, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig19_12 <- (length(size.select(dig19, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig20 <- insilico.digest(genome, enzyme$forward[20], enzyme$reverse[20], enzyme$forward2[20], enzyme$reverse2[20], verbose = F)
    dig20 <- adapt.select(dig20, type = "AB+BA", enzyme$forward[20], as.character(enzyme$reverse[20]), enzyme$forward2[20], as.character(enzyme$reverse2[20]))
    dig20_1 <- (length(dig20)-1)*ratio
    dig20_2 <- (length(size.select(dig20, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig20_3 <- (length(size.select(dig20, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig20_4 <- (length(size.select(dig20, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig20_5 <- (length(size.select(dig20, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig20_6 <- (length(size.select(dig20, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig20_7 <- (length(size.select(dig20, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig20_8 <- (length(size.select(dig20, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig20_9 <- (length(size.select(dig20, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig20_10 <- (length(size.select(dig20, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig20_11 <- (length(size.select(dig20, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig20_12 <- (length(size.select(dig20, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    dig21 <- insilico.digest(genome, enzyme$forward[21], enzyme$reverse[21], enzyme$forward2[21], enzyme$reverse2[21], verbose = F)
    dig21 <- adapt.select(dig21, type = "AB+BA", enzyme$forward[21], as.character(enzyme$reverse[21]), enzyme$forward2[21], as.character(enzyme$reverse2[21]))
    dig21_1 <- (length(dig21)-1)*ratio
    dig21_2 <- (length(size.select(dig21, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig21_3 <- (length(size.select(dig21, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig21_4 <- (length(size.select(dig21, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig21_5 <- (length(size.select(dig21, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig21_6 <- (length(size.select(dig21, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig21_7 <- (length(size.select(dig21, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig21_8 <- (length(size.select(dig21, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig21_9 <- (length(size.select(dig21, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig21_10 <- (length(size.select(dig21, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig21_11 <- (length(size.select(dig21, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig21_12 <- (length(size.select(dig21, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio

    fragments_all <- c(dig1_1, dig2_1, dig3_1, dig4_1, dig5_1, dig6_1, dig7_1, dig8_1, dig9_1, dig10_1, dig11_1, dig12_1, dig13_1, dig14_1, dig15_1, dig16_1, dig17_1, dig18_1, dig19_1, dig20_1, dig21_1)
    fragments_2 <- c(dig1_2, dig2_2, dig3_2, dig4_2, dig5_2, dig6_2, dig7_2, dig8_2, dig9_2, dig10_2, dig11_2, dig12_2, dig13_2, dig14_2, dig15_2, dig16_2, dig17_2, dig18_2, dig19_2, dig20_2, dig21_2)
    fragments_3 <- c(dig1_3, dig2_3, dig3_3, dig4_3, dig5_3, dig6_3, dig7_3, dig8_3, dig9_3, dig10_3, dig11_3, dig12_3, dig13_3, dig14_3, dig15_3, dig16_3, dig17_3, dig18_3, dig19_3, dig20_3, dig21_3)
    fragments_4 <- c(dig1_4, dig2_4, dig3_4, dig4_4, dig5_4, dig6_4, dig7_4, dig8_4, dig9_4, dig10_4, dig11_4, dig12_4, dig13_4, dig14_4, dig15_4, dig16_4, dig17_4, dig18_4, dig19_4, dig20_4, dig21_4)
    fragments_5 <- c(dig1_5, dig2_5, dig3_5, dig4_5, dig5_5, dig6_5, dig7_5, dig8_5, dig9_5, dig10_5, dig11_5, dig12_5, dig13_5, dig14_5, dig15_5, dig16_5, dig17_5, dig18_5, dig19_5, dig20_5, dig21_5)
    fragments_6 <- c(dig1_6, dig2_6, dig3_6, dig4_6, dig5_6, dig6_6, dig7_6, dig8_6, dig9_6, dig10_6, dig11_6, dig12_6, dig13_6, dig14_6, dig15_6, dig16_6, dig17_6, dig18_6, dig19_6, dig20_6, dig21_6)
    fragments_7 <- c(dig1_7, dig2_7, dig3_7, dig4_7, dig5_7, dig6_7, dig7_7, dig8_7, dig9_7, dig10_7, dig11_7, dig12_7, dig13_7, dig14_7, dig15_7, dig16_7, dig17_7, dig18_7, dig19_7, dig20_7, dig21_7)
    fragments_8 <- c(dig1_8, dig2_8, dig3_8, dig4_8, dig5_8, dig6_8, dig7_8, dig8_8, dig9_8, dig10_8, dig11_8, dig12_8, dig13_8, dig14_8, dig15_8, dig16_8, dig17_8, dig18_8, dig19_8, dig20_8, dig21_8)
    fragments_9 <- c(dig1_9, dig2_9, dig3_9, dig4_9, dig5_9, dig6_9, dig7_9, dig8_9, dig9_9, dig10_9, dig11_9, dig12_9, dig13_9, dig14_9, dig15_9, dig16_9, dig17_9, dig18_9, dig19_9, dig20_9, dig21_9)
    fragments_10 <- c(dig1_10, dig2_10, dig3_10, dig4_10, dig5_10, dig6_10, dig7_10, dig8_10, dig9_10, dig10_10, dig11_10, dig12_10, dig13_10, dig14_10, dig15_10, dig16_10, dig17_10, dig18_10, dig19_10, dig20_10, dig21_10)
    fragments_11 <- c(dig1_11, dig2_11, dig3_11, dig4_11, dig5_11, dig6_11, dig7_11, dig8_11, dig9_11, dig10_11, dig11_11, dig12_11, dig13_11, dig14_11, dig15_11, dig16_11, dig17_11, dig18_11, dig19_11, dig20_11, dig21_11)
    fragments_12 <- c(dig1_12, dig2_12, dig3_12, dig4_12, dig5_12, dig6_12, dig7_12, dig8_12, dig9_12, dig10_12, dig11_12, dig12_12, dig13_12, dig14_12, dig15_12, dig16_12, dig17_12, dig18_12, dig19_12, dig20_12, dig21_12)
    
    x <- data.frame(enzyme$run, enzyme$names, enzyme$forward, enzyme$reverse, enzyme$forward2, enzyme$reverse2,
                    fragments_all, fragments_2, fragments_3, fragments_4, fragments_5, fragments_6, fragments_7,
                    fragments_8, fragments_9, fragments_10, fragments_11, fragments_12)
    names(x)[names(x) == "fragments_2"] <- c(paste("fragments", minsize[1], maxsize[1], sep = "_"))
    names(x)[names(x) == "fragments_3"] <- c(paste("fragments", minsize[2], maxsize[2], sep = "_"))
    names(x)[names(x) == "fragments_4"] <- c(paste("fragments", minsize[3], maxsize[3], sep = "_"))
    names(x)[names(x) == "fragments_5"] <- c(paste("fragments", minsize[4], maxsize[4], sep = "_"))
    names(x)[names(x) == "fragments_6"] <- c(paste("fragments", minsize[5], maxsize[5], sep = "_"))
    names(x)[names(x) == "fragments_7"] <- c(paste("fragments", minsize[6], maxsize[6], sep = "_"))
    names(x)[names(x) == "fragments_8"] <- c(paste("fragments", minsize[7], maxsize[7], sep = "_"))
    names(x)[names(x) == "fragments_9"] <- c(paste("fragments", minsize[8], maxsize[8], sep = "_"))
    names(x)[names(x) == "fragments_10"] <- c(paste("fragments", minsize[9], maxsize[9], sep = "_"))
    names(x)[names(x) == "fragments_11"] <- c(paste("fragments", minsize[10], maxsize[10], sep = "_"))
    names(x)[names(x) == "fragments_12"] <- c(paste("fragments", minsize[11], maxsize[11], sep = "_"))
    names(x)[names(x) == "enzyme.run"] <- "run"
    return(x)
}
#####

## Load packages
library("dartR")
library("ggplot2")
library("dplyr")
setwd("~/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/outlier_analyses/outFLANK")
# ----------------- #
#
# OutFLANK
#
# ----------------- #

source("outflank.rdata")

##### Outflank and Bayescan results ####
#22 loci
bayescan_idx <- c(165, 179, 373, 428, 481, 654, 747, 799, 880, 907, 1081, 1086, 1117, 1134, 1170, 1402, 1433, 1789, 1864, 1905, 2129, 2130)
## PR = 100 (7 loci)
bayescan_idx <- c(373,  654,  799,  907, 1864, 2129, 2130)
#53 loci
outflank_idx <- c(165, 179, 295, 315, 335, 373, 390, 465, 481, 515, 548, 643, 654, 747, 799, 839, 880, 907, 967, 1061, 1081, 1086, 1117, 1134, 1170, 1179, 1188, 1402, 1433, 1446, 1578, 1789, 1864, 1882, 1905, 1910, 1991, 2001, 2029, 2050, 2081, 2129, 2130)

#### Running outflank #### 
infile <- "../../data/kinless.gen"
pop_names <- c("BRZ","KEY", "MRT", "PNS", "PR", "SCA", "TX", "VZ")
## Import data
full_dataset <- import2genind(infile, ncode = 3L)

## Change group labels
popNames(full_dataset) <- pop_names
  
## Run OutFLANK using dartR wrapper script
full_outflnk <- gl.outflank(
                  full_dataset, 
                  Hmin= 0.01, 
                  qthreshold = 0.05, 
                  LeftTrimFraction = 0.05,
                  RightTrimFraction = 0.05
                )

#### Processing results ####
## Outliers
full_outflnk <- full_outflnk$outflank$results

## Remove duplicated rows for each SNP
toRemove <- seq(1, nrow(full_outflnk), by=2)
full_outflnk <- full_outflnk[-toRemove, ]
#summary(full_outflnk)

## Get indexes for outliers
out_index <- which(full_outflnk$OutlierFlag==TRUE)
outflank_names <- locNames(full_dataset)[out_index]
#outlier_sub <- full_outflnk[out_index,]
#bayes_sub <- full_outflnk[bayescan_idx,]

#### Write output files ####
write.csv(full_outflnk, file = "outFLANK_kinless_FST.csv", row.names = FALSE)

full_out.outflnk <- locNames(full_dataset)[out_index]
full_out.outflnk
write.csv(full_out.outflnk, file="outlier_kinless_OutFLANK.csv")

#### Plots ####
## qvalue hist ##
hist(full_outflnk$qvalues, breaks = 100, 
     xlab = "Q values", main = "Q Value distribution of SNPs")
plot(x = full_outflnk$qvalues)


## fancy ggplot ##
full_outflnk$outlier <- "Neither"
outlier_both <- intersect(bayescan_idx, out_index)
full_outflnk$outlier[out_index] <- "outFLANK"
bayescan_safe <- bayescan_idx[bayescan_idx <= length(full_outflnk$LocusName)]
full_outflnk$outlier[bayescan_safe] <- "BayeScan"
full_outflnk$outlier[outlier_both] <- "Both"
full_outflnk$outlier <- factor(full_outflnk$outlier, levels = c("Neither", "outFLANK", "BayeScan", "Both"), ordered = TRUE)

mycolors <- c("#bbbbbb", "#4095b5", "#8a556e", "#f4cf30")
ggplot(data = full_outflnk, x = He, y = FST) + 
  theme_classic() +
  geom_point(aes(x = He, y = FST, col = outlier), alpha = 0.8,  size = 2.1) + 
  geom_point(data = subset(full_outflnk, outlier %in% "BayeScan"), aes(x = He, y = FST), color = "#8a556e") +
  geom_vline(xintercept = 0.1, alpha = 0.8, linetype = "dashed", size = 0.3) +
  labs(x = "Heterozygosity", y = "FST", color = "Detection") +
  scale_color_manual(values = mycolors) +
  ggtitle("Outlier Loci") +
  theme(plot.title = element_text(hjust = 0.5))


bayescan_out <- locNames(full_dataset)[bayescan_idx]
outflank_out <- locNames(full_dataset)[out_index]

write.table(outflank_out, file = "outflank_outliers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("outflank.rdata")

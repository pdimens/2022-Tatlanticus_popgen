## Load packages
library("ggplot2")
library("dplyr")
setwd("/mnt/Win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/extramiss/outlier_analyses/")

### outflank dataframe ####
all_loci <- read.csv("outFLANK/outFLANK_kinless_FST.csv")%>% 
  select(LocusName, He, FST, FSTNoCorr, OutlierFlag) %>% 
  mutate(LocusName = gsub("\\.[0-9]+", "", LocusName))

all_outflank <- filter(all_loci, OutlierFlag == TRUE)
all_outflank_putatives <- all_outflank$LocusName

##### Outflank and Bayescan results ####
# snps
out_bayescan <- readLines("bayescan/bayescan_outliers.txt")
bayescan_idx <- c()
# snps
all_outflank <- read.csv("outFLANK/outFLANK_kinless_FST.csv")
all_outflank$LocusName <- gsub("\\.[0-9]+", "", all_outflank$LocusName)
out_outflank <- all_outflank[all_outflank$OutlierFlag == T,]
locinames <- all_outflank$LocusName
write.table(out_outflank$LocusName, file = "outflank.outliers", row.names = F, col.names = F, quote = F)
outflank_idx <- c()


full_outflnk <- all_outflank %>% 
  mutate(putative = 
           case_when(
             (He >= 0.1 & OutlierFlag == TRUE) ~ TRUE,
             (He >= 0.1 & OutlierFlag == FALSE) ~ FALSE,
             (He < 0.1) ~ FALSE
           )
  )

putatives <- full_outflnk %>% filter(He >= 0.1 & OutlierFlag == TRUE)

# all shared putative outliers regardless of He <- used in graphic
all_shared <- intersect(out_bayescan, all_outflank_putatives)
all_onlybayes <- setdiff(out_bayescan, all_shared)
all_onlyoutflnk <- setdiff(all_outflank_putatives, all_shared)
all_nonout <- setdiff(all_loci$LocusName, c(all_onlybayes, all_onlyoutflnk, all_shared))

# only shared putative outliers (He threshold) <- used in publication
out_shared <- intersect(out_bayescan, out_outflank)
out_only_baye <- setdiff(out_bayescan, out_shared)
out_only_outf <- setdiff(putatives$LocusName, out_shared)
non_out <- setdiff(all_loci$LocusName, c(out_bayescan, out_outflank, out_shared))

## fancy ggplot ##
all_loci <-all_loci %>% 
  mutate(
    outlier = case_when(
      (LocusName %in% all_onlybayes) ~ "bayescan",
      (LocusName %in% all_onlyoutflnk) ~ "outflank",
      (LocusName %in% all_shared) ~ "bayescan + outflank",
      (LocusName %in% all_nonout ) ~ "neutral"
    )
  )


mycolors <- c("#8a556e", "#bbbbbb", "#f4cf30")

all_loci %>% 
  ggplot(x = He, y = FST) +
  geom_point(aes(x = He, y = FST, col = outlier), alpha = 0.8,  size = 2.1) + 
  geom_vline(xintercept = 0.1, alpha = 0.8, linetype = "dashed", size = 0.3) +
  labs(title = "Putative Outlier Loci", x = "Heterozygosity", y = "FST", color = "Outlier Detection") +
  scale_color_manual(values = mycolors)
ggsave("outlierplots.png", height = 6, width = 12, units = "in")
save.image(file = "outliers.rdata")
write.table(all_loci, file = "bscan_outflank_outliers.txt", row.names = FALSE, quote = FALSE)

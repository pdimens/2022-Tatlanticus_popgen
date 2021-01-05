library(adegenet)


#### Neutral ####
# load neutral env
load("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/DAPC/dapc_maf_kinless_neutral.rdata")
load("dapc_maf_kinless_neutral.rdata")

# load strata 
bft_strata <- read.delim("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/AMOVA/bft_strata.txt")
bft_strata <- read.delim("../AMOVA/bft_strata.txt")

# generate locality_year names
prefix <- unlist(lapply(strsplit(bft_strata$Locale, "_"), function(x) {x[1]}))
suffix <- bft_strata$Year
loc_year <- paste(prefix, suffix, sep = "_")    

# rename neutral dapc
dapc_neutral <- pramx

# generate alternative-pop DAPC elements
neutral_locality <- dapc_neutral

neutral_year <- dapc_neutral
neutral_year$DAPC$grp <- as.factor(bft_strata$Year)

neutral_loc_year <- dapc_neutral
neutral_loc_year$DAPC$grp <- as.factor(loc_year)


#### the plots ####
#params
mycolors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#263141","#a65628","#f781bf","#999999")
par(mfrow = c(3,1))

scatter(
  neutral_locality$DAPC,
  cex = 1.1, 
  legend = TRUE,
  txt.leg=c("Brazil","Brazil St. Peter","FL Keys", "Martinique","Pensacola","Puerto Rico", "South Carolina", "Texas", "Venezuela"),
  main = "Descriminant Analysis of Principal Components",
  cstar = FALSE,
  cellipse = TRUE,
  posi.leg = "topleft", 
  scree.pca = TRUE, 
  posi.pca = "topright", 
  posi.da="bottomright",
  cleg = 0.75, 
  xax = 1, 
  yax = 2, 
  solid = 0.8,
  pch= 1:9,
  bg="#EEEEEE",
  col = mycolors
)

scatter(neutral_year$DAPC, posi.da="bottomright", legend = TRUE, cstar = FALSE, col = mycolors, bg="#EEEEEE", pch=16:19, main = "Neutral By Year")
scatter(neutral_loc_year$DAPC, posi.da="bottomright", legend = TRUE, cstar = FALSE, clabel = FALSE, col = mycolors, bg="#EEEEEE", pch=1:14, main = "Neutral By Locality-Year")


#### Outliers ####
load("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/DAPC/dapc_maf_kinless_outlier.rdata")
load("dapc_maf_kinless_outlier.rdata")
# rename neutral dapc
dapc_outlier <- pramx

# generate alternative-pop DAPC elements
outlier_locality <- dapc_outlier

outlier_year <- dapc_outlier
outlier_year$DAPC$grp <- as.factor(bft_strata$Year)

outlier_loc_year <- dapc_outlier
outlier_loc_year$DAPC$grp <- as.factor(loc_year)

par(mfrow = c(3,1))
scatter(
  outlier_locality$DAPC,
  cex = 1.1, 
  legend = TRUE,
  txt.leg=c("Brazil","Brazil St. Peter","FL Keys", "Martinique","Pensacola","Puerto Rico", "South Carolina", "Texas", "Venezuela"),
  main = "Descriminant Analysis of Principal Components",
  cstar = FALSE,
  cellipse = TRUE,
  posi.leg = "topleft", 
  scree.pca = TRUE, 
  posi.pca = "topright", 
  posi.da="bottomright",
  cleg = 0.75, 
  xax = 1, 
  yax = 2, 
  solid = 0.8,
  pch= 1:9,
  bg="#EEEEEE",
  col = mycolors
)

#scatter(outlier_locality$DAPC, posi.da="bottomright", legend = TRUE, clabel = FALSE, bg="white", pch=17:25, main = "Outlier By Locality")
scatter(outlier_year$DAPC, posi.da="bottomright", legend = TRUE, clabel = TRUE, cstar = FALSE, col = mycolors, bg="#EEEEEE", pch=16:19, main = "Outlier By Year")
scatter(outlier_loc_year$DAPC, posi.da="bottomright", legend = TRUE, clabel = FALSE, cstar = FALSE, col = mycolors, bg="#EEEEEE", pch=1:14,, main = "Outlier By Locality-Year")

save.image("dapc_maf_kinless_both.rdata")

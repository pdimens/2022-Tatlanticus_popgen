### AMOVA ####

library(adegenet)
library(poppr)
library(pegas)
library(ade4)


#setwd("~/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/AMOVA")
setwd(getwd())
infile_outlier <- "../data_maf_01/bft_maf01_kinless_outlier.gen"
infile_neutral <- "../data_maf_01/bft_maf01_kinless_neutral.gen"

bft_neut <- read.genepop(infile_neutral, ncode = 3L)
  
  
#num.loci <- length(locNames(bft))
skipval <- nLoc(bft_neut) + 2
ind.tab <- read.table(infile_neutral, skip=skipval, fill=TRUE) # Read in igenepop file as table, skipping the header, the list of loci, and the first 'POP'
ind.tab <- subset(ind.tab, V1 !='POP') # Remove remaining 'POP' lines (rows) in table
inds <- ind.tab$V1 # Return the column of individual names as list of characters
inds <- gsub(',', '', inds) # Remove the trailing comma
indNames(bft_neut) <- inds # assign new list of individual (genotype) names to genind object
PopNames <- c("BRZ", "BRZ", "KEY","MRT","PNS","PR", "SCA","TX", "VZ")
popNames(bft_neut) <- PopNames

bft_out <- read.genepop(infile_outlier, ncode = 3L)
skipval_out <- nLoc(bft_out) + 2
ind.tab <- read.table(infile_outlier, skip=skipval_out, fill=TRUE) # Read in igenepop file as table, skipping the header, the list of loci, and the first 'POP'
ind.tab <- subset(ind.tab, V1 !='POP') # Remove remaining 'POP' lines (rows) in table
inds <- ind.tab$V1 # Return the column of individual names as list of characters
inds <- gsub(',', '', inds) # Remove the trailing comma
indNames(bft_out) <- inds # assign new list of individual (genotype) names to genind object
PopNames <- c("BRZ", "BRZ_SP", "KEY","MRT","PNS","PR", "SCA","TX", "VZ")
popNames(bft_out) <- PopNames
write.table(inds, "samplenames.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Add Strata and perform AMOVA ####

bft_strata <- read.delim("bft_strata.txt", stringsAsFactors = TRUE)
bft_strata$Locale <- factor(bft_strata$Locale)
bft_strata$Year <- factor(bft_strata$Year)
bft_strata$SampleYearWithin <- factor(bft_strata$SampleYearWithin)
bft_strata$SampleYearBtwn <- factor(bft_strata$SampleYearBtwn)

strata(bft_neut) <- bft_strata
strata(bft_out) <- bft_strata

amova_results_neut <- poppr.amova(
  bft_neut,
  hier = ~Locale/Year,
  clonecorrect = FALSE,
  within = FALSE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 4,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = "pegas",
  nperm = 10000
)

amova_results_out <- poppr.amova(
  bft_out,
  hier = ~Locale/Year,
  clonecorrect = FALSE,
  within = FALSE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 4,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = "pegas",
  nperm = 10000
)

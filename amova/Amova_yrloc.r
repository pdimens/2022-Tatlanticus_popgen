### AMOVA ####

library(adegenet) |> suppressPackageStartupMessages()
library(poppr)
library(pegas)
library(ade4)


#setwd("~/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/AMOVA")
setwd(getwd())
infile_outlier <- "../data/kinless.outlier.gen"
infile_neutral <- "../data/kinless.neutral.gen"
PopNames <- c("BRZ","KEY","MRT","PNS","PR", "SCA","TX", "VZ")

bft_neut <- read.genepop(infile_neutral, ncode = 3L)
popNames(bft_neut) <- PopNames

bft_out <- read.genepop(infile_outlier, ncode = 3L)
popNames(bft_out) <- PopNames

#### Add Strata and perform AMOVA ####
bft_strata <- read.table("bft.strata", stringsAsFactors = TRUE, header = T)
bft_strata$year <- factor(bft_strata$year)

strata(bft_neut) <- bft_strata
strata(bft_out) <- bft_strata

amova_results_neut <- poppr.amova(
  bft_neut,
  hier = ~year/population,
  clonecorrect = FALSE,
  within = TRUE,
  squared = TRUE,
  correction = "quasieuclid",
  algorithm = "farthest_neighbor",
  threads = 4,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = "pegas",
  nperm = 50000
)

amova_results_out <- poppr.amova(
  bft_out,
  hier = ~year/population,
  clonecorrect = FALSE,
  within = TRUE,
  squared = TRUE,
  correction = "quasieuclid",
  algorithm = "farthest_neighbor",
  threads = 20,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = "pegas",
  nperm = 50000
)


save.image("amova_yrloc.rdata")

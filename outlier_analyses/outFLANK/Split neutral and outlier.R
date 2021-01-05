library(adegenet)
library(radiator)
library(graph4lg)
library(hierfstat)
library(plyr)
infile <- "BFT_haplo.gen"

bft <- read.genepop(infile, ncode = 3L)
## snp ##
outlier_idx <- c(204,  441,  472,  842,  884,  933,  988, 1096, 1213, 1253, 1335, 1400, 1616, 1645, 1669, 2009, 2033,
               2149, 2157, 2222, 2246, 2304, 2442, 2731, 2747, 2828, 2867, 2943, 2973, 2990, 3322, 3460, 3508, 3588,
               3620, 3893, 4420, 4422, 4575, 4607, 4660, 4832, 4882, 4905, 5027, 5225, 5226)
## haplo ##
outlier_idx <- c(172,  206,  235,  324,  409,  571,  579,  642,  895,  944, 1003, 1036, 1038)

bft_neutral <- full_dataset[loc=-outlier_idx]
bft_outlier <- bft[loc=outlier_idx]
#genomic_converter(bft_neutral,filename = "BFT_nohaplo_neutral.gen", output = "genepop")
genind_to_genepop(bft_outlier, output = "./BFT_haplo_outlier.txt")
genind_to_genepop(bft_neutral, output = "./BFT_haplo_neutral.txt")

HE_summ <- basic.stats(bft)
he_pop <- HE_summ$Ho
he_pop <- as.data.frame(he_pop)
col.names(he_pop) <- c("BRZ", "BRZ_SP", "KEY", "MRT", "PNS","PR","SCA", "TX","VZ")
He <- colwise(mean)(as.data.frame(he_pop))

genomic_converter(bft_neutral, output = "structure", filename = "BFT_haplo_neutral", 
                  filter.common.markers = FALSE, filter.monomorphic = FALSE)

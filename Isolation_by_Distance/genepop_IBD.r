library(genepop)

setwd(getwd())

infile <- "../data_maf_01/bft_maf01_kinless_neutral_IDB.gen"

ibd(
    infile,
    outputFile = "neurtral.txt",
    #settingsFile = "",
    dataType = "Diploid",
    #statistic = "F/(1-F)", #for group data
    statistic = "a",      # for individual data
    geographicScale = "2D",
    CIcoverage = 0.95,
    testPoint = 0,
    minimalDistance = 1e-04,
    maximalDistance = 1e+09,
    mantelPermutations = 1000,
    mantelRankTest = FALSE,
    verbose = interactive()
)
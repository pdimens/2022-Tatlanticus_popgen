library(coda)
library(adegenet)
setwd("/mnt/Win10/Users/pdime/Documents/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/outlier_analyses/bayescan")
infile <- "kinless.bscan_fst.txt"
bft_bscan_fst <- read.csv(infile, sep="")
genfile <- "../../data/kinless.gen"
bft <- read.genepop(genfile, ncode = 3L)

# convergence checks
chain <- read.table("kinless.bscan.sel", header = TRUE)
chain <- chain[-c(1)]
mc_chain <- mcmc(chain, thin = 10)
plot(mc_chain)
summary(mc_chain)
autocorr.diag(mc_chain)
effectiveSize(mc_chain)
geweke.diag(mc_chain, frac1 = 0.1, frac2 = 0.5)
heidel.diag(mc_chain, eps=0.1, pvalue = 0.05)


plt <- plot_bayescan(infile,0,FDR=0.05)
out_ <- bft_bscan_fst[bft_bscan_fst$qval<=0.05,]

## haplo ##
bayescan_idx <- plt$outliers

bayescan_out <- locNames(bft)[bayescan_idx]
write.csv(bayescan_out, "identified_outliers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("bayescan.convergence.rdata")

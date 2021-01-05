library(adegenet)
library(ade4)
setwd(getwd())
infile <- "../data_maf_01/bft_maf01_kinless_neutral.gen"
#bft <- read.genepop("../data_maf_01/bft_maf01_kinless_neutral_IDB.gen", ncode = 3L)
bft <- read.genepop(infile, ncode = 3L)
names(coords) 
## TODO need to reformat file to have POP after each indivset
bft_jitt_loc <- read.table("sans_markers.txt")[1:2]
#bft_jitt_loc <- read.table("~/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/Isolation by Distance/bft_jitt_loc.txt", quote="\"", comment.char="")
bft_jitt_loc <- as.matrix(bft_jitt_loc)
rownames(bft_jitt_loc) <- indNames(bft)
colnames(bft_jitt_loc) <- c("x", "y")
bft$other$xy <- as.matrix(bft_jitt_loc)
## Replace popnames with indnames to make each individual its own population
pop(bft) <- indNames(bft)
popNames(bft)

### Try doing the IBD
toto <- genind2genpop(bft)
Dgen <- dist.genpop(toto, method = 2)
Dgeo <- dist(bft$other$xy)
ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 10000)
ibd

dim(as.matrix(Dgen))
dim(as.matrix(Dgeo))

write.csv(as.matrix(Dgen), file = "distancematrix.csv")
write.csv(as.matrix(Dgeo), file = "distancematrix_dgeo.csv")


##### By population
#bft_pop <- read.genepop("~/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/Isolation by Distance/BFT.maxmiss80.gen", ncode = 3L)
bft_pop <- read.genepop(infile, ncode = 3L)

longitudes <- c(-6.438, 0.0917, 24.767, 14.586, 30.066, 18.915, 32.554, 27.686, 11.491)
latitudes <- c(-34.225, -29.346, -81.005, -61.172, -87.214, -66.536, -78.995, -95.576, -66.229)
xy <- as.matrix(data.frame(x = longitudes, y = latitudes))
popnames = c("BRZ", "BRZ_SP", "KEY", "MRT", "PNS", "PR", "SCA", "TX", "VZ")
popNames(bft_pop) <- popnames
rownames(xy) <- popNames(bft_pop)
colnames(xy) <- c("x", "y")
bft_pop$other$xy <- as.matrix(xy)

### Try doing the IBD
toto_pop <- genind2genpop(bft_pop)
Dgen_pop <- dist.genpop(toto_pop, method = 2)
Dgeo_pop <- dist(bft_pop$other$xy)
ibd_pop <- mantel.randtest(Dgen_pop, Dgeo_pop, nrepet = 10000)
ibd_pop

write.csv(as.matrix(Dgen_pop), file = "distancematrix_pop.csv")
write.csv(as.matrix(Dgeo_pop), file = "distancematrix_pop_dgeo.csv")

par(mfrow= c(1,2))
pdf("IBD_plots.pdf")
plot(ibd, main = "IBD by individual")
plot(ibd_pop, main = "IBD by locality")
dev.off()

save.image("IBD.rdata")
q()
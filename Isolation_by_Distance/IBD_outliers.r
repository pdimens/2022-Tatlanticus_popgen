library(adegenet)
library(ade4)

setwd(getwd())

infile <- "../data_maf_01/bft_maf01_kinless_outlier.gen"
bft <- read.genepop(infile, ncode = 3L)

bft_jitt_loc <- read.table("sans_markers.txt")[1:2]

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

write.csv(as.matrix(Dgen), file = "distancematrix_outlier.csv")
write.csv(as.matrix(Dgeo), file = "distancematrix_dgeo_outlier.csv")

##### By population
bft_pop <- read.genepop(infile, ncode = 3L)

longitudes <- c(-6.438, 24.767, 28.565, 14.586, 30.066, 18.915, 32.554, 27.686, 11.491)
latitudes <- c(-34.225, -81.005, -89.186, -61.172, -87.214, -66.536, -78.995, -95.576, -66.229)
xy <- as.matrix(data.frame(x = longitudes, y = latitudes))
rownames(xy) <- popNames(bft_pop)
colnames(xy) <- c("x", "y")
bft_pop$other$xy <- as.matrix(xy)
popnames = c("BRZ","KEY", "LA", "MRT", "PNS", "PR", "SCA", "TX", "VZ")
popNames(bft_pop) <- popnames

### Try doing the IBD
toto_pop <- genind2genpop(bft_pop)
Dgen_pop <- dist.genpop(toto_pop, method = 2)
Dgeo_pop <- dist(bft_pop$other$xy)
ibd_pop <- mantel.randtest(Dgen_pop, Dgeo_pop, nrepet = 10000)
ibd_pop

write.csv(as.matrix(Dgen_pop), file = "distancematrix_outlier_pop.csv")
write.csv(as.matrix(Dgeo_pop), file = "distancematrix_outlier_pop_dgeo.csv")

par(mfrow= c(1,2))
pdf("IBD_outlier_plots.pdf")
plot(ibd, main = "IBD by individual")
plot(ibd_pop, main = "IBD by locality")
dev.off()

save.image("IBD_outlier.rdata")
q()
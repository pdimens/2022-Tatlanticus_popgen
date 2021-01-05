#! /usr/bin/env Rscript

library("adegenet")
setwd(getwd())

##### Load in data and create backup object

infile <- "BFT_snps_kinless.gen"
bft <- read.genepop(infile, ncode = 3L)
#backup <- bft

#### Add pop names to formatting
num.loci <- length(locNames(bft))
ind.tab <- read.table(infile, skip=num.loci+2, fill=TRUE) # Read in igenepop file as table, skipping the header, the list of loci, and the first 'POP'
ind.tab <- subset(ind.tab, V1 !='POP') # Remove remaining 'POP' lines (rows) in table
inds <- ind.tab$V1 # Return the column of individual names as list of characters
inds <- gsub(',', '', inds) # Remove the trailing comma
indNames(bft) <- inds # assign new list of individual (genotype) names to genind object
PopNames <- c("BRZ","BRZ_SP", "KEY","MRT","PNS","PR", "SCA","TX", "VZ")
popNames(bft) <- PopNames
#pop(bft)

##### Cross validation
set.seed(6969)
pdf("outlier_cross_validation.pdf") 
pramx <- xvalDapc(tab(bft, NA.method="mean"), pop(bft), n.pca=1:100, n.rep=200)  #45 pc + 8 DC
scatter(pramx$DAPC, posi.da="bottomright", bg="white", pch=17:25)

pramx2 <- xvalDapc(tab(bft, NA.method="mean"), pop(bft), n.pca=150:200, n.rep=200) #164 + 8
scatter(pramx2$DAPC, posi.da="bottomright", bg="white", pch=17:25)
save.image("dapc_maf_kinless_outlier.rdata")
q()

#### Find K Clusters
bic<-find.clusters(
        bft, 
        clust=NULL, 
        n.pca=45,
        n.clust=NULL, 
        stat=c("BIC"),
        choose.n.clust=FALSE, 
        criterion=c("diffNgroup"),
        max.n.clust=12, 
        n.iter=1000000,
        n.start=1,
        scale=FALSE, 
        pca.select=c("nbEig"),
        perc.pca=NULL,
        glPca=NULL
)
  
bic.2<-as.matrix(bic$Kstat)
bic.2
save.image("dapc_snp.rdata")
q()
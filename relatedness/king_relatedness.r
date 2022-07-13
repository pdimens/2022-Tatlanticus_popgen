library(SNPRelate) |> suppressPackageStartupMessages()
library(GWASTools) |> suppressPackageStartupMessages()
library(GENESIS) |> suppressPackageStartupMessages()
library(ggplot2) |> suppressPackageStartupMessages()
library(dplyr) |> suppressPackageStartupMessages()
library(SeqVarTools) |> suppressPackageStartupMessages()

setwd("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/Relatedness/")
#setwd(getwd())

vcf_file <- "../data_maf_01/extramissfilter/extramiss.recode.vcf"
gds_file <- "../data_maf_01/extramissfilter/extramiss.recode.gds"
SNPRelate::snpgdsVCF2GDS(vcf.fn = vcf_file, 
                         out.fn = gds_file,
                         method = 'biallelic.only',
                         verbose = TRUE)
#

gdsfile <- "/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/Relatedness/KING/bft_kingpcrelate.gds"
gdsfile <- "../data/extramiss.recode.gds"

gdsobj <- snpgdsOpen(gdsfile)

# pruning
#snpset <- snpgdsLDpruning(gdsobj, method="corr", autosome.only = FALSE, slide.max.bp=10e6, ld.threshold=sqrt(0.1))
#pruned <- unlist(snpset)

# KING PC-based kinship
king <- snpgdsIBDKING(gdsobj, autosome.only = FALSE) #, snp.id=pruned)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
kinship <- snpgdsIBDSelection(king)

ggplot(kinship, aes(IBS0, kinship)) +
    geom_point(alpha=0.5, color = "dodgerblue") +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    ylab("kinship estimate")

ggsave("king_kinship.png")

kin_filter <- kinship %>% 
                filter(kinship > 0)


# GENESIS PC-AIR
# based on KING values
# thresholds set for "unrelated is less than first cousins"
sampset <- pcairPartition(
    kinobj=kingMat, kin.thresh=2^(-9/2),
    divobj=kingMat, div.thresh=-2^(-9/2)
)

# pca on unrelated individuals
pca.unrel <- snpgdsPCA(gdsobj, sample.id=sampset$unrels, autosome.only = FALSE) #, snp.id = pruned)
# project values for relatives
snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj=gdsobj)
samp.load <- snpgdsPCASampLoading(snp.load, gdsobj=gdsobj, sample.id=sampset$rels)

# combine unrelated and related PCs and order as in GDS file
pcs <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(pcs) <- c(pca.unrel$sample.id, samp.load$sample.id)
sample.id <- rownames(kingMat)
samp.ord <- match(sample.id, rownames(pcs))
pcs <- pcs[samp.ord,]

# add population information
annot <- read.table("bft_strata.txt", header = TRUE)
names(annot) <- c("sample.id", "population", "year")
#annot$population <- gsub("\\_\\d+", "", annot$sample.id)
pc.df <- as.data.frame(pcs)
names(pc.df) <- 1:ncol(pcs)
pc.df$sample.id <- row.names(pcs)
pc.df <- left_join(pc.df, annot, by="sample.id")

library(GGally)
library(RColorBrewer)
popn <- length(unique(annot$population))
pop.cols <- RColorBrewer::brewer.pal(n = 4, name = "Set2")
pop.cols <- setNames(brewer.pal(popn, "Paired"), unique(annot$population))
ggparcoord(pc.df, columns=1:12, groupColumn="population", scale="uniminmax") +
    scale_color_manual(values=pop.cols) +
    xlab("PC") + ylab("")

ggsave("KING/extramiss/parallel_coords.png")


# Finally, pcrelate
showfile.gds(closeall=TRUE)
geno <- GdsGenotypeReader(filename = gdsfile)
genodata <- GenotypeData(geno)
genodata <- GenotypeBlockIterator(genodata)

pcrel <- pcrelate(genodata, pcs=pcs[,1:4], training.set=sampset$unrels, 
                  sample.include=sample.id) # deprecated, snp.include=pruned)

pcrelMat <- pcrelateToMatrix(pcrel, scaleKin = 1)

pca <- pcair(genodata,
             kinobj=pcrelMat, kin.thresh=2^(-9/2),
             divobj=kingMat, div.thresh=-2^(-9/2),
             sample.include=sample.id, 
             autosome.only = FALSE)


pcs <- pca$vectors
pc.df <- as.data.frame(pcs)
names(pc.df) <- paste0("PC", 1:ncol(pcs))
pc.df$sample.id <- row.names(pcs)
pc.df <- left_join(pc.df, annot, by="sample.id")

ggplot(pc.df, aes(PC1, PC2, color=population)) + geom_point() +
    scale_color_manual(values=pop.cols)

ggsave("post_pcrelate.png")

# perform pcrelate again
pcrel <- pcrelate(genodata, pcs=pcs[,1:4], training.set=pca$unrels, 
                  sample.include=sample.id)

kinship <- pcrel$kinBtwn

pcrelate_filter <- kinship %>% 
    mutate(
        relationship = case_when(
            (kin <= 0.1767) ~ "unrelated",
            (kin >=0.1767 & kin < 0.3535) ~ "halfsib",
            (kin >= 0.3535) ~ "fullsib"
        )
    )

write.table(pcrelate_filter, file = "pcrelate_output.txt", row.names = F, quote = F)

mycolors <- c("#8a556e", "#f4cf30", "#bbbbbb")
ggplot(pcrelate_filter, aes(k0, kin)) +
    geom_point(alpha=0.5, aes(color = relationship)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    scale_color_manual(values = mycolors) +
    labs(title = "PCRelate kinship estimates sample pairs") +
    ylab("kinship estimate") +
    xlab("probability of sharing 0 alleles")

ggsave("KING/extramiss/kinship_pcarelate.png")


kintable <- pcrelate_filter %>% 
                filter(relationship != "unrelated") %>% 
                arrange(desc(kin))

write.table(kintable, file = "KING/extramiss/kin_pcrelate.txt", row.names = FALSE, quote = FALSE)

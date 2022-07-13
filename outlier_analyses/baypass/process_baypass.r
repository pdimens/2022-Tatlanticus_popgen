library(corrplot)
library(ape)
source("baypass_utils.R") 
RColorBrewer::brewer.pal(3, "Set2")

# Read omega matrix BayPass output
bft.omega=as.matrix(read.table("bft_mat_omega.out"))

popnames <- c("BRZ", "KEY","MRT","PNS","PR","SCA","TX","VZ")
# Identify pool names for, plotting
colnames(bft.omega) <- popnames
rownames(bft.omega) <- popnames

# Create a correlation matrix of the omega values, which can be used to assess genomic differentiation between pools
cor.mat <- cov2cor(bft.omega)
corrplot(
  cor.mat,method="color",mar=c(2,1,2,2)+0.1,
  main=expression("Correlation map based on"~hat(Omega))
)

#plot as hierarchical tree
bft.tree <- as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(
  bft.tree,type="p",
  main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")")
)

# Read the xtx BayPass output
snp.res <- read.table("bft_summary_pi_xtx.out", h=T)

# Get the Pi Beta distribution for POD generation
pi.beta.coef <- read.table("bft_summary_beta_params.out", h=T)$Mean

# Upload original data to get read counts
bft.data <- geno2YN("../../data/bft.kinless.baypass")

# create simulated data
bft.sims <- simulate.baypass(
  omega.mat=bft.omega,nsnp=10000,sample.size=bft.data$NN,
  beta.pi=pi.beta.coef,pi.maf=0,suffix="bft.BP.sim"
)


# Read the new omega matrix and compare to original. Here you want similar values between the two
pod.omega=as.matrix(read.table("sims_mat_omega.out"))
plot(pod.omega,bft.omega) ; abline(a=0,b=1)

# Get the Forstner and Moonen Distance (FMD) between simulated and original posterior estimates (here a smaller value is better) 
fmd.dist(pod.omega,bft.omega)


# Look at POD xtx values, and identify SNPs where the xtx values are above the 99% significance threshold from the POD. So in the plot, it is those loci (dots) which are above the abline
pod.xtx=read.table("sims_summary_pi_xtx.out",h=T)$M_XtX
pod.thresh=quantile(pod.xtx,probs=0.99)
plot(snp.res$M_XtX)
abline(h=pod.thresh,lty=2)

# Obtain the xtx threshold. Any SNPs with an xtx value greater than thjuis are identified as outlier SNPs
pod.thresh

snp.scores <- data.frame(snp_idx = 1:length(snp.res$M_XtX), M_XtX = snp.res$M_XtX, outlier = snp.res$M_XtX > pod.thresh)
snp.scores$outlier <- snp.scores$`snp.res$M_XtX` > pod.thresh

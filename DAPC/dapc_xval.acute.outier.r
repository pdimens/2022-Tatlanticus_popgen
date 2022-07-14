#! /usr/bin/env Rscript

library(adegenet) 
library(dplyr)
library(tidyr)
library(ggpubr)

#setwd("/mnt/Win10/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/data_maf_01/extramissfilter")
setwd(getwd())
##### Load in data and create backup object


infile <- "../../data/kinless.outlier.gen"

bft <- read.genepop(infile, ncode = 3L)
metadata <- read.csv("bft.strata.csv", sep = ",", header = T)
#backup <- bftS

#### Add pop names to formatting
PopNames <- c("BRZ","KEY","MRT","PNS","PR", "SCA","TX", "VZ")
popNames(bft) <- PopNames
original_pops <- metadata$pop4
pop(bft) <- metadata$pop3

##### Cross validation
set.seed(6969)
#pramx <- xvalDapc(tab(bft, NA.method="mean"), pop(bft), n.pca=95:110, n.rep=200,  parallel = "multicore", ncpus = 25)  #15 pc + 8 DC
#pramx <- xvalDapc(tab(bft, NA.method="mean"), pop(bft), n.pca.max=200, n.rep=100,  parallel = "multicore", ncpus = 4)  #15 pc + 8 DC
#pramx <- dapc(bft, NA.method="mean", pop(bft), n.pca=90, n.rep=200,  parallel = "multicore", ncpus = 4)  #15 pc + 8 DC

#save.image("neutral.xvaldapc.rdata")
#q()



pcval <- 15
# outlier
#pcval <- 40
daval <- 3

maxK <- 9
myMat <- matrix(nrow=50, ncol=maxK)
colnames(myMat) <- 1:maxK
for(i in 1:nrow(myMat)){
  grp <- find.clusters(bft, n.pca = pcval, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

neutral <- dapc(bft, pop = bft$pop, n.pca = pcval, n.da = daval)
#outlier
#neutral <- dapc(bft, pop = bft$pop, n.pca = 40, n.da = 2)

#my_pal <- RColorBrewer::brewer.pal(n=9, name = "Set1")
my_pal <- c("#f1c232", "#d08120", "#54438a", "#655595", "#2a2145","#66b2b2", "#008080", "#006666", "#ff8181")


dapc_l <- neutral
my_df <- as.data.frame(dapc_l$ind.coord)
my_df$Group <- dapc_l$grp

df <- as.data.frame(neutral$ind.coord)
df <- cbind(rownames(df), bft$pop, df)
names(df) <- c("name", "population", "LD1", "LD2", "LD3") #, "LD4")
df$reclass <- metadata$pop3
df$reclass <- as.factor(df$reclass)
df <- df %>% mutate(population = factor(original_pops, levels = c("BRZ", "BRZSP", "MRT", "VZ", "PR", "TX", "PNS", "KEY", "SCA"), ordered = TRUE))
df <- df %>% group_by(population) %>% 
  mutate(centroid1 = mean(LD1), centroid2 = mean(LD2))

tmp <- as.data.frame(dapc_l$posterior)
tmp <- cbind(indNames(bft), bft$pop, tmp) #drop redundant col
names(tmp) <- c("name", "origin", "BRZ", "GULF", "CRB", "ATL")
tmp$origin <- factor(original_pops, levels = c("BRZ", "BRZSP", "MRT", "VZ", "PR", "TX", "PNS", "KEY", "SCA"), ordered = TRUE)
assignments <- pivot_longer(tmp, c(-name, -origin), names_to = "population", values_to = "posterior")
#assignments$origin <- factor(gsub("_\\d+", "", assignments$name), levels = c("BRZ", "MRT", "VZ", "PR", "TX", "PNS", "KEY", "SCA"), ordered = TRUE)


#### the plots ####

k_df <- reshape2::melt(myMat)
colnames(k_df)[1:3] <- c("Group", "K", "BIC")
k_df$K <- as.factor(k_df$K)

k_plot <- k_df %>% 
  ggplot(aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Bayesian Information Criterion (BIC)") +
  xlab("Number of groups (K)")

dapc_plot <- df %>% 
  ggplot(aes(x = LD1, y = LD2, shape = reclass)) +
  geom_point(size = 2, alpha = 0.7, aes(color = population, fill = population)) +
  stat_ellipse(linetype = 2, size = 0.5, show.legend = FALSE, type = "norm", level = 0.75, segments = 100) +
  geom_segment(data = df, aes(x = centroid1, y = centroid2, xend = LD1, yend = LD2, color = population), alpha = 0.75) +
  xlab("Linear Discriminant 1") +
  ylab("Linear Discriminant 2") +
  scale_color_manual(values=c(my_pal), name = "Location") +
  scale_shape_manual(values = c(16,17,15,18)) +
  scale_fill_manual(values=c(paste(my_pal, "96", sep = "")), name = "Location") +
  guides(shape=guide_legend(title="Population")) +
  theme_bw()

posterior_plot <- assignments %>% 
  ggplot(aes(x = name, y = posterior, fill = population)) +
  geom_bar(stat = "identity", width = 1.0, alpha = 0.7) +
  scale_fill_manual(values = my_pal[c(9,1,4,6)]) +
  #scale_fill_manual(values=c(paste(my_pal, "98", sep = ""))) +
  ylab("Posterior Membership Probability") +
  xlab("Samples") +
  facet_grid(~origin, scales = "free_x", space = "free" ) + 
  theme_classic() +
  guides(fill=guide_legend(title="Membership")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), axis.line = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off")

#### combining the plots ####
ggarrange(
  ggarrange(
    k_plot,
    dapc_plot,
    ncol = 2, labels = c("A", "B")
  ),
  posterior_plot,
  nrow = 2,
  labels = c("", "C"),
  heights = c(2, 1.3)
)

ggsave("DAPC_missingrm.4pop.pc15.png", height = 7, width = 11.5, units = "in")
ggsave("DAPC_missingrm.pc40.outliers.png", height = 7, width = 11, units = "in")


library(plotly)


locdf <- as.data.frame(neutral$ind.coord)
locdf$population <- as.factor(bft$pop)  

fig <- plot_ly(locdf, x = ~LD1, y = ~LD2, z = ~LD3, color = ~population, colors = my_pal[c(1,5,3,8)], size = 15)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'LD1'),
                                   yaxis = list(title = 'LD2'),
                                   zaxis = list(title = 'LD3')))

fig
htmlwidgets::saveWidget(as_widget(fig), paste0("dapc_", "outlier", ".html"))


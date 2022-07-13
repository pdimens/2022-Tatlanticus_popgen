library(network)
library(sna)
library(igraph)
library(intergraph)
library(ggnet)
library(ggplot2)

tmp <- read.table("C:/Users/Pavel/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/extramiss/relatedness/kin_pcrelate.txt", header = T)
kingraph <- graph.edgelist(as.matrix(tmp[,c(1,2)], directed=FALSE))
relatn <- tmp$relationship
relatn[relatn == "fullsib"] <- 1
relatn[relatn == "halfsib"] <- 2
relatn <- as.numeric(relatn)
E(kingraph)$weight <- relatn
kingraph <- as.undirected(kingraph)
mtx <- as_adjacency_matrix(kingraph, attr = "weight") %>% as.matrix
net <- as.network.matrix(mtx, matrix.type = "adjacency", directed = FALSE)
nodes <- as.edgelist(net)
net %v% "relationship" <- tmp$relationship
pops <- c()
for(i in 1:23){
  pops <- c(pops, net$val[i][[1]]$vertex.names)  
}
pops <- gsub("_\\d+", "", pops)
pops <- gsub("TXL", "TX", pops)
net %v% "population" <- pops
pal <- RColorBrewer::brewer.pal(n=3, "Dark2")
edgepal <- c("grey60", "dodgerblue")

ggnet2(net, 
       color = "population", 
       palette = "Set2", 
       edge.color = edgepal[relatn],
       legend.size = 12, 
       size = 6, 
       edge.size = 0.5, 
       legend.position = "right",
       mode = "kamadakawai",
       color.legend = "location") +
  geom_point(aes(color = color), size = 6, color = "white") +
  geom_point(aes(color = color), size = 6, alpha = 0.5) +
  geom_point(aes(color = color), size = 4)

ggsave("kinship.network.png", height = 1000, width = 1000, units = "px")

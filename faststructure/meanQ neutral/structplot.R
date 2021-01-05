library(ggplot2)
library(tidyr)
library(gridExtra)

setwd(getwd())
k2 <- read.table("bft_maf01_kinless_neutral.faststructure_out.2.meanQ")
k3 <- read.table("bft_maf01_kinless_neutral.faststructure_out.3.meanQ")
k4 <- read.table("bft_maf01_kinless_neutral.faststructure_out.4.meanQ")
k5 <- read.table("bft_maf01_kinless_neutral.faststructure_out.5.meanQ")
k6 <- read.table("bft_maf01_kinless_neutral.faststructure_out.6.meanQ")
k7 <- read.table("bft_maf01_kinless_neutral.faststructure_out.7.meanQ")
k8 <- read.table("bft_maf01_kinless_neutral.faststructure_out.8.meanQ")
k9 <- read.table("bft_maf01_kinless_neutral.faststructure_out.9.meanQ")
bft_names <- as.vector(unlist(read.table("../../../data_maf_01/sample_names.txt")))

# add names column to each
for (i in ls(pattern = "k")) {      
    assign(i,transform(get(i),sample = bft_names))
}

popnames <- as.factor(
    unlist(
        lapply(strsplit(bft_names, "_"), function(x) {x[1]})
    )
)

popcounts <- as.data.frame(summary(popnames))
colnames(popcounts) <- "freq"
# make the values cumulative
for(i in 2:8){
    popcounts$freq[i] <- popcounts$freq[i] + popcounts$freq[i-1]
}

poptext <- c(0,popcounts$freq)
for(i in 1:8){
    poptext[i] <- poptext[i] + (poptext[i+1] -  poptext[i])/2
}
poptext <- poptext[1:8]



# stack each into long format
k2 <- pivot_longer(k2, cols = !(dim(k2)[2]))
colnames(k2) <- c("sample", "k", "probability")
k3 <- pivot_longer(k3, cols = !(dim(k3)[2]))
colnames(k3) <- c("sample", "k", "probability")
k4 <- pivot_longer(k4, cols = !(dim(k4)[2]))
colnames(k4) <- c("sample", "k", "probability")
k5 <- pivot_longer(k5, cols = !(dim(k5)[2]))
colnames(k5) <- c("sample", "k", "probability")
k6 <- pivot_longer(k6, cols = !(dim(k6)[2]))
colnames(k6) <- c("sample", "k", "probability")
k7 <- pivot_longer(k7, cols = !(dim(k7)[2]))
colnames(k7) <- c("sample", "k", "probability")
k8 <- pivot_longer(k8, cols = !(dim(k8)[2]))
colnames(k8) <- c("sample", "k", "probability")
k9 <- pivot_longer(k9, cols = !(dim(k9)[2]))
colnames(k9) <- c("sample", "k", "probability")

plots <- list()
dfs <- c("k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9")
dfs <- dfs[1]
mycolors <- c("#4C566A", "#81afc1")

for (i in 1:length(dfs)) {
   str_plot <- ggplot(get(dfs[i]), aes_string(fill = "k", y = "probability", x = "sample")) +
            scale_fill_manual(values = mycolors) +
            ylab("Probability of Membership") + 
            geom_bar(position = "fill", stat = "identity", width = 1.0) +
            geom_vline(xintercept = popcounts$freq, size = 0.4, color = "#000000") +
            coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off") +
            annotate("text", x = poptext[1], y = 0.1, label = "BRZ") +
            annotate("text", x = poptext[2], y = 0.1, label = "KEY") +
            annotate("text", x = poptext[3], y = 0.1, label = "MRT") +
            annotate("text", x = poptext[4], y = 0.1, label = "PNS") +
            annotate("text", x = poptext[5], y = 0.1, label = "PR") +
            annotate("text", x = poptext[6], y = 0.1, label = "SCA") +
            annotate("text", x = poptext[7], y = 0.1, label = "TX") +
            annotate("text", x = poptext[8], y = 0.1, label = "VZ") +
            theme(
                axis.title = element_text(size=rel(1)),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y =  element_text(size=rel(1)),
                axis.ticks.x = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank(),
                panel.grid.minor = element_blank()
            )
   plots[[i]] <- str_plot
}
do.call(grid.arrange,plots)


#### individual plots ####
k2$k <- as.factor(k2$k)
levels(k2$k) <- c("1","2")
names(k2)[2] <- "population"

ggplot(k2, aes_string(fill = "population", y = "probability", x = "sample")) +
    scale_fill_manual(values = mycolors) +
    ylab("Probability of Membership") + 
    geom_bar(position = "fill", stat = "identity", width = rel(1.0)) +
    geom_vline(xintercept = popcounts$freq, size = rel(0.4), color = "#D8DEE9") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off") +
    annotate("text", x = poptext[1], y = 0.1, label = "BRZ") +
    annotate("text", x = poptext[2], y = 0.1, label = "KEY") +
    annotate("text", x = poptext[3], y = 0.1, label = "MRT") +
    annotate("text", x = poptext[4], y = 0.1, label = "PNS") +
    annotate("text", x = poptext[5], y = 0.1, label = "PR") +
    annotate("text", x = poptext[6], y = 0.1, label = "SCA") +
    annotate("text", x = poptext[7], y = 0.1, label = "TX") +
    annotate("text", x = poptext[8], y = 0.1, label = "VZ") +
    theme(
        axis.title = element_text(size=rel(1.1)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=rel(1.1)),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=rel(1)),
        legend.text = element_text(size=rel(1))
    )
